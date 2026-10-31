#include "families.h"
#include "group_agg.h"
#include "threads.h"
#include <RcppParallel.h>
#include <map>
using namespace Rcpp;
using namespace RcppParallel;
using panglm::compute_group_agg;

namespace {

// Per-observation terms of the beta-negative-binomial likelihood:
//   t   = lgamma(l + y) - lgamma(l) - lgamma(y + 1)
//   dlg = digamma(l + y) - digamma(l)
//   tdf = trigamma(l + y) - trigamma(l)
// For a count y the differences are finite sums over j = 0, ..., y - 1 of
// log(l + j), 1 / (l + j) and -1 / (l + j)^2 (the recurrences of the gamma,
// digamma and trigamma functions), which are much cheaper than the special
// functions for the small counts of typical data and are used up to y = 50.
// With y = 0 all three are zero. The terms are stored per observation and
// summed serially by the caller, so the result does not depend on the
// number of threads.
struct NBObsWorker : public Worker {
  const double *l, *y, *lgy;
  double *t, *dlg, *tdf;
  NBObsWorker(const double* l, const double* y, const double* lgy,
              double* t, double* dlg, double* tdf)
    : l(l), y(y), lgy(lgy), t(t), dlg(dlg), tdf(tdf) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t i = begin; i < end; ++i) {
      const double yi = y[i];
      if (yi == 0.0) {
        t[i] = 0.0 - lgy[i];
        dlg[i] = 0.0;
        tdf[i] = 0.0;
      } else if (yi <= 50.0 && yi == std::floor(yi)) {
        const int m = (int) yi;
        double sl = 0.0, s1 = 0.0, s2 = 0.0;
        for (int j = 0; j < m; ++j) {
          const double x = l[i] + j, r = 1.0 / x;
          sl += std::log(x);
          s1 += r;
          s2 += r * r;
        }
        t[i] = sl - lgy[i];
        dlg[i] = s1;
        tdf[i] = -s2;
      } else {
        const double ly = l[i] + y[i];
        t[i] = R::lgammafn(ly) - R::lgammafn(l[i]) - lgy[i];
        dlg[i] = R::digamma(ly) - R::digamma(l[i]);
        tdf[i] = R::trigamma(ly) - R::trigamma(l[i]);
      }
    }
  }
};

// Per-group special functions of a + L_i and a + b + L_i + Y_i
struct NBGroupWorker : public Worker {
  const double *Li, *Yi;
  double a, b;
  double *lg_aL, *lg_abLY, *dg_aL, *dg_abLY, *tg_aL, *tg_abLY;
  NBGroupWorker(const double* Li, const double* Yi, double a, double b,
                double* lg_aL, double* lg_abLY, double* dg_aL, double* dg_abLY,
                double* tg_aL, double* tg_abLY)
    : Li(Li), Yi(Yi), a(a), b(b), lg_aL(lg_aL), lg_abLY(lg_abLY), dg_aL(dg_aL),
      dg_abLY(dg_abLY), tg_aL(tg_aL), tg_abLY(tg_abLY) {}
  void operator()(std::size_t begin, std::size_t end) {
    for (std::size_t g = begin; g < end; ++g) {
      const double aL = a + Li[g], abLY = a + b + Li[g] + Yi[g];
      lg_aL[g] = R::lgammafn(aL);
      lg_abLY[g] = R::lgammafn(abLY);
      dg_aL[g] = R::digamma(aL);
      dg_abLY[g] = R::digamma(abLY);
      tg_aL[g] = R::trigamma(aL);
      tg_abLY[g] = R::trigamma(abLY);
    }
  }
};

} // namespace

// Random-effects negative binomial (beta-negative-binomial marginal),
// matching pglm's lnl.negbin.R "random" branch exactly: two free dispersion
// parameters a, b (vs. the single Gamma-shape parameter of the Poisson-Gamma
// random model), a genuinely different mixing distribution from
// random_poisson_fit_cpp -- not an alias.
//
// a and b are optimized via kappa_a = log(a), kappa_b = log(b), which keeps
// them positive automatically and avoids the ad hoc boundary clipping that
// makes plain-scale Newton steps fragile on sparse/skewed panels (the same
// fix applied to random_poisson_fit_cpp).
//
// [[Rcpp::export]]
List random_negbin_fit_cpp(const arma::mat& X, const arma::vec& y,
                           IntegerVector group_start, IntegerVector group_size,
                           int maxit = 100, double tol = 1e-10) {
  int n = X.n_rows, k = X.n_cols, G = group_start.size();

  arma::vec group_id(n);
  for (int g = 0; g < G; ++g)
    for (int t = 0; t < group_size[g]; ++t)
      group_id[group_start[g] + t] = g;

  arma::vec Yi(G, arma::fill::zeros);
  for (int i = 0; i < n; ++i) Yi[(int)group_id[i]] += y[i];

  // lgamma(y + 1) does not depend on the parameters
  arma::vec lgy(n);
  for (int i = 0; i < n; ++i) lgy[i] = R::lgammafn(y[i] + 1.0);

  // The group totals Y_i take few distinct values: the functions of b + Y_i
  // are computed once per distinct value.
  std::vector<double> yu;
  std::vector<int> yidx(G);
  {
    std::map<double, int> pos;
    for (int g = 0; g < G; ++g) {
      auto it = pos.find(Yi[g]);
      if (it == pos.end()) {
        it = pos.emplace(Yi[g], (int) yu.size()).first;
        yu.push_back(Yi[g]);
      }
      yidx[g] = it->second;
    }
  }
  const int U = yu.size();

  // theta = (beta_1..beta_k, kappa_a, kappa_b); a = exp(kappa_a), b = exp(kappa_b)
  arma::vec theta(k + 2, arma::fill::zeros);
  theta[k] = 0.0;   // a = 1
  theta[k + 1] = 0.0; // b = 1

  auto eval = [&](const arma::vec& th, arma::vec& grad, arma::mat& hess, double& ll) {
    arma::vec b_coef = th.subvec(0, k - 1);
    double a = std::exp(th[k]), b = std::exp(th[k + 1]);

    arma::vec eta = arma::clamp(X * b_coef, -30.0, 30.0);
    arma::vec l = arma::exp(eta);

    arma::vec Li; arma::mat lXi;
    compute_group_agg(X, l, group_start, group_size, Li, lXi);

    arma::vec tA(n), dlg(n), tdf(n);
    NBObsWorker ow(l.memptr(), y.memptr(), lgy.memptr(), tA.memptr(), dlg.memptr(), tdf.memptr());
    parallelFor(0, n, ow, panglm::grain(n));
    double lnA = 0.0;
    for (int i = 0; i < n; ++i) lnA += tA[i];

    arma::vec lg_aL(G), lg_abLY(G), dg_apLi(G), dg_apbpLiYi(G), tg_apLi(G), tg_apbpLiYi(G);
    NBGroupWorker gw(Li.memptr(), Yi.memptr(), a, b, lg_aL.memptr(), lg_abLY.memptr(),
                     dg_apLi.memptr(), dg_apbpLiYi.memptr(), tg_apLi.memptr(),
                     tg_apbpLiYi.memptr());
    parallelFor(0, G, gw, panglm::grain(G));
    std::vector<double> lg_bYu(U), dg_bYu(U), tg_bYu(U);
    for (int u = 0; u < U; ++u) {
      lg_bYu[u] = R::lgammafn(b + yu[u]);
      dg_bYu[u] = R::digamma(b + yu[u]);
      tg_bYu[u] = R::trigamma(b + yu[u]);
    }
    const double lg_ab = R::lgammafn(a + b), lg_a = R::lgammafn(a), lg_b = R::lgammafn(b);
    double lnC = 0.0;
    arma::vec dg_bpYi(G), tg_bpYi(G);
    for (int g = 0; g < G; ++g) {
      lnC += lg_ab + lg_aL[g] + lg_bYu[yidx[g]] - lg_a - lg_b - lg_abLY[g];
      dg_bpYi[g] = dg_bYu[yidx[g]];
      tg_bpYi[g] = tg_bYu[yidx[g]];
    }
    ll = lnA + lnC;

    // beta gradient
    arma::vec gradi(n);
    for (int i = 0; i < n; ++i) {
      int g = (int)group_id[i];
      gradi[i] = dlg[i] * l[i] + (dg_apLi[g] - dg_apbpLiYi[g]) * l[i];
    }
    arma::vec grad_beta = X.t() * gradi;

    double dg_apb = R::digamma(a + b), dg_a = R::digamma(a), dg_b = R::digamma(b);
    double grad_a = 0.0, grad_b = 0.0;
    for (int g = 0; g < G; ++g) {
      grad_a += dg_apb + dg_apLi[g] - dg_a - dg_apbpLiYi[g];
      grad_b += dg_apb + dg_bpYi[g] - dg_b - dg_apbpLiYi[g];
    }

    grad.set_size(k + 2);
    grad.subvec(0, k - 1) = grad_beta;
    grad[k] = grad_a * a;       // chain rule: d/dkappa_a = d/da * a
    grad[k + 1] = grad_b * b;

    // beta-beta Hessian
    arma::vec lnA_bb(n);
    for (int i = 0; i < n; ++i) {
      double v = dlg[i] * l[i] + tdf[i] * l[i] * l[i];
      lnA_bb[i] = (std::isfinite(v) && v > 0) ? v : 0.0;
    }
    arma::mat Xw = X.each_col() % arma::sqrt(lnA_bb);
    arma::mat H_A = Xw.t() * Xw;

    arma::vec weight1(n);
    for (int i = 0; i < n; ++i) {
      int g = (int)group_id[i];
      weight1[i] = -(dg_apLi[g] - dg_apbpLiYi[g]) * l[i];
      if (weight1[i] < 0) weight1[i] = 0;
    }
    arma::mat Xw1 = X.each_col() % arma::sqrt(weight1);
    arma::mat H_C1 = -(Xw1.t() * Xw1);

    arma::vec scale2 = arma::sqrt(arma::clamp(tg_apLi - tg_apbpLiYi, 0.0, arma::datum::inf));
    arma::mat lXi_scaled = lXi.each_col() % scale2;
    arma::mat H_C2 = lXi_scaled.t() * lXi_scaled;

    arma::mat Hbb = H_A + H_C1 + H_C2;

    // beta-a, beta-b cross terms (plain scale, chain-ruled below)
    arma::vec coef_a = tg_apLi - tg_apbpLiYi;      // per group
    arma::vec coef_b = -tg_apbpLiYi;               // per group
    arma::mat lXi_a = lXi.each_col() % coef_a;
    arma::mat lXi_b = lXi.each_col() % coef_b;
    arma::vec Hba = arma::sum(lXi_a, 0).t();
    arma::vec Hbb_ = arma::sum(lXi_b, 0).t();

    double tg_apb = R::trigamma(a + b), tg_a = R::trigamma(a), tg_b = R::trigamma(b);
    double Haa = 0.0, Hab = 0.0, Hbb2 = 0.0;
    for (int g = 0; g < G; ++g) {
      Haa += tg_apb + tg_apLi[g] - tg_a - tg_apbpLiYi[g];
      Hab += tg_apb - tg_apbpLiYi[g];
      Hbb2 += tg_apb + tg_bpYi[g] - tg_b - tg_apbpLiYi[g];
    }

    // Chain rule to (kappa_a, kappa_b) scale:
    // d2L/dka2 = Haa*a^2 + grad_a*a ; d2L/dka dkb = Hab*a*b ; etc.
    // d2L/dka dbeta = Hba * a
    hess.set_size(k + 2, k + 2);
    hess.submat(0, 0, k - 1, k - 1) = Hbb;
    hess.submat(0, k, k - 1, k) = Hba * a;
    hess.submat(0, k + 1, k - 1, k + 1) = Hbb_ * b;
    hess.submat(k, 0, k, k - 1) = (Hba * a).t();
    hess.submat(k + 1, 0, k + 1, k - 1) = (Hbb_ * b).t();
    hess(k, k) = Haa * a * a + grad_a * a;
    hess(k + 1, k + 1) = Hbb2 * b * b + grad_b * b;
    hess(k, k + 1) = Hab * a * b;
    hess(k + 1, k) = Hab * a * b;
  };

  arma::vec grad; arma::mat hess; double ll, ll_old = -arma::datum::inf;
  int iter = 0;
  // The gradient and Hessian of the accepted point are carried into the
  // next Newton step, so each step costs one evaluation (plus any
  // step-halvings) instead of re-evaluating the point just accepted.
  eval(theta, grad, hess, ll);
  for (iter = 0; iter < maxit; ++iter) {
    arma::vec step = arma::solve(hess, grad, arma::solve_opts::likely_sympd);

    double lambda = 1.0;
    arma::vec theta_new = theta - lambda * step;
    arma::vec g2; arma::mat h2; double ll_new;
    eval(theta_new, g2, h2, ll_new);
    int halvings = 0;
    while (!std::isfinite(ll_new) || (ll_new < ll && halvings < 30)) {
      lambda *= 0.5;
      theta_new = theta - lambda * step;
      eval(theta_new, g2, h2, ll_new);
      halvings++;
    }

    bool converged = std::fabs(ll_new - ll) < tol * (std::fabs(ll) + 1.0);
    theta = theta_new;
    grad = g2;
    hess = h2;
    ll = ll_new;
    ll_old = ll_new;
    if (converged) { iter++; break; }
  }

  arma::mat vcov_kappa = arma::inv_sympd(-hess);

  // Delta method: Var(a) = a^2 * Var(kappa_a), etc.
  arma::vec scale_vec(k + 2, arma::fill::ones);
  scale_vec[k] = std::exp(theta[k]);
  scale_vec[k + 1] = std::exp(theta[k + 1]);
  arma::mat vcov = vcov_kappa % (scale_vec * scale_vec.t());

  return List::create(
    Named("coefficients") = theta.subvec(0, k - 1),
    Named("a") = std::exp(theta[k]),
    Named("b") = std::exp(theta[k + 1]),
    Named("vcov_unscaled") = vcov,
    Named("loglik") = ll_old,
    Named("iterations") = iter
  );
}
