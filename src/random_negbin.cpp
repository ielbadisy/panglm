#include "families.h"
#include "group_agg.h"
using namespace Rcpp;
using panglm::compute_group_agg;

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

    double lnA = 0.0;
    arma::vec dgl_ly(n), dgl_l(n);
    for (int i = 0; i < n; ++i) {
      lnA += R::lgammafn(l[i] + y[i]) - R::lgammafn(l[i]) - R::lgammafn(y[i] + 1.0);
      dgl_ly[i] = R::digamma(l[i] + y[i]);
      dgl_l[i]  = R::digamma(l[i]);
    }

    double lnC = 0.0;
    arma::vec dg_apLi(G), dg_apbpLiYi(G), dg_bpYi(G);
    for (int g = 0; g < G; ++g) {
      lnC += R::lgammafn(a + b) + R::lgammafn(a + Li[g]) + R::lgammafn(b + Yi[g])
           - R::lgammafn(a) - R::lgammafn(b) - R::lgammafn(a + b + Li[g] + Yi[g]);
      dg_apLi[g] = R::digamma(a + Li[g]);
      dg_apbpLiYi[g] = R::digamma(a + b + Li[g] + Yi[g]);
      dg_bpYi[g] = R::digamma(b + Yi[g]);
    }
    ll = lnA + lnC;

    // beta gradient
    arma::vec gradi(n);
    for (int i = 0; i < n; ++i) {
      int g = (int)group_id[i];
      gradi[i] = (dgl_ly[i] - dgl_l[i]) * l[i] + (dg_apLi[g] - dg_apbpLiYi[g]) * l[i];
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
      double v = (dgl_ly[i] - dgl_l[i]) * l[i] +
                 (R::trigamma(l[i] + y[i]) - R::trigamma(l[i])) * l[i] * l[i];
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

    arma::vec tg_apLi(G), tg_apbpLiYi(G), tg_bpYi(G);
    for (int g = 0; g < G; ++g) {
      tg_apLi[g] = R::trigamma(a + Li[g]);
      tg_apbpLiYi[g] = R::trigamma(a + b + Li[g] + Yi[g]);
      tg_bpYi[g] = R::trigamma(b + Yi[g]);
    }
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

    double tg_apb = R::trigamma(a + b);
    double Haa = 0.0, Hab = 0.0, Hbb2 = 0.0;
    for (int g = 0; g < G; ++g) {
      Haa += tg_apb + tg_apLi[g] - R::trigamma(a) - tg_apbpLiYi[g];
      Hab += tg_apb - tg_apbpLiYi[g];
      Hbb2 += tg_apb + tg_bpYi[g] - R::trigamma(b) - tg_apbpLiYi[g];
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
  for (iter = 0; iter < maxit; ++iter) {
    eval(theta, grad, hess, ll);
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
    ll_old = ll_new;
    if (converged) { iter++; break; }
  }

  eval(theta, grad, hess, ll_old);
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
