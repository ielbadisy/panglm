#include "families.h"
#include "group_agg.h"
using namespace Rcpp;
using panglm::compute_group_agg;

// Poisson random-effects (Poisson-Gamma / negative-binomial marginal)
// panel model, ported from pglm's lnl.poisson "random" branch. Jointly
// estimates the K regression coefficients beta and the Gamma dispersion
// parameter d > 0 (frailty variance is 1/d) by Newton-Raphson with step
// halving; per-group sufficient statistics computed in parallel.
//
// [[Rcpp::export]]
List random_poisson_fit_cpp(const arma::mat& X, const arma::vec& y,
                             IntegerVector group_start, IntegerVector group_size,
                             int maxit = 100, double tol = 1e-10) {
  int n = X.n_rows, k = X.n_cols, G = group_start.size();

  arma::vec group_id(n);
  for (int g = 0; g < G; ++g)
    for (int t = 0; t < group_size[g]; ++t)
      group_id[group_start[g] + t] = g;

  arma::vec Yi(G, arma::fill::zeros);
  for (int i = 0; i < n; ++i) Yi[(int)group_id[i]] += y[i];

  arma::vec beta = arma::zeros<arma::vec>(k);
  double d = 1.0; // dispersion, param k

  auto eval = [&](const arma::vec& b, double dd,
                  arma::vec& grad, arma::mat& hess, double& ll) {
    arma::vec eta = X * b;
    eta = arma::clamp(eta, -30.0, 30.0);
    arma::vec lit = arma::exp(eta);

    arma::vec Li; arma::mat lXi;
    compute_group_agg(X, lit, group_start, group_size, Li, lXi);

    arma::vec dYiLi = (dd + Yi) / (dd + Li);      // per-group
    arma::vec w_obs(n), weight1(n);
    for (int i = 0; i < n; ++i) {
      w_obs[i] = dYiLi[(int)group_id[i]];
      weight1[i] = w_obs[i] * lit[i];
    }
    arma::vec gradi = y - w_obs % lit;
    arma::vec grad_beta = X.t() * gradi;

    double grad_d = 0.0;
    for (int g = 0; g < G; ++g) {
      grad_d += 1.0 + std::log(dd) - std::log(dd + Li[g]) - (dd + Yi[g]) / (dd + Li[g])
              + R::digamma(dd + Yi[g]) - R::digamma(dd);
    }

    grad.set_size(k + 1);
    grad.subvec(0, k - 1) = grad_beta;
    grad[k] = grad_d;

    // Hessian beta-beta
    arma::mat Xw = X.each_col() % arma::sqrt(weight1);
    arma::mat H1 = -(Xw.t() * Xw);
    arma::vec scale2 = arma::sqrt((dd + Yi) / arma::square(dd + Li));
    arma::mat lXi_scaled = lXi.each_col() % scale2;
    arma::mat H2 = lXi_scaled.t() * lXi_scaled;
    arma::mat Hbb = H1 + H2;

    // Hessian beta-d
    arma::vec coef_bd = -(Li - Yi) / arma::square(dd + Li); // per group
    arma::mat lXi_bd = lXi.each_col() % coef_bd;
    arma::vec Hbd = arma::sum(lXi_bd, 0).t();

    double Hdd = 0.0;
    for (int g = 0; g < G; ++g) {
      Hdd += 1.0 / dd - 1.0 / (dd + Li[g]) - (Li[g] - Yi[g]) / std::pow(dd + Li[g], 2)
           + R::trigamma(dd + Yi[g]) - R::trigamma(dd);
    }

    hess.set_size(k + 1, k + 1);
    hess.submat(0, 0, k - 1, k - 1) = Hbb;
    hess.submat(0, k, k - 1, k) = Hbd;
    hess.submat(k, 0, k, k - 1) = Hbd.t();
    hess(k, k) = Hdd;

    double lnA = arma::sum(y % eta) - arma::sum(arma::lgamma(y + 1.0));
    double lnC = 0.0;
    for (int g = 0; g < G; ++g) {
      lnC += dd * std::log(dd) - (dd + Yi[g]) * std::log(Li[g] + dd)
           + R::lgammafn(dd + Yi[g]) - R::lgammafn(dd);
    }
    ll = lnA + lnC;
  };

  arma::vec grad; arma::mat hess; double ll, ll_old = -arma::datum::inf;
  arma::vec theta(k + 1);
  theta.subvec(0, k - 1) = beta; theta[k] = d;

  int iter = 0;
  for (iter = 0; iter < maxit; ++iter) {
    eval(theta.subvec(0, k - 1), theta[k], grad, hess, ll);
    arma::vec step = arma::solve(hess, grad, arma::solve_opts::likely_sympd);

    double lambda = 1.0;
    arma::vec theta_new = theta - lambda * step;
    if (theta_new[k] <= 1e-6) theta_new[k] = 1e-6;
    arma::vec g2; arma::mat h2; double ll_new;
    eval(theta_new.subvec(0, k - 1), theta_new[k], g2, h2, ll_new);
    int halvings = 0;
    while (!std::isfinite(ll_new) || (ll_new < ll && halvings < 30)) {
      lambda *= 0.5;
      theta_new = theta - lambda * step;
      if (theta_new[k] <= 1e-6) theta_new[k] = 1e-6;
      eval(theta_new.subvec(0, k - 1), theta_new[k], g2, h2, ll_new);
      halvings++;
    }

    bool converged = std::fabs(ll_new - ll) < tol * (std::fabs(ll) + 1.0);
    theta = theta_new;
    ll_old = ll_new;
    if (converged) { iter++; break; }
  }

  eval(theta.subvec(0, k - 1), theta[k], grad, hess, ll_old); // Hessian at converged theta
  arma::mat vcov = arma::inv_sympd(-hess); // hess is negative-definite (concave loglik)

  return List::create(
    Named("coefficients") = theta.subvec(0, k - 1),
    Named("dispersion_param") = theta[k],
    Named("vcov_unscaled") = vcov,
    Named("loglik") = ll_old,
    Named("iterations") = iter
  );
}
