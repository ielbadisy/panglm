#include "families.h"
#include "group_agg.h"
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using panglm::compute_group_agg;

// Concentrated (conditional) Poisson fixed-effects likelihood, matching
// pglm's lnl.poisson "within" branch. Newton-Raphson with step-halving;
// the per-group sufficient statistics (Li, lXi) are computed in parallel
// each iteration.
//
// [[Rcpp::export]]
List within_poisson_fit_cpp(const arma::mat& X, const arma::vec& y,
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

  auto eval = [&](const arma::vec& b, arma::vec& grad, arma::mat& hess, double& ll) {
    arma::vec eta = X * b;
    eta = arma::clamp(eta, -30.0, 30.0);
    arma::vec lit = arma::exp(eta);

    arma::vec Li; arma::mat lXi;
    compute_group_agg(X, lit, group_start, group_size, Li, lXi);

    arma::vec w = Yi / Li;               // per-group ratio
    arma::vec w_obs(n), weight1(n);
    for (int i = 0; i < n; ++i) {
      w_obs[i] = w[(int)group_id[i]];
      weight1[i] = w_obs[i] * lit[i];
    }
    arma::vec gradi = y - w_obs % lit;
    grad = X.t() * gradi;

    arma::mat Xw = X.each_col() % arma::sqrt(weight1);
    arma::mat H1 = -(Xw.t() * Xw);

    arma::vec scale2 = arma::sqrt(Yi / arma::square(Li));
    arma::mat lXi_scaled = lXi.each_col() % scale2;
    arma::mat H2 = lXi_scaled.t() * lXi_scaled;

    hess = H1 + H2;

    double lnA = arma::sum(y % eta) - arma::sum(arma::lgamma(y + 1.0));
    double lnB = arma::sum(arma::lgamma(Yi + 1.0) - arma::log(Li) % Yi);
    ll = lnA + lnB;
  };

  arma::vec grad; arma::mat hess; double ll, ll_old = -arma::datum::inf;
  int iter = 0;
  for (iter = 0; iter < maxit; ++iter) {
    eval(beta, grad, hess, ll);
    arma::vec step = arma::solve(hess, grad, arma::solve_opts::likely_sympd);

    double lambda = 1.0;
    arma::vec beta_new = beta - lambda * step;
    arma::vec g2; arma::mat h2; double ll_new;
    eval(beta_new, g2, h2, ll_new);
    int halvings = 0;
    while (!std::isfinite(ll_new) || (ll_new < ll && halvings < 30)) {
      lambda *= 0.5;
      beta_new = beta - lambda * step;
      eval(beta_new, g2, h2, ll_new);
      halvings++;
    }

    bool converged = std::fabs(ll_new - ll) < tol * (std::fabs(ll) + 1.0);
    beta = beta_new;
    ll_old = ll_new;
    if (converged) { iter++; break; }
  }

  eval(beta, grad, hess, ll_old); // Hessian at the converged beta, for the vcov below
  arma::mat vcov = arma::inv_sympd(-hess); // hess is negative-definite (concave loglik)

  return List::create(
    Named("coefficients") = beta,
    Named("vcov_unscaled") = vcov,
    Named("loglik") = ll_old,
    Named("iterations") = iter
  );
}
