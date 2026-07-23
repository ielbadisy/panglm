#include "families.h"
#include "group_agg.h"
using namespace Rcpp;
using panglm::compute_group_agg;

namespace {

// Levenberg-Marquardt-style damped solve: try solving (H + lambda*I) x = g
// with escalating ridge lambda until it succeeds, rather than throwing on a
// singular/near-singular Hessian (common on sparse/skewed panels where the
// dispersion parameter is weakly identified). Returns the damping actually
// used (0 if none was needed) so callers can warn when it engages.
double damped_solve(const arma::mat& H, const arma::vec& g, arma::vec& x) {
  double diag_scale = std::max(1e-8, arma::mean(arma::abs(H.diag())));
  double lambda = 0.0;
  arma::mat I = arma::eye(H.n_rows, H.n_cols);
  for (int attempt = 0; attempt < 20; ++attempt) {
    bool ok = arma::solve(x, H - lambda * I, g, arma::solve_opts::no_approx + arma::solve_opts::likely_sympd);
    if (ok && x.is_finite()) return lambda;
    lambda = (lambda == 0.0) ? (1e-6 * diag_scale) : (lambda * 10.0);
  }
  x = arma::solve(H - lambda * I, g, arma::solve_opts::likely_sympd); // last resort, may throw
  return lambda;
}

// Same idea for the final covariance inversion: damp (-Hessian) toward
// positive-definiteness rather than erroring, and fall back to a
// pseudo-inverse if damping alone doesn't suffice.
double damped_inv_sympd(const arma::mat& negH, arma::mat& out) {
  double diag_scale = std::max(1e-8, arma::mean(arma::abs(negH.diag())));
  double lambda = 0.0;
  arma::mat I = arma::eye(negH.n_rows, negH.n_cols);
  for (int attempt = 0; attempt < 20; ++attempt) {
    bool ok = arma::inv_sympd(out, negH + lambda * I);
    if (ok) return lambda;
    lambda = (lambda == 0.0) ? (1e-6 * diag_scale) : (lambda * 10.0);
  }
  out = arma::pinv(negH); // last resort
  return lambda;
}

} // namespace

// Poisson random-effects (Poisson-Gamma / negative-binomial marginal)
// panel model, ported from pglm's lnl.poisson "random" branch. Jointly
// estimates the K regression coefficients beta and the Gamma dispersion
// parameter d > 0 (frailty variance is 1/d) by Newton-Raphson with step
// halving; per-group sufficient statistics computed in parallel.
//
// beta_start/d_start let the caller seed the optimizer with the pooled
// Poisson MLE and a method-of-moments dispersion estimate rather than
// zeros/1, which matters on sparse or skewed panels where the plain
// zero-start Newton path can wander into a near-singular Hessian region.
// Hessian solves are Levenberg-Marquardt damped instead of throwing on
// (near-)singularity; `damping_used` reports the largest ridge applied,
// so R can warn the user rather than silently degrading precision.
//
// [[Rcpp::export]]
List random_poisson_fit_cpp(const arma::mat& X, const arma::vec& y,
                             IntegerVector group_start, IntegerVector group_size,
                             int maxit = 100, double tol = 1e-10,
                             Rcpp::Nullable<Rcpp::NumericVector> beta_start = R_NilValue,
                             double d_start = 1.0) {
  int n = X.n_rows, k = X.n_cols, G = group_start.size();

  arma::vec group_id(n);
  for (int g = 0; g < G; ++g)
    for (int t = 0; t < group_size[g]; ++t)
      group_id[group_start[g] + t] = g;

  arma::vec Yi(G, arma::fill::zeros);
  for (int i = 0; i < n; ++i) Yi[(int)group_id[i]] += y[i];

  arma::vec beta = arma::zeros<arma::vec>(k);
  if (beta_start.isNotNull()) {
    NumericVector bs(beta_start);
    if ((int)bs.size() == k) for (int j = 0; j < k; ++j) beta[j] = bs[j];
  }
  double d = (d_start > 0 && std::isfinite(d_start)) ? d_start : 1.0;

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
  double max_damping = 0.0;

  int iter = 0;
  for (iter = 0; iter < maxit; ++iter) {
    eval(theta.subvec(0, k - 1), theta[k], grad, hess, ll);
    arma::vec step;
    double damping = damped_solve(-hess, -grad, step); // Newton step on -hess (neg-definite -> pos-definite)
    max_damping = std::max(max_damping, damping);

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
  arma::mat vcov;
  double vcov_damping = damped_inv_sympd(-hess, vcov); // hess is negative-definite (concave loglik)
  max_damping = std::max(max_damping, vcov_damping);

  return List::create(
    Named("coefficients") = theta.subvec(0, k - 1),
    Named("dispersion_param") = theta[k],
    Named("vcov_unscaled") = vcov,
    Named("loglik") = ll_old,
    Named("iterations") = iter,
    Named("damping_used") = max_damping
  );
}
