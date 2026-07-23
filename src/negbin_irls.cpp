#include <RcppArmadillo.h>
using namespace Rcpp;

namespace {

// One Newton step on the theta profile log-likelihood (mu held fixed),
// the same update used by MASS::glm.nb's theta.ml. theta is parametrized
// on the log scale internally to keep it positive without ad hoc clipping.
double update_theta(const arma::vec& y, const arma::vec& mu, double theta, int newton_steps = 4) {
  double kappa = std::log(theta);
  for (int s = 0; s < newton_steps; ++s) {
    double th = std::exp(kappa);
    double g = 0.0, h = 0.0;
    for (arma::uword i = 0; i < y.n_elem; ++i) {
      double ymu = th + mu[i];
      g += R::digamma(y[i] + th) - R::digamma(th) + std::log(th) + 1.0
         - std::log(ymu) - (y[i] + th) / ymu;
      h += R::trigamma(y[i] + th) - R::trigamma(th) + 1.0 / th
         - 2.0 / ymu + (y[i] + th) / (ymu * ymu);
    }
    // chain rule to kappa = log(theta): d/dkappa = g*th, d2/dkappa2 = h*th^2 + g*th
    double g_k = g * th;
    double h_k = h * th * th + g * th;
    if (!std::isfinite(g_k) || !std::isfinite(h_k) || h_k == 0.0) break;
    double step = g_k / h_k;
    if (!std::isfinite(step)) break;
    kappa -= step;
    kappa = std::max(-15.0, std::min(15.0, kappa));
  }
  return std::exp(kappa);
}

double negbin_loglik(const arma::vec& y, const arma::vec& mu, double theta) {
  double ll = 0.0;
  for (arma::uword i = 0; i < y.n_elem; ++i) {
    ll += R::lgammafn(y[i] + theta) - R::lgammafn(theta) - R::lgammafn(y[i] + 1.0)
        + theta * std::log(theta) - theta * std::log(theta + mu[i])
        + y[i] * std::log(mu[i]) - y[i] * std::log(theta + mu[i]);
  }
  return ll;
}

} // namespace

// Genuine NB2 regression (log link) via alternating IRLS-for-beta /
// Newton-for-theta (the standard MASS::glm.nb algorithm), implemented
// directly rather than through the generic canonical-link IRLS since the
// NB2 variance V(mu) = mu + mu^2/theta depends on the jointly-estimated
// dispersion theta.
//
// Used both for pooled negbin (X = covariates) and for the
// Allison-Waterman unconditional fixed-effects negative binomial
// (X = covariates augmented with individual dummy columns, feasible for
// small-to-moderate N since it just adds columns to a dense IRLS problem).
//
// [[Rcpp::export]]
List negbin_irls_fit_cpp(const arma::mat& X, const arma::vec& y,
                         int maxit = 100, double tol = 1e-10) {
  int n = X.n_rows, k = X.n_cols;

  // Starting values: OLS on log(y + 0.1) for beta; method-of-moments for theta.
  arma::vec ystart = arma::log(y + 0.1);
  arma::vec beta = arma::solve(X.t() * X, X.t() * ystart, arma::solve_opts::likely_sympd);
  double ybar = arma::mean(y), yvar = arma::var(y);
  double theta = (yvar > ybar && ybar > 0) ? (ybar * ybar) / (yvar - ybar) : 100.0;
  if (!std::isfinite(theta) || theta <= 0) theta = 100.0;

  arma::mat XtWX_inv(k, k, arma::fill::eye);
  double ll_old = -arma::datum::inf;
  int iter = 0;

  for (iter = 0; iter < maxit; ++iter) {
    // --- IRLS for beta given theta ---
    arma::vec beta_inner = beta;
    for (int inner = 0; inner < 50; ++inner) {
      arma::vec eta = arma::clamp(X * beta_inner, -30.0, 30.0);
      arma::vec mu = arma::clamp(arma::exp(eta), 1e-8, arma::datum::inf);
      arma::vec w = (mu * theta) / (theta + mu);
      arma::vec z = eta + (y - mu) / mu;

      arma::mat Xw = X.each_col() % arma::sqrt(w);
      arma::vec zw = z % arma::sqrt(w);
      arma::mat XtWX = Xw.t() * Xw;
      XtWX_inv = arma::inv_sympd(XtWX + 1e-10 * arma::eye(k, k));
      arma::vec beta_new = XtWX_inv * (Xw.t() * zw);

      double delta = arma::norm(beta_new - beta_inner, 2) / (arma::norm(beta_inner, 2) + 1e-8);
      beta_inner = beta_new;
      if (delta < tol) break;
    }
    beta = beta_inner;

    // --- Newton for theta given beta ---
    arma::vec eta = arma::clamp(X * beta, -30.0, 30.0);
    arma::vec mu = arma::clamp(arma::exp(eta), 1e-8, arma::datum::inf);
    theta = update_theta(y, mu, theta);

    double ll_new = negbin_loglik(y, mu, theta);
    bool converged = std::fabs(ll_new - ll_old) < tol * (std::fabs(ll_old) + 1.0);
    ll_old = ll_new;
    if (converged) { iter++; break; }
  }

  return List::create(
    Named("coefficients") = beta,
    Named("theta") = theta,
    Named("vcov_unscaled") = XtWX_inv,
    Named("loglik") = ll_old,
    Named("iterations") = iter
  );
}
