#include "families.h"
using namespace Rcpp;
using namespace panglm;

// Generic IRLS for pooled GLM (gaussian/poisson/binomial), any of the
// canonical links defined in families.h. Returns coefficients, the
// unscaled covariance (XtWX)^-1, fitted values, log-likelihood and the
// number of iterations used.
//
// [[Rcpp::export]]
List irls_fit_cpp(const arma::mat& X, const arma::vec& y,
                   int family_id, int link_id,
                   int maxit = 100, double tol = 1e-10) {
  FamilyType family = family_from_int(family_id);
  LinkType   link   = link_from_int(link_id);

  int n = X.n_rows, k = X.n_cols;
  arma::vec beta = arma::zeros<arma::vec>(k);

  // sane starting values
  {
    arma::vec mu0;
    if (family == BINOMIAL) mu0 = (y + 0.5) / 2.0;
    else if (family == POISSON) mu0 = y + 0.1;
    else mu0 = y;
    arma::vec eta0;
    if (link == LOG) eta0 = arma::log(mu0);
    else if (link == LOGIT) eta0 = arma::log(mu0 / (1.0 - mu0));
    else if (link == PROBIT) {
      eta0 = arma::zeros<arma::vec>(n);
      for (int i = 0; i < n; ++i) eta0[i] = R::qnorm(mu0[i], 0.0, 1.0, 1, 0);
    } else eta0 = mu0;
    beta = arma::solve(X.t() * X, X.t() * eta0, arma::solve_opts::likely_sympd);
  }

  arma::mat XtWX_inv;
  double dev_old = 0.0;
  int iter = 0;
  for (iter = 0; iter < maxit; ++iter) {
    arma::vec eta = X * beta;
    arma::vec mu  = linkinv(eta, link);
    // clip to keep variance/derivatives well-defined
    if (family == BINOMIAL) mu = arma::clamp(mu, 1e-8, 1.0 - 1e-8);
    if (family == POISSON)  mu = arma::clamp(mu, 1e-8, arma::datum::inf);

    arma::vec deta = mu_eta(eta, link);
    arma::vec varmu = variance(mu, family);
    varmu = arma::clamp(varmu, 1e-10, arma::datum::inf);

    arma::vec w = (deta % deta) / varmu;
    arma::vec z = eta + (y - mu) / deta;

    arma::mat Xw = X.each_col() % arma::sqrt(w);
    arma::vec zw = z % arma::sqrt(w);

    arma::mat XtWX = Xw.t() * Xw;
    arma::vec XtWz = Xw.t() * zw;
    arma::mat XtWX_i = arma::inv_sympd(XtWX + 1e-12 * arma::eye(k, k));
    arma::vec beta_new = XtWX_i * XtWz;

    double dev_new = 0.0;
    arma::vec mu_new = linkinv(X * beta_new, link);
    if (family == BINOMIAL) mu_new = arma::clamp(mu_new, 1e-8, 1.0 - 1e-8);
    if (family == POISSON)  mu_new = arma::clamp(mu_new, 1e-8, arma::datum::inf);
    for (int i = 0; i < n; ++i) dev_new += -2.0 * loglik_obs(y[i], mu_new[i], family);

    XtWX_inv = XtWX_i;
    beta = beta_new;

    if (iter > 0 && std::fabs(dev_new - dev_old) / (std::fabs(dev_old) + 0.1) < tol) {
      dev_old = dev_new;
      iter++;
      break;
    }
    dev_old = dev_new;
  }

  arma::vec eta_final = X * beta;
  arma::vec mu_final = linkinv(eta_final, link);
  if (family == BINOMIAL) mu_final = arma::clamp(mu_final, 1e-8, 1.0 - 1e-8);
  if (family == POISSON)  mu_final = arma::clamp(mu_final, 1e-8, arma::datum::inf);

  double loglik = 0.0;
  for (int i = 0; i < n; ++i) loglik += loglik_obs(y[i], mu_final[i], family);

  double phi = 1.0;
  if (family == GAUSSIAN) {
    arma::vec resid = y - mu_final;
    phi = arma::dot(resid, resid) / std::max(1, n - k);
    loglik = 0.0;
    for (int i = 0; i < n; ++i)
      loglik += -0.5 * std::log(2.0 * M_PI * phi) - 0.5 * (y[i] - mu_final[i]) * (y[i] - mu_final[i]) / phi;
  }

  return List::create(
    Named("coefficients") = beta,
    Named("vcov_unscaled") = XtWX_inv,
    Named("fitted.values") = mu_final,
    Named("loglik") = loglik,
    Named("dispersion") = phi,
    Named("iterations") = iter
  );
}
