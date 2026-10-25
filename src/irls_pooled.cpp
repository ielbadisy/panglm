#include "families.h"
using namespace Rcpp;
using namespace panglm;

namespace {

arma::vec clamp_mean(arma::vec mu, FamilyType family) {
  if (family == BINOMIAL) mu = arma::clamp(mu, 1e-8, 1.0 - 1e-8);
  if (family == POISSON || family == GAMMA) mu = arma::clamp(mu, 1e-10, arma::datum::inf);
  return mu;
}

} // namespace

// Generic IRLS for pooled GLM (gaussian/poisson/binomial/gamma), any of
// the links defined in families.h, with optional prior (observation)
// weights. Returns coefficients, the unscaled covariance (X'WX)^-1, fitted
// values, the (weighted) log-likelihood, the dispersion and the number of
// iterations used. Observations with zero prior weight do not contribute
// to the fit, the log-likelihood or the residual degrees of freedom.
//
// [[Rcpp::export]]
List irls_fit_cpp(const arma::mat& X, const arma::vec& y,
                   int family_id, int link_id,
                   int maxit = 100, double tol = 1e-10,
                   Rcpp::Nullable<Rcpp::NumericVector> weights = R_NilValue) {
  FamilyType family = family_from_int(family_id);
  LinkType   link   = link_from_int(link_id);

  int n = X.n_rows, k = X.n_cols;
  arma::vec prior = arma::ones<arma::vec>(n);
  if (weights.isNotNull()) {
    prior = Rcpp::as<arma::vec>(weights.get());
    if ((int) prior.n_elem != n) stop("weights must have one value per row of X");
  }
  double n_pos = (double) arma::accu(prior > 0.0);
  arma::vec beta = arma::zeros<arma::vec>(k);

  // sane starting values
  {
    arma::vec mu0;
    if (family == BINOMIAL) mu0 = (y + 0.5) / 2.0;
    else if (family == POISSON) mu0 = y + 0.1;
    else if (family == GAMMA) mu0 = arma::clamp(y, 1e-10, arma::datum::inf);
    else mu0 = y;
    arma::vec eta0;
    if (link == LOG) eta0 = arma::log(mu0);
    else if (link == LOGIT) eta0 = arma::log(mu0 / (1.0 - mu0));
    else if (link == PROBIT) {
      eta0 = arma::zeros<arma::vec>(n);
      for (int i = 0; i < n; ++i) eta0[i] = R::qnorm(mu0[i], 0.0, 1.0, 1, 0);
    } else eta0 = mu0;
    arma::mat Xp = X.each_col() % prior;
    beta = arma::solve(Xp.t() * X + 1e-12 * arma::eye(k, k), Xp.t() * eta0,
                       arma::solve_opts::likely_sympd);
  }

  arma::mat XtWX_inv;
  double dev_old = 0.0;
  int iter = 0;
  for (iter = 0; iter < maxit; ++iter) {
    arma::vec eta = X * beta;
    arma::vec mu  = clamp_mean(linkinv(eta, link), family);

    arma::vec deta = mu_eta(eta, link);
    arma::vec varmu = variance(mu, family);
    varmu = arma::clamp(varmu, 1e-10, arma::datum::inf);

    arma::vec w = prior % (deta % deta) / varmu;
    arma::vec z = eta + (y - mu) / deta;

    arma::mat Xw = X.each_col() % arma::sqrt(w);
    arma::vec zw = z % arma::sqrt(w);

    arma::mat XtWX = Xw.t() * Xw;
    arma::vec XtWz = Xw.t() * zw;
    arma::mat XtWX_i = arma::inv_sympd(XtWX + 1e-12 * arma::eye(k, k));
    arma::vec beta_new = XtWX_i * XtWz;

    double dev_new = 0.0;
    arma::vec mu_new = clamp_mean(linkinv(X * beta_new, link), family);
    for (int i = 0; i < n; ++i) {
      if (prior[i] > 0.0) dev_new += -2.0 * prior[i] * loglik_obs(y[i], mu_new[i], family);
    }

    XtWX_inv = XtWX_i;
    beta = beta_new;

    if (iter > 0 && std::fabs(dev_new - dev_old) / (std::fabs(dev_old) + 0.1) < tol) {
      dev_old = dev_new;
      iter++;
      break;
    }
    dev_old = dev_new;
  }

  // Covariance evaluated at the final coefficients.
  arma::vec eta_final = X * beta;
  arma::vec mu_final = clamp_mean(linkinv(eta_final, link), family);
  {
    arma::vec deta = mu_eta(eta_final, link);
    arma::vec varmu = arma::clamp(variance(mu_final, family), 1e-10, arma::datum::inf);
    arma::vec w = prior % (deta % deta) / varmu;
    arma::mat Xw = X.each_col() % arma::sqrt(w);
    XtWX_inv = arma::inv_sympd(Xw.t() * Xw + 1e-12 * arma::eye(k, k));
  }

  double loglik = 0.0;
  for (int i = 0; i < n; ++i) {
    if (prior[i] > 0.0) loglik += prior[i] * loglik_obs(y[i], mu_final[i], family);
  }

  double df = std::max(1.0, n_pos - k);
  double phi = 1.0;
  if (family == GAUSSIAN) {
    arma::vec resid = y - mu_final;
    double rss = arma::dot(prior % resid, resid);
    phi = rss / df;
    // Same convention as stats::glm(): -2 logLik = n (log(2 pi RSS / n) + 1)
    // - sum(log w), with n the number of observations with positive weight.
    double slogw = 0.0;
    for (int i = 0; i < n; ++i) if (prior[i] > 0.0) slogw += std::log(prior[i]);
    loglik = -0.5 * n_pos * (std::log(2.0 * M_PI * rss / n_pos) + 1.0) + 0.5 * slogw;
  } else if (family == GAMMA) {
    arma::vec r = (y - mu_final) / mu_final;
    phi = arma::dot(prior % r, r) / df;
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
