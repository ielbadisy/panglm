// Weighted logistic regression by iteratively reweighted least squares,
// following stats::glm.fit(family = binomial()) step by step: same starting
// values, same clamping of the linear predictor, same deviance and
// convergence rule. It only skips the R-level overhead of glm.fit. Each
// least squares step is solved by a QR decomposition, as in glm.fit.
//
// ok = false asks the caller to fall back on glm.fit (near-singular design,
// non-finite values, no convergence or fitted probabilities numerically 0
// or 1), so unusual cases keep the exact glm.fit behaviour.

#include <RcppArmadillo.h>
#include <cfloat>
#include <cmath>

// [[Rcpp::depends(RcppArmadillo)]]

namespace {

// binomial()$linkinv and binomial()$mu.eta for the logit link, as in R's
// C code (src/library/stats/src/family.c: THRESH = 30)
const double THRESH = 30.0;

inline double logit_linkinv(double eta) {
  if (eta < -THRESH) return DBL_EPSILON;
  if (eta > THRESH) return 1.0 / DBL_EPSILON;
  return std::exp(eta);
}

inline double logit_mu(double eta) {
  double t = logit_linkinv(eta);
  return t / (1.0 + t);
}

inline double logit_mu_eta(double eta) {
  if (eta > THRESH || eta < -THRESH) return DBL_EPSILON;
  double opexp = 1.0 + std::exp(eta);
  return std::exp(eta) / (opexp * opexp);
}

// binomial()$dev.resids for y in [0, 1]
inline double binom_dev(double y, double mu, double w) {
  double r = 0.0;
  if (y > 0.0) r += y * std::log(y / mu);
  if (y < 1.0) r += (1.0 - y) * std::log((1.0 - y) / (1.0 - mu));
  return 2.0 * w * r;
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List logit_irls_cpp(const arma::mat& X, const arma::vec& y, const arma::vec& w,
                          double epsilon, int maxit) {
  const arma::uword n = X.n_rows, p = X.n_cols;
  arma::uvec good = arma::find(w > 0);
  const arma::uword ng = good.n_elem;
  arma::vec mu(n), eta(n), beta(p, arma::fill::zeros);
  // binomial()$initialize with one trial per row
  for (arma::uword i = 0; i < n; ++i) {
    double m = (w[i] * y[i] + 0.5) / (w[i] + 1.0);
    eta[i] = std::log(m / (1.0 - m));
    mu[i] = logit_mu(eta[i]);
  }
  // deviances are summed in long double, as R's sum() does
  long double acc = 0.0L;
  for (arma::uword i = 0; i < n; ++i) acc += binom_dev(y[i], mu[i], w[i]);
  double devold = (double) acc;
  bool conv = false, ok = (ng >= p) && std::isfinite(devold);
  int iter = 0;
  arma::mat Xw(ng, p);
  arma::vec zw(ng);
  for (iter = 1; ok && iter <= maxit; ++iter) {
    for (arma::uword k = 0; k < ng; ++k) {
      arma::uword i = good[k];
      double me = logit_mu_eta(eta[i]);
      double var = mu[i] * (1.0 - mu[i]);
      double z = eta[i] + (y[i] - mu[i]) / me;
      double ww = std::sqrt(w[i] * me * me / var);
      zw[k] = z * ww;
      for (arma::uword j = 0; j < p; ++j) Xw(k, j) = X(i, j) * ww;
    }
    arma::mat Q, R;
    if (!arma::qr_econ(Q, R, Xw)) { ok = false; break; }
    // a diagonal of R that is small next to its column norm means a nearly
    // collinear design; glm.fit handles that case
    for (arma::uword j = 0; j < p; ++j) {
      double cn = arma::norm(Xw.col(j));
      if (!(std::fabs(R(j, j)) > 1e-7 * cn)) { ok = false; break; }
    }
    if (!ok) break;
    beta = arma::solve(arma::trimatu(R), Q.t() * zw, arma::solve_opts::no_approx);
    if (!beta.is_finite()) { ok = false; break; }
    eta = X * beta;
    acc = 0.0L;
    for (arma::uword i = 0; i < n; ++i) {
      mu[i] = logit_mu(eta[i]);
      acc += binom_dev(y[i], mu[i], w[i]);
    }
    double dev = (double) acc;
    if (!std::isfinite(dev)) { ok = false; break; }
    if (std::fabs(dev - devold) / (std::fabs(dev) + 0.1) < epsilon) {
      conv = true;
      break;
    }
    devold = dev;
  }
  if (!conv) ok = false;
  // fitted probabilities numerically 0 or 1 (the case glm.fit warns about)
  // come from separated data, where the coefficients are not determined;
  // glm.fit handles them
  if (ok) {
    const double eps = 10.0 * DBL_EPSILON;
    for (arma::uword i = 0; i < n; ++i) {
      if (mu[i] > 1.0 - eps || mu[i] < eps) { ok = false; break; }
    }
  }
  return Rcpp::List::create(Rcpp::Named("coefficients") = Rcpp::NumericVector(beta.begin(), beta.end()),
                            Rcpp::Named("fitted.values") = Rcpp::NumericVector(mu.begin(), mu.end()),
                            Rcpp::Named("iter") = iter,
                            Rcpp::Named("converged") = conv,
                            Rcpp::Named("ok") = ok);
}
