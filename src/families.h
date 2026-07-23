#ifndef PANGLM_FAMILIES_H
#define PANGLM_FAMILIES_H

#include <RcppArmadillo.h>

namespace panglm {

enum FamilyType { GAUSSIAN = 0, POISSON = 1, BINOMIAL = 2 };
enum LinkType    { IDENTITY = 0, LOG = 1, LOGIT = 2, PROBIT = 3 };

inline arma::vec linkinv(const arma::vec& eta, LinkType link) {
  switch (link) {
    case IDENTITY:
      return eta;
    case LOG:
      return arma::exp(eta);
    case LOGIT:
      return 1.0 / (1.0 + arma::exp(-eta));
    case PROBIT: {
      arma::vec out(eta.n_elem);
      for (arma::uword i = 0; i < eta.n_elem; ++i)
        out[i] = R::pnorm(eta[i], 0.0, 1.0, 1, 0);
      return out;
    }
  }
  return eta;
}

inline arma::vec mu_eta(const arma::vec& eta, LinkType link) {
  switch (link) {
    case IDENTITY:
      return arma::ones<arma::vec>(eta.n_elem);
    case LOG:
      return arma::exp(eta);
    case LOGIT: {
      arma::vec p = linkinv(eta, LOGIT);
      return p % (1.0 - p);
    }
    case PROBIT: {
      arma::vec out(eta.n_elem);
      for (arma::uword i = 0; i < eta.n_elem; ++i)
        out[i] = R::dnorm(eta[i], 0.0, 1.0, 0);
      return out;
    }
  }
  return arma::ones<arma::vec>(eta.n_elem);
}

inline arma::vec variance(const arma::vec& mu, FamilyType family) {
  switch (family) {
    case GAUSSIAN:
      return arma::ones<arma::vec>(mu.n_elem);
    case POISSON:
      return mu;
    case BINOMIAL:
      return mu % (1.0 - mu);
  }
  return arma::ones<arma::vec>(mu.n_elem);
}

// log-likelihood contribution per observation (up to family-specific
// normalizing terms that do not depend on the parameters, dropped where
// they would only add a constant, e.g. binomial's log-choose(1,y) = 0).
inline double loglik_obs(double y, double mu, FamilyType family) {
  switch (family) {
    case GAUSSIAN:
      return -0.5 * (y - mu) * (y - mu);
    case POISSON:
      return y * std::log(mu) - mu - R::lgammafn(y + 1.0);
    case BINOMIAL:
      return y * std::log(mu) + (1.0 - y) * std::log(1.0 - mu);
  }
  return 0.0;
}

inline FamilyType family_from_int(int f) { return static_cast<FamilyType>(f); }
inline LinkType   link_from_int(int l)   { return static_cast<LinkType>(l); }

} // namespace panglm

#endif
