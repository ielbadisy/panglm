#include "families.h"
#include <RcppParallel.h>
using namespace Rcpp;
using namespace RcppParallel;
using namespace panglm;

namespace {

// Gauss-Hermite quadrature log-likelihood + gradient for a binomial (logit
// or probit) random-intercept panel model, integrating out the individual
// effect u_i ~ N(0, sigma^2) numerically. Reduced in parallel over
// individuals (panel units) - each unit's quadrature sum is independent.
struct GHReduceWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  const arma::vec& beta;
  double sigma;
  const RVector<double> z; // Gauss-Hermite nodes
  const RVector<double> w; // weights, already normalized (w_r/sqrt(pi))
  LinkType link;

  double loglik;
  arma::vec grad; // size k+1: beta grad, sigma grad

  GHReduceWorker(const NumericMatrix& X, const NumericVector& y,
                 const IntegerVector& group_start, const IntegerVector& group_size,
                 const arma::vec& beta, double sigma,
                 const NumericVector& z, const NumericVector& w, LinkType link)
    : X(X), y(y), group_start(group_start), group_size(group_size),
      beta(beta), sigma(sigma), z(z), w(w), link(link),
      loglik(0.0), grad(arma::zeros<arma::vec>(beta.n_elem + 1)) {}

  GHReduceWorker(const GHReduceWorker& other, RcppParallel::Split)
    : X(other.X), y(other.y), group_start(other.group_start), group_size(other.group_size),
      beta(other.beta), sigma(other.sigma), z(other.z), w(other.w), link(other.link),
      loglik(0.0), grad(arma::zeros<arma::vec>(other.beta.n_elem + 1)) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t k = X.ncol();
    std::size_t R = z.length();
    std::vector<double> node_ll(R), node_score_sum(R);
    std::vector<std::vector<double>> node_grad(R, std::vector<double>(k));
    const double sqrt2 = std::sqrt(2.0);

    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g];
      int size  = group_size[g];

      for (std::size_t r = 0; r < R; ++r) {
        double b = sqrt2 * sigma * z[r];
        double ll = 0.0, sscore = 0.0;
        std::vector<double> gk(k, 0.0);
        for (int t = 0; t < size; ++t) {
          double eta = b;
          for (std::size_t j = 0; j < k; ++j) eta += X(start + t, j) * beta[j];
          double yt = y[start + t];

          double mu, score, logp;
          if (link == LOGIT) {
            mu = 1.0 / (1.0 + std::exp(-eta));
            score = yt - mu;
            logp = yt * std::log(mu) + (1.0 - yt) * std::log(1.0 - mu);
          } else { // PROBIT
            double Phi = R::pnorm(eta, 0.0, 1.0, 1, 0);
            double phi = R::dnorm(eta, 0.0, 1.0, 0);
            Phi = std::min(std::max(Phi, 1e-12), 1.0 - 1e-12);
            score = yt * (phi / Phi) - (1.0 - yt) * (phi / (1.0 - Phi));
            logp = yt * std::log(Phi) + (1.0 - yt) * std::log(1.0 - Phi);
          }
          ll += logp;
          sscore += score;
          for (std::size_t j = 0; j < k; ++j) gk[j] += score * X(start + t, j);
        }
        node_ll[r] = ll;
        node_score_sum[r] = sscore;
        node_grad[r] = gk;
      }

      double m = *std::max_element(node_ll.begin(), node_ll.end());
      double denom = 0.0;
      std::vector<double> wexp(R);
      for (std::size_t r = 0; r < R; ++r) {
        wexp[r] = w[r] * std::exp(node_ll[r] - m);
        denom += wexp[r];
      }
      double loglik_i = m + std::log(denom);
      loglik += loglik_i;

      for (std::size_t r = 0; r < R; ++r) {
        double post = wexp[r] / denom;
        for (std::size_t j = 0; j < k; ++j) grad[j] += post * node_grad[r][j];
        grad[k] += post * node_score_sum[r] * sqrt2 * z[r];
      }
    }
  }

  void join(const GHReduceWorker& rhs) {
    loglik += rhs.loglik;
    grad += rhs.grad;
  }
};

} // namespace

// [[Rcpp::export]]
List random_binomial_loglik_grad_cpp(const arma::vec& beta, double sigma,
                                      const arma::mat& X, const arma::vec& y,
                                      IntegerVector group_start, IntegerVector group_size,
                                      NumericVector nodes, NumericVector weights,
                                      int link_id) {
  int n = X.n_rows, k = X.n_cols;
  NumericMatrix Xr(n, k);
  for (int i = 0; i < n; ++i) for (int j = 0; j < k; ++j) Xr(i, j) = X(i, j);
  NumericVector yr(y.begin(), y.end());

  LinkType link = link_from_int(link_id);
  GHReduceWorker worker(Xr, yr, group_start, group_size, beta, sigma, nodes, weights, link);
  parallelReduce(0, group_start.size(), worker);

  return List::create(
    Named("loglik") = worker.loglik,
    Named("gradient") = worker.grad
  );
}
