#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include "threads.h"
using namespace Rcpp;
using namespace RcppParallel;

namespace {

// Conditional logistic regression (Chamberlain 1980 fixed-effects logit):
// for a group of T observations with linear predictors eta_0..eta_{T-1}
// and total successes S = sum(y), the fixed effect is profiled out exactly
// by conditioning on S. The conditional likelihood needs the elementary
// symmetric polynomial (ESP) of degree S in the weights w_t = exp(eta_t),
// and the score needs each observation's "leave-one-out" marginal
// probability p_t = P(y_t = 1 | sum = S), obtained via a forward/backward
// DP (the same recursion used by exact conditional logistic regression,
// e.g. survival::clogit(method = "exact")) in O(T^2) per group.
//
// Groups with S = 0 or S = T carry no information (every 0/1 sequence is
// equally consistent with them) and are skipped, matching clogit/fixest
// convention.
struct GroupResult { double ll; std::vector<double> grad_eta; bool used; };

GroupResult conditional_logit_group(const std::vector<double>& eta, const std::vector<double>& y) {
  int T = eta.size();
  int S = 0;
  for (int t = 0; t < T; ++t) S += (int) std::round(y[t]);

  GroupResult out;
  out.grad_eta.assign(T, 0.0);
  out.ll = 0.0;
  out.used = false;
  if (S == 0 || S == T) return out; // no information, drop
  out.used = true;

  double maxeta = eta[0];
  for (int t = 1; t < T; ++t) if (eta[t] > maxeta) maxeta = eta[t];
  std::vector<double> w(T);
  for (int t = 0; t < T; ++t) w[t] = std::exp(eta[t] - maxeta);

  // Forward: A[t][s] = ESP of degree s using w_0..w_{t-1}
  std::vector<std::vector<double>> A(T + 1, std::vector<double>(T + 1, 0.0));
  A[0][0] = 1.0;
  for (int t = 1; t <= T; ++t) {
    A[t][0] = 1.0;
    for (int s = 1; s <= t; ++s) A[t][s] = A[t - 1][s] + w[t - 1] * A[t - 1][s - 1];
  }

  // Backward: B[t][s] = ESP of degree s using w_t..w_{T-1}
  std::vector<std::vector<double>> B(T + 1, std::vector<double>(T + 1, 0.0));
  B[T][0] = 1.0;
  for (int t = T - 1; t >= 0; --t) {
    B[t][0] = 1.0;
    int maxs = T - t;
    for (int s = 1; s <= maxs; ++s) B[t][s] = B[t + 1][s] + w[t] * B[t + 1][s - 1];
  }

  double CS = A[T][S]; // total ESP of degree S (the conditional normalizing constant)

  double lnA = 0.0; // sum_t y_t * eta_t (centering cancels: y-weighted sum uses original eta)
  for (int t = 0; t < T; ++t) lnA += y[t] * eta[t];
  out.ll = lnA - (std::log(CS) + S * maxeta); // undo the max-eta centering on log C(S)

  for (int t = 0; t < T; ++t) {
    double Dt = 0.0; // leave-one-out ESP of degree S-1, excluding index t
    int lo = std::max(0, S - 1 - (T - t - 1));
    int hi = std::min(S - 1, t);
    for (int j = lo; j <= hi; ++j) Dt += A[t][j] * B[t + 1][S - 1 - j];
    double p_t = (CS > 0) ? (w[t] * Dt / CS) : 0.0;
    out.grad_eta[t] = y[t] - p_t;
  }
  return out;
}

struct ConditionalLogitWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> eta;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  double ll;
  std::vector<double> grad;
  std::vector<double> score_eta;
  int n_used_groups;

  ConditionalLogitWorker(const NumericMatrix& X, const NumericVector& eta, const NumericVector& y,
                         const IntegerVector& group_start, const IntegerVector& group_size)
    : X(X), eta(eta), y(y), group_start(group_start), group_size(group_size),
      ll(0.0), grad(X.ncol(), 0.0), score_eta(X.nrow(), 0.0), n_used_groups(0) {}

  ConditionalLogitWorker(const ConditionalLogitWorker& other, Split)
    : X(other.X), eta(other.eta), y(other.y), group_start(other.group_start), group_size(other.group_size),
      ll(0.0), grad(other.X.ncol(), 0.0), score_eta(other.X.nrow(), 0.0), n_used_groups(0) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::size_t k = X.ncol();
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g], size = group_size[g];
      std::vector<double> eta_g(size), y_g(size);
      for (int t = 0; t < size; ++t) { eta_g[t] = eta[start + t]; y_g[t] = y[start + t]; }

      GroupResult res = conditional_logit_group(eta_g, y_g);
      if (!res.used) continue;
      n_used_groups++;
      ll += res.ll;
      for (int t = 0; t < size; ++t) {
        score_eta[start + t] = res.grad_eta[t];
        for (std::size_t j = 0; j < k; ++j)
          grad[j] += res.grad_eta[t] * X(start + t, j);
      }
    }
  }

  void join(const ConditionalLogitWorker& rhs) {
    ll += rhs.ll;
    n_used_groups += rhs.n_used_groups;
    for (std::size_t j = 0; j < grad.size(); ++j) grad[j] += rhs.grad[j];
    for (std::size_t i = 0; i < score_eta.size(); ++i) score_eta[i] += rhs.score_eta[i];
  }
};

} // namespace

// Conditional (fixed-effects) logistic regression log-likelihood and
// gradient, exact (Chamberlain 1980), parallelized over panel units.
//
// [[Rcpp::export]]
List conditional_logit_loglik_grad_cpp(const arma::vec& beta, const arma::mat& X, const arma::vec& y,
                                       IntegerVector group_start, IntegerVector group_size) {
  int n = X.n_rows, k = X.n_cols;
  arma::vec eta_arma = X * beta;
  NumericMatrix Xr(n, k);
  for (int i = 0; i < n; ++i) for (int j = 0; j < k; ++j) Xr(i, j) = X(i, j);
  NumericVector eta(eta_arma.begin(), eta_arma.end());
  NumericVector yr(y.begin(), y.end());

  ConditionalLogitWorker worker(Xr, eta, yr, group_start, group_size);
  parallelReduce(0, group_start.size(), worker, panglm::grain(group_start.size()));

  return List::create(
    Named("loglik") = worker.ll,
    Named("gradient") = worker.grad,
    Named("score_eta") = worker.score_eta,
    Named("n_used_groups") = worker.n_used_groups
  );
}

namespace {

// Log-likelihood, gradient and Hessian of the conditional logit with respect
// to beta, by the recursion of exact conditional logistic regression (Gail,
// Lubin and Rubinstein 1981; survival::clogit(method = "exact")): the
// elementary symmetric polynomial C(S) of the weights w_t = exp(x_t'beta)
// and its first and second derivatives in beta are built one observation
// at a time, C_t(s) = C_{t-1}(s) + w_t C_{t-1}(s - 1), so the derivatives
// follow from d w_t = w_t x_t. Then ll_i = y'eta - log C(S),
// grad_i = X'y - C'/C and Hess_i = -(C''/C - C' C'^T / C^2).
struct ConditionalLogitHessWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> eta;
  const RVector<double> y;
  const RVector<int> group_start;
  const RVector<int> group_size;
  std::size_t k;
  double ll;
  std::vector<double> grad;
  std::vector<double> hess;
  int n_used_groups;

  ConditionalLogitHessWorker(const NumericMatrix& X, const NumericVector& eta, const NumericVector& y,
                             const IntegerVector& group_start, const IntegerVector& group_size)
    : X(X), eta(eta), y(y), group_start(group_start), group_size(group_size), k(X.ncol()),
      ll(0.0), grad(X.ncol(), 0.0), hess(X.ncol() * X.ncol(), 0.0), n_used_groups(0) {}

  ConditionalLogitHessWorker(const ConditionalLogitHessWorker& o, Split)
    : X(o.X), eta(o.eta), y(o.y), group_start(o.group_start), group_size(o.group_size), k(o.k),
      ll(0.0), grad(o.k, 0.0), hess(o.k * o.k, 0.0), n_used_groups(0) {}

  void operator()(std::size_t begin, std::size_t end) {
    std::vector<double> A, dA, d2A;
    for (std::size_t g = begin; g < end; ++g) {
      int start = group_start[g], T = group_size[g];
      int S = 0;
      for (int t = 0; t < T; ++t) S += (int) std::round(y[start + t]);
      if (S == 0 || S == T) continue;
      n_used_groups++;
      // Use the smaller of S and T - S: conditioning on S successes is the
      // same as conditioning on T - S failures with eta replaced by -eta.
      bool flip = 2 * S > T;
      int Sx = flip ? T - S : S;
      double sgn = flip ? -1.0 : 1.0;
      double maxeta = -INFINITY;
      for (int t = 0; t < T; ++t) maxeta = std::max(maxeta, sgn * eta[start + t]);
      A.assign(Sx + 1, 0.0);
      dA.assign((Sx + 1) * k, 0.0);
      d2A.assign((Sx + 1) * k * k, 0.0);
      A[0] = 1.0;
      for (int t = 0; t < T; ++t) {
        double w = std::exp(sgn * eta[start + t] - maxeta);
        int top = std::min(t + 1, Sx);
        for (int s = top; s >= 1; --s) {
          double a = A[s - 1];
          const double* da = &dA[(s - 1) * k];
          const double* d2a = &d2A[(s - 1) * k * k];
          double* dn = &dA[s * k];
          double* d2n = &d2A[s * k * k];
          for (std::size_t i = 0; i < k; ++i) {
            double xi = sgn * X(start + t, i);
            for (std::size_t j = 0; j <= i; ++j) {
              double xj = sgn * X(start + t, j);
              double v = d2a[i * k + j] + xi * da[j] + da[i] * xj + xi * xj * a;
              d2n[i * k + j] += w * v;
            }
          }
          for (std::size_t i = 0; i < k; ++i) dn[i] += w * (da[i] + sgn * X(start + t, i) * a);
          A[s] += w * a;
        }
      }
      double C = A[Sx];
      const double* dC = &dA[Sx * k];
      const double* d2C = &d2A[Sx * k * k];
      double lnum = 0.0;
      for (int t = 0; t < T; ++t) {
        double yt = flip ? 1.0 - y[start + t] : y[start + t];
        lnum += yt * sgn * eta[start + t];
      }
      ll += lnum - (std::log(C) + Sx * maxeta);
      for (std::size_t i = 0; i < k; ++i) {
        double xy = 0.0;
        for (int t = 0; t < T; ++t) {
          double yt = flip ? 1.0 - y[start + t] : y[start + t];
          xy += yt * sgn * X(start + t, i);
        }
        // the flip changes the sign of beta's effect: chain rule factor sgn
        // is already inside the x's, so the gradient in beta is direct
        grad[i] += xy - dC[i] / C;
        for (std::size_t j = 0; j <= i; ++j) {
          double h = -(d2C[i * k + j] / C - dC[i] * dC[j] / (C * C));
          hess[i * k + j] += h;
        }
      }
    }
  }

  void join(const ConditionalLogitHessWorker& rhs) {
    ll += rhs.ll;
    n_used_groups += rhs.n_used_groups;
    for (std::size_t j = 0; j < grad.size(); ++j) grad[j] += rhs.grad[j];
    for (std::size_t j = 0; j < hess.size(); ++j) hess[j] += rhs.hess[j];
  }
};

} // namespace

// Conditional logit log-likelihood, gradient and Hessian in beta (exact
// recursion), parallelized over panel units.
//
// [[Rcpp::export]]
List conditional_logit_hessian_cpp(const arma::vec& beta, const arma::mat& X, const arma::vec& y,
                                   IntegerVector group_start, IntegerVector group_size) {
  int n = X.n_rows, k = X.n_cols;
  arma::vec eta_arma = X * beta;
  NumericMatrix Xr(n, k);
  std::copy(X.begin(), X.end(), Xr.begin());
  NumericVector eta(eta_arma.begin(), eta_arma.end());
  NumericVector yr(y.begin(), y.end());
  ConditionalLogitHessWorker worker(Xr, eta, yr, group_start, group_size);
  parallelReduce(0, group_start.size(), worker, panglm::grain(group_start.size()));
  NumericMatrix H(k, k);
  for (int i = 0; i < k; ++i) for (int j = 0; j <= i; ++j) {
    H(i, j) = worker.hess[i * k + j];
    H(j, i) = H(i, j);
  }
  return List::create(Named("loglik") = worker.ll, Named("gradient") = worker.grad,
                      Named("hessian") = H, Named("n_used_groups") = worker.n_used_groups);
}
