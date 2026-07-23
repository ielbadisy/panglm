#include <RcppArmadillo.h>
#include <RcppParallel.h>
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
  int n_used_groups;

  ConditionalLogitWorker(const NumericMatrix& X, const NumericVector& eta, const NumericVector& y,
                         const IntegerVector& group_start, const IntegerVector& group_size)
    : X(X), eta(eta), y(y), group_start(group_start), group_size(group_size),
      ll(0.0), grad(X.ncol(), 0.0), n_used_groups(0) {}

  ConditionalLogitWorker(const ConditionalLogitWorker& other, Split)
    : X(other.X), eta(other.eta), y(other.y), group_start(other.group_start), group_size(other.group_size),
      ll(0.0), grad(other.X.ncol(), 0.0), n_used_groups(0) {}

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
      for (int t = 0; t < size; ++t)
        for (std::size_t j = 0; j < k; ++j)
          grad[j] += res.grad_eta[t] * X(start + t, j);
    }
  }

  void join(const ConditionalLogitWorker& rhs) {
    ll += rhs.ll;
    n_used_groups += rhs.n_used_groups;
    for (std::size_t j = 0; j < grad.size(); ++j) grad[j] += rhs.grad[j];
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
  parallelReduce(0, group_start.size(), worker);

  return List::create(
    Named("loglik") = worker.ll,
    Named("gradient") = worker.grad,
    Named("n_used_groups") = worker.n_used_groups
  );
}
