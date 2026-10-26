#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include "threads.h"
using namespace Rcpp;
using namespace RcppParallel;

// Adaptive Gauss-Hermite quadrature for random-intercept models whose
// observations depend on the individual effect u_i ~ N(0, sigma^2) only
// through e_it = x_it'beta + u_i: binomial (logit, probit) and ordered
// (logit, probit) outcomes.
//
// For individual i let
//   h(u) = sum_t log p(y_it | x_it'beta + u) + log phi(u; 0, sigma).
// The nodes are centred at the mode u* of h and scaled by
// tau = (-h''(u*))^(-1/2), the observed curvature (as in lme4, ordinal and
// GLMMadaptive):
//   L_i ~= sqrt(2) tau sum_q w_q exp(z_q^2) exp(h(u* + sqrt(2) tau z_q)),
// with (z_q, w_q) the Gauss-Hermite rule. One node gives the Laplace
// approximation. The gradient returned is the exact gradient of this
// approximation, including the dependence of u* and tau on the parameters
// (implicit function theorem), so that a Newton or quasi-Newton search
// converges to the maximizer of the approximate likelihood.
//
// Parameters: theta = (beta [k], cuts [m], log sigma); m = 0 for the
// binomial model. Output: one log-likelihood and one score row per
// individual.

namespace {

enum Kind { BIN_LOGIT = 0, BIN_PROBIT = 1, ORD_LOGIT = 2, ORD_PROBIT = 3 };

struct ObsDeriv {
  double lp, d1, d2, d3;   // log p and its derivatives in e
  int nc;                  // number of cutpoints involved (0, 1 or 2)
  int ci[2];               // their (0-based) indices
  double lc[2], lce[2], lcee[2]; // d/dc, d2/(dc de), d3/(dc de2) of log p
};

// log Phi(q) and its first three derivatives in q
inline void log_phi_derivs(double q, double& l, double& l1, double& l2, double& l3) {
  l = R::pnorm(q, 0.0, 1.0, 1, 1);
  double lam = std::exp(R::dnorm(q, 0.0, 1.0, 1) - l);
  l1 = lam;
  l2 = -lam * (q + lam);
  l3 = -(l2 * (q + lam) + lam * (1.0 + l2));
}

// density f and its first two derivatives at x (zero at +-Inf)
inline void dens_derivs(double x, bool logit, double& f, double& f1, double& f2) {
  if (!std::isfinite(x)) { f = f1 = f2 = 0.0; return; }
  if (logit) {
    double F = 1.0 / (1.0 + std::exp(-x));
    f = F * (1.0 - F);
    f1 = f * (1.0 - 2.0 * F);
    f2 = f * (1.0 - 2.0 * F) * (1.0 - 2.0 * F) - 2.0 * f * f;
  } else {
    f = R::dnorm(x, 0.0, 1.0, 0);
    f1 = -x * f;
    f2 = (x * x - 1.0) * f;
  }
}

inline double cdf(double x, bool logit) {
  if (x == R_PosInf) return 1.0;
  if (x == R_NegInf) return 0.0;
  return logit ? 1.0 / (1.0 + std::exp(-x)) : R::pnorm(x, 0.0, 1.0, 1, 0);
}

inline ObsDeriv obs_deriv(int kind, double y, double e, const std::vector<double>& cuts) {
  ObsDeriv o;
  o.nc = 0;
  if (kind == BIN_LOGIT) {
    // log mu = -log(1 + exp(-e)), log(1 - mu) = -log(1 + exp(e))
    double lmu = (e > 0) ? -std::log1p(std::exp(-e)) : e - std::log1p(std::exp(e));
    double l1mu = (e > 0) ? -e - std::log1p(std::exp(-e)) : -std::log1p(std::exp(e));
    double mu = std::exp(lmu);
    o.lp = 0.0;
    if (y != 0.0) o.lp += y * lmu;
    if (y != 1.0) o.lp += (1.0 - y) * l1mu;
    double v = mu * (1.0 - mu);
    o.d1 = y - mu;
    o.d2 = -v;
    o.d3 = -v * (1.0 - 2.0 * mu);
  } else if (kind == BIN_PROBIT) {
    o.lp = o.d1 = o.d2 = o.d3 = 0.0;
    double l, l1, l2, l3;
    if (y != 0.0) {
      log_phi_derivs(e, l, l1, l2, l3);
      o.lp += y * l; o.d1 += y * l1; o.d2 += y * l2; o.d3 += y * l3;
    }
    if (y != 1.0) {
      log_phi_derivs(-e, l, l1, l2, l3);
      double w = 1.0 - y;
      o.lp += w * l; o.d1 -= w * l1; o.d2 += w * l2; o.d3 -= w * l3;
    }
  } else {
    bool logit = (kind == ORD_LOGIT);
    int J = (int) cuts.size() + 1;
    int yc = (int) y; // category code 1..J
    double a = (yc <= J - 1) ? cuts[yc - 1] - e : R_PosInf;
    double b = (yc >= 2) ? cuts[yc - 2] - e : R_NegInf;
    double P = (b > 0) ? cdf(-b, logit) - cdf(-a, logit) : cdf(a, logit) - cdf(b, logit);
    P = std::max(P, 1e-300);
    double fa, fa1, fa2, fb, fb1, fb2;
    dens_derivs(a, logit, fa, fa1, fa2);
    dens_derivs(b, logit, fb, fb1, fb2);
    double Pe = -(fa - fb), Pee = fa1 - fb1, Peee = -(fa2 - fb2);
    o.lp = std::log(P);
    o.d1 = Pe / P;
    o.d2 = Pee / P - o.d1 * o.d1;
    o.d3 = Peee / P - 3.0 * Pee * Pe / (P * P) + 2.0 * o.d1 * o.d1 * o.d1;
    auto add_cut = [&](int idx, double Pc, double Pce, double Pcee) {
      int s = o.nc++;
      o.ci[s] = idx;
      o.lc[s] = Pc / P;
      o.lce[s] = Pce / P - Pc * Pe / (P * P);
      o.lcee[s] = Pcee / P - 2.0 * Pce * Pe / (P * P) - Pc * Pee / (P * P) +
        2.0 * Pc * Pe * Pe / (P * P * P);
    };
    if (yc <= J - 1) add_cut(yc - 1, fa, -fa1, fa2);
    if (yc >= 2) add_cut(yc - 2, -fb, fb1, -fb2);
  }
  return o;
}

struct AGQWorker : public Worker {
  const RMatrix<double> X;
  const RVector<double> y;
  const RVector<int> gstart;
  const RVector<int> gsize;
  const std::vector<double>& beta;
  const std::vector<double>& cuts;
  double sigma;
  const std::vector<double>& z;
  const std::vector<double>& logc; // log(w_q sqrt(pi)) + z_q^2
  int kind;
  RVector<double> ll;
  RMatrix<double> score;

  AGQWorker(const NumericMatrix& X, const NumericVector& y, const IntegerVector& gstart,
            const IntegerVector& gsize, const std::vector<double>& beta,
            const std::vector<double>& cuts, double sigma, const std::vector<double>& z,
            const std::vector<double>& logc, int kind, NumericVector ll, NumericMatrix score)
    : X(X), y(y), gstart(gstart), gsize(gsize), beta(beta), cuts(cuts), sigma(sigma),
      z(z), logc(logc), kind(kind), ll(ll), score(score) {}

  void operator()(std::size_t begin, std::size_t end) {
    const std::size_t k = X.ncol();
    const std::size_t m = cuts.size();
    const std::size_t p = k + m + 1;
    const std::size_t R = z.size();
    const double s2 = sigma * sigma;
    const double sqrt2 = std::sqrt(2.0);
    const double lnorm = -0.5 * std::log(2.0 * M_PI) - std::log(sigma);
    std::vector<double> eta, hq(R), huq(R), dmode(p), dlt(p), grad(p);
    std::vector<std::vector<double>> htq(R, std::vector<double>(p));
    std::vector<double> Hu_t(p), Huu_t(p);

    for (std::size_t g = begin; g < end; ++g) {
      const int st = gstart[g], T = gsize[g];
      eta.assign(T, 0.0);
      for (int t = 0; t < T; ++t) {
        double s = 0.0;
        for (std::size_t j = 0; j < k; ++j) s += X(st + t, j) * beta[j];
        eta[t] = s;
      }
      auto hval = [&](double u) {
        double s = lnorm - 0.5 * u * u / s2;
        for (int t = 0; t < T; ++t) s += obs_deriv(kind, y[st + t], eta[t] + u, cuts).lp;
        return s;
      };
      // mode of the strictly concave h by damped Newton
      double u = 0.0, hu_cur = hval(u);
      for (int it = 0; it < 200; ++it) {
        double g1 = -u / s2, g2 = -1.0 / s2;
        for (int t = 0; t < T; ++t) {
          ObsDeriv o = obs_deriv(kind, y[st + t], eta[t] + u, cuts);
          g1 += o.d1; g2 += o.d2;
        }
        double step = g1 / g2;
        if (!std::isfinite(step)) break;
        double lim = 10.0 * std::max(1.0, sigma);
        if (step > lim) step = lim;
        if (step < -lim) step = -lim;
        double un = u - step, hn = hval(un);
        int halv = 0;
        while ((!std::isfinite(hn) || hn < hu_cur - 1e-12 * std::fabs(hu_cur)) && halv < 40) {
          step /= 2.0; un = u - step; hn = hval(un); ++halv;
        }
        u = un; hu_cur = hn;
        if (std::fabs(step) < 1e-11 * (1.0 + std::fabs(u))) break;
      }
      // derivatives at the mode
      double H2 = -1.0 / s2, H3 = 0.0;
      std::fill(Hu_t.begin(), Hu_t.end(), 0.0);
      std::fill(Huu_t.begin(), Huu_t.end(), 0.0);
      for (int t = 0; t < T; ++t) {
        ObsDeriv o = obs_deriv(kind, y[st + t], eta[t] + u, cuts);
        H2 += o.d2; H3 += o.d3;
        for (std::size_t j = 0; j < k; ++j) {
          double x = X(st + t, j);
          Hu_t[j] += o.d2 * x; Huu_t[j] += o.d3 * x;
        }
        for (int c = 0; c < o.nc; ++c) {
          Hu_t[k + o.ci[c]] += o.lce[c]; Huu_t[k + o.ci[c]] += o.lcee[c];
        }
      }
      Hu_t[p - 1] = 2.0 * u / s2;
      Huu_t[p - 1] = 2.0 / s2;
      const double tau = 1.0 / std::sqrt(-H2);
      for (std::size_t j = 0; j < p; ++j) {
        dmode[j] = -Hu_t[j] / H2;
        dlt[j] = -0.5 * (Huu_t[j] + H3 * dmode[j]) / H2;
      }
      // the nodes
      double top = R_NegInf;
      for (std::size_t q = 0; q < R; ++q) {
        double b = u + sqrt2 * tau * z[q];
        double h = lnorm - 0.5 * b * b / s2, hu = -b / s2;
        std::vector<double>& ht = htq[q];
        std::fill(ht.begin(), ht.end(), 0.0);
        for (int t = 0; t < T; ++t) {
          ObsDeriv o = obs_deriv(kind, y[st + t], eta[t] + b, cuts);
          h += o.lp; hu += o.d1;
          for (std::size_t j = 0; j < k; ++j) ht[j] += o.d1 * X(st + t, j);
          for (int c = 0; c < o.nc; ++c) ht[k + o.ci[c]] += o.lc[c];
        }
        ht[p - 1] = b * b / s2 - 1.0;
        hq[q] = logc[q] + h;
        huq[q] = hu;
        if (hq[q] > top) top = hq[q];
      }
      double sum = 0.0;
      for (std::size_t q = 0; q < R; ++q) sum += std::exp(hq[q] - top);
      const double lse = top + std::log(sum);
      ll[g] = std::log(sqrt2 * tau) + lse;
      std::fill(grad.begin(), grad.end(), 0.0);
      for (std::size_t q = 0; q < R; ++q) {
        double pq = std::exp(hq[q] - lse);
        if (pq == 0.0) continue;
        for (std::size_t j = 0; j < p; ++j) {
          grad[j] += pq * (htq[q][j] + huq[q] * (dmode[j] + sqrt2 * z[q] * tau * dlt[j]));
        }
      }
      for (std::size_t j = 0; j < p; ++j) score(g, j) = dlt[j] + grad[j];
    }
  }
};

} // namespace

// [[Rcpp::export]]
List random_agq_cpp(NumericVector beta, NumericVector cuts, double sigma,
                    NumericMatrix X, NumericVector y,
                    IntegerVector group_start, IntegerVector group_size,
                    NumericVector nodes, NumericVector weights, int kind) {
  const int G = group_start.size();
  const int p = X.ncol() + cuts.size() + 1;
  std::vector<double> b(beta.begin(), beta.end()), c(cuts.begin(), cuts.end());
  std::vector<double> z(nodes.begin(), nodes.end()), logc(nodes.size());
  for (int q = 0; q < nodes.size(); ++q) {
    logc[q] = std::log(weights[q] * std::sqrt(M_PI)) + z[q] * z[q];
  }
  NumericVector ll(G);
  NumericMatrix score(G, p);
  AGQWorker worker(X, y, group_start, group_size, b, c, sigma, z, logc, kind, ll, score);
  parallelFor(0, G, worker, panglm::grain(G));
  return List::create(Named("ll") = ll, Named("score") = score);
}
