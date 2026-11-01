#include <Rcpp.h>
#include <cmath>
#include <vector>
#include <algorithm>
using namespace Rcpp;

// Correlated random-effects two-part model (see R/twopartcre.R).
//
// Individual g has binary indicators d_t (all observations) and log costs
// ly_t (positive observations). Given the binary-part intercept u,
//   h(u) = sum_t log F(s_t (a_t + u)) + P(u) + log phi(u; 0, sigma1),
// where s_t = 2 d_t - 1, a = Xb gamma, and P(u) is the closed-form Gaussian
// random-intercept log density of the positive part with residuals
// e_t = r_t - c u (r = ly - Xp beta), intercept variance tv and residual
// variance se2, including the Jacobian -sum(ly). With m positive
// observations, S = sum e_t, D = se2 + m tv, k = tv / D and
// Q = sum e_t^2 - k S^2:
//   P(u) = -(m log(2 pi) + (m - 1) log se2 + log D) / 2 - Q / (2 se2) - sum ly.
// The integral over u uses adaptive Gauss-Hermite quadrature at the mode u*
// of h with scale tau = (-h''(u*))^(-1/2):
//   ll = log(sqrt(2) tau) + log sum_q exp(log_w_q + h(u* + sqrt(2) tau z_q)).
// The score is the exact gradient of this approximation, including the
// dependence of u* and tau on the parameters (implicit function theorem).
//
// Derivatives are first taken in the working parameters
//   (gamma, beta, c, tv, se2, l1), with c = rho sigma2 / sigma1,
//   tv = sigma2^2 (1 - rho^2), se2 = sigma_e^2, l1 = log sigma1 (direct),
// and then mapped to theta = (gamma, beta, log sigma1, log sigma2,
// log sigma_e, atanh rho).

namespace {

// log F(z) and its first three derivatives in z
inline void logF_derivs(double z, bool probit, double& l0, double& l1, double& l2,
                        double& l3) {
  if (probit) {
    l0 = R::pnorm(z, 0.0, 1.0, 1, 1);
    double lam = std::exp(R::dnorm(z, 0.0, 1.0, 1) - l0);
    l1 = lam;
    l2 = -lam * (z + lam);
    l3 = -(l2 * (z + lam) + lam * (1.0 + l2));
  } else {
    l0 = (z > 0) ? -std::log1p(std::exp(-z)) : z - std::log1p(std::exp(z));
    double p = std::exp(l0);
    double q = (z > 0) ? std::exp(-z) * p : 1.0 / (1.0 + std::exp(z));  // 1 - p
    l1 = q;
    l2 = -p * q;
    l3 = -p * q * (q - p);
  }
}

} // namespace

// [[Rcpp::export]]
List twopart_cre_cpp(const NumericMatrix& Xb, const NumericVector& sgn,
                     const IntegerVector& group_b, const NumericMatrix& Xp,
                     const NumericVector& ly, const IntegerVector& group_p, int G,
                     const NumericVector& theta, bool probit,
                     const NumericVector& nodes, const NumericVector& log_w,
                     bool score) {
  const int nb = Xb.nrow(), kb = Xb.ncol(), np = Xp.nrow(), kp = Xp.ncol();
  const int nq = nodes.size();
  const int p = kb + kp + 4;
  const int iv = kb + kp;  // offset of the variance parameters
  NumericVector ll(G);
  NumericMatrix sc(score ? G : 0, score ? p : 0);

  const double s1 = std::exp(theta[iv]), s2 = std::exp(theta[iv + 1]);
  const double se2 = std::exp(2.0 * theta[iv + 2]), rho = std::tanh(theta[iv + 3]);
  const double c = rho * s2 / s1, tv = s2 * s2 * (1.0 - rho * rho);
  const double is1sq = 1.0 / (s1 * s1);
  const double log2pi = std::log(2.0 * M_PI);
  const double sqrt2 = std::sqrt(2.0);

  // linear predictors and residuals
  std::vector<double> a(nb, 0.0), r(np, 0.0);
  for (int j = 0; j < kb; ++j) {
    const double g = theta[j];
    for (int i = 0; i < nb; ++i) a[i] += Xb(i, j) * g;
  }
  for (int i = 0; i < np; ++i) r[i] = ly[i];
  for (int j = 0; j < kp; ++j) {
    const double b = theta[kb + j];
    for (int i = 0; i < np; ++i) r[i] -= Xp(i, j) * b;
  }

  // observations by individual (1-based group ids)
  std::vector<int> cb(G + 1, 0), cp(G + 1, 0);
  for (int i = 0; i < nb; ++i) cb[group_b[i]]++;
  for (int i = 0; i < np; ++i) cp[group_p[i]]++;
  for (int g = 0; g < G; ++g) { cb[g + 1] += cb[g]; cp[g + 1] += cp[g]; }
  std::vector<int> ob(nb), op(np);
  {
    std::vector<int> pos(cb.begin(), cb.end() - 1);
    for (int i = 0; i < nb; ++i) ob[pos[group_b[i] - 1]++] = i;
    std::vector<int> pp(cp.begin(), cp.end() - 1);
    for (int i = 0; i < np; ++i) op[pp[group_p[i] - 1]++] = i;
  }

  std::vector<double> hq(nq), uq(nq), pq(nq);
  std::vector<double> fix(p), hut(p), huut(p), du(p);
  std::vector<double> Xs(kp), Rx(kp);

  for (int g = 0; g < G; ++g) {
    const int b0 = cb[g], b1 = cb[g + 1], p0 = cp[g], p1 = cp[g + 1];
    const double m = p1 - p0;
    double R1 = 0.0, R2 = 0.0, sly = 0.0;
    for (int t = p0; t < p1; ++t) {
      const int i = op[t];
      R1 += r[i]; R2 += r[i] * r[i]; sly += ly[i];
    }
    const double D = se2 + m * tv;
    const double k = tv / D;

    // h and its u-derivatives
    auto h_at = [&](double u, double& h1, double& h2) {
      double h = 0.0; h1 = 0.0; h2 = 0.0;
      for (int t = b0; t < b1; ++t) {
        const int i = ob[t];
        const double z = sgn[i] * (a[i] + u);
        double l0, l1, l2, l3;
        logF_derivs(z, probit, l0, l1, l2, l3);
        h += l0; h1 += sgn[i] * l1; h2 += l2;
      }
      if (m > 0) {
        const double S = R1 - c * m * u;
        const double Q = R2 - 2.0 * c * u * R1 + c * c * u * u * m - k * S * S;
        h += -0.5 * (m * log2pi + (m - 1.0) * std::log(se2) + std::log(D)) -
          0.5 * Q / se2 - sly;
        h1 += c * S / D;
        h2 -= c * c * m / D;
      }
      h += -0.5 * log2pi - theta[iv] - 0.5 * u * u * is1sq;
      h1 -= u * is1sq;
      h2 -= is1sq;
      return h;
    };

    // mode by Newton steps (clipped, as for the other quadrature fits)
    double u = 0.0, h1, h2;
    bool ok = true;
    for (int step = 0; step < 100; ++step) {
      h_at(u, h1, h2);
      double delta = h1 / h2;
      if (!std::isfinite(delta)) { ok = false; break; }
      delta = std::max(-3.0, std::min(3.0, delta));
      u -= delta;
      if (std::fabs(delta) < 1e-10) break;
    }
    if (!ok) {
      std::fill(ll.begin(), ll.end(), R_NegInf);
      return List::create(Named("ll") = ll, Named("score") = R_NilValue);
    }
    const double umode = u;
    h_at(umode, h1, h2);
    const double tau = 1.0 / std::sqrt(-h2);

    double top = R_NegInf;
    for (int q = 0; q < nq; ++q) {
      double d1, d2;
      uq[q] = umode + sqrt2 * tau * nodes[q];
      hq[q] = log_w[q] + h_at(uq[q], d1, d2);
      pq[q] = d1;  // h'(u_q), kept for the score
      if (hq[q] > top) top = hq[q];
    }
    double tot = 0.0;
    for (int q = 0; q < nq; ++q) tot += std::exp(hq[q] - top);
    ll[g] = std::log(sqrt2 * tau) + top + std::log(tot);
    if (!score || !std::isfinite(ll[g])) continue;

    // node weights pi_q, and J_u = sum pi h'(u_q), J_tau
    double Ju = 0.0, Jt = 1.0 / tau, ubar = 0.0;
    for (int q = 0; q < nq; ++q) {
      const double w = std::exp(hq[q] - top) / tot;
      Ju += w * pq[q];
      Jt += w * pq[q] * sqrt2 * nodes[q];
      ubar += w * uq[q];
      hq[q] = w;  // reuse as the weight
    }

    std::fill(fix.begin(), fix.end(), 0.0);
    std::fill(hut.begin(), hut.end(), 0.0);
    std::fill(huut.begin(), huut.end(), 0.0);

    // binary part
    double huuu = 0.0;
    for (int t = b0; t < b1; ++t) {
      const int i = ob[t];
      double l0, l1, l2, l3;
      double fsum = 0.0;
      for (int q = 0; q < nq; ++q) {
        logF_derivs(sgn[i] * (a[i] + uq[q]), probit, l0, l1, l2, l3);
        fsum += hq[q] * l1;
      }
      logF_derivs(sgn[i] * (a[i] + umode), probit, l0, l1, l2, l3);
      huuu += sgn[i] * l3;
      for (int j = 0; j < kb; ++j) {
        const double x = Xb(i, j);
        fix[j] += sgn[i] * fsum * x;
        hut[j] += l2 * x;
        huut[j] += sgn[i] * l3 * x;
      }
    }

    // positive part
    if (m > 0) {
      std::fill(Xs.begin(), Xs.end(), 0.0);
      std::fill(Rx.begin(), Rx.end(), 0.0);
      for (int t = p0; t < p1; ++t) {
        const int i = op[t];
        for (int j = 0; j < kp; ++j) { Xs[j] += Xp(i, j); Rx[j] += r[i] * Xp(i, j); }
      }
      for (int j = 0; j < kp; ++j) {
        fix[kb + j] = (Rx[j] - k * R1 * Xs[j]) / se2 - ubar * c * Xs[j] / D;
        hut[kb + j] = -c * Xs[j] / D;
      }
      double fc = 0.0, ftv = 0.0, fse = 0.0;
      for (int q = 0; q < nq; ++q) {
        const double uu = uq[q];
        const double S = R1 - c * m * uu;
        const double Q = R2 - 2.0 * c * uu * R1 + c * c * uu * uu * m - k * S * S;
        fc += hq[q] * uu * S / D;
        ftv += hq[q] * (-0.5 * m / D + 0.5 * S * S / (D * D));
        fse += hq[q] * (-0.5 * ((m - 1.0) / se2 + 1.0 / D) + 0.5 * Q / (se2 * se2) -
          0.5 * tv * S * S / (D * D * se2));
      }
      fix[iv] = fc; fix[iv + 1] = ftv; fix[iv + 2] = fse;
      const double S = R1 - c * m * umode;
      hut[iv] = (S - c * m * umode) / D;
      hut[iv + 1] = -c * S * m / (D * D);
      hut[iv + 2] = -c * S / (D * D);
      huut[iv] = -2.0 * c * m / D;
      huut[iv + 1] = c * c * m * m / (D * D);
      huut[iv + 2] = c * c * m / (D * D);
    }

    // log phi(u; sigma1), direct dependence on l1
    double fl1 = 0.0;
    for (int q = 0; q < nq; ++q) fl1 += hq[q] * (-1.0 + uq[q] * uq[q] * is1sq);
    fix[iv + 3] = fl1;
    hut[iv + 3] = 2.0 * umode * is1sq;
    huut[iv + 3] = 2.0 * is1sq;

    // total derivative in the working parameters
    const double tau3 = tau * tau * tau;
    for (int j = 0; j < p; ++j) {
      const double dum = -hut[j] / h2;
      const double dtau = 0.5 * tau3 * (huut[j] + huuu * dum);
      fix[j] += Ju * dum + Jt * dtau;
    }

    // map (c, tv, se2, l1) to (log sigma1, log sigma2, log sigma_e, atanh rho)
    const double dc = fix[iv], dtv = fix[iv + 1], dse = fix[iv + 2], dl1 = fix[iv + 3];
    for (int j = 0; j < iv; ++j) sc(g, j) = fix[j];
    sc(g, iv) = -c * dc + dl1;
    sc(g, iv + 1) = c * dc + 2.0 * tv * dtv;
    sc(g, iv + 2) = 2.0 * se2 * dse;
    sc(g, iv + 3) = (s2 / s1) * (1.0 - rho * rho) * dc -
      2.0 * rho * s2 * s2 * (1.0 - rho * rho) * dtv;
  }
  if (!score) return List::create(Named("ll") = ll, Named("score") = R_NilValue);
  return List::create(Named("ll") = ll, Named("score") = sc);
}
