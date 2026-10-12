// Dormand-Prince 5(4) integrator (D-021).
//
// Explicit embedded Runge-Kutta pair of Dormand & Prince (1980) with the step
// size control, starting step and stiffness detection of Hairer's DOPRI5
// (Hairer, Norsett & Wanner, Solving Ordinary Differential Equations I,
// 2nd ed., Springer 1993, section II.4 and the DOPRI5 code):
// * 7 stages, "first same as last": 6 new derivative evaluations per step;
// * error norm: root mean square over components of err_j / sk_j with
//   sk_j = atol + rtol * max(|y_j|, |y_new_j|);
// * PI step size control (beta = 0.04), safety factor 0.9, step ratio
//   limited to [0.2, 10];
// * stiffness test: h * |lambda| estimated from the last two stages; more than
//   15 steps with an estimate above 3.25 (the stability boundary) flag the
//   problem as stiff (6 consecutive non-stiff steps reset the count).
// A fixed-step mode (fixed_h > 0, no error control) exists for convergence
// tests.
//
// The right-hand side is any callable f(t, y, dydt) with y and dydt of
// length n. The integrator allocates its work vectors once and can be reused
// for many integrations of the same size.

#ifndef SGDRIVE_DOPRI5_H
#define SGDRIVE_DOPRI5_H

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

namespace sgdrive {

struct Dopri5Options {
  double rtol = 1e-8;
  double atol = 1e-10;
  long max_steps = 100000;  // per call to integrate()
  double hmax = 0.0;        // 0: the length of the interval
  double fixed_h = 0.0;     // > 0: fixed steps of this size, no error control
  bool detect_stiff = true;
  int stiff_steps = 15;     // Hairer's threshold
  bool stop_when_stiff = false;  // return (status 3) as soon as stiffness is detected
};

struct Dopri5Stats {
  long steps = 0;        // attempted steps
  long accepted = 0;
  long rejected = 0;
  long evaluations = 0;  // calls of the right-hand side
  bool stiff = false;    // stiffness detected during this call
  double h_last = 0.0;   // last accepted step size
  int status = 0;        // 0 ok, 1 too many steps, 2 step size too small, 3 stopped: stiff
};

class Dopri5 {
 public:
  explicit Dopri5(int n)
      : n_(n), k1_(n), k2_(n), k3_(n), k4_(n), k5_(n), k6_(n), k7_(n), ytmp_(n),
        ynew_(n), ysti_(n) {}

  // Integrate y from t0 to t1 in place.
  template <class F>
  Dopri5Stats integrate(F& f, double t0, double t1, double* y, const Dopri5Options& opt) {
    Dopri5Stats st;
    const double span = t1 - t0;
    if (span <= 0.0) return st;
    const double uround = std::numeric_limits<double>::epsilon();
    const double hmax = opt.hmax > 0.0 ? std::min(opt.hmax, span) : span;

    f(t0, y, k1_.data());
    st.evaluations++;
    double h;
    if (opt.fixed_h > 0.0) {
      h = std::min(opt.fixed_h, span);
    } else {
      h = initial_step(f, t0, y, hmax, opt, st);
    }

    double t = t0;
    double facold = 1e-4;
    bool last = false, reject = false;
    int iasti = 0, nonsti = 0;

    while (true) {
      if (st.steps >= opt.max_steps) {
        st.status = 1;
        return st;
      }
      if (0.1 * std::fabs(h) <= std::fabs(t) * uround) {
        st.status = 2;
        return st;
      }
      if ((t + 1.01 * h - t1) > 0.0) {
        h = t1 - t;
        last = true;
      }
      st.steps++;

      step(f, t, h, y, st);  // fills ynew_, ysti_, k7_ (= f(t + h, ynew))

      if (opt.fixed_h > 0.0) {
        st.accepted++;
        accept(t, h, y, st);
        if (last) return st;
        continue;
      }

      const double err = error_norm(h, y, opt);
      // PI step size control (Hairer's DOPRI5: beta = 0.04, expo1 = 0.2 - 0.75 beta)
      const double beta = 0.04, expo1 = 0.2 - beta * 0.75, safe = 0.9;
      const double facc1 = 1.0 / 0.2, facc2 = 1.0 / 10.0;
      const double fac11 = std::pow(err, expo1);
      double fac = fac11 / std::pow(facold, beta);
      fac = std::max(facc2, std::min(facc1, fac / safe));
      double hnew = h / fac;

      if (err <= 1.0) {
        facold = std::max(err, 1e-4);
        st.accepted++;
        if (opt.detect_stiff) {
          double stnum = 0.0, stden = 0.0;
          for (int i = 0; i < n_; ++i) {
            const double a = k7_[i] - k6_[i];
            const double b = ynew_[i] - ysti_[i];
            stnum += a * a;
            stden += b * b;
          }
          if (stden > 0.0) {
            const double hlamb = h * std::sqrt(stnum / stden);
            if (hlamb > 3.25) {
              nonsti = 0;
              if (++iasti >= opt.stiff_steps) {
                st.stiff = true;
                if (opt.stop_when_stiff) {
                  st.status = 3;
                  return st;
                }
              }
            } else if (++nonsti >= 6) {
              iasti = 0;
            }
          }
        }
        accept(t, h, y, st);
        if (last) return st;
        hnew = std::min(std::fabs(hnew), hmax);
        if (reject) hnew = std::min(hnew, h);
        reject = false;
      } else {
        hnew = h / std::min(facc1, fac11 / safe);
        reject = true;
        if (st.accepted >= 1) st.rejected++;
        last = false;
      }
      h = hnew;
    }
  }

 private:
  int n_;
  std::vector<double> k1_, k2_, k3_, k4_, k5_, k6_, k7_, ytmp_, ynew_, ysti_;

  // Dormand-Prince coefficients
  static constexpr double c2 = 1.0 / 5.0, c3 = 3.0 / 10.0, c4 = 4.0 / 5.0, c5 = 8.0 / 9.0;
  static constexpr double a21 = 1.0 / 5.0;
  static constexpr double a31 = 3.0 / 40.0, a32 = 9.0 / 40.0;
  static constexpr double a41 = 44.0 / 45.0, a42 = -56.0 / 15.0, a43 = 32.0 / 9.0;
  static constexpr double a51 = 19372.0 / 6561.0, a52 = -25360.0 / 2187.0,
                          a53 = 64448.0 / 6561.0, a54 = -212.0 / 729.0;
  static constexpr double a61 = 9017.0 / 3168.0, a62 = -355.0 / 33.0, a63 = 46732.0 / 5247.0,
                          a64 = 49.0 / 176.0, a65 = -5103.0 / 18656.0;
  static constexpr double a71 = 35.0 / 384.0, a73 = 500.0 / 1113.0, a74 = 125.0 / 192.0,
                          a75 = -2187.0 / 6784.0, a76 = 11.0 / 84.0;
  // error coefficients (5th minus embedded 4th order weights)
  static constexpr double e1 = 71.0 / 57600.0, e3 = -71.0 / 16695.0, e4 = 71.0 / 1920.0,
                          e5 = -17253.0 / 339200.0, e6 = 22.0 / 525.0, e7 = -1.0 / 40.0;

  template <class F>
  void step(F& f, double t, double h, const double* y, Dopri5Stats& st) {
    for (int i = 0; i < n_; ++i) ytmp_[i] = y[i] + h * a21 * k1_[i];
    f(t + c2 * h, ytmp_.data(), k2_.data());
    for (int i = 0; i < n_; ++i) ytmp_[i] = y[i] + h * (a31 * k1_[i] + a32 * k2_[i]);
    f(t + c3 * h, ytmp_.data(), k3_.data());
    for (int i = 0; i < n_; ++i) {
      ytmp_[i] = y[i] + h * (a41 * k1_[i] + a42 * k2_[i] + a43 * k3_[i]);
    }
    f(t + c4 * h, ytmp_.data(), k4_.data());
    for (int i = 0; i < n_; ++i) {
      ytmp_[i] = y[i] + h * (a51 * k1_[i] + a52 * k2_[i] + a53 * k3_[i] + a54 * k4_[i]);
    }
    f(t + c5 * h, ytmp_.data(), k5_.data());
    for (int i = 0; i < n_; ++i) {
      ysti_[i] = y[i] + h * (a61 * k1_[i] + a62 * k2_[i] + a63 * k3_[i] + a64 * k4_[i] +
                             a65 * k5_[i]);
    }
    f(t + h, ysti_.data(), k6_.data());
    for (int i = 0; i < n_; ++i) {
      ynew_[i] = y[i] + h * (a71 * k1_[i] + a73 * k3_[i] + a74 * k4_[i] + a75 * k5_[i] +
                             a76 * k6_[i]);
    }
    f(t + h, ynew_.data(), k7_.data());
    st.evaluations += 6;
  }

  double error_norm(double h, const double* y, const Dopri5Options& opt) const {
    double sum = 0.0;
    for (int i = 0; i < n_; ++i) {
      const double e = h * (e1 * k1_[i] + e3 * k3_[i] + e4 * k4_[i] + e5 * k5_[i] +
                            e6 * k6_[i] + e7 * k7_[i]);
      const double sk = opt.atol + opt.rtol * std::max(std::fabs(y[i]), std::fabs(ynew_[i]));
      const double r = e / sk;
      sum += r * r;
    }
    return std::sqrt(sum / n_);
  }

  void accept(double& t, double h, double* y, Dopri5Stats& st) {
    for (int i = 0; i < n_; ++i) {
      y[i] = ynew_[i];
      k1_[i] = k7_[i];  // first same as last
    }
    t += h;
    st.h_last = h;
  }

  // Hairer's starting step (DOPRI5, HINIT), using k1_ = f(t0, y).
  template <class F>
  double initial_step(F& f, double t0, const double* y, double hmax, const Dopri5Options& opt,
                      Dopri5Stats& st) {
    double dnf = 0.0, dny = 0.0;
    for (int i = 0; i < n_; ++i) {
      const double sk = opt.atol + opt.rtol * std::fabs(y[i]);
      dnf += (k1_[i] / sk) * (k1_[i] / sk);
      dny += (y[i] / sk) * (y[i] / sk);
    }
    double h = (dnf <= 1e-10 || dny <= 1e-10) ? 1e-6 : std::sqrt(dny / dnf) * 0.01;
    h = std::min(h, hmax);
    for (int i = 0; i < n_; ++i) ytmp_[i] = y[i] + h * k1_[i];
    f(t0 + h, ytmp_.data(), k2_.data());
    st.evaluations++;
    double der2 = 0.0;
    for (int i = 0; i < n_; ++i) {
      const double sk = opt.atol + opt.rtol * std::fabs(y[i]);
      const double r = (k2_[i] - k1_[i]) / sk;
      der2 += r * r;
    }
    der2 = std::sqrt(der2) / h;
    const double der12 = std::max(std::fabs(der2), std::sqrt(dnf));
    const double h1 =
        der12 <= 1e-15 ? std::max(1e-6, std::fabs(h) * 1e-3) : std::pow(0.01 / der12, 0.2);
    return std::min({100.0 * std::fabs(h), h1, hmax});
  }
};

}  // namespace sgdrive

#endif  // SGDRIVE_DOPRI5_H
