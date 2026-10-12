// Exact advancement along a chain of sub-stages over one time step (D-022,
// model specification S16), for the chain-binomial solver.
//
// A chain has classes 0..n-1. Individuals in class i < last move to class
// i + 1 at rate q[i]; from class last - 1 they leave the chain ("exit": to
// the next life stage, or death by old age for adults). Classes >= last are
// never reached (D-018). Every class of the chain has the same hazard mu;
// death acts only while an individual is in the chain, so an individual that
// leaves it during the step no longer dies at this chain's rate.
//
// advance() fills P, row-major with n rows and n + 1 columns: P[i][j] is the
// probability that an individual in class i at the start of the step is
// alive in class j (i <= j < last) at its end, and P[i][n] that it has left
// the chain alive; 1 - sum_j P[i][j] is the probability that it died in the
// chain. Rows i >= last are zero. P = A(dt) = exp(G dt) restricted to these
// rows, for the generator G of advancement and death (exit absorbing). With
// mu = 0, rows sum to 1.
//
// advance_average() also fills Pbar = (1/dt) int_0^dt A(u) du: the same for
// an individual that entered the chain at a time uniformly distributed within
// the step (D-022: arrivals are placed at their expected position, having
// faced this chain's hazard only since their arrival).
//
// Both are computed by uniformization: with Lambda = max q + mu and
// S = I + G / Lambda, A(h) = sum_k Pois(k; Lambda h) S^k and
// Abar(h) = sum_k P(N > k) / (Lambda h) S^k with N ~ Pois(Lambda h). All
// terms are non-negative, so there is no cancellation. If Lambda dt > 300 the
// step is halved until it is not, and the results are doubled back with
// A(2h) = A(h)^2 and Abar(2h) = (Abar(h) + A(h) Abar(h)) / 2. A row whose sum
// exceeds 1 by round-off is divided by its sum.
//
// Work arrays are kept between calls; repeated calls with the same chain
// length allocate nothing.

#ifndef SGDRIVE_CHAINS_H
#define SGDRIVE_CHAINS_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgdrive {

class ChainAdvance {
 public:
  explicit ChainAdvance(int n)
      : n_(n),
        v_(n + 1),
        vn_(n + 1),
        a_((n + 1) * (n + 1)),
        abar_((n + 1) * (n + 1)),
        tmp_((n + 1) * (n + 1)) {}

  int size() const { return n_; }

  void advance(const double* q, int last, double mu, double dt, double* P) {
    compute(q, last, mu, dt, P, nullptr);
  }

  void advance_average(const double* q, int last, double mu, double dt, double* P,
                       double* Pbar) {
    compute(q, last, mu, dt, P, Pbar);
  }

 private:
  static constexpr double kMaxLambdaDt = 300.0;
  int n_;
  std::vector<double> v_, vn_, a_, abar_, tmp_, w_, c_;

  // Next class after j (the exit state n after the last reached class).
  int next(int j, int last) const { return j + 1 < last ? j + 1 : n_; }

  void compute(const double* q, int last, double mu, double dt, double* P, double* Pbar) {
    const int n = n_, m = n + 1;
    const bool avg = Pbar != nullptr;
    std::fill(P, P + static_cast<std::size_t>(n) * m, 0.0);
    if (avg) std::fill(Pbar, Pbar + static_cast<std::size_t>(n) * m, 0.0);
    if (last <= 0) return;
    double lambda = 0.0;
    for (int i = 0; i < last; ++i) lambda = std::max(lambda, q[i]);
    lambda += mu;
    if (lambda <= 0.0 || dt <= 0.0) {
      for (int i = 0; i < last; ++i) {
        P[i * m + i] = 1.0;
        if (avg) Pbar[i * m + i] = 1.0;
      }
      return;
    }
    int halvings = 0;
    double h = dt;
    while (lambda * h > kMaxLambdaDt) {
      h *= 0.5;
      ++halvings;
    }
    uniformize(q, last, mu, lambda, h, avg);
    for (int s = 0; s < halvings; ++s) {
      if (avg) {  // Abar(2h) = (Abar(h) + A(h) Abar(h)) / 2, with A(h) before squaring
        multiply(a_.data(), abar_.data(), tmp_.data());
        for (int k = 0; k < m * m; ++k) abar_[k] = 0.5 * (abar_[k] + tmp_[k]);
      }
      multiply(a_.data(), a_.data(), tmp_.data());
      std::swap(a_, tmp_);
    }
    copy_rows(a_, last, P);
    if (avg) copy_rows(abar_, last, Pbar);
  }

  // Rows 0..last-1 of a full (n + 1) x (n + 1) matrix into P (a row sum
  // above 1 from round-off is scaled back to 1).
  void copy_rows(const std::vector<double>& full, int last, double* P) const {
    const int m = n_ + 1;
    for (int i = 0; i < last; ++i) {
      double sum = 0.0;
      for (int j = i; j < m; ++j) sum += full[i * m + j];
      const double scale = sum > 1.0 ? 1.0 / sum : 1.0;
      for (int j = i; j < m; ++j) P[i * m + j] = full[i * m + j] * scale;
    }
  }

  // a_ = A(h) and, if avg, abar_ = Abar(h), on all n + 1 states (the exit
  // state is absorbing; rows of classes >= last stay zero; mass lost to death
  // is not tracked).
  void uniformize(const double* q, int last, double mu, double lambda, double h, bool avg) {
    const int m = n_ + 1;
    std::fill(a_.begin(), a_.end(), 0.0);
    a_[n_ * m + n_] = 1.0;
    if (avg) {
      std::fill(abar_.begin(), abar_.end(), 0.0);
      abar_[n_ * m + n_] = 1.0;
    }
    // Poisson weights w_k = Pois(k; a), enough for a tail below double
    // precision, and c_k = P(N > k) / a (they sum to 1)
    const double a = lambda * h;
    const int k_max = static_cast<int>(std::ceil(a + 12.0 * std::sqrt(a) + 30.0));
    w_.resize(k_max + 1);
    c_.resize(k_max + 1);
    w_[0] = std::exp(-a);
    int k_end = k_max;
    for (int k = 1; k <= k_max; ++k) {
      w_[k] = w_[k - 1] * a / k;
      if (k > a && w_[k] < 1e-20) {
        k_end = k;
        break;
      }
    }
    double tail = 0.0;
    for (int k = k_end; k >= 0; --k) {
      c_[k] = tail / a;  // P(N > k) / a
      tail += w_[k];
    }
    for (int i = 0; i < last; ++i) {
      std::fill(v_.begin(), v_.end(), 0.0);
      v_[i] = 1.0;
      for (int k = 0; k <= k_end; ++k) {
        if (k > 0) {
          // v <- v S: advance with probability q / lambda, die with
          // probability mu / lambda, otherwise stay
          std::fill(vn_.begin(), vn_.end(), 0.0);
          vn_[n_] = v_[n_];
          for (int j = i; j < last; ++j) {
            if (v_[j] == 0.0) continue;
            const double p = q[j] / lambda;
            vn_[j] += v_[j] * (1.0 - p - mu / lambda);
            vn_[next(j, last)] += v_[j] * p;
          }
          std::swap(v_, vn_);
        }
        for (int j = i; j < m; ++j) {
          a_[i * m + j] += w_[k] * v_[j];
          if (avg) abar_[i * m + j] += c_[k] * v_[j];
        }
      }
    }
  }

  // C = A B for upper-triangular (n + 1) x (n + 1) matrices.
  void multiply(const double* A, const double* B, double* C) const {
    const int m = n_ + 1;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        double s = 0.0;
        for (int k = i; k <= j; ++k) s += A[i * m + k] * B[k * m + j];
        C[i * m + j] = s;
      }
    }
  }
};

}  // namespace sgdrive

#endif  // SGDRIVE_CHAINS_H
