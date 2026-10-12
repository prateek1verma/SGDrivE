// Exact advancement along a chain of sub-stages over one time step (D-022,
// model specification S16), for the chain-binomial solver.
//
// A chain has classes 0..n-1. Individuals in class i < last move to class
// i + 1 at rate q[i]; from class last - 1 they leave the chain ("exit": to
// the next life stage, or death by old age for adults). Classes >= last are
// never reached (D-018). Death is not included: the hazard is the same in
// every class of a chain, so survival over the step, exp(-mu dt), is
// independent of how far an individual advances.
//
// advance() fills P, row-major with n rows and n + 1 columns: P[i][j] is the
// probability that an individual in class i at the start of the step is in
// class j (i <= j < last) at its end, and P[i][n] that it has left the
// chain. Rows i >= last are zero. P = exp(G dt) restricted to these rows, for
// the pure-advancement generator G, computed by uniformization: with
// Lambda = max q, P = sum_k Pois(k; Lambda dt) S^k with S = I + G / Lambda.
// All terms are non-negative, so there is no cancellation. If Lambda dt > 300
// the step is halved until it is not, and the result squared back. Each row
// is finally divided by its sum (1 up to about 1e-15) so that it is an exact
// probability vector.
//
// Work arrays are allocated once per chain length; advance() allocates
// nothing.

#ifndef SGDRIVE_CHAINS_H
#define SGDRIVE_CHAINS_H

#include <algorithm>
#include <cmath>
#include <vector>

namespace sgdrive {

class ChainAdvance {
 public:
  explicit ChainAdvance(int n)
      : n_(n), v_(n + 1), vn_(n + 1), a_((n + 1) * (n + 1)), b_((n + 1) * (n + 1)) {}

  int size() const { return n_; }

  void advance(const double* q, int last, double dt, double* P) {
    const int n = n_, m = n + 1;
    std::fill(P, P + static_cast<std::size_t>(n) * m, 0.0);
    double lambda = 0.0;
    for (int i = 0; i < last; ++i) lambda = std::max(lambda, q[i]);
    if (last <= 0) return;
    if (lambda <= 0.0 || dt <= 0.0) {
      for (int i = 0; i < last; ++i) P[i * m + i] = 1.0;
      return;
    }
    int halvings = 0;
    double h = dt;
    while (lambda * h > kMaxLambdaDt) {
      h *= 0.5;
      ++halvings;
    }
    // full (n + 1) x (n + 1) matrix including the absorbing exit state
    uniformize(q, last, lambda, h, a_.data());
    for (int s = 0; s < halvings; ++s) {
      square(a_.data(), b_.data());
      std::swap(a_, b_);
    }
    // rows sum to 1 up to round-off accumulated over the terms; normalize so
    // that they are exact probability vectors for multinomial draws
    for (int i = 0; i < last; ++i) {
      double sum = 0.0;
      for (int j = i; j < m; ++j) sum += a_[i * m + j];
      for (int j = i; j < m; ++j) P[i * m + j] = a_[i * m + j] / sum;
    }
  }

 private:
  static constexpr double kMaxLambdaDt = 300.0;
  int n_;
  std::vector<double> v_, vn_, a_, b_;

  // Next class after j (the exit state n after the last reached class).
  int next(int j, int last) const { return j + 1 < last ? j + 1 : n_; }

  // A = exp(G h) on all n + 1 states (rows >= last other than exit stay 0).
  void uniformize(const double* q, int last, double lambda, double h, double* A) {
    const int m = n_ + 1;
    std::fill(A, A + static_cast<std::size_t>(m) * m, 0.0);
    A[n_ * m + n_] = 1.0;  // exit is absorbing
    const double a = lambda * h;
    // enough Poisson terms for a tail below double precision
    const int k_max = static_cast<int>(std::ceil(a + 12.0 * std::sqrt(a) + 30.0));
    for (int i = 0; i < last; ++i) {
      std::fill(v_.begin(), v_.end(), 0.0);
      v_[i] = 1.0;
      double w = std::exp(-a);  // Pois(0; a)
      for (int j = 0; j < m; ++j) A[i * m + j] = w * v_[j];
      for (int k = 1; k <= k_max; ++k) {
        // v <- S v (row vector times S): advance with probability q / lambda
        std::fill(vn_.begin(), vn_.end(), 0.0);
        vn_[n_] = v_[n_];
        for (int j = i; j < last; ++j) {
          if (v_[j] == 0.0) continue;
          const double p = q[j] / lambda;
          vn_[j] += v_[j] * (1.0 - p);
          vn_[next(j, last)] += v_[j] * p;
        }
        std::swap(v_, vn_);
        w *= a / k;
        for (int j = i; j < m; ++j) A[i * m + j] += w * v_[j];
        if (k > a && w < 1e-20) break;
      }
    }
  }

  // B = A A for upper-triangular (n + 1) x (n + 1) matrices.
  void square(const double* A, double* B) const {
    const int m = n_ + 1;
    for (int i = 0; i < m; ++i) {
      for (int j = 0; j < m; ++j) {
        double s = 0.0;
        for (int k = i; k <= j; ++k) s += A[i * m + k] * A[k * m + j];
        B[i * m + j] = s;
      }
    }
  }
};

}  // namespace sgdrive

#endif  // SGDRIVE_CHAINS_H
