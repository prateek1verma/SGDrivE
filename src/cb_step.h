// Chain-binomial step of the SGDrivE v0.1 model (D-022, D-024; model
// specification S15-S19).
//
// A block-by-block translation of the R reference step .sg_cb_step_r()
// (R/cb-reference.R). The random draws are made with R's own generators
// (R::rbinom, R::rpois, the functions behind stats::rbinom and stats::rpois)
// in the canonical order of the specification, so with the same seed the two
// give identical results; with ExpectedDraws every draw is replaced by its
// expectation (the deterministic map). The R version is the readable
// reference and the test oracle.
//
// prepare_day() computes the chains of the day (death and destination
// probabilities, arrival distributions; S16) once per day; step() advances
// a state by one step of length dt. Working arrays are allocated in the
// constructor; neither function allocates.

#ifndef SGDRIVE_CB_STEP_H
#define SGDRIVE_CB_STEP_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "chains.h"
#include "model_view.h"
#include "rates.h"

namespace sgdrive {

// Random draws from R's generator (call within an Rcpp::RNGScope).
struct RandomDraws {
  double binom(double n, double p) { return (n == 0.0 || p == 0.0) ? 0.0 : R::rbinom(n, p); }
  double pois(double lambda) { return lambda <= 0.0 ? 0.0 : R::rpois(lambda); }
  // Total of n zero-truncated Poisson(lambda) numbers (eggs of n batches with
  // at least one egg; D-026): one Poisson draw above lambda = 30, where the
  // truncation is below 1e-13, otherwise each batch by inversion (as in R).
  double ztp_sum(double n, double lambda) {
    if (lambda > 30.0) return pois(n * lambda);
    const double k_max = std::ceil(lambda + 12.0 * std::sqrt(lambda) + 30.0);
    double total = 0.0;
    for (double i = 0.0; i < n; i += 1.0) {
      const double u = R::runif(0.0, 1.0) * -std::expm1(-lambda);
      double p = std::exp(-lambda) * lambda;
      double cum = p;
      double k = 1.0;
      while (u > cum && k < k_max) {
        k += 1.0;
        p = p * lambda / k;
        cum = cum + p;
      }
      total += k;
    }
    return total;
  }
  // Multinomial as sequential conditional binomials in category order; the
  // remaining probabilities are summed from the last category (as in R).
  void multinom(double n, const double* prob, int k, double* out, double* rest) {
    std::fill(out, out + k, 0.0);
    rest[k - 1] = prob[k - 1];
    for (int j = k - 2; j >= 0; --j) rest[j] = prob[j] + rest[j + 1];
    double left = n;
    for (int j = 0; j < k - 1; ++j) {
      if (left == 0.0) break;
      const double pj = rest[j] > 0.0 ? std::min(1.0, prob[j] / rest[j]) : 0.0;
      out[j] = binom(left, pj);
      left -= out[j];
    }
    out[k - 1] = left;
  }
};

// Expectations instead of draws (the deterministic map).
struct ExpectedDraws {
  double binom(double n, double p) { return n * p; }
  double pois(double lambda) { return lambda; }
  double ztp_sum(double n, double lambda) { return n * lambda / -std::expm1(-lambda); }
  void multinom(double n, const double* prob, int k, double* out, double*) {
    for (int j = 0; j < k; ++j) out[j] = n * prob[j];
  }
};

// One chain prepared for steps of length dt (see .sg_cb_chain() in R).
struct PreparedChain {
  int n = 0, last = 0;
  std::vector<double> die;   // [n]: probability of dying in the chain
  std::vector<double> dest;  // [n x (n + 1)], row-major: survivors' destinations
  std::vector<double> arr;   // [n]: placement of arrivals that survive
  double arr_die = 0.0;      // probability that an arrival dies in the step
  void resize(int size) {
    n = size;
    die.assign(n, 0.0);
    dest.assign(static_cast<std::size_t>(n) * (n + 1), 0.0);
    arr.assign(n, 0.0);
  }
};

class CbStepper {
 public:
  CbStepper(const ModelView& mv, double dt);

  // Forcing-dependent chains of simulation day `day` (1-based).
  void prepare_day(int day);

  // Advance state x (length n_states) by one step in place; flows (length
  // n_counters, order of .sg_counter_names()) receive the step's flows.
  template <class D>
  void step(double* x, double* flows, D& draws);

 private:
  const ModelView& mv_;
  double dt_;
  RateKernel kernel_;
  DayRates r_;
  int day_ = 0;
  ChainAdvance adv_e_, adv_l_, adv_p_, adv_a_;
  PreparedChain ch_e_, ch_l_, ch_p_;
  std::vector<PreparedChain> ch_m_, ch_f_;
  // work arrays
  std::vector<double> q_, p_, pbar_;           // chain preparation
  std::vector<double> cur_, nx_;               // state at step start (after mating), new state
  std::vector<double> eggs_, out_e_, out_l_, out_p_, em_m_, em_f_, matings_, ym_;
  std::vector<double> buf_, rest_;             // multinomial output and work

  void prepare(PreparedChain& ch, ChainAdvance& adv, const double* q, int n, int last,
               double mu);
  template <class D>
  void move(const double* x, double* nx, const PreparedChain& ch, D& draws, double* out,
            double* deaths);
  template <class D>
  void place(double n, const PreparedChain& ch, double* nx, D& draws, double* deaths);
};

// ---- template members --------------------------------------------------------

// Chain moves of one chain (S16): death in the chain, then destination.
template <class D>
void CbStepper::move(const double* x, double* nx, const PreparedChain& ch, D& draws,
                     double* out, double* deaths) {
  const int n = ch.n;
  for (int i = 0; i < n; ++i) {
    if (x[i] == 0.0) continue;
    const double d = draws.binom(x[i], ch.die[i]);
    *deaths += d;
    const int k = n + 1 - i;  // categories: sub-stages i..n-1, then left the chain
    draws.multinom(x[i] - d, ch.dest.data() + static_cast<std::size_t>(i) * (n + 1) + i, k,
                   buf_.data(), rest_.data());
    for (int j = 0; j < k - 1; ++j) nx[i + j] += buf_[j];
    *out += buf_[k - 1];
  }
}

// Arrivals in one chain (S16): death for the rest of the step, then placement.
template <class D>
void CbStepper::place(double n, const PreparedChain& ch, double* nx, D& draws,
                      double* deaths) {
  if (n == 0.0) return;
  const double d = draws.binom(n, ch.arr_die);
  *deaths += d;
  draws.multinom(n - d, ch.arr.data(), ch.n, buf_.data(), rest_.data());
  for (int j = 0; j < ch.n; ++j) nx[j] += buf_[j];
}

template <class D>
void CbStepper::step(double* x, double* flows, D& draws) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const double dt = dt_;

  // ---- aggregates (S15): from the state at the start of the step ----------
  kernel_.aggregate(r_, x);
  const Aggregates& a = kernel_.agg();
  std::copy(x, x + v.n_states, cur_.begin());
  std::fill(nx_.begin(), nx_.end(), 0.0);
  double* egg = cur_.data() + v.offsets[kE];
  double* lar = cur_.data() + v.offsets[kL];
  double* pup = cur_.data() + v.offsets[kP];
  double* mal = cur_.data() + v.offsets[kM];
  double* unm = cur_.data() + v.offsets[kU];
  double* fem = cur_.data() + v.offsets[kF];
  double* n_e = nx_.data() + v.offsets[kE];
  double* n_l = nx_.data() + v.offsets[kL];
  double* n_p = nx_.data() + v.offsets[kP];
  double* n_m = nx_.data() + v.offsets[kM];
  double* n_u = nx_.data() + v.offsets[kU];
  double* n_f = nx_.data() + v.offsets[kF];

  // ---- eggs (S17): oviposition events, those with at least one viable egg
  // (fertile, the case proxy; D-026) and their viable eggs, or single eggs --
  const bool batches = v.oviposition_batches;
  const double factor = r_.batch_factor, clutch = r_.clutch;
  std::fill(eggs_.begin(), eggs_.end(), 0.0);
  double n_batches = 0.0, n_fertile = 0.0, expected_fertile = 0.0;
  const int n_cross = ng * ng;
  for (int g = 0; g < ng; ++g) {
    for (int m = 0; m < ng; ++m) {
      const int pair = g * ng + m;
      const double tot = kernel_.pair_viability()[pair];
      const double share = -std::expm1(-clutch * tot);
      double n_eggs;
      if (batches) {
        const double nb = draws.pois(a.psi * factor * a.w_lay[pair] * dt);
        n_batches += nb;
        if (nb == 0.0 || share <= 0.0) continue;
        const double nf = draws.binom(nb, share);
        n_fertile += nf;
        if (nf == 0.0) continue;
        n_eggs = draws.ztp_sum(nf, clutch * tot);
      } else {
        expected_fertile += a.w_lay[pair] * share;
        if (tot <= 0.0) continue;
        n_eggs = draws.pois(v.R * r_.phi_ov * a.psi * a.w_lay[pair] * tot * dt);
      }
      for (int k = 0; k < ng; ++k) {
        ym_[k] = v.egg_tensor[static_cast<std::size_t>(k) * n_cross + pair] / tot;
      }
      draws.multinom(n_eggs, ym_.data(), ng, buf_.data(), rest_.data());
      for (int k = 0; k < ng; ++k) eggs_[k] += buf_[k];
    }
  }

  // ---- mating (S18), before ageing ------------------------------------------
  std::fill(matings_.begin(), matings_.end(), 0.0);
  const double nu_g = v.nu * a.g_find;
  for (int g = 0; g < ng; ++g) {
    if (!a.mate_ok[g]) continue;
    const double* pi_g = a.pi.data() + static_cast<std::size_t>(g) * ng;
    for (int s = 0; s < nA; ++s) {
      const double nu_s = v.mating_instant ? nu_g : v.nu * v.receptive[s] * a.g_find;
      const double y = draws.binom(unm[g * nA + s], 1.0 - std::exp(-nu_s * dt));
      if (y == 0.0) continue;
      draws.multinom(y, pi_g, ng, buf_.data(), rest_.data());
      unm[g * nA + s] -= y;
      for (int m = 0; m < ng; ++m) {
        fem[(g * ng + m) * nA + s] += buf_[m];
        matings_[g * ng + m] += buf_[m];
      }
    }
  }

  // ---- chain moves (S16), block by block in state order ------------------------
  if (v.larval_per_wound && v.H > 0) {  // the larval hazard depends on crowding
    for (int s = 0; s < nL; ++s) q_[s] = r_.q_l;
    prepare(ch_l_, adv_l_, q_.data(), nL, nL, a.mu_l);
  }
  double dth[6] = {0, 0, 0, 0, 0, 0};  // deaths E, L, P, M, U, F
  double sen_m = 0.0, sen_u = 0.0, sen_f = 0.0;
  for (int g = 0; g < ng; ++g) {
    out_e_[g] = 0.0;
    move(egg + g * nE, n_e + g * nE, ch_e_, draws, &out_e_[g], &dth[0]);
  }
  for (int g = 0; g < ng; ++g) {
    out_l_[g] = 0.0;
    move(lar + g * nL, n_l + g * nL, ch_l_, draws, &out_l_[g], &dth[1]);
  }
  for (int g = 0; g < ng; ++g) {
    out_p_[g] = 0.0;
    move(pup + g * nP, n_p + g * nP, ch_p_, draws, &out_p_[g], &dth[2]);
  }
  for (int g = 0; g < ng; ++g) move(mal + g * nA, n_m + g * nA, ch_m_[g], draws, &sen_m, &dth[3]);
  for (int g = 0; g < ng; ++g) move(unm + g * nA, n_u + g * nA, ch_f_[g], draws, &sen_u, &dth[4]);
  for (int g = 0; g < ng; ++g) {
    for (int m = 0; m < ng; ++m) {
      const int off = (g * ng + m) * nA;
      move(fem + off, n_f + off, ch_f_[g], draws, &sen_f, &dth[5]);
    }
  }

  // ---- emergence (S19): sex, then survival ------------------------------------
  double sum_em_m = 0.0, sum_em_f = 0.0, sum_z = 0.0;
  for (int g = 0; g < ng; ++g) {
    const double z = out_p_[g];
    const double zf = draws.binom(z, v.phi[g]);
    em_f_[g] = draws.binom(zf, v.xiF[g]);
    em_m_[g] = draws.binom(z - zf, v.xiM[g]);
    sum_z += z;
    sum_em_m += em_m_[g];
    sum_em_f += em_f_[g];
  }

  // ---- arrivals (S16), in state order -----------------------------------------
  for (int k = 0; k < ng; ++k) place(eggs_[k], ch_e_, n_e + k * nE, draws, &dth[0]);
  for (int g = 0; g < ng; ++g) place(out_e_[g], ch_l_, n_l + g * nL, draws, &dth[1]);
  for (int g = 0; g < ng; ++g) place(out_l_[g], ch_p_, n_p + g * nP, draws, &dth[2]);
  for (int g = 0; g < ng; ++g) place(em_m_[g], ch_m_[g], n_m + g * nA, draws, &dth[3]);
  // emerged females (S18, S19): death, then mating (at emergence with
  // "instant", otherwise for the rest of the step), then placement
  const double x_arr = (v.mating_instant ? nu_g : v.nu * v.receptive[0] * a.g_find) * dt;
  const double p_arr_mate = x_arr > 0.0 ? 1.0 + std::expm1(-x_arr) / x_arr : 0.0;
  for (int g = 0; g < ng; ++g) {
    const PreparedChain& ch = ch_f_[g];
    std::fill(ym_.begin(), ym_.end(), 0.0);
    double unmated = 0.0;
    if (em_f_[g] > 0.0) {
      const bool mated_now = v.mating_instant && a.mate_ok[g];
      const double d = draws.binom(em_f_[g], ch.arr_die);
      dth[mated_now ? 5 : 4] += d;
      const double surv = em_f_[g] - d;
      double y = 0.0;
      if (mated_now) {
        y = surv;
      } else if (!v.mating_instant && a.mate_ok[g]) {
        y = draws.binom(surv, p_arr_mate);
      }
      if (y > 0.0) {
        draws.multinom(y, a.pi.data() + static_cast<std::size_t>(g) * ng, ng, ym_.data(),
                       rest_.data());
      }
      unmated = surv - y;
      for (int m = 0; m < ng; ++m) matings_[g * ng + m] += ym_[m];
    }
    if (unmated > 0.0) {
      draws.multinom(unmated, ch.arr.data(), nA, buf_.data(), rest_.data());
      for (int s = 0; s < nA; ++s) n_u[g * nA + s] += buf_[s];
    }
    for (int m = 0; m < ng; ++m) {
      if (ym_[m] > 0.0) {
        draws.multinom(ym_[m], ch.arr.data(), nA, buf_.data(), rest_.data());
        for (int s = 0; s < nA; ++s) n_f[(g * ng + m) * nA + s] += buf_[s];
      }
    }
  }

  // ---- flows (order of .sg_counter_names()) -----------------------------------
  int c = 0;
  for (int k = 0; k < ng; ++k) flows[c++] = eggs_[k];
  flows[c++] = batches ? n_batches : a.psi * a.sum_w * factor * dt;  // egg_batches
  flows[c++] = batches ? n_fertile : a.psi * factor * expected_fertile * dt;  // fertile_batches
  double sum_out_e = 0.0, sum_out_l = 0.0;
  for (int g = 0; g < ng; ++g) {
    sum_out_e += out_e_[g];
    sum_out_l += out_l_[g];
  }
  flows[c++] = sum_out_e;  // larvae_in
  flows[c++] = sum_out_l;  // pupae_in
  flows[c++] = sum_em_m;
  flows[c++] = sum_em_f;
  flows[c++] = sum_z - sum_em_m - sum_em_f;  // emergence_deaths
  for (int i = 0; i < ng * ng; ++i) flows[c++] = matings_[i];
  for (int b = 0; b < 6; ++b) flows[c++] = dth[b];
  flows[c++] = sen_m;
  flows[c++] = sen_u;
  flows[c++] = sen_f;

  std::copy(nx_.begin(), nx_.end(), x);
}

}  // namespace sgdrive

#endif  // SGDRIVE_CB_STEP_H
