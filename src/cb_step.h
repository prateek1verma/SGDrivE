// Chain-binomial step of the ScrewDrivE model (D-022, D-024; model
// specification S15-S19; several nodes: section 10 and N4, D-032).
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

#ifndef SCREWDRIVE_CB_STEP_H
#define SCREWDRIVE_CB_STEP_H

#include <Rcpp.h>

#include <algorithm>
#include <cmath>
#include <vector>

#include "chains.h"
#include "model_view.h"
#include "rates.h"

namespace screwdrive {

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
  // Distinct hosts hit when n batches fall at random on w wounded hosts
  // (D-030): batch by batch, a new host with probability (w - k) / w; stops
  // once every host is hit (as in R).
  double occupancy(double n, double w) {
    double k = 0.0;
    for (double i = 0.0; i < n; i += 1.0) {
      if (k >= w) break;
      if (R::runif(0.0, 1.0) < (w - k) / w) k += 1.0;
    }
    return k;
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
  double occupancy(double n, double w) { return w > 0.0 ? w * -std::expm1(-n / w) : 0.0; }
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

// The day's prepared chains of one forcing series (shared by the nodes on
// that series, D-032): day rates, eggs, larvae without crowding, pupae, and
// the adult chains of each genotype.
struct SeriesChains {
  DayRates r;
  PreparedChain e, l, p;
  std::vector<PreparedChain> m, f;  // [n_g]
};

// What part 1 of a node's step hands to part 3 (MODEL_SPEC section 10).
struct NodeWork {
  Aggregates agg;   // aggregates at the start of the step (S15)
  PreparedChain l;  // larval chain with crowding, prepared this step
  bool crowd = false, empty = false;
  std::vector<double> eggs, out_e, out_l, out_p, em_m, em_f, matings, host_i;
  double dth[6] = {0, 0, 0, 0, 0, 0};  // deaths E, L, P, M, U, F
  double sen_m = 0, sen_u = 0, sen_f = 0;
  double n_batches = 0, n_fertile = 0, expected_fertile = 0;
  double h_s = 0, h_w = 0, infested = 0, healed = 0, ended = 0, treated = 0, wounded = 0;
};

class CbStepper {
 public:
  CbStepper(const ModelView& mv, double dt);

  // Forcing-dependent chains of simulation day `day` (1-based), once per
  // distinct forcing series.
  void prepare_day(int day);

  // Advance state x (length n_states) by one step in place; flows (length
  // n_counters, order of .sg_model_counter_names()) receive the step's flows.
  // With several nodes: part 1 of every node, movement, part 3 of every node
  // (MODEL_SPEC section 10).
  template <class D>
  void step(double* x, double* flows, D& draws);

 private:
  const ModelView& mv_;
  double dt_;
  RateKernel kernel_;
  int day_ = 0;
  ChainAdvance adv_e_, adv_l_, adv_p_, adv_a_, adv_i_;
  std::vector<SeriesChains> series_;  // [n_series]
  PreparedChain ch_i_;                // infestation chain (D-030), the same everywhere
  std::vector<NodeWork> work_;        // [n_node]
  // work arrays
  std::vector<double> q_, p_, pbar_;  // chain preparation
  std::vector<double> cur_, nx_;      // state at step start (after mating), new state
  std::vector<double> ym_;
  std::vector<double> buf_, rest_;    // multinomial output and work
  // movement (N4): destinations of each node, the union of the groups' edges
  // in increasing node order (CSR), and the union position of every group
  // edge; per day the probability of leaving [3, n_node], the destination
  // probabilities [3, union edges] and the probability of being lost
  // [3, n_node]; immigrants of the step; one row's probabilities
  std::vector<int> move_row_ptr_, move_col_, move_map_;
  std::vector<double> move_pe_, move_prob_, move_lost_, imm_;
  std::vector<double> move_p_, move_buf_, move_rest_;

  void prepare(PreparedChain& ch, ChainAdvance& adv, const double* q, int n, int last,
               double mu);
  template <class D>
  void move(const double* x, double* nx, const PreparedChain& ch, D& draws, double* out,
            double* deaths);
  template <class D>
  void place(double n, const PreparedChain& ch, double* nx, D& draws, double* deaths);
  template <class D>
  void part1(int node, D& draws);
  template <class D>
  void part3(int node, double* flows, D& draws);
  void prepare_movement(int day);
  template <class D>
  void movement(double* move_counters, D& draws);
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
  std::copy(x, x + v.n_states, cur_.begin());
  std::fill(nx_.begin(), nx_.end(), 0.0);
  for (int i = 0; i < v.n_node; ++i) part1(i, draws);
  if (v.n_groups > 0) movement(flows + static_cast<std::size_t>(v.n_node) * v.n_counters_node, draws);
  for (int i = 0; i < v.n_node; ++i) {
    part3(i, flows + static_cast<std::size_t>(i) * v.n_counters_node, draws);
  }
  std::copy(nx_.begin(), nx_.end(), x);
}

// (N4) Movement of all nodes, between part 1 and part 3: the adults after
// the chain moves (in nx_) leave with the day's probability, compartment by
// compartment in state order (node slowest; M, U, F; genotype, mate genotype,
// sub-stage); emigrants are split over the destinations in increasing node
// order and "lost" (sequential conditional binomials); immigrants are added
// after all draws, so no fly moves twice in a step. Counters: emigrants,
// immigrants and lost by class for each node (9 per node).
template <class D>
void CbStepper::movement(double* mc, D& draws) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nA = v.n_A;
  const int block_off[3] = {v.offsets[kM], v.offsets[kU], v.offsets[kF]};
  const int block_len[3] = {ng * nA, ng * nA, ng * ng * nA};
  std::fill(mc, mc + ModelView::kMoveCounters * v.n_node, 0.0);
  std::fill(imm_.begin(), imm_.end(), 0.0);
  for (int i = 0; i < v.n_node; ++i) {
    const std::size_t base = static_cast<std::size_t>(i) * v.node_stride;
    const int e0 = move_row_ptr_[i], deg = move_row_ptr_[i + 1] - move_row_ptr_[i];
    for (int s = 0; s < 3; ++s) {
      const double pe = move_pe_[static_cast<std::size_t>(i) * 3 + s];
      if (pe == 0.0) continue;
      for (int d = 0; d < deg; ++d) move_p_[d] = move_prob_[static_cast<std::size_t>(e0 + d) * 3 + s];
      move_p_[deg] = move_lost_[static_cast<std::size_t>(i) * 3 + s];
      double* x = nx_.data() + base + block_off[s];
      for (int k = 0; k < block_len[s]; ++k) {
        if (x[k] == 0.0) continue;
        const double e = draws.binom(x[k], pe);
        if (e == 0.0) continue;
        x[k] -= e;
        mc[i * ModelView::kMoveCounters + s] += e;  // emigrants
        draws.multinom(e, move_p_.data(), deg + 1, move_buf_.data(), move_rest_.data());
        for (int d = 0; d < deg; ++d) {
          const int j = move_col_[e0 + d];
          imm_[static_cast<std::size_t>(j) * v.node_stride + block_off[s] + k] += move_buf_[d];
          mc[j * ModelView::kMoveCounters + 3 + s] += move_buf_[d];  // immigrants
        }
        mc[i * ModelView::kMoveCounters + 6 + s] += move_buf_[deg];  // lost
      }
    }
  }
  for (std::size_t q = 0; q < nx_.size(); ++q) nx_[q] += imm_[q];
}

// Part 1 of node i's step: aggregates (S15), eggs (S17), hosts (H3), mating
// (S18) and chain moves (S16), from the node's counts at the start of the
// step. A node without flies skips its fly draws (draws of zero individuals
// use no random numbers, so the result is the same).
template <class D>
void CbStepper::part1(int node, D& draws) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const double dt = dt_;
  const std::size_t base = static_cast<std::size_t>(node) * v.node_stride;
  const SeriesChains& sc = series_[v.forcing_index[node]];
  const DayRates& r = sc.r;
  NodeWork& w = work_[node];

  // ---- aggregates (S15): from the state at the start of the step ----------
  kernel_.set_node(node);
  kernel_.aggregate(r, cur_.data() + base);
  w.agg = kernel_.agg();
  const Aggregates& a = w.agg;
  double* egg = cur_.data() + base + v.offsets[kE];
  double* lar = cur_.data() + base + v.offsets[kL];
  double* pup = cur_.data() + base + v.offsets[kP];
  double* mal = cur_.data() + base + v.offsets[kM];
  double* unm = cur_.data() + base + v.offsets[kU];
  double* fem = cur_.data() + base + v.offsets[kF];
  double* n_e = nx_.data() + base + v.offsets[kE];
  double* n_l = nx_.data() + base + v.offsets[kL];
  double* n_p = nx_.data() + base + v.offsets[kP];
  double* n_m = nx_.data() + base + v.offsets[kM];
  double* n_u = nx_.data() + base + v.offsets[kU];
  double* n_f = nx_.data() + base + v.offsets[kF];
  const int n_fly = v.offsets[kF] + ng * ng * nA;  // fly states come first in a node
  w.empty = std::all_of(egg, egg + n_fly, [](double z) { return z == 0.0; });

  std::fill(w.eggs.begin(), w.eggs.end(), 0.0);
  std::fill(w.out_e.begin(), w.out_e.end(), 0.0);
  std::fill(w.out_l.begin(), w.out_l.end(), 0.0);
  std::fill(w.out_p.begin(), w.out_p.end(), 0.0);
  std::fill(w.matings.begin(), w.matings.end(), 0.0);
  std::fill(w.dth, w.dth + 6, 0.0);
  w.sen_m = w.sen_u = w.sen_f = 0.0;
  w.n_batches = w.n_fertile = w.expected_fertile = 0.0;

  // ---- eggs (S17): oviposition events, those with at least one viable egg
  // (fertile, the case proxy; D-026) and their viable eggs, or single eggs --
  const bool batches = v.oviposition_batches;
  const double factor = r.batch_factor, clutch = r.clutch;
  const int n_cross = ng * ng;
  if (!w.empty) {
    for (int g = 0; g < ng; ++g) {
      for (int m = 0; m < ng; ++m) {
        const int pair = g * ng + m;
        const double tot = kernel_.pair_egg_total()[pair];
        // females lay batches (and search) in proportion to s_g; viable eggs
        // per egg of a batch without s_g (D-031)
        const double w_search = kernel_.pair_search()[pair] * a.w_lay[pair];
        const double viab = kernel_.pair_batch_viability()[pair];
        const double share = -std::expm1(-clutch * viab);
        double n_eggs;
        if (batches) {
          const double nb = draws.pois(a.psi * factor * w_search * dt);
          w.n_batches += nb;
          if (nb == 0.0 || share <= 0.0) continue;
          const double nf = draws.binom(nb, share);
          w.n_fertile += nf;
          if (nf == 0.0) continue;
          n_eggs = draws.ztp_sum(nf, clutch * viab);
        } else {
          w.expected_fertile += w_search * share;
          if (tot <= 0.0) continue;
          n_eggs = draws.pois(v.R * r.phi_ov * a.psi * a.w_lay[pair] * tot * dt);
        }
        for (int k = 0; k < ng; ++k) {
          ym_[k] = v.egg_tensor[static_cast<std::size_t>(k) * n_cross + pair] / tot;
        }
        draws.multinom(n_eggs, ym_.data(), ng, buf_.data(), rest_.data());
        for (int k = 0; k < ng; ++k) w.eggs[k] += buf_[k];
      }
    }
  }

  // ---- hosts (H3, D-030): fertile batches on uninfested wounds infest the
  // distinct hosts they hit (cases); the rest are superinfestations; then
  // healing, the infestation chain (end or treatment) and wounding -------------
  w.h_s = w.h_w = w.infested = w.healed = w.ended = w.treated = w.wounded = 0.0;
  if (v.host_model) {
    if (!batches) {  // single eggs: fertile oviposition events drawn for the hosts
      w.n_fertile = draws.pois(a.psi * factor * w.expected_fertile * dt);
    }
    w.h_s = cur_[base + v.host_offsets[0]];
    w.h_w = cur_[base + v.host_offsets[1]];
    const double sites = a.sites;
    const double n_on_w =
        (sites > 0.0 && w.n_fertile > 0.0) ? draws.binom(w.n_fertile, w.h_w / sites) : 0.0;
    w.infested = draws.occupancy(n_on_w, w.h_w);
    w.healed = draws.binom(w.h_w - w.infested, 1.0 - std::exp(-v.host_healing * dt));
    std::fill(w.host_i.begin(), w.host_i.end(), 0.0);
    move(cur_.data() + base + v.host_offsets[2], w.host_i.data(), ch_i_, draws, &w.ended,
         &w.treated);
    w.wounded = draws.binom(w.h_s, 1.0 - std::exp(-v.host_wounding * dt));
  }
  if (w.empty) return;

  // ---- mating (S18), before ageing ------------------------------------------
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
        w.matings[g * ng + m] += buf_[m];
      }
    }
  }

  // ---- chain moves (S16), block by block in state order ------------------------
  // the larval hazard depends on crowding
  w.crowd = (v.larval_per_wound && v.node_H[node] > 0) || v.larval_per_infested_host;
  if (w.crowd) {
    for (int s = 0; s < nL; ++s) q_[s] = r.q_l;
    prepare(w.l, adv_l_, q_.data(), nL, nL, a.mu_l);
  }
  const PreparedChain& ch_l = w.crowd ? w.l : sc.l;
  for (int g = 0; g < ng; ++g) {
    move(egg + g * nE, n_e + g * nE, sc.e, draws, &w.out_e[g], &w.dth[0]);
  }
  for (int g = 0; g < ng; ++g) {
    move(lar + g * nL, n_l + g * nL, ch_l, draws, &w.out_l[g], &w.dth[1]);
  }
  for (int g = 0; g < ng; ++g) {
    move(pup + g * nP, n_p + g * nP, sc.p, draws, &w.out_p[g], &w.dth[2]);
  }
  for (int g = 0; g < ng; ++g) {
    move(mal + g * nA, n_m + g * nA, sc.m[g], draws, &w.sen_m, &w.dth[3]);
  }
  for (int g = 0; g < ng; ++g) {
    move(unm + g * nA, n_u + g * nA, sc.f[g], draws, &w.sen_u, &w.dth[4]);
  }
  for (int g = 0; g < ng; ++g) {
    for (int m = 0; m < ng; ++m) {
      const int off = (g * ng + m) * nA;
      move(fem + off, n_f + off, sc.f[g], draws, &w.sen_f, &w.dth[5]);
    }
  }
}

// Part 3 of node i's step: emergence (S19), arrivals, and the node's flows
// (written to `flows`, the node's block of counters) and host states.
template <class D>
void CbStepper::part3(int node, double* flows, D& draws) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const double dt = dt_;
  const std::size_t base = static_cast<std::size_t>(node) * v.node_stride;
  const SeriesChains& sc = series_[v.forcing_index[node]];
  const DayRates& r = sc.r;
  NodeWork& w = work_[node];
  const Aggregates& a = w.agg;
  double* n_e = nx_.data() + base + v.offsets[kE];
  double* n_l = nx_.data() + base + v.offsets[kL];
  double* n_p = nx_.data() + base + v.offsets[kP];
  double* n_m = nx_.data() + base + v.offsets[kM];
  double* n_u = nx_.data() + base + v.offsets[kU];
  double* n_f = nx_.data() + base + v.offsets[kF];
  const PreparedChain& ch_l = w.crowd ? w.l : sc.l;
  const bool batches = v.oviposition_batches;
  const double factor = r.batch_factor;
  const double nu_g = v.nu * a.g_find;

  // ---- emergence (S19): sex, then survival ------------------------------------
  double sum_em_m = 0.0, sum_em_f = 0.0, sum_z = 0.0;
  std::fill(w.em_m.begin(), w.em_m.end(), 0.0);
  std::fill(w.em_f.begin(), w.em_f.end(), 0.0);
  if (!w.empty) {
    for (int g = 0; g < ng; ++g) {
      const double z = w.out_p[g];
      const double zf = draws.binom(z, v.phi[g]);
      w.em_f[g] = draws.binom(zf, v.xiF[g]);
      w.em_m[g] = draws.binom(z - zf, v.xiM[g]);
      sum_z += z;
      sum_em_m += w.em_m[g];
      sum_em_f += w.em_f[g];
    }

    // ---- arrivals (S16), in state order ---------------------------------------
    for (int k = 0; k < ng; ++k) place(w.eggs[k], sc.e, n_e + k * nE, draws, &w.dth[0]);
    for (int g = 0; g < ng; ++g) place(w.out_e[g], ch_l, n_l + g * nL, draws, &w.dth[1]);
    for (int g = 0; g < ng; ++g) place(w.out_l[g], sc.p, n_p + g * nP, draws, &w.dth[2]);
    for (int g = 0; g < ng; ++g) place(w.em_m[g], sc.m[g], n_m + g * nA, draws, &w.dth[3]);
    // emerged females (S18, S19): death, then mating (at emergence with
    // "instant", otherwise for the rest of the step), then placement
    const double x_arr = (v.mating_instant ? nu_g : v.nu * v.receptive[0] * a.g_find) * dt;
    const double p_arr_mate = x_arr > 0.0 ? 1.0 + std::expm1(-x_arr) / x_arr : 0.0;
    for (int g = 0; g < ng; ++g) {
      const PreparedChain& ch = sc.f[g];
      std::fill(ym_.begin(), ym_.end(), 0.0);
      double unmated = 0.0;
      if (w.em_f[g] > 0.0) {
        const bool mated_now = v.mating_instant && a.mate_ok[g];
        const double d = draws.binom(w.em_f[g], ch.arr_die);
        w.dth[mated_now ? 5 : 4] += d;
        const double surv = w.em_f[g] - d;
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
        for (int m = 0; m < ng; ++m) w.matings[g * ng + m] += ym_[m];
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
  }

  // ---- flows (order of .sg_counter_names()) -----------------------------------
  int c = 0;
  for (int k = 0; k < ng; ++k) flows[c++] = w.eggs[k];
  flows[c++] = batches ? w.n_batches : a.psi * a.sum_search * factor * dt;  // egg_batches
  flows[c++] = (batches || v.host_model) ? w.n_fertile
                                          : a.psi * factor * w.expected_fertile * dt;  // fertile
  double sum_out_e = 0.0, sum_out_l = 0.0;
  for (int g = 0; g < ng; ++g) {
    sum_out_e += w.out_e[g];
    sum_out_l += w.out_l[g];
  }
  flows[c++] = sum_out_e;  // larvae_in
  flows[c++] = sum_out_l;  // pupae_in
  flows[c++] = sum_em_m;
  flows[c++] = sum_em_f;
  flows[c++] = sum_z - sum_em_m - sum_em_f;  // emergence_deaths
  for (int i = 0; i < ng * ng; ++i) flows[c++] = w.matings[i];
  for (int b = 0; b < 6; ++b) flows[c++] = w.dth[b];
  flows[c++] = w.sen_m;
  flows[c++] = w.sen_u;
  flows[c++] = w.sen_f;
  if (v.host_model) {  // host flows and new host states (D-030)
    flows[c++] = w.wounded;
    flows[c++] = w.healed;
    flows[c++] = w.infested;
    flows[c++] = w.ended;
    flows[c++] = w.treated;
    flows[c++] = w.n_fertile - w.infested;  // superinfestations
    nx_[base + v.host_offsets[0]] = w.h_s - w.wounded + w.healed + w.ended + w.treated;
    nx_[base + v.host_offsets[1]] = w.h_w - w.infested - w.healed + w.wounded;
    for (int j = 0; j < v.n_I; ++j) nx_[base + v.host_offsets[2] + j] = w.host_i[j];
    nx_[base + v.host_offsets[2]] += w.infested;
  }
}

}  // namespace screwdrive

#endif  // SCREWDRIVE_CB_STEP_H
