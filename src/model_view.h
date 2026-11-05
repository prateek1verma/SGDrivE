// Model view: the flat, read-only description of a model used by the C++
// backend (D-021). Built once per simulation from the R list produced by
// .sg_model_view() (R/model-view.R), whose layout it checks.
//
// Conventions (see R/model-view.R and dev/CONVENTIONS.md):
// * all indices and offsets are 0-based;
// * matrices are column-major std::vector<double>; field comments give their
//   dimensions as [rows, cols];
// * the state vector holds blocks E, L, P, M, U, F starting at offsets[0..5];
//   within a block the sub-stage varies fastest, then the mate genotype (F
//   only), then the genotype. With the host-wound model (D-030) the host
//   blocks HS, HW, HI (sub-stages 1..n_I) follow, at host_offsets[0..2].
// * several nodes (MODEL_SPEC 8a, D-032): the state vector holds n_node
//   per-node vectors of node_stride states (node slowest); all offsets are
//   within a node. Daily forcing is [n_days, n_series], with the series of
//   node i at forcing_index[i]. Counters: n_counters_node per node, node after
//   node, then 9 movement counters per node when there is movement.

#ifndef SCREWDRIVE_MODEL_VIEW_H
#define SCREWDRIVE_MODEL_VIEW_H

#include <Rcpp.h>

#include <string>
#include <vector>

namespace screwdrive {

// Positions of the state blocks in ModelView::offsets.
enum Block { kE = 0, kL = 1, kP = 2, kM = 3, kU = 4, kF = 5 };

struct ModelView {
  // dimensions; n_states and n_counters are totals over all nodes
  int n_g, n_E, n_L, n_P, n_A, n_states, n_days, n_counters;
  std::vector<int> offsets;  // [6], 0-based starts of E, L, P, M, U, F within a node

  // nodes (D-032)
  int n_node, node_stride, n_counters_node, n_series;
  std::vector<int> forcing_index;            // [n_node]: forcing series of each node
  std::vector<double> node_H, node_h_M, node_N_H;  // [n_node]: wounds, h_M A_i, hosts

  // movement groups (MODEL_SPEC 8a, N1-N3); n_groups = 0 without movement
  int n_groups;
  std::vector<int> group_clock;       // [n_groups]: 0 physiological, 1 calendar
  std::vector<double> group_rate;     // [3, n_node, n_groups]: delta^G (classes M, U, F)
  std::vector<int> group_row_ptr;     // [(n_node + 1) * n_groups], local to each group
  std::vector<int> group_edge_start;  // [n_groups]: first edge of each group
  std::vector<int> group_col;         // [n_edges]: destination node
  std::vector<double> group_prob;     // [3, n_edges]: pi^G per class
  std::vector<double> group_loss;     // [3, n_node, n_groups]: share leaving the landscape
  double dx_ref, kappa_max;           // dx(T_ref); cap of the clock (may be +Inf)
  static constexpr int kMoveCounters = 9;  // emigrants, immigrants, lost x (M, U, F)

  // host-wound model (D-030): host count of node 1 (N_H) and per-host rates per day
  bool host_model;
  int n_I;                        // infestation sub-stages (0 without hosts)
  std::vector<int> host_offsets;  // [3], 0-based starts of HS, HW, HI; -1 without hosts
  double N_H, host_wounding, host_healing, host_advance, host_treatment,
      host_attractiveness;

  // options
  bool mating_instant, search_all_adults, larval_per_wound;
  bool larval_per_infested_host;  // crowding per infested host (host model, D-030)
  bool oviposition_batches;  // chain-binomial egg laying in batches (D-024)

  // scalar parameters (counts of node 1 where relevant; node_* hold all nodes)
  double X_E, tau_L, X_P, R, a_s, sigma0, nu, mu_L0, gamma_L, mu_W, mu_b_F, mu_b_M;
  double H, h_M, h_M_param;
  double clutch_size;  // mean eggs per batch; 0 = R * phi_ov (D-024)

  // adult chains (S6, D-018)
  std::vector<double> inv_w_f, inv_w_m;  // [n_A, n_g]; 0 for classes never reached
  std::vector<int> last_f, last_m;       // [n_g]; number of classes reached (1-based)
  std::vector<int> laying, receptive, mature;  // [n_A]; 0/1 flags per class

  // cube
  std::vector<double> egg_tensor;  // [n_g * n_g, n_g]; rows (g, m) with m fastest
  std::vector<double> search_weight;    // [n_g]: female fecundity s_g (D-031)
  std::vector<double> batch_viability;  // [n_g * n_g]: sum_k ih tau per cross (D-031)
  std::vector<double> eta;         // [n_g, n_g]; eta[g, m]
  std::vector<double> phi, xiF, xiM;  // [n_g]
  std::vector<int> wild;              // [n_g]: 1 for the cube's wild type(s)

  // daily forcing (S1-S4, S7)
  std::vector<double> dx, dx_egg, mu_T, mu_E, phi_ov, moisture;  // [n_days, n_series]

  // releases: added at the start of day rel_day (1-based) to state rel_index
  std::vector<int> rel_day, rel_index;
  std::vector<double> rel_count;

  explicit ModelView(const Rcpp::List& view);

  // Position of a state in the full state vector (0-based).
  int idx(Block b, int g, int i) const {  // blocks E, L, P, M, U
    return offsets[b] + g * substages(b) + i;
  }
  int idx_f(int g, int m, int i) const {  // block F
    return offsets[kF] + (g * n_g + m) * n_A + i;
  }
  int substages(Block b) const {
    switch (b) {
      case kE:
        return n_E;
      case kL:
        return n_L;
      case kP:
        return n_P;
      default:
        return n_A;
    }
  }
};

}  // namespace screwdrive

#endif  // SCREWDRIVE_MODEL_VIEW_H
