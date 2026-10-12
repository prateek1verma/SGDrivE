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
//   only), then the genotype.

#ifndef SGDRIVE_MODEL_VIEW_H
#define SGDRIVE_MODEL_VIEW_H

#include <Rcpp.h>

#include <string>
#include <vector>

namespace sgdrive {

// Positions of the state blocks in ModelView::offsets.
enum Block { kE = 0, kL = 1, kP = 2, kM = 3, kU = 4, kF = 5 };

struct ModelView {
  // dimensions
  int n_g, n_E, n_L, n_P, n_A, n_states, n_days, n_counters;
  std::vector<int> offsets;  // [6], 0-based starts of E, L, P, M, U, F

  // options
  bool mating_instant, search_all_adults, larval_per_wound;

  // scalar parameters (counts per node where relevant)
  double X_E, tau_L, X_P, R, a_s, sigma0, nu, mu_L0, gamma_L, mu_W, mu_b_F, mu_b_M;
  double H, h_M, h_M_param;

  // adult chains (S6, D-018)
  std::vector<double> inv_w_f, inv_w_m;  // [n_A, n_g]; 0 for classes never reached
  std::vector<int> last_f, last_m;       // [n_g]; number of classes reached (1-based)
  std::vector<int> laying, receptive, mature;  // [n_A]; 0/1 flags per class

  // cube
  std::vector<double> egg_tensor;  // [n_g * n_g, n_g]; rows (g, m) with m fastest
  std::vector<double> eta;         // [n_g, n_g]; eta[g, m]
  std::vector<double> phi, xiF, xiM;  // [n_g]

  // daily forcing (S1-S4, S7)
  std::vector<double> dx, dx_egg, mu_T, mu_E, phi_ov, moisture;  // [n_days]

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

}  // namespace sgdrive

#endif  // SGDRIVE_MODEL_VIEW_H
