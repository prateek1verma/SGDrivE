// C++ rate kernel of the SGDrivE v0.1 model (D-021).
//
// A block-by-block translation of the R reference rate function
// .sg_rates_r() (R/reference-rates.R), which implements dev/MODEL_SPEC.md
// S7-S14. The two must agree to floating-point round-off; the R version is
// the readable reference and the test oracle.
//
// Working arrays are allocated once in the constructor; eval() allocates
// nothing.

#ifndef SGDRIVE_RATES_H
#define SGDRIVE_RATES_H

#include <vector>

#include "model_view.h"

namespace sgdrive {

class RateKernel {
 public:
  explicit RateKernel(const ModelView& mv);

  // Derivatives of all states (dydt, length n_states) and rates of the
  // cumulative counters (counters, length n_counters; order of
  // .sg_counter_names()) on simulation day `day` (1-based), whose forcing is
  // held constant. `counters` may be nullptr when not needed.
  void eval(int day, const double* y, double* dydt, double* counters);

  // Largest total rate (transfer + mortality) of any compartment in use on
  // `day` for state `y` (solver diagnostics; .sg_max_rate() in R).
  double max_rate(int day, const double* y);

  const ModelView& view() const { return mv_; }

 private:
  const ModelView& mv_;
  // working arrays (sizes in comments)
  std::vector<double> w_lay_;        // [n_g * n_g]: egg-laying females, (g, m) with m fastest
  std::vector<double> births_;       // [n_g]
  std::vector<double> m_bar_;        // [n_g]: mature males
  std::vector<double> pi_;           // [n_g * n_g]: pi[g, m] at g * n_g + m
  std::vector<int> mate_ok_;         // [n_g]
  std::vector<double> in_m_, in_f_, in_u_, in_f_mated_;  // [n_g]
  std::vector<double> mating_flux_;  // [n_A * n_g]: [substage, genotype]
};

}  // namespace sgdrive

#endif  // SGDRIVE_RATES_H
