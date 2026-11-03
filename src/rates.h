// C++ rate kernel of the ScrewDrivE v0.1 model (D-021).
//
// A block-by-block translation of the R reference rate function
// .sg_rates_r() (R/reference-rates.R), which implements dev/MODEL_SPEC.md
// S7-S14. The two must agree to floating-point round-off; the R version is
// the readable reference and the test oracle.
//
// The day's forcing-dependent rates (day_rates) and the state-dependent
// aggregates (aggregate) are separate steps so that the chain-binomial solver
// reuses them (D-010, D-022). Working arrays are allocated once in the
// constructor; eval() and aggregate() allocate nothing.

#ifndef SCREWDRIVE_RATES_H
#define SCREWDRIVE_RATES_H

#include <vector>

#include "model_view.h"

namespace screwdrive {

// Quantities that depend only on the day's forcing (S1-S5, S7 without
// density dependence). Shared by the ODE and chain-binomial solvers.
struct DayRates {
  double dx = 0.0;      // degree-days per day, adults and pupae
  double mu_e = 0.0;    // egg hazard
  double phi_ov = 0.0;  // oviposition scalar
  double q_e = 0.0, q_l = 0.0, q_p = 0.0;  // juvenile transfer rates per sub-stage
  double mu_p = 0.0, mu_f = 0.0, mu_m = 0.0;  // pupal, female and male hazards
  double clutch = 0.0;        // mean eggs per oviposition event (D-024)
  double batch_factor = 1.0;  // oviposition events per wound found (D-024)
};

// Quantities that depend on the state (S7 larval crowding, S8, S9, S12; host
// sums H1-H2 of the host-wound model, D-030).
// Shared by the ODE and chain-binomial solvers.
struct Aggregates {
  double sum_e = 0.0, sum_l = 0.0, sum_p = 0.0, sum_m = 0.0, sum_u = 0.0, sum_f = 0.0;
  double host_s = 0.0, host_w = 0.0, sum_hi = 0.0;  // hosts: no wound, wounded, infested
  double sites = 0.0;   // attractive sites A_s = HW + alpha * sum(HI) (host model)
  double mu_l = 0.0;    // larval hazard, with crowding
  double sum_w = 0.0;   // egg-laying mated females
  double sum_search = 0.0;  // egg-laying females weighted by fecundity s_g (D-031)
  double n_s = 0.0;     // wound-search competitors
  double psi = 0.0;     // wound-finding rate
  double m_tot = 0.0;   // mature males
  double g_find = 0.0;  // probability scale of finding a male
  std::vector<double> w_lay;  // [n_g * n_g]: egg-laying females (g, m), m fastest
  std::vector<double> births; // [n_g]: egg production by offspring genotype
  std::vector<double> m_bar;  // [n_g]: mature males by genotype
  std::vector<double> pi;     // [n_g * n_g]: mate shares pi[g, m] at g * n_g + m
  std::vector<int> mate_ok;   // [n_g]: compatible mature males present
};

class RateKernel {
 public:
  explicit RateKernel(const ModelView& mv);

  // Forcing-dependent rates of simulation day `day` (1-based).
  DayRates day_rates(int day) const;

  // State-dependent aggregates for state `y` (stored in agg()).
  void aggregate(const DayRates& r, const double* y);
  const Aggregates& agg() const { return agg_; }

  // Derivatives of all states (dydt, length n_states) and rates of the
  // cumulative counters (counters, length n_counters; order of
  // .sg_counter_names()) on simulation day `day` (1-based), whose forcing is
  // held constant. `counters` may be nullptr when not needed.
  void eval(int day, const double* y, double* dydt, double* counters);

  // Largest total rate (transfer + mortality) of any compartment in use on
  // `day` for state `y` (solver diagnostics; .sg_max_rate() in R).
  double max_rate(int day, const double* y);

  const ModelView& view() const { return mv_; }

  // Per female x mate pair ([n_g * n_g], m fastest): viable eggs per egg
  // laid, sum_k s ih tau; search and batch weight, the female's s_g; and
  // viable eggs per egg of a batch, sum_k ih tau (D-031).
  const std::vector<double>& pair_egg_total() const { return pair_egg_total_; }
  const std::vector<double>& pair_search() const { return pair_search_; }
  const std::vector<double>& pair_batch_viability() const { return v_batch_; }

 private:
  const ModelView& mv_;
  Aggregates agg_;
  std::vector<double> pair_egg_total_, pair_search_, v_batch_;
  // working arrays of eval() (sizes in comments)
  std::vector<double> in_m_, in_f_, in_u_, in_f_mated_;  // [n_g]
  std::vector<double> mating_flux_;  // [n_A * n_g]: [substage, genotype]

  double larval_hazard(double sum_l, double sum_hi) const;
  // fertile oviposition events per day (D-026): psi * batch factor * sum of
  // searching females (s_g-weighted, D-031) times the share of batches with a
  // viable egg
  double fertile_rate(const DayRates& r) const;
};

}  // namespace screwdrive

#endif  // SCREWDRIVE_RATES_H
