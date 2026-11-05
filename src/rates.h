// C++ rate kernel of the ScrewDrivE model (D-021).
//
// A block-by-block translation of the R reference rate function
// .sg_rates_r() (R/reference-rates.R), which implements dev/MODEL_SPEC.md
// S7-S14 and, for several nodes, section 8a N1-N3 (D-032). The two must agree
// to floating-point round-off; the R version is the readable reference and
// the test oracle.
//
// Several nodes: eval() runs the one-node body on each node's part of the
// state (with the node's forcing, wounds, mate finding and hosts, selected by
// set_node()), then adds the movement terms. A model with one node and no
// movement runs the one-node body exactly as before.
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

  // Select the node used by day_rates() and aggregate() (default node 0):
  // its forcing series, wounds H, mate-finding h_M and host count N_H.
  void set_node(int node);

  // Forcing-dependent rates of simulation day `day` (1-based) for the
  // selected node.
  DayRates day_rates(int day) const;

  // State-dependent aggregates for the selected node's state `y` (pointer to
  // the node's first state; stored in agg()).
  void aggregate(const DayRates& r, const double* y);
  const Aggregates& agg() const { return agg_; }

  // Derivatives of all states (dydt, length n_states) and rates of the
  // cumulative counters (counters, length n_counters; order of
  // .sg_model_counter_names()) on simulation day `day` (1-based), whose
  // forcing is held constant. `counters` may be nullptr when not needed.
  void eval(int day, const double* y, double* dydt, double* counters);

  // Largest total rate (transfer + mortality) of any compartment in use on
  // `day` for state `y` (solver diagnostics; .sg_max_rate() in R); with
  // several nodes the largest over nodes of the node's rate plus its total
  // leave rate.
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
  // selected node (set_node): forcing offset and node quantities
  int node_ = 0, forcing_off_ = 0;
  double H_ = 0.0, h_M_ = 0.0, N_H_ = 0.0;
  // [3]: total leave rate of the selected node by class, for max_rate()
  double lambda_[3] = {0.0, 0.0, 0.0};
  std::vector<double> pair_egg_total_, pair_search_, v_batch_;
  // working arrays of eval() (sizes in comments)
  std::vector<double> in_m_, in_f_, in_u_, in_f_mated_;  // [n_g]
  std::vector<double> mating_flux_;  // [n_A * n_g]: [substage, genotype]

  // one-node body of eval() and max_rate() for the selected node; y, dydt and
  // counters point to the node's first state and counter
  void eval_node(int day, const double* y, double* dydt, double* counters);
  double max_rate_node(int day, const double* y);
  // movement clock kappa (N3) of group g for node i on day `day`
  double kappa(int g, int i, int day) const;
  // movement terms (N1, N2) of all nodes; move_counters (9 per node) may be
  // nullptr
  void movement(int day, const double* y, double* dydt, double* move_counters);
  // total leave rate of node i by class, summed over groups (lambda_)
  void node_leave_rates(int i, int day);
  double larval_hazard(double sum_l, double sum_hi) const;
  // fertile oviposition events per day (D-026): psi * batch factor * sum of
  // searching females (s_g-weighted, D-031) times the share of batches with a
  // viable egg
  double fertile_rate(const DayRates& r) const;
};

}  // namespace screwdrive

#endif  // SCREWDRIVE_RATES_H
