// C++ rate kernel (see rates.h). Each block mirrors the corresponding block of
// .sg_rates_r() in R/reference-rates.R and cites the specification equation.

#include "rates.h"

#include <algorithm>
#include <cmath>

namespace screwdrive {

RateKernel::RateKernel(const ModelView& mv)
    : mv_(mv),
      in_m_(mv.n_g),
      in_f_(mv.n_g),
      in_u_(mv.n_g),
      in_f_mated_(mv.n_g),
      mating_flux_(static_cast<std::size_t>(mv.n_A) * mv.n_g) {
  const std::size_t ng = mv.n_g;
  agg_.w_lay.assign(ng * ng, 0.0);
  agg_.births.assign(ng, 0.0);
  agg_.m_bar.assign(ng, 0.0);
  agg_.pi.assign(ng * ng, 0.0);
  agg_.mate_ok.assign(ng, 0);
  const std::size_t n_cross = ng * ng;
  pair_egg_total_.assign(n_cross, 0.0);
  pair_search_.assign(n_cross, 0.0);
  v_batch_ = mv.batch_viability;
  for (std::size_t pair = 0; pair < n_cross; ++pair) {
    for (std::size_t k = 0; k < ng; ++k) pair_egg_total_[pair] += mv.egg_tensor[k * n_cross + pair];
    pair_search_[pair] = mv.search_weight[pair / ng];  // female genotype g = pair / n_g
  }
}

// Larval hazard (S7); with the host-wound model plus treatment and crowding
// per infested host (H2, floor 1e-6 * N_H as in .sg_crowd_floor).
double RateKernel::larval_hazard(double sum_l, double sum_hi) const {
  double mu_l = mv_.mu_L0;
  if (mv_.larval_per_wound && mv_.H > 0) mu_l += mv_.gamma_L * sum_l / mv_.H;
  if (mv_.host_model) {
    mu_l += mv_.host_treatment;
    if (mv_.larval_per_infested_host) {
      mu_l += mv_.gamma_L * sum_l / (sum_hi + 1e-6 * mv_.N_H);
    }
  }
  return mu_l;
}

double RateKernel::fertile_rate(const DayRates& r) const {
  const Aggregates& a = agg_;
  double fertile = 0.0;
  for (int pair = 0; pair < mv_.n_g * mv_.n_g; ++pair) {
    fertile += a.w_lay[pair] * pair_search_[pair] * -std::expm1(-r.clutch * v_batch_[pair]);
  }
  return a.psi * r.batch_factor * fertile;
}

DayRates RateKernel::day_rates(int day) const {
  const ModelView& v = mv_;
  const int d = day - 1;  // 0-based index into the daily forcing
  DayRates r;
  // daily forcing (S1-S4, S7)
  r.dx = v.dx[d];
  r.mu_e = v.mu_E[d];
  r.phi_ov = v.phi_ov[d];
  // transfer rates (S5, S6); adult rates per sub-stage are dx * inv_w
  r.q_e = v.n_E * v.dx_egg[d] / v.X_E;
  r.q_l = v.n_L / v.tau_L;
  r.q_p = v.n_P * r.dx / v.X_P;
  // hazards (S7)
  r.mu_p = v.mu_T[d] + v.mu_W * (1.0 - v.moisture[d]);
  r.mu_f = v.mu_T[d] + v.mu_b_F;
  r.mu_m = v.mu_T[d] + v.mu_b_M;
  // oviposition events (D-024): a batch is the day's eggs unless clutch_size is set
  r.clutch = v.clutch_size > 0 ? v.clutch_size : v.R * r.phi_ov;
  r.batch_factor = v.clutch_size > 0 ? v.R * r.phi_ov / v.clutch_size : 1.0;
  return r;
}

void RateKernel::aggregate(const DayRates& r, const double* y) {
  const ModelView& v = mv_;
  Aggregates& a = agg_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const double* egg = y + v.offsets[kE];
  const double* lar = y + v.offsets[kL];
  const double* pup = y + v.offsets[kP];
  const double* mal = y + v.offsets[kM];
  const double* unm = y + v.offsets[kU];
  const double* fem = y + v.offsets[kF];

  // ---- block totals ----------------------------------------------------------
  a.sum_e = a.sum_l = a.sum_p = a.sum_m = a.sum_u = a.sum_f = 0.0;
  for (int i = 0; i < ng * nE; ++i) a.sum_e += egg[i];
  for (int i = 0; i < ng * nL; ++i) a.sum_l += lar[i];
  for (int i = 0; i < ng * nP; ++i) a.sum_p += pup[i];
  for (int i = 0; i < ng * nA; ++i) a.sum_m += mal[i];
  for (int i = 0; i < ng * nA; ++i) a.sum_u += unm[i];
  for (int i = 0; i < ng * ng * nA; ++i) a.sum_f += fem[i];

  // ---- hosts (H1, H2) ------------------------------------------------------------
  a.host_s = a.host_w = a.sum_hi = a.sites = 0.0;
  if (v.host_model) {
    a.host_s = y[v.host_offsets[0]];
    a.host_w = y[v.host_offsets[1]];
    const double* hi = y + v.host_offsets[2];
    for (int j = 0; j < v.n_I; ++j) a.sum_hi += hi[j];
    a.sites = a.host_w + v.host_attractiveness * a.sum_hi;
  }

  // ---- larval crowding (S7) ------------------------------------------------------
  a.mu_l = larval_hazard(a.sum_l, a.sum_hi);

  // ---- wound search and egg production (S8, S9) ----------------------------
  a.sum_w = a.sum_search = 0.0;
  for (int g = 0; g < ng; ++g) {
    for (int m = 0; m < ng; ++m) {
      const double* x = fem + (g * ng + m) * nA;
      double w = 0.0;
      for (int s = 0; s < nA; ++s) {
        if (v.laying[s]) w += x[s];
      }
      a.w_lay[g * ng + m] = w;
      a.sum_w += w;
      a.sum_search += v.search_weight[g] * w;
    }
  }
  // searching females are weighted by fecundity: females that produce no
  // eggs do not compete for wounds (D-031)
  a.n_s = v.search_all_adults ? (a.sum_m + a.sum_u + a.sum_f) : a.sum_search;
  const double wounds = v.host_model ? a.sites : v.H;  // fixed H, or A_s (H2)
  a.psi = a.n_s > 0.0 ? v.sigma0 * (1.0 - std::exp(-v.a_s * r.dx * wounds / a.n_s)) : v.sigma0;
  const int n_cross = ng * ng;
  for (int k = 0; k < ng; ++k) {
    const double* col = v.egg_tensor.data() + static_cast<std::size_t>(k) * n_cross;
    double b = 0.0;
    for (int c = 0; c < n_cross; ++c) b += col[c] * a.w_lay[c];
    a.births[k] = v.R * r.phi_ov * a.psi * b;
  }

  // ---- male availability and mating (S12) ----------------------------------
  a.m_tot = 0.0;
  for (int m = 0; m < ng; ++m) {
    const double* x = mal + m * nA;
    double mb = 0.0;
    for (int s = 0; s < nA; ++s) {
      if (v.mature[s]) mb += x[s];
    }
    a.m_bar[m] = mb;
    a.m_tot += mb;
  }
  if (v.h_M_param > 0) {
    a.g_find = a.m_tot / (v.h_M + a.m_tot);
  } else {
    a.g_find = a.m_tot > 0.0 ? 1.0 : 0.0;
  }
  for (int g = 0; g < ng; ++g) {
    double pref_sum = 0.0;
    for (int m = 0; m < ng; ++m) {
      const double pref = v.eta[g + m * ng] * a.m_bar[m];  // eta[g, m] * M_bar[m]
      a.pi[g * ng + m] = pref;
      pref_sum += pref;
    }
    a.mate_ok[g] = pref_sum > 0.0;
    for (int m = 0; m < ng; ++m) {
      a.pi[g * ng + m] = a.mate_ok[g] ? a.pi[g * ng + m] / pref_sum : 0.0;
    }
  }
}

void RateKernel::eval(int day, const double* y, double* dydt, double* counters) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const double* egg = y + v.offsets[kE];
  const double* lar = y + v.offsets[kL];
  const double* pup = y + v.offsets[kP];
  const double* mal = y + v.offsets[kM];
  const double* unm = y + v.offsets[kU];
  const double* fem = y + v.offsets[kF];
  double* d_e = dydt + v.offsets[kE];
  double* d_l = dydt + v.offsets[kL];
  double* d_p = dydt + v.offsets[kP];
  double* d_m = dydt + v.offsets[kM];
  double* d_u = dydt + v.offsets[kU];
  double* d_f = dydt + v.offsets[kF];

  const DayRates r = day_rates(day);
  aggregate(r, y);
  const Aggregates& a = agg_;
  const double dx = r.dx, mu_e = r.mu_e, q_e = r.q_e, q_l = r.q_l, q_p = r.q_p;
  const double mu_l = a.mu_l, mu_p = r.mu_p, mu_f = r.mu_f, mu_m = r.mu_m;
  const double g_find = a.g_find, psi = a.psi;
  const std::vector<int>& mate_ok_ = a.mate_ok;
  const std::vector<double>& births_ = a.births;
  const std::vector<double>& pi_ = a.pi;

  // ---- emergence (S10, S11) ------------------------------------------------------
  double sum_emerge = 0.0, sum_in_m = 0.0, sum_in_f = 0.0;
  for (int g = 0; g < ng; ++g) {
    const double emerge = q_p * pup[g * nP + (nP - 1)];
    in_m_[g] = (1.0 - v.phi[g]) * v.xiM[g] * emerge;
    in_f_[g] = v.phi[g] * v.xiF[g] * emerge;
    sum_emerge += emerge;
    sum_in_m += in_m_[g];
    sum_in_f += in_f_[g];
  }
  const double emergence_deaths = sum_emerge - sum_in_m - sum_in_f;

  // Mating flux per virgin [substage, genotype] (S12, S14).
  for (int g = 0; g < ng; ++g) {
    if (v.mating_instant) {
      in_u_[g] = mate_ok_[g] ? 0.0 : in_f_[g];
      in_f_mated_[g] = mate_ok_[g] ? in_f_[g] : 0.0;
    } else {
      in_u_[g] = in_f_[g];
      in_f_mated_[g] = 0.0;
    }
    for (int s = 0; s < nA; ++s) {
      double rate = 0.0;
      if (mate_ok_[g]) {
        rate = v.mating_instant ? v.nu * g_find : v.nu * v.receptive[s] * g_find;
      }
      mating_flux_[g * nA + s] = rate * unm[g * nA + s];
    }
  }

  // ---- derivatives: juvenile Erlang chains (S10) -------------------------------
  auto chain_const = [](const double* x, double* dx_out, int n, double inflow, double q,
                        double mu) {
    for (int s = 0; s < n; ++s) dx_out[s] = -(q + mu) * x[s];
    dx_out[0] += inflow;
    for (int s = 1; s < n; ++s) dx_out[s] += q * x[s - 1];
  };
  for (int g = 0; g < ng; ++g) {
    chain_const(egg + g * nE, d_e + g * nE, nE, births_[g], q_e, mu_e);
    chain_const(lar + g * nL, d_l + g * nL, nL, q_e * egg[g * nE + (nE - 1)], q_l, mu_l);
    chain_const(pup + g * nP, d_p + g * nP, nP, q_l * lar[g * nL + (nL - 1)], q_p, mu_p);
  }

  // ---- derivatives: segmented adult chains (S6, S13) ---------------------------
  // Outflow of class s goes to s + 1 while s < last - 1 (0-based), otherwise
  // it is death by old age.
  double sen_m = 0.0, sen_u = 0.0, sen_f = 0.0;
  for (int g = 0; g < ng; ++g) {
    const double* iw_m = v.inv_w_m.data() + static_cast<std::size_t>(g) * nA;
    const double* iw_f = v.inv_w_f.data() + static_cast<std::size_t>(g) * nA;
    const int last_m = v.last_m[g], last_f = v.last_f[g];

    // males
    const double* xm = mal + g * nA;
    double* dm = d_m + g * nA;
    for (int s = 0; s < nA; ++s) dm[s] = -mu_m * xm[s];
    for (int s = 0; s < nA; ++s) {
      const double out = dx * iw_m[s] * xm[s];
      dm[s] -= out;
      if (s < last_m - 1) {
        dm[s + 1] += out;
      } else {
        sen_m += out;
      }
    }
    dm[0] += in_m_[g];

    // unmated females
    const double* xu = unm + g * nA;
    double* du = d_u + g * nA;
    for (int s = 0; s < nA; ++s) du[s] = -mu_f * xu[s] - mating_flux_[g * nA + s];
    for (int s = 0; s < nA; ++s) {
      const double out = dx * iw_f[s] * xu[s];
      du[s] -= out;
      if (s < last_f - 1) {
        du[s + 1] += out;
      } else {
        sen_u += out;
      }
    }
    du[0] += in_u_[g];

    // mated females of genotype g with each mate genotype m
    for (int m = 0; m < ng; ++m) {
      const double* xf = fem + (g * ng + m) * nA;
      double* df = d_f + (g * ng + m) * nA;
      const double pi_gm = pi_[g * ng + m];
      for (int s = 0; s < nA; ++s) df[s] = -mu_f * xf[s] + mating_flux_[g * nA + s] * pi_gm;
      for (int s = 0; s < nA; ++s) {
        const double out = dx * iw_f[s] * xf[s];
        df[s] -= out;
        if (s < last_f - 1) {
          df[s + 1] += out;
        } else {
          sen_f += out;
        }
      }
      df[0] += in_f_mated_[g] * pi_gm;
    }
  }

  // ---- hosts (H1, H2) ------------------------------------------------------------
  // Fertile batches are shared among wounds by attractiveness: on an uninfested
  // wound a batch starts an infestation (a case), on an infested one it only
  // adds larvae (superinfestation).
  const double fertile = fertile_rate(r);
  double new_inf = 0.0, superinf = 0.0, wounded = 0.0, healed = 0.0, ended = 0.0;
  double treated = 0.0;
  if (v.host_model) {
    const int nI = v.n_I;
    const double* hi = y + v.host_offsets[2];
    double* d_hi = dydt + v.host_offsets[2];
    if (a.sites > 0.0) {
      new_inf = fertile * a.host_w / a.sites;
      superinf = fertile * v.host_attractiveness * a.sum_hi / a.sites;
    }
    wounded = v.host_wounding * a.host_s;
    healed = v.host_healing * a.host_w;
    ended = v.host_advance * hi[nI - 1];
    treated = v.host_treatment * a.sum_hi;
    for (int j = 0; j < nI; ++j) d_hi[j] = -(v.host_advance + v.host_treatment) * hi[j];
    d_hi[0] += new_inf;
    for (int j = 1; j < nI; ++j) d_hi[j] += v.host_advance * hi[j - 1];
    dydt[v.host_offsets[0]] = -wounded + healed + ended + treated;
    dydt[v.host_offsets[1]] = wounded - healed - new_inf;
  }

  if (counters == nullptr) return;

  // ---- counters (rates per day; order of .sg_counter_names()) ---------------------
  int c = 0;
  for (int k = 0; k < ng; ++k) counters[c++] = births_[k];
  double last_e = 0.0, last_l = 0.0;
  for (int g = 0; g < ng; ++g) {
    last_e += egg[g * nE + (nE - 1)];
    last_l += lar[g * nL + (nL - 1)];
  }
  // egg_batches: oviposition events per day (D-024); fertile_batches: those
  // with at least one viable egg, the case proxy (D-026)
  counters[c++] = psi * a.sum_search * r.batch_factor;
  counters[c++] = fertile;
  counters[c++] = q_e * last_e;    // larvae_in
  counters[c++] = q_l * last_l;    // pupae_in
  counters[c++] = sum_in_m;        // emerged_male
  counters[c++] = sum_in_f;        // emerged_female
  counters[c++] = emergence_deaths;
  for (int g = 0; g < ng; ++g) {  // matings[g|m], m fastest
    double flux_g = in_f_mated_[g];
    for (int s = 0; s < nA; ++s) flux_g += mating_flux_[g * nA + s];
    for (int m = 0; m < ng; ++m) counters[c++] = flux_g * pi_[g * ng + m];
  }
  counters[c++] = mu_e * a.sum_e;
  counters[c++] = mu_l * a.sum_l;
  counters[c++] = mu_p * a.sum_p;
  counters[c++] = mu_m * a.sum_m;
  counters[c++] = mu_f * a.sum_u;
  counters[c++] = mu_f * a.sum_f;
  counters[c++] = sen_m;
  counters[c++] = sen_u;
  counters[c++] = sen_f;
  if (v.host_model) {  // host flows (D-030), order of .sg_host_counter_names
    counters[c++] = wounded;
    counters[c++] = healed;
    counters[c++] = new_inf;
    counters[c++] = ended;
    counters[c++] = treated;
    counters[c++] = superinf;
  }
}

double RateKernel::max_rate(int day, const double* y) {
  const ModelView& v = mv_;
  const DayRates r = day_rates(day);
  aggregate(r, y);  // larval crowding, male availability (S7, S12)
  const Aggregates& a = agg_;
  const double dx = r.dx, mu_f = r.mu_f, mu_m = r.mu_m;
  double best = std::max({r.q_e + r.mu_e, r.q_l + a.mu_l, r.q_p + r.mu_p});
  // unmated females: ageing, death and mating (D-026); mated females age and
  // die at the same rates without mating, so they never exceed this
  for (int g = 0; g < v.n_g; ++g) {
    for (int s = 0; s < v.n_A; ++s) {
      const double iw = v.inv_w_f[static_cast<std::size_t>(g) * v.n_A + s];
      if (iw <= 0) continue;
      double nu_s = 0.0;
      if (a.mate_ok[g]) {
        nu_s = v.mating_instant ? v.nu * a.g_find : v.nu * v.receptive[s] * a.g_find;
      }
      best = std::max(best, dx * iw + mu_f + nu_s);
    }
  }
  double max_iw_m = -1.0;
  for (double iw : v.inv_w_m) {
    if (iw > 0) max_iw_m = std::max(max_iw_m, iw);
  }
  if (max_iw_m > 0) best = std::max(best, dx * max_iw_m + mu_m);
  if (v.host_model) {  // host rates (D-030): wounding, healing + infestation, chain
    const double hazard = a.host_w > 0.0 && a.sites > 0.0 ? fertile_rate(r) / a.sites : 0.0;
    best = std::max({best, v.host_wounding, v.host_healing + hazard,
                     v.host_advance + v.host_treatment});
  }
  return best;
}

}  // namespace screwdrive

// Internal test entry points: evaluate the kernel for one state.
// [[Rcpp::export]]
Rcpp::List sg_rates_cpp(const Rcpp::List& view, int day, const Rcpp::NumericVector& state) {
  const screwdrive::ModelView mv(view);
  if (state.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  screwdrive::RateKernel kernel(mv);
  Rcpp::NumericVector deriv(mv.n_states), counters(mv.n_counters);
  kernel.eval(day, state.begin(), deriv.begin(), counters.begin());
  return Rcpp::List::create(Rcpp::Named("deriv") = deriv, Rcpp::Named("counters") = counters);
}

// [[Rcpp::export]]
double sg_max_rate_cpp(const Rcpp::List& view, int day, const Rcpp::NumericVector& state) {
  const screwdrive::ModelView mv(view);
  if (state.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  screwdrive::RateKernel kernel(mv);
  return kernel.max_rate(day, state.begin());
}

// Day rates and state aggregates for one state (tests of the shared parts
// used by both solvers).
// [[Rcpp::export]]
Rcpp::List sg_aggregates_cpp(const Rcpp::List& view, int day, const Rcpp::NumericVector& state) {
  const screwdrive::ModelView mv(view);
  if (state.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  screwdrive::RateKernel kernel(mv);
  const screwdrive::DayRates r = kernel.day_rates(day);
  kernel.aggregate(r, state.begin());
  const screwdrive::Aggregates& a = kernel.agg();
  return Rcpp::List::create(
      Rcpp::Named("q") = Rcpp::NumericVector::create(r.q_e, r.q_l, r.q_p),
      Rcpp::Named("mu") = Rcpp::NumericVector::create(r.mu_e, a.mu_l, r.mu_p, r.mu_f, r.mu_m),
      Rcpp::Named("dx") = r.dx, Rcpp::Named("phi_ov") = r.phi_ov,
      Rcpp::Named("w_lay") = Rcpp::wrap(a.w_lay), Rcpp::Named("n_s") = a.n_s,
      Rcpp::Named("psi") = a.psi, Rcpp::Named("births") = Rcpp::wrap(a.births),
      Rcpp::Named("m_bar") = Rcpp::wrap(a.m_bar), Rcpp::Named("g_find") = a.g_find,
      Rcpp::Named("pi") = Rcpp::wrap(a.pi), Rcpp::Named("mate_ok") = Rcpp::wrap(a.mate_ok),
      Rcpp::Named("sites") = a.sites, Rcpp::Named("sum_hi") = a.sum_hi);
}
