// C++ rate kernel (see rates.h). Each block mirrors the corresponding block of
// .sg_rates_r() in R/reference-rates.R and cites the specification equation.

#include "rates.h"

#include <algorithm>
#include <cmath>

namespace sgdrive {

RateKernel::RateKernel(const ModelView& mv)
    : mv_(mv),
      w_lay_(static_cast<std::size_t>(mv.n_g) * mv.n_g),
      births_(mv.n_g),
      m_bar_(mv.n_g),
      pi_(static_cast<std::size_t>(mv.n_g) * mv.n_g),
      mate_ok_(mv.n_g),
      in_m_(mv.n_g),
      in_f_(mv.n_g),
      in_u_(mv.n_g),
      in_f_mated_(mv.n_g),
      mating_flux_(static_cast<std::size_t>(mv.n_A) * mv.n_g) {}

void RateKernel::eval(int day, const double* y, double* dydt, double* counters) {
  const ModelView& v = mv_;
  const int ng = v.n_g, nE = v.n_E, nL = v.n_L, nP = v.n_P, nA = v.n_A;
  const int d = day - 1;  // 0-based index into the daily forcing
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

  // ---- daily forcing (S1-S4, S7) -------------------------------------------
  const double dx = v.dx[d];
  const double mu_e = v.mu_E[d];
  const double phi_ov = v.phi_ov[d];

  // ---- transfer rates (S5, S6) ---------------------------------------------
  const double q_e = nE * v.dx_egg[d] / v.X_E;
  const double q_l = nL / v.tau_L;
  const double q_p = nP * dx / v.X_P;
  // adult rates per sub-stage: q[s, g] = dx * inv_w[s + g * nA] (0 if never reached)

  // ---- block totals ----------------------------------------------------------
  double sum_e = 0.0, sum_l = 0.0, sum_p = 0.0, sum_m = 0.0, sum_u = 0.0, sum_f = 0.0;
  for (int i = 0; i < ng * nE; ++i) sum_e += egg[i];
  for (int i = 0; i < ng * nL; ++i) sum_l += lar[i];
  for (int i = 0; i < ng * nP; ++i) sum_p += pup[i];
  for (int i = 0; i < ng * nA; ++i) sum_m += mal[i];
  for (int i = 0; i < ng * nA; ++i) sum_u += unm[i];
  for (int i = 0; i < ng * ng * nA; ++i) sum_f += fem[i];

  // ---- hazards (S7) ------------------------------------------------------------
  double mu_l = v.mu_L0;
  if (v.larval_per_wound && v.H > 0) mu_l += v.gamma_L * sum_l / v.H;
  const double mu_p = v.mu_T[d] + v.mu_W * (1.0 - v.moisture[d]);
  const double mu_f = v.mu_T[d] + v.mu_b_F;
  const double mu_m = v.mu_T[d] + v.mu_b_M;

  // ---- wound search and egg production (S8, S9) ----------------------------
  double sum_w = 0.0;
  for (int g = 0; g < ng; ++g) {
    for (int m = 0; m < ng; ++m) {
      const double* x = fem + (g * ng + m) * nA;
      double w = 0.0;
      for (int s = 0; s < nA; ++s) {
        if (v.laying[s]) w += x[s];
      }
      w_lay_[g * ng + m] = w;
      sum_w += w;
    }
  }
  const double n_s = v.search_all_adults ? (sum_m + sum_u + sum_f) : sum_w;
  const double psi =
      n_s > 0.0 ? v.sigma0 * (1.0 - std::exp(-v.a_s * dx * v.H / n_s)) : v.sigma0;
  const int n_cross = ng * ng;
  for (int k = 0; k < ng; ++k) {
    const double* col = v.egg_tensor.data() + static_cast<std::size_t>(k) * n_cross;
    double b = 0.0;
    for (int r = 0; r < n_cross; ++r) b += col[r] * w_lay_[r];
    births_[k] = v.R * phi_ov * psi * b;
  }

  // ---- male availability and mating (S12) ----------------------------------
  double m_tot = 0.0;
  for (int m = 0; m < ng; ++m) {
    const double* x = mal + m * nA;
    double mb = 0.0;
    for (int s = 0; s < nA; ++s) {
      if (v.mature[s]) mb += x[s];
    }
    m_bar_[m] = mb;
    m_tot += mb;
  }
  double g_find;
  if (v.h_M_param > 0) {
    g_find = m_tot / (v.h_M + m_tot);
  } else {
    g_find = m_tot > 0.0 ? 1.0 : 0.0;
  }
  for (int g = 0; g < ng; ++g) {
    double pref_sum = 0.0;
    for (int m = 0; m < ng; ++m) {
      const double pref = v.eta[g + m * ng] * m_bar_[m];  // eta[g, m] * M_bar[m]
      pi_[g * ng + m] = pref;
      pref_sum += pref;
    }
    mate_ok_[g] = pref_sum > 0.0;
    for (int m = 0; m < ng; ++m) {
      pi_[g * ng + m] = mate_ok_[g] ? pi_[g * ng + m] / pref_sum : 0.0;
    }
  }

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

  if (counters == nullptr) return;

  // ---- counters (rates per day; order of .sg_counter_names()) ---------------------
  int c = 0;
  for (int k = 0; k < ng; ++k) counters[c++] = births_[k];
  double last_e = 0.0, last_l = 0.0;
  for (int g = 0; g < ng; ++g) {
    last_e += egg[g * nE + (nE - 1)];
    last_l += lar[g * nL + (nL - 1)];
  }
  counters[c++] = psi * sum_w;     // egg_batches
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
  counters[c++] = mu_e * sum_e;
  counters[c++] = mu_l * sum_l;
  counters[c++] = mu_p * sum_p;
  counters[c++] = mu_m * sum_m;
  counters[c++] = mu_f * sum_u;
  counters[c++] = mu_f * sum_f;
  counters[c++] = sen_m;
  counters[c++] = sen_u;
  counters[c++] = sen_f;
}

double RateKernel::max_rate(int day, const double* y) {
  const ModelView& v = mv_;
  const int d = day - 1;
  const double dx = v.dx[d];
  double sum_l = 0.0;
  const double* lar = y + v.offsets[kL];
  for (int i = 0; i < v.n_g * v.n_L; ++i) sum_l += lar[i];
  double mu_l = v.mu_L0;
  if (v.larval_per_wound && v.H > 0) mu_l += v.gamma_L * sum_l / v.H;
  const double mu_p = v.mu_T[d] + v.mu_W * (1.0 - v.moisture[d]);
  const double mu_f = v.mu_T[d] + v.mu_b_F;
  const double mu_m = v.mu_T[d] + v.mu_b_M;
  double best = std::max({v.n_E * v.dx_egg[d] / v.X_E + v.mu_E[d], v.n_L / v.tau_L + mu_l,
                          v.n_P * dx / v.X_P + mu_p});
  double max_iw_f = -1.0, max_iw_m = -1.0;
  for (double iw : v.inv_w_f) {
    if (iw > 0) max_iw_f = std::max(max_iw_f, iw);
  }
  for (double iw : v.inv_w_m) {
    if (iw > 0) max_iw_m = std::max(max_iw_m, iw);
  }
  if (max_iw_f > 0) best = std::max(best, dx * max_iw_f + mu_f);
  if (max_iw_m > 0) best = std::max(best, dx * max_iw_m + mu_m);
  return best;
}

}  // namespace sgdrive

// Internal test entry points: evaluate the kernel for one state.
// [[Rcpp::export]]
Rcpp::List sg_rates_cpp(const Rcpp::List& view, int day, const Rcpp::NumericVector& state) {
  const sgdrive::ModelView mv(view);
  if (state.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  sgdrive::RateKernel kernel(mv);
  Rcpp::NumericVector deriv(mv.n_states), counters(mv.n_counters);
  kernel.eval(day, state.begin(), deriv.begin(), counters.begin());
  return Rcpp::List::create(Rcpp::Named("deriv") = deriv, Rcpp::Named("counters") = counters);
}

// [[Rcpp::export]]
double sg_max_rate_cpp(const Rcpp::List& view, int day, const Rcpp::NumericVector& state) {
  const sgdrive::ModelView mv(view);
  if (state.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  sgdrive::RateKernel kernel(mv);
  return kernel.max_rate(day, state.begin());
}
