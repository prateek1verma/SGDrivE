// Construction and validation of the C++ model view (see model_view.h).

#include "model_view.h"

#include <cmath>
#include <string>

namespace sgdrive {

namespace {

// Typed field accessors with informative errors.
SEXP field(const Rcpp::List& view, const char* name) {
  if (!view.containsElementNamed(name)) {
    Rcpp::stop("model view: missing field '%s'", name);
  }
  return view[name];
}

int get_int(const Rcpp::List& view, const char* name) {
  Rcpp::IntegerVector x(field(view, name));
  if (x.size() != 1 || x[0] == NA_INTEGER) {
    Rcpp::stop("model view: field '%s' must be a single integer", name);
  }
  return x[0];
}

double get_double(const Rcpp::List& view, const char* name) {
  Rcpp::NumericVector x(field(view, name));
  if (x.size() != 1 || !std::isfinite(x[0])) {
    Rcpp::stop("model view: field '%s' must be a single finite number", name);
  }
  return x[0];
}

std::vector<int> get_ints(const Rcpp::List& view, const char* name, R_xlen_t n) {
  Rcpp::IntegerVector x(field(view, name));
  if (n >= 0 && x.size() != n) {
    Rcpp::stop("model view: field '%s' has length %d, expected %d", name,
               static_cast<int>(x.size()), static_cast<int>(n));
  }
  return Rcpp::as<std::vector<int>>(x);
}

std::vector<double> get_doubles(const Rcpp::List& view, const char* name, R_xlen_t n) {
  Rcpp::NumericVector x(field(view, name));
  if (n >= 0 && x.size() != n) {
    Rcpp::stop("model view: field '%s' has length %d, expected %d", name,
               static_cast<int>(x.size()), static_cast<int>(n));
  }
  for (R_xlen_t i = 0; i < x.size(); ++i) {
    if (!std::isfinite(x[i])) {
      Rcpp::stop("model view: field '%s' contains a non-finite value", name);
    }
  }
  return Rcpp::as<std::vector<double>>(x);
}

void check_flags(const std::vector<int>& x, const char* name) {
  for (int v : x) {
    if (v != 0 && v != 1) Rcpp::stop("model view: field '%s' must contain only 0/1", name);
  }
}

}  // namespace

ModelView::ModelView(const Rcpp::List& view) {
  n_g = get_int(view, "n_g");
  n_E = get_int(view, "n_E");
  n_L = get_int(view, "n_L");
  n_P = get_int(view, "n_P");
  n_A = get_int(view, "n_A");
  n_states = get_int(view, "n_states");
  n_days = get_int(view, "n_days");
  n_counters = get_int(view, "n_counters");
  if (n_g < 1 || n_E < 1 || n_L < 1 || n_P < 1 || n_A < 1 || n_days < 1) {
    Rcpp::stop("model view: dimensions must be >= 1");
  }
  host_model = get_int(view, "host_model") == 1;
  n_I = get_int(view, "n_I");
  if (host_model && (n_I < 1 || n_I > 10)) {
    Rcpp::stop("model view: n_I = %d must be between 1 and 10", n_I);
  }
  if (!host_model && n_I != 0) Rcpp::stop("model view: n_I must be 0 without hosts");
  const int fly_states = n_g * (n_E + n_L + n_P + 2 * n_A) + n_g * n_g * n_A;
  const int expected_states = fly_states + (host_model ? 2 + n_I : 0);
  if (n_states != expected_states) {
    Rcpp::stop("model view: n_states = %d, expected %d", n_states, expected_states);
  }
  // counters: eggs_laid[k] (n_g), 6 flow totals, matings[g|m] (n_g^2),
  // 6 deaths by block, 3 deaths by old age, and 6 host flows with the
  // host-wound model (see .sg_counter_names())
  if (n_counters != n_g + n_g * n_g + 16 + (host_model ? 6 : 0)) {
    Rcpp::stop("model view: n_counters = %d does not match n_g = %d", n_counters, n_g);
  }

  offsets = get_ints(view, "offsets", 6);
  const int block_len[6] = {n_g * n_E, n_g * n_L, n_g * n_P, n_g * n_A, n_g * n_A,
                            n_g * n_g * n_A};
  int expected_offset = 0;
  for (int b = 0; b < 6; ++b) {
    if (offsets[b] != expected_offset) {
      Rcpp::stop("model view: block %d starts at %d, expected %d", b, offsets[b],
                 expected_offset);
    }
    expected_offset += block_len[b];
  }
  host_offsets = get_ints(view, "host_offsets", 3);
  const int host_expected[3] = {fly_states, fly_states + 1, fly_states + 2};
  for (int b = 0; b < 3; ++b) {
    if (host_offsets[b] != (host_model ? host_expected[b] : -1)) {
      Rcpp::stop("model view: host block %d starts at %d", b, host_offsets[b]);
    }
  }
  N_H = get_double(view, "N_H");
  host_wounding = get_double(view, "host_wounding");
  host_healing = get_double(view, "host_healing");
  host_advance = get_double(view, "host_advance");
  host_treatment = get_double(view, "host_treatment");
  host_attractiveness = get_double(view, "host_attractiveness");
  if (N_H < 0 || host_wounding < 0 || host_healing < 0 || host_advance < 0 ||
      host_treatment < 0 || host_attractiveness < 0) {
    Rcpp::stop("model view: host count and rates must be >= 0");
  }
  if (host_model && (N_H <= 0 || host_healing <= 0 || host_advance <= 0)) {
    Rcpp::stop("model view: host count, healing and infestation rates must be > 0");
  }

  mating_instant = get_int(view, "mating_instant") == 1;
  search_all_adults = get_int(view, "search_all_adults") == 1;
  larval_per_wound = get_int(view, "larval_per_wound") == 1;
  larval_per_infested_host = get_int(view, "larval_per_infested_host") == 1;
  if (larval_per_infested_host && !host_model) {
    Rcpp::stop("model view: crowding per infested host needs the host model");
  }
  oviposition_batches = get_int(view, "oviposition_batches") == 1;

  X_E = get_double(view, "X_E");
  tau_L = get_double(view, "tau_L");
  X_P = get_double(view, "X_P");
  R = get_double(view, "R");
  a_s = get_double(view, "a_s");
  sigma0 = get_double(view, "sigma0");
  nu = get_double(view, "nu");
  mu_L0 = get_double(view, "mu_L0");
  gamma_L = get_double(view, "gamma_L");
  mu_W = get_double(view, "mu_W");
  mu_b_F = get_double(view, "mu_b_F");
  mu_b_M = get_double(view, "mu_b_M");
  H = get_double(view, "H");
  h_M = get_double(view, "h_M");
  h_M_param = get_double(view, "h_M_param");
  clutch_size = get_double(view, "clutch_size");
  if (clutch_size < 0) Rcpp::stop("model view: clutch_size must be >= 0");

  inv_w_f = get_doubles(view, "inv_w_f", static_cast<R_xlen_t>(n_A) * n_g);
  inv_w_m = get_doubles(view, "inv_w_m", static_cast<R_xlen_t>(n_A) * n_g);
  last_f = get_ints(view, "last_f", n_g);
  last_m = get_ints(view, "last_m", n_g);
  for (int g = 0; g < n_g; ++g) {
    if (last_f[g] < 1 || last_f[g] > n_A || last_m[g] < 1 || last_m[g] > n_A) {
      Rcpp::stop("model view: last reached adult class out of range for genotype %d", g + 1);
    }
  }
  laying = get_ints(view, "laying", n_A);
  receptive = get_ints(view, "receptive", n_A);
  mature = get_ints(view, "mature", n_A);
  check_flags(laying, "laying");
  check_flags(receptive, "receptive");
  check_flags(mature, "mature");

  egg_tensor = get_doubles(view, "egg_tensor", static_cast<R_xlen_t>(n_g) * n_g * n_g);
  search_weight = get_doubles(view, "search_weight", n_g);
  batch_viability = get_doubles(view, "batch_viability", static_cast<R_xlen_t>(n_g) * n_g);
  eta = get_doubles(view, "eta", static_cast<R_xlen_t>(n_g) * n_g);
  phi = get_doubles(view, "phi", n_g);
  xiF = get_doubles(view, "xiF", n_g);
  xiM = get_doubles(view, "xiM", n_g);
  wild = get_ints(view, "wild", n_g);
  check_flags(wild, "wild");

  dx = get_doubles(view, "dx", n_days);
  dx_egg = get_doubles(view, "dx_egg", n_days);
  mu_T = get_doubles(view, "mu_T", n_days);
  mu_E = get_doubles(view, "mu_E", n_days);
  phi_ov = get_doubles(view, "phi_ov", n_days);
  moisture = get_doubles(view, "moisture", n_days);

  rel_day = get_ints(view, "rel_day", -1);
  rel_index = get_ints(view, "rel_index", static_cast<R_xlen_t>(rel_day.size()));
  rel_count = get_doubles(view, "rel_count", static_cast<R_xlen_t>(rel_day.size()));
  for (std::size_t j = 0; j < rel_day.size(); ++j) {
    if (rel_day[j] < 1 || rel_day[j] > n_days || rel_index[j] < 0 ||
        rel_index[j] >= n_states || rel_count[j] < 0) {
      Rcpp::stop("model view: release %d is out of range", static_cast<int>(j) + 1);
    }
  }
}

}  // namespace sgdrive

// Build a model view from its R description and report what was read
// (internal; used by tests to check that R and C++ agree on the layout).
// [[Rcpp::export]]
Rcpp::List sg_model_view_info_cpp(const Rcpp::List& view) {
  const sgdrive::ModelView mv(view);
  double sum_inv_w_f = 0.0, sum_egg = 0.0, sum_dx = 0.0;
  for (double v : mv.inv_w_f) sum_inv_w_f += v;
  for (double v : mv.egg_tensor) sum_egg += v;
  for (double v : mv.dx) sum_dx += v;
  return Rcpp::List::create(
      Rcpp::Named("dims") = Rcpp::IntegerVector::create(mv.n_g, mv.n_E, mv.n_L, mv.n_P,
                                                        mv.n_A, mv.n_states, mv.n_days,
                                                        mv.n_counters),
      Rcpp::Named("offsets") = mv.offsets,
      Rcpp::Named("options") = Rcpp::LogicalVector::create(
          mv.mating_instant, mv.search_all_adults, mv.larval_per_wound),
      Rcpp::Named("last_f") = mv.last_f, Rcpp::Named("last_m") = mv.last_m,
      Rcpp::Named("n_releases") = static_cast<int>(mv.rel_day.size()),
      Rcpp::Named("sum_inv_w_f") = sum_inv_w_f, Rcpp::Named("sum_egg_tensor") = sum_egg,
      Rcpp::Named("sum_dx") = sum_dx,
      Rcpp::Named("idx_F_last") = mv.idx_f(mv.n_g - 1, mv.n_g - 1, mv.n_A - 1),
      Rcpp::Named("idx_M_first") = mv.idx(sgdrive::kM, 0, 0),
      Rcpp::Named("host_model") = mv.host_model, Rcpp::Named("n_I") = mv.n_I,
      Rcpp::Named("host_offsets") = mv.host_offsets,
      Rcpp::Named("host_rates") = Rcpp::NumericVector::create(
          mv.N_H, mv.host_wounding, mv.host_healing, mv.host_advance, mv.host_treatment,
          mv.host_attractiveness));
}
