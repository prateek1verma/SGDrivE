// C++ simulation driver for the deterministic model (D-021).
//
// A KernelHandle (model view + rate kernel + releases indexed by day) is
// created once per simulation and kept alive as an R external pointer, so
// R can drive it repeatedly: the C++ day loop (sg_ode_run_cpp), single rate
// evaluations for deSolve (sg_kernel_rates_cpp) and the largest-rate
// diagnostic (sg_kernel_max_rate_cpp).
//
// Day loop (as in the R backend, D-019): for each day, add the day's releases
// to the state, integrate from day - 1 to day with the day's forcing held
// constant (Dormand-Prince 5(4), dopri5.h), integrating the flow counters
// alongside the states (reset each day), then record. If a day fails or is
// detected as stiff (when stop_on_stiff is set), the loop stops and returns
// the state at the start of that day, after its releases, so R can integrate
// that day with a stiff solver and resume.
//
// Node extinction rule (MODEL_SPEC N5, D-032; optional): at the end of each
// day, every node whose flies (E, L, P, M, U, F) total less than the
// extinction threshold has its fly states set to zero; with the host-wound
// model its infested hosts return to HS and are counted in the node's
// infestations_ended counter. The removed flies by stage and node are
// returned as extra counters after the integrated ones.
//
// Recorder (phase 1 WP9; optional; recorder.h): on output days, weighted
// sums y = W_x x of the state are stored (and the full state only if asked
// for); every day, weighted sums W_c c of the counters (and the full counters
// only if asked for). W_x and W_c come from R in CSR form (R/record.R).

#include <Rcpp.h>

#include <algorithm>
#include <vector>

#include "dopri5.h"
#include "model_view.h"
#include "rates.h"
#include "recorder.h"

namespace screwdrive {

struct KernelHandle {
  ModelView mv;
  RateKernel kernel;
  std::vector<std::vector<int>> releases_by_day;  // [n_days][j]: release rows on that day

  explicit KernelHandle(const Rcpp::List& view) : mv(view), kernel(mv), releases_by_day(mv.n_days) {
    for (std::size_t j = 0; j < mv.rel_day.size(); ++j) {
      releases_by_day[mv.rel_day[j] - 1].push_back(static_cast<int>(j));
    }
  }
};

}  // namespace screwdrive

namespace {

// (N5) Apply the node extinction rule to state x at the end of a day; add the
// removed flies (6 per node: E, L, P, M, U, F) to removed and the infested
// hosts returned to HS to the node's infestations_ended counter in counters
// (index ended_index within a node's counters; -1 without hosts).
void node_extinction(const screwdrive::ModelView& mv, double threshold, int ended_index,
                     double* x, double* counters, double* removed) {
  const int ng = mv.n_g;
  const int len[6] = {ng * mv.n_E, ng * mv.n_L, ng * mv.n_P, ng * mv.n_A, ng * mv.n_A,
                      ng * ng * mv.n_A};
  for (int i = 0; i < mv.n_node; ++i) {
    double* xi = x + static_cast<std::size_t>(i) * mv.node_stride;
    double block_sum[6];
    double total = 0.0;
    for (int b = 0; b < 6; ++b) {
      block_sum[b] = 0.0;
      const double* xb = xi + mv.offsets[b];
      for (int k = 0; k < len[b]; ++k) block_sum[b] += xb[k];
      total += block_sum[b];
    }
    if (!(total < threshold)) continue;
    for (int b = 0; b < 6; ++b) {
      std::fill(xi + mv.offsets[b], xi + mv.offsets[b] + len[b], 0.0);
      removed[i * 6 + b] += block_sum[b];
    }
    if (mv.host_model) {
      double* hi = xi + mv.host_offsets[2];
      double infested = 0.0;
      for (int j = 0; j < mv.n_I; ++j) {
        infested += hi[j];
        hi[j] = 0.0;
      }
      xi[mv.host_offsets[0]] += infested;
      counters[static_cast<std::size_t>(i) * mv.n_counters_node + ended_index] += infested;
    }
  }
}

screwdrive::KernelHandle* handle(SEXP ptr) {
  Rcpp::XPtr<screwdrive::KernelHandle> h(ptr);
  if (h.get() == nullptr) Rcpp::stop("invalid kernel handle (create it with sg_kernel_new_cpp)");
  return h.get();
}

void check_day(const screwdrive::ModelView& mv, int day) {
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
}

}  // namespace

// Create the kernel handle for a model view.
// [[Rcpp::export]]
SEXP sg_kernel_new_cpp(const Rcpp::List& view) {
  Rcpp::XPtr<screwdrive::KernelHandle> ptr(new screwdrive::KernelHandle(view), true);
  return ptr;
}

// Derivatives and counter rates for deSolve: c(dydt, counters), with y of
// length n_states + n_counters (the counter part of y is not used).
// [[Rcpp::export]]
Rcpp::NumericVector sg_kernel_rates_cpp(SEXP ptr, int day, const Rcpp::NumericVector& y) {
  screwdrive::KernelHandle* h = handle(ptr);
  const screwdrive::ModelView& mv = h->mv;
  check_day(mv, day);
  if (y.size() < mv.n_states) Rcpp::stop("state has the wrong length");
  Rcpp::NumericVector out(mv.n_states + mv.n_counters);
  h->kernel.eval(day, y.begin(), out.begin(), out.begin() + mv.n_states);
  return out;
}

// [[Rcpp::export]]
double sg_kernel_max_rate_cpp(SEXP ptr, int day, const Rcpp::NumericVector& x) {
  screwdrive::KernelHandle* h = handle(ptr);
  check_day(h->mv, day);
  if (x.size() < h->mv.n_states) Rcpp::stop("state has the wrong length");
  return h->kernel.max_rate(day, x.begin());
}

// Run days day_from..day_to (1-based) with the Dormand-Prince integrator.
//   x0: state at the start of day_from, before its releases;
//   record: logical, length n_days: store the state at the end of that day.
// Returns the states recorded (with their days), daily counters and
// diagnostics for the days completed (day_from to completed), the state after
// the last completed day (or, if the loop stopped on day d, the state at the
// start of day d after its releases) and stopped_day (0 if all days were
// completed) with stop_status (1 too many steps, 2 step size too small,
// 3 stiff, 4 non-finite state). With more than one node also limiting_node:
// for each completed day the node (1-based) whose state most often had the
// largest scaled error in the accepted Dormand-Prince steps (0 if none).
// With extinction = true the node extinction rule (N5) is applied at the end
// of each day, with threshold `ext_threshold`; the counters then have 6 more
// columns per node (removed E, L, P, M, U, F), and ended_index is the
// position of infestations_ended within a node's counters (-1 without hosts).
// `recorder` (NULL: full states and counters) is a list with logical `state`
// and `counters` and CSR matrices `wx` [sums x states] and `wc` [sums x
// counters, including the removed ones]; the result then also has
// `recorded` (W_x x on recorded days) and `counter_sums` (W_c c per day), and
// `states` or `counters` with no columns when they are not kept.
// [[Rcpp::export]]
Rcpp::List sg_ode_run_cpp(SEXP ptr, const Rcpp::NumericVector& x0, int day_from, int day_to,
                          double rtol, double atol, const Rcpp::LogicalVector& record,
                          bool stop_on_stiff, double max_steps, bool extinction = false,
                          double ext_threshold = 1.0, int ended_index = -1,
                          SEXP recorder = R_NilValue) {
  screwdrive::KernelHandle* h = handle(ptr);
  const screwdrive::ModelView& mv = h->mv;
  const int n = mv.n_states, nc = mv.n_counters;
  if (extinction && mv.host_model && (ended_index < 0 || ended_index >= mv.n_counters_node)) {
    Rcpp::stop("ended_index is out of range");
  }
  const int n_removed = extinction ? 6 * mv.n_node : 0;
  std::vector<double> removed(n_removed);
  const screwdrive::Recorder rec(recorder, n, nc + n_removed);
  std::vector<double> counter_row(static_cast<std::size_t>(nc) + n_removed);
  if (x0.size() != n) Rcpp::stop("initial state has the wrong length");
  if (record.size() != mv.n_days) Rcpp::stop("`record` must have one value per day");
  check_day(mv, day_from);
  check_day(mv, day_to);
  if (day_to < day_from) Rcpp::stop("day_to is before day_from");

  screwdrive::Dopri5Options opt;
  opt.rtol = rtol;
  opt.atol = atol;
  opt.max_steps = static_cast<long>(max_steps);
  opt.detect_stiff = true;
  opt.stop_when_stiff = stop_on_stiff;

  // results are accumulated in growable vectors and copied once at the end, so
  // a run that stops early (stiff fallback) costs no more than the days it did
  std::vector<double> st_states, st_counters, st_max_rate, st_min_state, st_recorded, st_sums;
  std::vector<int> st_days, st_steps, st_rejected, st_evals, st_stiff, st_limiting;
  // limiting node (D-032): counts per node of the steps it limited
  const bool track_nodes = mv.n_node > 1;
  std::vector<int> node_counts(mv.n_node, 0);
  if (track_nodes) {
    opt.argmax_counts = node_counts.data();
    opt.argmax_group = mv.node_stride;
    opt.argmax_limit = n;  // states only, not the counters integrated alongside
  }

  std::vector<double> x(x0.begin(), x0.end());
  std::vector<double> y(static_cast<std::size_t>(n) + nc);
  screwdrive::Dopri5 solver(n + nc);
  int stopped_day = 0, stop_status = 0, completed = 0;

  for (int d = day_from; d <= day_to; ++d) {
    Rcpp::checkUserInterrupt();
    // releases: instantaneous jumps at the start of the day
    for (int j : h->releases_by_day[d - 1]) x[mv.rel_index[j]] += mv.rel_count[j];

    std::copy(x.begin(), x.end(), y.begin());
    std::fill(y.begin() + n, y.end(), 0.0);
    const double mr = h->kernel.max_rate(d, x.data());
    std::fill(node_counts.begin(), node_counts.end(), 0);

    auto rhs = [&](double, const double* yy, double* dy) {
      h->kernel.eval(d, yy, dy, dy + n);
    };
    const screwdrive::Dopri5Stats st = solver.integrate(rhs, d - 1.0, static_cast<double>(d),
                                                     y.data(), opt);
    // a solution that is not finite stops the run (status 4; defensive: the
    // error control rejects such steps)
    const bool finite = std::all_of(y.begin(), y.end(), [](double v) { return std::isfinite(v); });
    if (st.status != 0 || !finite) {
      stopped_day = d;  // x holds the state at the start of day d, after releases
      stop_status = st.status != 0 ? st.status : 4;
      break;
    }
    std::copy(y.begin(), y.begin() + n, x.begin());
    if (extinction) {
      std::fill(removed.begin(), removed.end(), 0.0);
      node_extinction(mv, ext_threshold, ended_index, x.data(), y.data() + n, removed.data());
    }
    std::copy(y.begin() + n, y.end(), counter_row.begin());
    std::copy(removed.begin(), removed.end(), counter_row.begin() + nc);
    if (rec.counters) st_counters.insert(st_counters.end(), counter_row.begin(), counter_row.end());
    rec.wc.multiply_append(counter_row.data(), st_sums);
    st_steps.push_back(static_cast<int>(st.steps));
    st_rejected.push_back(static_cast<int>(st.rejected));
    st_evals.push_back(static_cast<int>(st.evaluations));
    st_stiff.push_back(st.stiff ? 1 : 0);
    st_max_rate.push_back(mr);
    st_min_state.push_back(*std::min_element(x.begin(), x.end()));
    if (track_nodes) {
      const auto top = std::max_element(node_counts.begin(), node_counts.end());
      st_limiting.push_back(*top > 0 ? static_cast<int>(top - node_counts.begin()) + 1 : 0);
    }
    if (record[d - 1]) {
      if (rec.state) st_states.insert(st_states.end(), x.begin(), x.end());
      rec.wx.multiply_append(x.data(), st_recorded);
      st_days.push_back(d);
    }
    completed = d;
  }

  const int n_done = static_cast<int>(st_steps.size());
  const int n_rec = static_cast<int>(st_days.size());
  using screwdrive::as_matrix;
  const int nc_all = nc + n_removed;
  Rcpp::NumericMatrix states = as_matrix(st_states, n_rec, rec.state ? n : 0);
  Rcpp::NumericMatrix counters = as_matrix(st_counters, n_done, rec.counters ? nc_all : 0);
  Rcpp::NumericMatrix recorded = as_matrix(st_recorded, n_rec, rec.wx.rows());
  Rcpp::NumericMatrix counter_sums = as_matrix(st_sums, n_done, rec.wc.rows());
  Rcpp::LogicalVector stiff(n_done);
  for (int r = 0; r < n_done; ++r) stiff[r] = st_stiff[r];

  return Rcpp::List::create(
      Rcpp::Named("x") = Rcpp::NumericVector(x.begin(), x.end()),
      Rcpp::Named("completed") = completed, Rcpp::Named("stopped_day") = stopped_day,
      Rcpp::Named("stop_status") = stop_status,
      Rcpp::Named("states") = states, Rcpp::Named("state_days") = Rcpp::wrap(st_days),
      Rcpp::Named("counters") = counters, Rcpp::Named("steps") = Rcpp::wrap(st_steps),
      Rcpp::Named("rejected") = Rcpp::wrap(st_rejected),
      Rcpp::Named("evaluations") = Rcpp::wrap(st_evals), Rcpp::Named("stiff") = stiff,
      Rcpp::Named("max_rate") = Rcpp::wrap(st_max_rate),
      Rcpp::Named("min_state") = Rcpp::wrap(st_min_state),
      Rcpp::Named("limiting_node") = Rcpp::wrap(st_limiting),
      Rcpp::Named("recorded") = recorded, Rcpp::Named("counter_sums") = counter_sums);
}
