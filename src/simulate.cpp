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

#include <Rcpp.h>

#include <algorithm>
#include <vector>

#include "dopri5.h"
#include "model_view.h"
#include "rates.h"

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
// [[Rcpp::export]]
Rcpp::List sg_ode_run_cpp(SEXP ptr, const Rcpp::NumericVector& x0, int day_from, int day_to,
                          double rtol, double atol, const Rcpp::LogicalVector& record,
                          bool stop_on_stiff, double max_steps) {
  screwdrive::KernelHandle* h = handle(ptr);
  const screwdrive::ModelView& mv = h->mv;
  const int n = mv.n_states, nc = mv.n_counters;
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
  std::vector<double> st_states, st_counters, st_max_rate, st_min_state;
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
    st_counters.insert(st_counters.end(), y.begin() + n, y.end());
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
      st_states.insert(st_states.end(), x.begin(), x.end());
      st_days.push_back(d);
    }
    completed = d;
  }

  const int n_done = static_cast<int>(st_steps.size());
  const int n_rec = static_cast<int>(st_days.size());
  Rcpp::NumericMatrix states(n_rec, n);
  for (int r = 0; r < n_rec; ++r) {
    for (int i = 0; i < n; ++i) states(r, i) = st_states[static_cast<std::size_t>(r) * n + i];
  }
  Rcpp::NumericMatrix counters(n_done, nc);
  for (int r = 0; r < n_done; ++r) {
    for (int c = 0; c < nc; ++c) counters(r, c) = st_counters[static_cast<std::size_t>(r) * nc + c];
  }
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
      Rcpp::Named("limiting_node") = Rcpp::wrap(st_limiting));
}
