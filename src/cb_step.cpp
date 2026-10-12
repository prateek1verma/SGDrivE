// Chain-binomial step (see cb_step.h): construction, daily preparation of the
// chains, and the R entry points (one step, and the day loop).

#include "cb_step.h"

#include <limits>

namespace sgdrive {

CbStepper::CbStepper(const ModelView& mv, double dt)
    : mv_(mv),
      dt_(dt),
      kernel_(mv),
      adv_e_(mv.n_E),
      adv_l_(mv.n_L),
      adv_p_(mv.n_P),
      adv_a_(mv.n_A),
      ch_m_(mv.n_g),
      ch_f_(mv.n_g) {
  const int n_max = std::max({mv.n_E, mv.n_L, mv.n_P, mv.n_A, mv.n_g});
  ch_e_.resize(mv.n_E);
  ch_l_.resize(mv.n_L);
  ch_p_.resize(mv.n_P);
  for (int g = 0; g < mv.n_g; ++g) {
    ch_m_[g].resize(mv.n_A);
    ch_f_[g].resize(mv.n_A);
  }
  q_.assign(n_max, 0.0);
  p_.assign(static_cast<std::size_t>(n_max) * (n_max + 1), 0.0);
  pbar_.assign(p_.size(), 0.0);
  cur_.assign(mv.n_states, 0.0);
  nx_.assign(mv.n_states, 0.0);
  for (auto* w : {&eggs_, &out_e_, &out_l_, &out_p_, &em_m_, &em_f_, &ym_}) w->assign(mv.n_g, 0.0);
  matings_.assign(static_cast<std::size_t>(mv.n_g) * mv.n_g, 0.0);
  buf_.assign(n_max + 1, 0.0);
  rest_.assign(n_max + 1, 0.0);
}

// Death and destination probabilities, and the arrival distribution, of one
// chain (S16); see .sg_cb_chain() in R/cb-reference.R.
void CbStepper::prepare(PreparedChain& ch, ChainAdvance& adv, const double* q, int n, int last,
                        double mu) {
  const double dt = dt_;
  const int m = n + 1;
  adv.advance_average(q, last, mu, dt, p_.data(), pbar_.data());
  for (int i = 0; i < n; ++i) {
    const double* row = p_.data() + static_cast<std::size_t>(i) * m;
    double* dest = ch.dest.data() + static_cast<std::size_t>(i) * m;
    if (i >= last) {  // never reached: individuals placed here stay and only die
      std::fill(dest, dest + m, 0.0);
      dest[i] = 1.0;
      ch.die[i] = 1.0 - std::exp(-mu * dt);
      continue;
    }
    double keep = 0.0;
    for (int j = 0; j < m; ++j) keep += row[j];
    ch.die[i] = mu > 0.0 ? std::min(std::max(1.0 - keep, 0.0), 1.0) : 0.0;
    const double scale = keep > 0.0 ? keep : 1.0;
    for (int j = 0; j < m; ++j) dest[j] = row[j] / scale;
  }
  // arrivals enter sub-stage 1; leaving again within the step is kept in the
  // last reached sub-stage
  double first_sum = 0.0;
  for (int j = 0; j < m; ++j) first_sum += pbar_[j];
  for (int j = 0; j < n; ++j) ch.arr[j] = pbar_[j];
  if (last >= 1) ch.arr[last - 1] += pbar_[n];
  double arr_sum = 0.0;
  for (int j = 0; j < n; ++j) arr_sum += ch.arr[j];
  if (arr_sum > 0.0) {
    for (int j = 0; j < n; ++j) ch.arr[j] /= arr_sum;
  }
  ch.arr_die = mu > 0.0 ? std::min(std::max(1.0 - first_sum, 0.0), 1.0) : 0.0;
  ch.last = last;
}

void CbStepper::prepare_day(int day) {
  const ModelView& v = mv_;
  day_ = day;
  r_ = kernel_.day_rates(day);
  for (int s = 0; s < v.n_E; ++s) q_[s] = r_.q_e;
  prepare(ch_e_, adv_e_, q_.data(), v.n_E, v.n_E, r_.mu_e);
  for (int s = 0; s < v.n_L; ++s) q_[s] = r_.q_l;
  prepare(ch_l_, adv_l_, q_.data(), v.n_L, v.n_L, v.mu_L0);  // redone per step if crowding
  for (int s = 0; s < v.n_P; ++s) q_[s] = r_.q_p;
  prepare(ch_p_, adv_p_, q_.data(), v.n_P, v.n_P, r_.mu_p);
  for (int g = 0; g < v.n_g; ++g) {
    const double* iw_m = v.inv_w_m.data() + static_cast<std::size_t>(g) * v.n_A;
    const double* iw_f = v.inv_w_f.data() + static_cast<std::size_t>(g) * v.n_A;
    for (int s = 0; s < v.n_A; ++s) q_[s] = r_.dx * iw_m[s];
    prepare(ch_m_[g], adv_a_, q_.data(), v.n_A, v.last_m[g], r_.mu_m);
    for (int s = 0; s < v.n_A; ++s) q_[s] = r_.dx * iw_f[s];
    prepare(ch_f_[g], adv_a_, q_.data(), v.n_A, v.last_f[g], r_.mu_f);
  }
}

}  // namespace sgdrive

namespace {

int steps_per_day(double dt) {
  const double steps = std::round(1.0 / dt);
  if (!(dt > 0.0) || steps < 1.0 || std::fabs(steps * dt - 1.0) > 1e-9) {
    Rcpp::stop("`dt` must divide one day into a whole number of steps");
  }
  return static_cast<int>(steps);
}

bool is_whole(const Rcpp::NumericVector& x) {
  for (double v : x) {
    if (v < 0.0 || v != std::floor(v) || !std::isfinite(v)) return false;
  }
  return true;
}

}  // namespace

// One step from state x on simulation day `day` (tests of the step against
// the R reference): list(x, flows).
// [[Rcpp::export]]
Rcpp::List sg_cb_step_cpp(const Rcpp::List& view, int day, double dt,
                          const Rcpp::NumericVector& x, bool stochastic) {
  const sgdrive::ModelView mv(view);
  if (x.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  steps_per_day(dt);
  if (stochastic && !is_whole(x)) Rcpp::stop("stochastic steps need whole, non-negative counts");
  sgdrive::CbStepper stepper(mv, dt);
  stepper.prepare_day(day);
  Rcpp::NumericVector y = Rcpp::clone(x), flows(mv.n_counters);
  if (stochastic) {
    sgdrive::RandomDraws draws;
    stepper.step(y.begin(), flows.begin(), draws);
  } else {
    sgdrive::ExpectedDraws draws;
    stepper.step(y.begin(), flows.begin(), draws);
  }
  return Rcpp::List::create(Rcpp::Named("x") = y, Rcpp::Named("flows") = flows);
}

// Day loop of the chain-binomial solver for days 1..n_days.
//   x0: initial state (whole numbers when stochastic);
//   rel_count: release counts in the order of the model's release events
//     (rounded in R when stochastic);
//   record: logical, length n_days: store the state at the end of that day.
// Releases are added at the start of their day, then 1/dt steps follow with
// the day's forcing. The run stops early once every compartment is zero and
// no release is left (later days stay zero). Returns the states recorded,
// daily counters, the first day ending with every compartment zero
// (extinct_day, 0 if none), the first day ending with no individual of the
// wild-type genotype(s) (wild_extinct_day; mated females by their own
// genotype; 0 if none) and the number of days simulated.
// [[Rcpp::export]]
Rcpp::List sg_cb_run_cpp(const Rcpp::List& view, const Rcpp::NumericVector& x0,
                         const Rcpp::NumericVector& rel_count, double dt, bool stochastic,
                         int n_days, const Rcpp::LogicalVector& record) {
  const sgdrive::ModelView mv(view);
  const int n = mv.n_states, nc = mv.n_counters;
  if (x0.size() != n) Rcpp::stop("initial state has the wrong length");
  if (n_days < 1 || n_days > mv.n_days) {
    Rcpp::stop("n_days must be between 1 and the forcing length");
  }
  if (record.size() != n_days) Rcpp::stop("`record` must have one value per day");
  if (rel_count.size() != static_cast<R_xlen_t>(mv.rel_day.size())) {
    Rcpp::stop("`rel_count` must have one value per release event");
  }
  const int steps = steps_per_day(dt);
  if (stochastic && (!is_whole(x0) || !is_whole(rel_count))) {
    Rcpp::stop("stochastic runs need whole, non-negative initial counts and releases");
  }

  std::vector<std::vector<int>> releases_by_day(n_days);
  int last_release_day = 0;
  for (std::size_t j = 0; j < mv.rel_day.size(); ++j) {
    const int d = mv.rel_day[j];
    if (d <= n_days) {
      releases_by_day[d - 1].push_back(static_cast<int>(j));
      last_release_day = std::max(last_release_day, d);
    }
  }
  int n_record = 0;
  for (int d = 0; d < n_days; ++d) n_record += record[d] ? 1 : 0;

  Rcpp::NumericMatrix states(n_record, n);
  Rcpp::IntegerVector state_days(n_record);
  Rcpp::NumericMatrix counters(n_days, nc);
  std::vector<double> x(x0.begin(), x0.end()), flows(nc);
  sgdrive::CbStepper stepper(mv, dt);
  sgdrive::RandomDraws random_draws;
  sgdrive::ExpectedDraws expected_draws;
  int row = 0, extinct_day = 0, wild_extinct_day = 0, days_run = 0;
  // individuals of the wild-type genotype(s): blocks E..U by genotype, F by
  // the female's genotype
  auto wild_total = [&](const std::vector<double>& y) {
    double tot = 0.0;
    const int ng = mv.n_g;
    const int len[5] = {mv.n_E, mv.n_L, mv.n_P, mv.n_A, mv.n_A};
    for (int b = 0; b < 5; ++b) {
      for (int g = 0; g < ng; ++g) {
        if (!mv.wild[g]) continue;
        const double* p = y.data() + mv.offsets[b] + g * len[b];
        for (int s = 0; s < len[b]; ++s) tot += p[s];
      }
    }
    for (int g = 0; g < ng; ++g) {
      if (!mv.wild[g]) continue;
      const double* p = y.data() + mv.offsets[sgdrive::kF] + g * ng * mv.n_A;
      for (int s = 0; s < ng * mv.n_A; ++s) tot += p[s];
    }
    return tot;
  };

  for (int d = 1; d <= n_days; ++d) {
    Rcpp::checkUserInterrupt();
    for (int j : releases_by_day[d - 1]) x[mv.rel_index[j]] += rel_count[j];
    stepper.prepare_day(d);
    for (int k = 0; k < steps; ++k) {
      if (stochastic) {
        stepper.step(x.data(), flows.data(), random_draws);
      } else {
        stepper.step(x.data(), flows.data(), expected_draws);
      }
      for (int c = 0; c < nc; ++c) counters(d - 1, c) += flows[c];
    }
    days_run = d;
    if (record[d - 1]) {
      for (int i = 0; i < n; ++i) states(row, i) = x[i];
      state_days[row++] = d;
    }
    const bool empty = std::all_of(x.begin(), x.end(), [](double v) { return v == 0.0; });
    if (empty && extinct_day == 0) extinct_day = d;
    if (wild_extinct_day == 0 && (empty || wild_total(x) == 0.0)) wild_extinct_day = d;
    if (empty && d >= last_release_day) break;  // nothing can happen any more
  }
  // days after an early stop stay zero; give them their rows
  for (int d = days_run + 1; d <= n_days; ++d) {
    if (record[d - 1]) state_days[row++] = d;
  }
  return Rcpp::List::create(Rcpp::Named("states") = states, Rcpp::Named("state_days") = state_days,
                            Rcpp::Named("counters") = counters,
                            Rcpp::Named("extinct_day") = extinct_day,
                            Rcpp::Named("wild_extinct_day") = wild_extinct_day,
                            Rcpp::Named("days_run") = days_run);
}

// n_rep independent stochastic steps from the same state x (statistical tests
// of the step): a matrix with one row per replicate, the new state followed by
// the step's flows.
// [[Rcpp::export]]
Rcpp::NumericMatrix sg_cb_step_rep_cpp(const Rcpp::List& view, int day, double dt,
                                       const Rcpp::NumericVector& x, int n_rep) {
  const sgdrive::ModelView mv(view);
  if (x.size() != mv.n_states) Rcpp::stop("state has the wrong length");
  if (day < 1 || day > mv.n_days) Rcpp::stop("day %d is outside the forcing", day);
  if (n_rep < 1) Rcpp::stop("n_rep must be >= 1");
  steps_per_day(dt);
  if (!is_whole(x)) Rcpp::stop("stochastic steps need whole, non-negative counts");
  sgdrive::CbStepper stepper(mv, dt);
  stepper.prepare_day(day);
  sgdrive::RandomDraws draws;
  const int n = mv.n_states, nc = mv.n_counters;
  Rcpp::NumericMatrix out(n_rep, n + nc);
  std::vector<double> y(n), flows(nc);
  for (int r = 0; r < n_rep; ++r) {
    std::copy(x.begin(), x.end(), y.begin());
    stepper.step(y.data(), flows.data(), draws);
    for (int i = 0; i < n; ++i) out(r, i) = y[i];
    for (int c = 0; c < nc; ++c) out(r, n + c) = flows[c];
  }
  return out;
}
