// Recorder of the simulation drivers (multi-node work phase 1 WP9, phase 2
// WP5; D-032): what a run stores. On output days weighted sums y = W_x x of
// the state (and the full state only if asked for); every day weighted sums
// W_c c of the counters (and the full counters only if asked for), added
// into the row of the day's period (day, week, month or year; phase 3 WP4).
// W_x and W_c come from R in compressed sparse row form (R/record.R). Used by the
// ODE day loop (simulate.cpp) and the chain-binomial day loop (cb_step.cpp);
// recording uses no random numbers.

#ifndef SCREWDRIVE_RECORDER_H
#define SCREWDRIVE_RECORDER_H

#include <Rcpp.h>

#include <cmath>
#include <vector>

namespace screwdrive {

// A sparse matrix in compressed sparse row form: row r has the entries
// row_ptr[r] .. row_ptr[r + 1] - 1 of col (0-based columns) and w.
struct Csr {
  std::vector<int> row_ptr, col;
  std::vector<double> w;
  int rows() const { return static_cast<int>(row_ptr.size()) - 1; }
  // append the rows of W v to out
  void multiply_append(const double* v, std::vector<double>& out) const {
    for (int r = 0; r < rows(); ++r) {
      double s = 0.0;
      for (int k = row_ptr[r]; k < row_ptr[r + 1]; ++k) s += w[k] * v[col[k]];
      out.push_back(s);
    }
  }
};

inline Csr read_csr(const Rcpp::List& rec, const char* name, int n_cols) {
  const Rcpp::List m = rec[name];
  Csr out;
  out.row_ptr = Rcpp::as<std::vector<int>>(m["row_ptr"]);
  out.col = Rcpp::as<std::vector<int>>(m["col"]);
  out.w = Rcpp::as<std::vector<double>>(m["w"]);
  if (out.row_ptr.empty() || out.row_ptr[0] != 0 ||
      out.row_ptr.back() != static_cast<int>(out.col.size()) || out.w.size() != out.col.size()) {
    Rcpp::stop("recorder: %s is not a valid CSR matrix", name);
  }
  for (std::size_t r = 1; r < out.row_ptr.size(); ++r) {
    if (out.row_ptr[r] < out.row_ptr[r - 1]) Rcpp::stop("recorder: %s row pointers decrease", name);
  }
  for (std::size_t k = 0; k < out.col.size(); ++k) {
    if (out.col[k] < 0 || out.col[k] >= n_cols || !std::isfinite(out.w[k])) {
      Rcpp::stop("recorder: %s has a column or weight out of range", name);
    }
  }
  return out;
}

// What a run stores: full states, weighted state sums W_x, full counters,
// weighted counter sums W_c by period. Without a recorder (NULL): full
// states and counters. `c_period` (optional, 0-based, one per day of the
// forcing, non-decreasing) gives each day's period; without it every day is
// its own period.
struct Recorder {
  bool state = true, counters = true;
  Csr wx, wc;
  std::vector<int> c_period;
  Recorder(SEXP rec, int n_states, int n_counters) {
    wx.row_ptr = {0};
    wc.row_ptr = {0};
    if (Rf_isNull(rec)) return;
    const Rcpp::List r(rec);
    state = Rcpp::as<bool>(r["state"]);
    counters = Rcpp::as<bool>(r["counters"]);
    wx = read_csr(r, "wx", n_states);
    wc = read_csr(r, "wc", n_counters);
    if (r.containsElementNamed("c_period") && !Rf_isNull(r["c_period"])) {
      c_period = Rcpp::as<std::vector<int>>(r["c_period"]);
      for (std::size_t k = 0; k < c_period.size(); ++k) {
        if (c_period[k] < 0 || (k > 0 && c_period[k] < c_period[k - 1])) {
          Rcpp::stop("recorder: c_period must be non-decreasing and >= 0");
        }
      }
    }
  }
  // Add the counter sums W_c c of day `day` (1-based) into the row of its
  // period: a new row when the period changes (its 0-based period in
  // `periods`), else the last row.
  void add_sums(int day, const double* c, std::vector<double>& sums,
                std::vector<int>& periods) const {
    const int p = c_period.empty() ? day - 1 : c_period.at(static_cast<std::size_t>(day - 1));
    const int nr = wc.rows();
    if (periods.empty() || periods.back() != p) {
      wc.multiply_append(c, sums);
      periods.push_back(p);
      return;
    }
    double* row = sums.data() + sums.size() - static_cast<std::size_t>(nr);
    for (int r = 0; r < nr; ++r) {
      double s = 0.0;
      for (int k = wc.row_ptr[r]; k < wc.row_ptr[r + 1]; ++k) s += wc.w[k] * c[wc.col[k]];
      row[r] += s;
    }
  }
};

// Row-major growable storage to a column-major R matrix.
inline Rcpp::NumericMatrix as_matrix(const std::vector<double>& v, int rows, int cols) {
  Rcpp::NumericMatrix m(rows, cols);
  for (int r = 0; r < rows; ++r) {
    for (int c = 0; c < cols; ++c) m(r, c) = v[static_cast<std::size_t>(r) * cols + c];
  }
  return m;
}

}  // namespace screwdrive

#endif  // SCREWDRIVE_RECORDER_H
