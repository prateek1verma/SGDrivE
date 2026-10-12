// Test entry point for the exact chain advancement (chains.h).

#include <Rcpp.h>

#include "chains.h"

// Advancement probabilities over dt for rates q (classes >= last never
// reached): a matrix with length(q) rows and length(q) + 1 columns (last
// column: left the chain).
// [[Rcpp::export]]
Rcpp::NumericMatrix sg_chain_advance_cpp(const Rcpp::NumericVector& q, int last, double dt) {
  const int n = q.size();
  if (n < 1) Rcpp::stop("the chain needs at least one class");
  if (last < 0 || last > n) Rcpp::stop("`last` must be between 0 and the number of classes");
  if (!(dt >= 0.0)) Rcpp::stop("`dt` must be >= 0");
  for (int i = 0; i < n; ++i) {
    if (!(q[i] >= 0.0) || !std::isfinite(q[i])) Rcpp::stop("rates must be finite and >= 0");
  }
  sgdrive::ChainAdvance chain(n);
  std::vector<double> p(static_cast<std::size_t>(n) * (n + 1));
  chain.advance(q.begin(), last, dt, p.data());
  Rcpp::NumericMatrix out(n, n + 1);
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j <= n; ++j) out(i, j) = p[static_cast<std::size_t>(i) * (n + 1) + j];
  }
  return out;
}
