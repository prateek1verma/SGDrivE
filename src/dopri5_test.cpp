// Test harness for the Dormand-Prince integrator (internal): integrates small
// problems with known behaviour so tests can check accuracy, convergence order
// and stiffness detection independently of the screwworm model.

#include <Rcpp.h>

#include <string>
#include <vector>

#include "dopri5.h"

// [[Rcpp::export]]
Rcpp::List sg_dopri5_test_cpp(std::string problem, Rcpp::NumericVector y0, double t0,
                              double t1, double rtol, double atol, double fixed_h,
                              double max_steps, double par) {
  const int n = y0.size();
  std::vector<double> y(y0.begin(), y0.end());
  screwdrive::Dopri5Options opt;
  opt.rtol = rtol;
  opt.atol = atol;
  opt.fixed_h = fixed_h;
  opt.max_steps = static_cast<long>(max_steps);

  auto run = [&](auto f) {
    screwdrive::Dopri5 solver(n);
    return solver.integrate(f, t0, t1, y.data(), opt);
  };
  screwdrive::Dopri5Stats st;
  if (problem == "decay") {  // y' = -par y
    st = run([&](double, const double* x, double* dx) {
      for (int i = 0; i < n; ++i) dx[i] = -par * x[i];
    });
  } else if (problem == "oscillator") {  // y0' = y1, y1' = -par^2 y0
    if (n != 2) Rcpp::stop("oscillator needs 2 states");
    st = run([&](double, const double* x, double* dx) {
      dx[0] = x[1];
      dx[1] = -par * par * x[0];
    });
  } else if (problem == "erlang") {  // chain with rate par in every compartment
    st = run([&](double, const double* x, double* dx) {
      dx[0] = -par * x[0];
      for (int i = 1; i < n; ++i) dx[i] = par * (x[i - 1] - x[i]);
    });
  } else if (problem == "prothero_robinson") {  // y' = -par (y - cos t) - sin t
    st = run([&](double t, const double* x, double* dx) {
      dx[0] = -par * (x[0] - std::cos(t)) - std::sin(t);
    });
  } else {
    Rcpp::stop("unknown test problem '%s'", problem);
  }
  return Rcpp::List::create(
      Rcpp::Named("y") = Rcpp::NumericVector(y.begin(), y.end()),
      Rcpp::Named("steps") = static_cast<double>(st.steps),
      Rcpp::Named("accepted") = static_cast<double>(st.accepted),
      Rcpp::Named("rejected") = static_cast<double>(st.rejected),
      Rcpp::Named("evaluations") = static_cast<double>(st.evaluations),
      Rcpp::Named("stiff") = st.stiff, Rcpp::Named("status") = st.status,
      Rcpp::Named("h_last") = st.h_last);
}
