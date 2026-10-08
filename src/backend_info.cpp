// Toolchain check for the SGDrivE C++ backend (milestone M0).
//
// Contains no model code. It only reports how the shared library was
// compiled, so tests and continuous integration can confirm that the C++
// toolchain works on every platform before the rate kernel is written.

#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::List sg_backend_info_cpp() {
#if defined(__clang__)
  const std::string compiler = std::string("clang ") + __clang_version__;
#elif defined(__GNUC__)
  const std::string compiler = std::string("gcc ") + __VERSION__;
#elif defined(_MSC_VER)
  const std::string compiler = "msvc " + std::to_string(_MSC_VER);
#else
  const std::string compiler = "unknown";
#endif

  return Rcpp::List::create(
      Rcpp::Named("cplusplus") = static_cast<double>(__cplusplus),
      Rcpp::Named("compiler") = compiler);
}
