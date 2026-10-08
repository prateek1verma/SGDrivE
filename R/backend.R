#' Information about the compiled C++ backend
#'
#' Reports the package version and how the C++ backend was compiled. This is
#' a toolchain check added in milestone M0: it confirms that the C++ code
#' compiles, links and can be called on the current platform. It contains no
#' model code.
#'
#' @return A named list with elements:
#'   \describe{
#'     \item{`package_version`}{Installed version of SGDrivE (character).}
#'     \item{`cplusplus`}{Value of the `__cplusplus` macro the backend was
#'       compiled with (numeric; e.g. `201703` for C++17).}
#'     \item{`compiler`}{Compiler identification string (character).}
#'   }
#' @export
#' @examples
#' sg_backend_info()
sg_backend_info <- function() {
  info <- sg_backend_info_cpp()
  c(list(package_version = as.character(utils::packageVersion("SGDrivE"))), info)
}
