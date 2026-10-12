# Run a test once for each simulation backend (milestone M3: every M2 check of
# the R reference must also hold for the C++ backend). Inside `code`,
# `backend` is the backend's name and `sim_ode()` is simulate_ode() with that
# backend; the test description gets the backend appended.
sg_test_backends <- c("R", "cpp")

sg_test_each_backend <- function(desc, code) {
  code <- substitute(code)
  env <- parent.frame()
  for (backend in sg_test_backends) {
    test_env <- new.env(parent = env)
    test_env$backend <- backend
    test_env$sim_ode <- local({
      b <- backend
      function(...) simulate_ode(..., backend = b)
    })
    testthat::test_that(paste0(desc, " [", backend, "]"), eval(code, test_env))
  }
}
