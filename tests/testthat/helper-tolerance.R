# Shared numerical tolerances for ScrewDrivE tests.
#
# Conventions (see dev/CONVENTIONS.md, "Testing"):
# * sg_tol_exact: identities that should hold to floating-point precision
#   (e.g. probabilities summing to one, flow accounting).
# * sg_tol_ode: agreement between two deterministic solvers of the same
#   model (e.g. C++ RK45 vs the pure-R deSolve reference), relative error.
# * sg_tol_reference: agreement with the legacy Gutierrez & Ponti (2014)
#   R translations, which use a different time discretisation.
sg_tol_exact <- 1e-12
sg_tol_ode <- 1e-6
sg_tol_reference <- 1e-2
