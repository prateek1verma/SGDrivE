# Minimal hand-built cubes for tests that must not depend on MGDrivE.

# Two genotypes: wild type W and a released sterile genotype S (SIT-like).
# Crosses with S males produce non-viable offspring (tau = residual).
sg_test_cube_sit <- function(residual = 0, omega_s = 1) {
  g <- c("W", "S")
  ih <- array(0, c(2, 2, 2), dimnames = list(g, g, g))
  ih[, , "W"] <- 1 # all offspring are W (sterility acts through tau)
  tau <- array(1, c(2, 2, 2), dimnames = list(g, g, g))
  tau[, "S", ] <- residual
  list(
    genotypesID = g, genotypesN = 2L, wildType = "W", releaseType = "S",
    ih = ih, tau = tau,
    eta = matrix(1, 2, 2, dimnames = list(g, g)),
    phi = c(W = 0.5, S = 0.5), xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1),
    s = c(W = 1, S = 1), omega = c(W = 1, S = omega_s)
  )
}

# Three-genotype Mendelian cube built by hand (same as MGDrivE::cubeMendelian()).
sg_test_cube_mendelian <- function() {
  g <- c("AA", "Aa", "aa")
  gam <- rbind(AA = c(1, 0), Aa = c(0.5, 0.5), aa = c(0, 1)) # P(A), P(a)
  ih <- array(0, c(3, 3, 3), dimnames = list(g, g, g))
  for (f in g) {
    for (m in g) {
      pa <- gam[f, ]
      pm <- gam[m, ]
      ih[f, m, ] <- c(pa[1] * pm[1], pa[1] * pm[2] + pa[2] * pm[1], pa[2] * pm[2])
    }
  }
  list(
    genotypesID = g, genotypesN = 3L, wildType = "AA", releaseType = "aa",
    ih = ih, tau = array(1, c(3, 3, 3), dimnames = list(g, g, g)),
    eta = matrix(1, 3, 3, dimnames = list(g, g)),
    phi = stats::setNames(rep(0.5, 3), g), xiF = stats::setNames(rep(1, 3), g),
    xiM = stats::setNames(rep(1, 3), g), s = stats::setNames(rep(1, 3), g),
    omega = stats::setNames(rep(1, 3), g)
  )
}

# One genotype only (wild type), the configuration of Gutierrez & Ponti (2014).
sg_test_cube_wild <- function() {
  g <- "W"
  list(
    genotypesID = g, genotypesN = 1L, wildType = g,
    ih = array(1, c(1, 1, 1), dimnames = list(g, g, g)),
    tau = array(1, c(1, 1, 1), dimnames = list(g, g, g)),
    eta = matrix(1, 1, 1, dimnames = list(g, g)),
    phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1), omega = c(W = 1)
  )
}
