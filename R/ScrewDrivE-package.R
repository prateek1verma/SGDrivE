#' ScrewDrivE: Screwworm Genetic-Control Dynamics Explorer
#'
#' Genotype-resolved population dynamics of the New World screwworm
#' (*Cochliomyia hominivorax*) under genetic control.
#'
#' The package is under development (milestone M0: foundations). The planned
#' model combines:
#' * the temperature-driven life history of Gutierrez & Ponti (2014):
#'   degree-day development, temperature-dependent mortality, adult ageing
#'   chains, a pre-oviposition period and wound-limited oviposition;
#' * genotype inheritance read unchanged from 'MGDrivE' inheritance cubes;
#' * the unmated/mated female mating structure of 'MGDrivE2'.
#'
#' Deterministic (ODE) and stochastic (chain-binomial) solvers will share one
#' C++ rate kernel. Populations are tracked as counts per node.
#'
#' @references
#' Gutierrez, A.P. & Ponti, L. (2014). The new world screwworm: prospective
#' distribution and role of weather in eradication. *Agricultural and Forest
#' Entomology* 16:158-173. \doi{10.1111/afe.12046}
#'
#' Sanchez C., H.M., Wu, S.L., Bennett, J.B. & Marshall, J.M. (2020). MGDrivE:
#' A modular simulation framework for the spread of gene drives through
#' spatially explicit mosquito populations. *Methods in Ecology and Evolution*
#' 11:229-239. \doi{10.1111/2041-210X.13318}
#'
#' Wu, S.L., Bennett, J.B., Sanchez C., H.M., Dolgert, A.J., Leon, T.M. &
#' Marshall, J.M. (2021). MGDrivE 2: A simulation framework for gene drive
#' systems incorporating seasonality and epidemiological dynamics. *PLoS
#' Computational Biology* 17(5):e1009030. \doi{10.1371/journal.pcbi.1009030}
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @useDynLib ScrewDrivE, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom stats update
## usethis namespace: end
NULL
