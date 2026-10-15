#' Initial state of a model
#'
#' Builds the state vector (counts per node) at the start of a simulation.
#'
#' @param model An `sg_model` object ([sg_model()]).
#' @param method How the initial state is defined:
#'   \describe{
#'     \item{`"constant"`}{The same count in every sub-stage of each block,
#'       given by `values`, split between genotypes by `genotypes`.}
#'     \item{`"hardy_weinberg"`}{As `"constant"`, with genotype weights from
#'       the allele frequencies `alleles` (Hardy-Weinberg proportions).}
#'     \item{`"state"`}{A full state vector given in `state`.}
#'     \item{`"burn_in"`}{The state reached by running the ODE (without
#'       releases) for `years` before the model's start, from a `"constant"`
#'       start (or `"hardy_weinberg"` if `alleles` is given). The burn-in
#'       weather is that of the years before the start date in `weather`, or
#'       without `weather` the model's first 365 days of forcing repeated.
#'       The ODE forgets its starting state within months, so the result
#'       depends on the weather and parameters, not on `values`.}
#'     \item{`"equilibrium"`}{The equilibrium of a wild-type population at the
#'       constant `temperature` (degC), in closed form (one-dimensional root
#'       search with larval crowding or an Allee effect). Needs a single
#'       wild-type genotype that breeds true. Gives an empty state, with a
#'       message, when the population cannot persist at that temperature.}
#'   }
#' @param values For `"constant"` and `"hardy_weinberg"`: named counts per
#'   sub-stage for blocks `E` (eggs), `L` (larvae), `P` (pupae), `M` (males),
#'   `U` (unmated females) and `F` (mated females), summed over genotypes.
#'   Missing blocks are 0. The default, 0.25 per sub-stage for E, L and P and
#'   0.125 for M and F, reproduces the starting state of Gutierrez & Ponti
#'   (2014): 0.25 per class, adults split equally between the sexes, all
#'   females mated.
#' @param genotypes For `"constant"`: named non-negative weights of genotypes
#'   (normalised to sum to 1). Default: equal weights for the cube's wild
#'   type(s).
#' @param alleles For `"hardy_weinberg"`: named allele frequencies (summing to
#'   1). Genotype names must consist of two allele letters, e.g. `"AA"`,
#'   `"Aa"`, as in single-locus MGDrivE cubes.
#' @param state For `"state"`: numeric vector of length `model$n_states`
#'   (>= 0), unnamed or named as [sg_state_names()].
#' @param per Unit of `values`: `"node"` (counts per node) or `"km2"`
#'   (per km^2, multiplied by the node area).
#' @param years For `"burn_in"`: number of years of burn-in (whole number).
#' @param weather For `"burn_in"`: optional daily weather (as for
#'   [sg_forcing()]) covering the `years` before the model's start date.
#' @param temperature For `"equilibrium"`: constant temperature (degC).
#'   Default: the temperature of the model's forcing, when it is constant
#'   (e.g. from [sg_constant_weather()]).
#'
#' @details
#' Weights are applied block by block, honouring each genotype's biology:
#' * Males (M) go only to genotypes that produce males (`phi < 1`) and females
#'   (U, F) only to genotypes that produce females (`phi > 0`); weights are
#'   renormalised among those genotypes. This matters for sex-linked cubes.
#' * Adults go only to adult classes that the genotype reaches (D-018).
#' * Mated females of genotype g are split across mate genotypes m in
#'   proportion to `eta[g, m]` times the initial male weight of m, i.e.
#'   random mating among the initial adults.
#'
#' @return An object of class `sg_state`: a named numeric vector of length
#'   `model$n_states` (names from [sg_state_names()]) with attribute
#'   `model_hash`, used by the solvers to check that the state belongs to the
#'   model. Burn-in and equilibrium states also carry an attribute `burn_in`
#'   (years, weather source, and `start_effect`: the relative difference to
#'   the same burn-in from a start ten times larger; a warning is given above
#'   1e-3) or
#'   `equilibrium` (temperature, eggs per egg at low density `R0`, wound
#'   search `psi`, eggs laid per day).
#' @seealso [sg_model()], [sg_state_names()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   mod <- sg_model(
#'     p, sg_cube(MGDrivE::cubeMendelian()),
#'     sg_forcing(sg_example_weather()[1:30, ], p)
#'   )
#'   sg_init(mod) # wild type only
#'   sg_init(mod, "hardy_weinberg", alleles = c(A = 0.7, a = 0.3))
#'   sg_init(mod, "equilibrium", temperature = 27) # wild type at 27 degC
#'   # with constant weather the temperature is taken from the forcing
#'   mod_20 <- sg_model(
#'     p, sg_cube(MGDrivE::cubeMendelian()),
#'     sg_forcing(sg_constant_weather(20, days = 30), p)
#'   )
#'   sg_init(mod_20, "equilibrium")
#'   # a year of burn-in over the weather before the start date
#'   wx <- sg_example_weather()
#'   mod_71 <- sg_model(
#'     p, sg_cube(MGDrivE::cubeMendelian()),
#'     sg_forcing(wx, p, start = "1971-01-01", end = "1971-03-31")
#'   )
#'   sg_init(mod_71, "burn_in", weather = wx)
#' }
sg_init <- function(model,
                    method = c("constant", "hardy_weinberg", "state", "burn_in", "equilibrium"),
                    values = c(E = 0.25, L = 0.25, P = 0.25, M = 0.125, U = 0, F = 0.125),
                    genotypes = NULL, alleles = NULL, state = NULL,
                    per = c("node", "km2"), years = 1L, weather = NULL, temperature = NULL) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  method <- .sg_match_choice(
    method, c("constant", "hardy_weinberg", "state", "burn_in", "equilibrium"), "method"
  )
  per <- .sg_match_choice(per, c("node", "km2"), "per")

  if (method %in% c("burn_in", "equilibrium")) {
    res <- if (method == "burn_in") {
      weights <- if (is.null(alleles)) {
        .sg_init_weights(model, genotypes)
      } else {
        .sg_hw_weights(model, alleles)
      }
      start <- .sg_init_constant(model, .sg_init_values(values, per, model$area), weights)
      .sg_init_burn_in(model, start, years, weather)
    } else {
      if (is.null(temperature)) {
        temperature <- .sg_forcing_constant_temp(model$forcing, model$params)
        if (is.null(temperature)) {
          .sg_stop(
            "method = \"equilibrium\" needs a `temperature` (degC): the model's forcing ",
            "does not have one constant temperature."
          )
        }
      }
      .sg_init_equilibrium(model, temperature)
    }
    x <- res$x
    names(x) <- sg_state_names(model)
    out <- structure(x, class = "sg_state", model_hash = model$meta$hash)
    attr(out, method) <- res$info
    return(out)
  }
  x <- switch(EXPR = method,
    state = .sg_init_state(model, state),
    constant = .sg_init_constant(
      model, .sg_init_values(values, per, model$area),
      .sg_init_weights(model, genotypes)
    ),
    hardy_weinberg = .sg_init_constant(
      model, .sg_init_values(values, per, model$area),
      .sg_hw_weights(model, alleles)
    )
  )
  names(x) <- sg_state_names(model)
  structure(x, class = "sg_state", model_hash = model$meta$hash)
}

#' @export
print.sg_state <- function(x, ...) {
  tot <- .sg_state_block_totals(unclass(x))
  cat("<sg_state> ", length(x), " states, total ", format(sum(x), digits = 6L), "\n", sep = "")
  b <- attr(x, "burn_in")
  if (!is.null(b)) {
    cat(sprintf(
      "  after a %d-year burn-in (%s); effect of the starting state %.2g\n",
      b$years, b$source, b$start_effect
    ))
  }
  e <- attr(x, "equilibrium")
  if (!is.null(e)) {
    cat(sprintf(
      "  equilibrium at %g degC: R0 = %.3g, wound search %.3g, %.4g eggs per day\n",
      e$temperature, e$R0, e$psi, e$eggs_per_day
    ))
  }
  print(round(tot, 6L))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# Totals by block (rows) and genotype (columns), from state names.
.sg_state_block_totals <- function(x) {
  nm <- names(x)
  block <- substr(nm, 1L, 1L)
  geno <- sub("^.\\[([^,|]+).*$", "\\1", nm)
  blocks <- c("E", "L", "P", "M", "U", "F")
  genos <- unique(geno)
  out <- matrix(0, length(blocks), length(genos), dimnames = list(blocks, genos))
  sums <- tapply(unname(x), list(block, geno), sum)
  out[rownames(sums), colnames(sums)] <- sums
  out[is.na(out)] <- 0
  out
}

.sg_init_values <- function(values, per, area) {
  blocks <- c("E", "L", "P", "M", "U", "F")
  ok <- is.numeric(values) && !is.null(names(values)) && all(nzchar(names(values)))
  if (!ok) {
    .sg_stop("`values` must be a named numeric vector, e.g. c(E = 0.25, M = 0.125).")
  }
  bad <- setdiff(names(values), blocks)
  if (length(bad) > 0L) {
    .sg_stop(
      "Unknown block(s) in `values`: ", paste(bad, collapse = ", "),
      ". Use E, L, P, M, U and F."
    )
  }
  if (any(!is.finite(values)) || any(values < 0)) {
    .sg_stop("`values` must be finite and >= 0.")
  }
  out <- stats::setNames(rep(0, length(blocks)), blocks)
  out[names(values)] <- values
  if (per == "km2") out <- out * area
  out
}

.sg_init_weights <- function(model, genotypes) {
  geno <- model$cube$genotypes
  if (is.null(genotypes)) {
    w <- stats::setNames(as.numeric(geno %in% model$cube$wild_type), geno)
    return(w / sum(w))
  }
  ok <- is.numeric(genotypes) && !is.null(names(genotypes)) && all(is.finite(genotypes)) &&
    all(genotypes >= 0) && sum(genotypes) > 0
  if (!ok) {
    .sg_stop("`genotypes` must be named non-negative weights with a positive sum.")
  }
  bad <- setdiff(names(genotypes), geno)
  if (length(bad) > 0L) {
    .sg_stop(
      "Genotype(s) in `genotypes` not in the cube: ", paste(bad, collapse = ", "), ".",
      .sg_suggest(bad[1L], geno)
    )
  }
  w <- stats::setNames(rep(0, length(geno)), geno)
  w[names(genotypes)] <- genotypes
  w / sum(w)
}

# Hardy-Weinberg genotype weights for two-letter single-locus genotype names.
.sg_hw_weights <- function(model, alleles) {
  geno <- model$cube$genotypes
  ok <- is.numeric(alleles) && !is.null(names(alleles)) && all(is.finite(alleles)) &&
    all(alleles >= 0) && abs(sum(alleles) - 1) < 1e-8
  if (!ok) {
    .sg_stop("`alleles` must be named allele frequencies (>= 0) summing to 1.")
  }
  if (any(nchar(geno) != 2L)) {
    .sg_stop(
      "method = \"hardy_weinberg\" needs genotype names of two allele letters ",
      "(e.g. \"AA\", \"Aa\"); this cube has ", .sg_wrap_names(geno),
      ". Use method = \"constant\" with `genotypes` instead."
    )
  }
  a1 <- substr(geno, 1L, 1L)
  a2 <- substr(geno, 2L, 2L)
  bad <- setdiff(unique(c(a1, a2)), names(alleles))
  if (length(bad) > 0L) {
    .sg_stop(
      "Allele(s) ", paste(bad, collapse = ", "),
      " of the cube's genotypes are missing from `alleles`."
    )
  }
  w <- alleles[a1] * alleles[a2] * ifelse(a1 == a2, 1, 2)
  w <- stats::setNames(as.numeric(w), geno)
  if (sum(w) <= 0) {
    .sg_stop("The allele frequencies give zero weight to every genotype.")
  }
  w / sum(w)
}

.sg_init_constant <- function(model, values, w) {
  cb <- model$cube
  geno <- cb$genotypes
  x <- numeric(model$n_states)
  renorm <- function(v) if (sum(v) > 0) v / sum(v) else v

  for (b in c("E", "L", "P")) {
    for (g in geno) {
      x[sg_state_index(model, b, genotype = g)] <- values[[b]] * w[[g]]
    }
  }
  w_m <- renorm(w * (cb$phi < 1))
  w_f <- renorm(w * (cb$phi > 0))
  for (g in geno) {
    reach_m <- seq_len(model$adult$male$last_class[[g]])
    reach_f <- seq_len(model$adult$female$last_class[[g]])
    x[sg_state_index(model, "M", genotype = g, substage = reach_m)] <- values[["M"]] * w_m[[g]]
    x[sg_state_index(model, "U", genotype = g, substage = reach_f)] <- values[["U"]] * w_f[[g]]
    if (values[["F"]] > 0 && w_f[[g]] > 0) {
      mate_share <- renorm(cb$eta[g, ] * w_m)
      if (sum(mate_share) == 0) {
        .sg_stop(
          "Mated females of genotype ", g, " cannot be initialised: no initial males ",
          "of a genotype it can mate with (eta > 0)."
        )
      }
      for (m in geno[mate_share > 0]) {
        idx <- sg_state_index(model, "F", genotype = g, mate = m, substage = reach_f)
        x[idx] <- values[["F"]] * w_f[[g]] * mate_share[[m]]
      }
    }
  }
  x
}

.sg_init_state <- function(model, state) {
  if (is.null(state)) {
    .sg_stop("method = \"state\" needs `state`, a numeric vector of length ", model$n_states, ".")
  }
  ok <- is.numeric(state) && length(state) == model$n_states && all(is.finite(state)) &&
    all(state >= 0)
  if (!ok) {
    .sg_stop(
      "`state` must be a numeric vector of length ", model$n_states,
      " (the model's number of states) with finite values >= 0."
    )
  }
  if (!is.null(names(state)) && !identical(names(state), sg_state_names(model))) {
    .sg_stop("Names of `state` do not match sg_state_names(model) (same names, same order).")
  }
  as.numeric(state)
}
