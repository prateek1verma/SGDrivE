# Movement layers from transition probabilities, and conversions between
# leave rates and leave probabilities (multi-node work phase 3, WP2; design
# note dev/MOVEMENT_DESIGN.md section 6.1, D-032).

#' Movement layer from daily transition probabilities
#'
#' Converts a matrix of daily transition probabilities between nodes into a
#' movement layer, as MGDrivE2's `movement_prob2rate()`: the leave rate of
#' node i is \eqn{\delta_i = -\ln P_{ii}}{delta[i] = -log(P[i, i])} and the
#' destination probabilities are \eqn{\pi_{ij} = P_{ij} / (1 - P_{ii})}{pi[i, j]
#' = P[i, j] / (1 - P[i, i])}.
#'
#' @param P Square matrix of daily probabilities: `P[i, j]` that a fly in
#'   node i is in node j a day later, `P[i, i]` that it stays (rows: from,
#'   columns: to; node ids as dimnames or in node order). Each row sums to at
#'   most 1; the remainder of a row leaves the landscape (open edges,
#'   counted as lost). `P[i, i]` must be positive.
#' @param clock Movement clock, as in [sg_movement_layer()]. MGDrivE2 moves
#'   flies in calendar time; use `"calendar"` to compare with it.
#'
#' @details The leave probability within a day, \eqn{1 - P_{ii}}{1 - P[i, i]},
#'   is reproduced exactly. The destination split ignores flies that move
#'   more than once within a day, as in MGDrivE2, so `P` is not the matrix
#'   exponential of the layer's generator when leave rates are large.
#'
#'   The same rates apply to males, unmated and mated females. The result is
#'   an ordinary layer, so class rates can be changed afterwards, e.g.
#'   `layer$rate$virgin <- 0` (MGDrivE2 does not move unmated females). The
#'   rates are taken as rates at `movement_t_ref` ([sg_params()]); with the
#'   physiological clock they change with temperature (D-032).
#'
#' @return An `sg_movement_layer` ([sg_movement_layer()]) with `n_nodes` set.
#' @seealso [sg_movement_layer()], [sg_movement_diffusion()], [sg_rate_from_fraction()]
#' @export
#' @examples
#' # three nodes: 90% stay each day, the rest go to the neighbours; node 3
#' # loses 2% a day to the outside
#' P <- rbind(
#'   c(0.90, 0.10, 0.00),
#'   c(0.05, 0.90, 0.05),
#'   c(0.00, 0.08, 0.90)
#' )
#' layer <- sg_movement_from_prob(P, clock = "calendar")
#' layer$rate$male # -log(0.9) in every node
#' layer$edges
sg_movement_from_prob <- function(P, # nolint: object_name_linter.
                                  clock = c("physiological", "calendar")) {
  clock <- .sg_match_choice(clock, c("physiological", "calendar"), "clock")
  if (!is.matrix(P) || !is.numeric(P) || nrow(P) != ncol(P) || nrow(P) < 1L) {
    .sg_stop("`P` must be a numeric square matrix of daily transition probabilities.")
  }
  if (any(!is.finite(P)) || any(P < 0) || any(P > 1)) {
    .sg_stop("Entries of `P` must be probabilities in [0, 1].")
  }
  if (!identical(rownames(P), colnames(P))) {
    .sg_stop("Row and column names of `P` must be the same node ids in the same order.")
  }
  stay <- diag(P)
  if (any(stay <= 0)) {
    .sg_stop(
      "Every node must keep some flies within a day (P[i, i] > 0); node ",
      which(stay <= 0)[1L], " has P[i, i] = 0."
    )
  }
  if (any(rowSums(P) > 1 + 1e-12)) {
    .sg_stop("Each row of `P` must sum to at most 1 (the remainder leaves the landscape).")
  }
  n <- nrow(P)
  delta <- -log(stay)
  dest <- P
  diag(dest) <- 0
  leave <- 1 - stay
  dest <- dest / ifelse(leave > 0, leave, 1)
  dest[leave == 0, ] <- 0
  dest <- pmin(dest, 1)
  layer <- sg_movement_layer(list(male = delta, virgin = delta, mated = delta), dest, clock = clock)
  layer$n_nodes <- n
  layer
}

#' Leave rates from leave probabilities
#'
#' Converts field measures of emigration into leave rates per day at the
#' reference temperature `movement_t_ref`, the rates [sg_movement_layer()]
#' takes (D-032):
#' * `sg_rate_from_fraction()`: the share of flies that leaves a node within
#'   a day, \eqn{\delta = -\ln(1 - f)}{delta = -log(1 - f)}; for example
#'   Matlock & Skoda (2009) estimated that about 16% of females and 10% of
#'   males leave a station array per day;
#' * `sg_rate_from_lifetime()`: the probability that a fly leaves its node at
#'   least once in its life, inverted exactly for the model's adult chain at
#'   `movement_t_ref`.
#'
#' @param f Daily leave fractions in \[0, 1).
#' @param p Lifetime leave probabilities in \[0, 1).
#' @param model An `sg_model` object ([sg_model()]): its adult chain
#'   (sub-stages, lifespan of the genotype), temperature mortality,
#'   background mortality and movement parameters.
#' @param class `"mated"`, `"virgin"` or `"male"`: the female chain for the
#'   first two, the male chain for males.
#' @param clock The layer's clock: with `"physiological"` (default) the rate
#'   at `movement_t_ref` is scaled by kappa = min(`movement_kappa_max`, 1)
#'   there; with `"calendar"` it is not.
#' @param genotype Genotype whose adult lifespan is used (default the first
#'   wild-type genotype).
#'
#' @details The lifetime leave probability of a fly with leave rate delta is
#'   \eqn{1 - E[e^{-\delta \tau}]}{1 - E[exp(-delta tau)]}, with tau the adult
#'   lifespan: the time until death by the adult hazard (temperature and
#'   background mortality) or by old age at the end of the adult chain,
#'   computed exactly from the chain's sub-stage rates at `movement_t_ref`
#'   (no simulation). Moving does not change a fly's age or survival. For
#'   unmated females it treats the whole adult life as unmated, so it is the
#'   probability for a female that never mates.
#'
#' @return Leave rates per day at `movement_t_ref`, one per value of `f` or
#'   `p`.
#' @seealso [sg_movement_layer()], [sg_movement_from_prob()], [sg_movement_diffusion()]
#' @export
#' @examples
#' sg_rate_from_fraction(c(female = 0.16, male = 0.10))
#'
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), sg_forcing(sg_constant_weather(27), p))
#'   # leave rate at which half of the mated females leave their node in their life
#'   sg_rate_from_lifetime(0.5, mod, class = "mated")
#' }
sg_rate_from_fraction <- function(f) {
  if (!is.numeric(f) || length(f) < 1L || any(!is.finite(f)) || any(f < 0 | f >= 1)) {
    .sg_stop("`f` must be daily leave fractions in [0, 1).")
  }
  -log1p(-f)
}

#' @rdname sg_rate_from_fraction
#' @export
sg_rate_from_lifetime <- function(p, model, class = c("mated", "virgin", "male"),
                                  clock = c("physiological", "calendar"), genotype = NULL) {
  if (!is.numeric(p) || length(p) < 1L || any(!is.finite(p)) || any(p < 0 | p >= 1)) {
    .sg_stop("`p` must be lifetime leave probabilities in [0, 1).")
  }
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  class <- .sg_match_choice(class, c("mated", "virgin", "male"), "class")
  clock <- .sg_match_choice(clock, c("physiological", "calendar"), "clock")
  t_ref <- model$params$movement_t_ref
  life <- .sg_adult_life(model, class, t_ref, genotype)
  kappa <- if (clock == "calendar") 1 else min(model$params$movement_kappa_max, 1)
  rate <- vapply(p, function(target) {
    if (target == 0) {
      return(0)
    }
    hi <- 1
    while (.sg_lifetime_leave(hi, life) < target) {
      hi <- hi * 2
      if (hi > 1e12) .sg_stop("No leave rate gives a lifetime leave probability of ", target, ".")
    }
    stats::uniroot(function(d) .sg_lifetime_leave(d, life) - target,
      lower = 0, upper = hi, tol = 1e-14 * hi, maxiter = 2000L
    )$root
  }, 0)
  stats::setNames(rate / kappa, names(p))
}

# ---- internal ------------------------------------------------------------------

# The adult chain of a class at constant temperature `temp` for one genotype
# (default the first wild-type genotype): transfer rates q of the classes it
# reaches (dd per day over sub-stage widths) and the hazard mu (temperature
# and background mortality, S2, S7).
.sg_adult_life <- function(model, class, temp, genotype = NULL) {
  sex <- if (class == "male") "male" else "female"
  ch <- model$adult[[sex]]
  if (is.null(genotype)) genotype <- model$cube$wild_type[1L]
  if (!genotype %in% rownames(ch$width)) {
    .sg_stop(
      "Unknown genotype \"", genotype, "\".", .sg_suggest(genotype, rownames(ch$width))
    )
  }
  p <- model$params
  w <- ch$width[genotype, seq_len(ch$last_class[[genotype]])]
  background <- if (sex == "male") p$male_background_mortality else p$female_background_mortality
  list(
    q = sg_degree_days(temp, p) / w,
    mu = sg_temp_mortality(temp, p) + background
  )
}

# Lifetime leave probability 1 - E[exp(-delta tau)] for the adult chain
# `life` (.sg_adult_life()), vectorised over delta. From class i the
# probability of dying (by the hazard or old age) before the first move is
# f_i = (mu + q_i f_{i+1}) / (delta + mu + q_i), with f = 1 after the last
# class (old age).
.sg_lifetime_leave <- function(delta, life) {
  f <- rep(1, length(delta))
  for (q in rev(life$q)) {
    den <- delta + life$mu + q
    f <- ifelse(den > 0, (life$mu + q * f) / den, 1)
  }
  1 - f
}

# Expected adult lifespan E[tau] (days) of the chain `life`:
# g_i = (1 + q_i g_{i+1}) / (mu + q_i), g = 0 after the last class.
.sg_expected_life <- function(life) {
  g <- 0
  for (q in rev(life$q)) {
    den <- life$mu + q
    if (den == 0) {
      return(Inf)
    }
    g <- (1 + q * g) / den
  }
  g
}
