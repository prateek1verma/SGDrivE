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
#'     \item{`"threshold"`}{The Allee threshold of a wild-type population at
#'       the constant `temperature`: the lower, unstable equilibrium that
#'       exists with an Allee effect (`male_density_half_mating > 0`), below
#'       which too few females find a mate and the population declines
#'       (D-029). Same requirements as `"equilibrium"`; gives an empty state,
#'       with a message, when there is no threshold.}
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
#'   (>= 0), unnamed or named as [sg_state_names()]; for a model with a
#'   landscape also a list with one vector of the states of a node per node,
#'   named by the node ids.
#' @param per Unit of `values`: `"node"` (counts per node) or `"km2"`
#'   (per km^2, multiplied by the node area).
#' @param years For `"burn_in"`: number of years of burn-in (whole number).
#' @param weather For `"burn_in"`: optional daily weather (as for
#'   [sg_forcing()]) covering the `years` before the model's start date. With
#'   a landscape: one data frame for all nodes or a list with one per node
#'   (named by node id, or in node order).
#' @param temperature For `"equilibrium"` and `"threshold"`: constant
#'   temperature (degC); with a landscape one value or one per node.
#'   Default: the temperature of the model's forcing (of each node's forcing),
#'   when it is constant (e.g. from [sg_constant_weather()]).
#' @param nodes For `"constant"`, `"hardy_weinberg"` and `"burn_in"` with a
#'   landscape ([sg_landscape()]): ids or numbers of the nodes that receive
#'   the starting population (default all nodes); the other nodes start
#'   without flies. Ignored for one-node models.
#' @param node_extinction For `"burn_in"`: apply the node extinction rule
#'   during the burn-in run (see [simulate_ode()]), so that a burn-in with
#'   movement does not leave fractions of flies in every node. Default
#'   `FALSE`.
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
#' **Several nodes** (D-032). Each node is initialised as a one-node model
#' with its own forcing, area, wounds or hosts. `"constant"` and
#' `"hardy_weinberg"` put `values` into each node of `nodes` (counts per node,
#' or per km^2 times the node's area); the default start therefore puts the
#' Gutierrez & Ponti starting state into every node, as for one node. Every
#' node gets its hosts at their steady state without flies. `"equilibrium"`
#' and `"threshold"` are computed node by node without movement (exact for a
#' node whose net movement is zero, otherwise a start for a burn-in); nodes
#' without an equilibrium start empty. `"burn_in"` runs all nodes with
#' movement.
#'
#' @return An object of class `sg_state`: a named numeric vector of length
#'   `model$n_states` (names from [sg_state_names()]) with attribute
#'   `model_hash`, used by the solvers to check that the state belongs to the
#'   model. Burn-in and equilibrium states also carry an attribute `burn_in`
#'   (years, weather source, and `start_effect`: the relative difference to
#'   the same burn-in from a start ten times larger; a warning is given above
#'   1e-3) or
#'   `equilibrium` or `threshold` (temperature, eggs per egg at low density
#'   `R0`, wound search `psi`, eggs laid per day, `adult_females`: unmated
#'   plus mated females in the node; with a landscape a data frame with one
#'   row per node).
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
                    method = c(
                      "constant", "hardy_weinberg", "state", "burn_in", "equilibrium",
                      "threshold"
                    ),
                    values = c(E = 0.25, L = 0.25, P = 0.25, M = 0.125, U = 0, F = 0.125),
                    genotypes = NULL, alleles = NULL, state = NULL,
                    per = c("node", "km2"), years = 1L, weather = NULL, temperature = NULL,
                    nodes = NULL, node_extinction = FALSE) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  method <- .sg_match_choice(
    method, c("constant", "hardy_weinberg", "state", "burn_in", "equilibrium", "threshold"),
    "method"
  )
  per <- .sg_match_choice(per, c("node", "km2"), "per")
  if (!is.logical(node_extinction) || length(node_extinction) != 1L || is.na(node_extinction)) {
    .sg_stop("`node_extinction` must be TRUE or FALSE.")
  }
  start_nodes <- .sg_node_numbers(model, nodes)
  multi <- !is.null(model$landscape)

  if (method %in% c("burn_in", "equilibrium", "threshold")) {
    res <- if (method == "burn_in") {
      weights <- if (is.null(alleles)) {
        .sg_init_weights(model, genotypes)
      } else {
        .sg_hw_weights(model, alleles)
      }
      start <- .sg_init_constant_nodes(model, values, per, weights, start_nodes)
      .sg_init_burn_in(model, start, years, weather, node_extinction = node_extinction)
    } else if (multi) {
      .sg_init_equilibrium_nodes(model, temperature, method)
    } else {
      if (is.null(temperature)) {
        temperature <- .sg_forcing_constant_temp(model$forcing, model$params)
        if (is.null(temperature)) {
          .sg_stop(
            "method = \"", method, "\" needs a `temperature` (degC): the model's forcing ",
            "does not have one constant temperature."
          )
        }
      }
      .sg_init_equilibrium(model, temperature,
        what = paste0("method = \"", method, "\""),
        root = if (method == "threshold") "lower" else "upper"
      )
    }
    x <- res$x
    names(x) <- sg_state_names(model)
    out <- structure(x, class = "sg_state", model_hash = model$meta$hash)
    attr(out, method) <- res$info
    return(out)
  }
  x <- switch(EXPR = method,
    state = .sg_init_state(model, .sg_state_from_list(model, state)),
    constant = .sg_init_constant_nodes(
      model, values, per, .sg_init_weights(model, genotypes), start_nodes
    ),
    hardy_weinberg = .sg_init_constant_nodes(
      model, values, per, .sg_hw_weights(model, alleles), start_nodes
    )
  )
  names(x) <- sg_state_names(model)
  structure(x, class = "sg_state", model_hash = model$meta$hash)
}

#' @export
print.sg_state <- function(x, ...) {
  if (any(grepl("@", names(x), fixed = TRUE))) {
    return(.sg_print_state_nodes(x))
  }
  tot <- .sg_state_block_totals(unclass(x))
  cat("<sg_state> ", length(x), " states, total ", format(sum(x), digits = 6L), "\n", sep = "")
  b <- attr(x, "burn_in")
  if (!is.null(b)) {
    cat(sprintf(
      "  after a %d-year burn-in (%s); effect of the starting state %.2g\n",
      b$years, b$source, b$start_effect
    ))
  }
  for (kind in c("equilibrium", "threshold")) {
    e <- attr(x, kind)
    if (is.null(e)) next
    cat(sprintf(
      "  %s at %g degC: R0 = %.3g, wound search %.3g, %.4g eggs per day\n",
      if (kind == "threshold") "Allee threshold (unstable equilibrium)" else kind,
      e$temperature, e$R0, e$psi, e$eggs_per_day
    ))
  }
  print(round(tot, 6L))
  host <- grepl("^H[SWI]\\[", names(x))
  if (any(host)) {
    v <- unclass(x)[host]
    cat(sprintf(
      "hosts: %s without a wound, %s wounded, %s infested (%s in all)\n",
      format(sum(v[startsWith(names(v), "HS")]), digits = 6L),
      format(sum(v[startsWith(names(v), "HW")]), digits = 6L),
      format(sum(v[startsWith(names(v), "HI")]), digits = 6L),
      format(sum(v), digits = 6L)
    ))
  }
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# print() of a state of several nodes: fly totals by stage and node, hosts by
# node, and how the state was made.
.sg_print_state_nodes <- function(x) {
  v <- unclass(x)
  nm <- names(v)
  node <- sub("^[^@]*@", "", nm)
  ids <- unique(node)
  host <- grepl("^H[SWI]\\[", nm)
  cat(
    "<sg_state> ", length(v), " states in ", length(ids), " nodes, total ",
    format(sum(v[!host]), digits = 6L), " flies\n",
    sep = ""
  )
  b <- attr(x, "burn_in")
  if (!is.null(b)) {
    cat(sprintf(
      "  after a %d-year burn-in (%s); effect of the starting state %.2g\n",
      b$years, b$source, b$start_effect
    ))
  }
  for (kind in c("equilibrium", "threshold")) {
    e <- attr(x, kind)
    if (is.null(e)) next
    cat(
      "  ", if (kind == "threshold") "Allee threshold" else "equilibrium",
      " per node (without movement):\n",
      sep = ""
    )
    print(e, row.names = FALSE, digits = 4L)
  }
  blocks <- c("E", "L", "P", "M", "U", "F")
  stage <- substr(nm, 1L, 1L)
  tot <- vapply(ids, function(id) {
    vapply(blocks, function(bk) sum(v[!host & node == id & stage == bk]), 0)
  }, numeric(length(blocks)))
  tot <- matrix(tot, length(blocks), length(ids), dimnames = list(blocks, ids))
  print(round(tot, 6L))
  if (any(host)) {
    hosts <- vapply(ids, function(id) {
      sel <- host & node == id
      c(
        without_wound = sum(v[sel & startsWith(nm, "HS")]),
        wounded = sum(v[sel & startsWith(nm, "HW")]),
        infested = sum(v[sel & startsWith(nm, "HI")])
      )
    }, numeric(3L))
    cat("hosts:\n")
    print(round(matrix(hosts, 3L, length(ids), dimnames = list(
      c("without wound", "wounded", "infested"), ids
    )), 6L))
  }
  invisible(x)
}

# Positions of all host states (D-030); empty without a host block.
.sg_host_index <- function(model) {
  if (is.null(model$hosts)) {
    return(integer(0))
  }
  unlist(lapply(.sg_host_blocks, function(b) sg_state_index(model, b)), use.names = FALSE)
}

# Totals by block (rows) and genotype (columns), from state names; host
# states are left out.
.sg_state_block_totals <- function(x) {
  x <- x[!grepl("^H[SWI]\\[", names(x))]
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
  .sg_init_hosts(model, x)
}

# Hosts at their steady state without flies (D-030): a fraction
# wounded_fraction wounded, none infested. No change without a host block.
.sg_init_hosts <- function(model, x) {
  h <- model$hosts
  if (is.null(h)) {
    return(x)
  }
  p0 <- model$params$wounded_fraction
  x[sg_state_index(model, "HS")] <- (1 - p0) * h$n
  x[sg_state_index(model, "HW")] <- p0 * h$n
  x[sg_state_index(model, "HI")] <- 0
  x
}

# "constant" (or "hardy_weinberg") start in the nodes `start_nodes` (all
# nodes of a one-node model): each such node gets the one-node start of its
# own one-node model (values per km^2 scaled by that node's area); every node
# gets its hosts at their fly-free steady state (D-032).
.sg_init_constant_nodes <- function(model, values, per, weights, start_nodes) {
  if (is.null(model$landscape)) {
    return(.sg_init_constant(model, .sg_init_values(values, per, model$area), weights))
  }
  x <- numeric(model$n_states)
  for (i in seq_len(model$dims$n_node)) {
    nm <- .sg_node_model(model, i)
    x[.sg_node_positions(model, i)] <- if (i %in% start_nodes) {
      .sg_init_constant(nm, .sg_init_values(values, per, nm$area), weights)
    } else {
      .sg_init_hosts(nm, numeric(nm$n_states))
    }
  }
  x
}

# method = "state" with a named list of per-node vectors (one per node, in any
# order) or a full vector.
.sg_state_from_list <- function(model, state) {
  if (!is.list(state)) {
    return(state)
  }
  ids <- .sg_node_ids(model)
  if (is.null(names(state)) || !setequal(names(state), ids) || anyDuplicated(names(state))) {
    .sg_stop(
      "A list `state` needs one element per node, named by the node ids (",
      .sg_wrap_names(ids), ")."
    )
  }
  n_s <- .sg_n_states_node(model)
  parts <- lapply(ids, function(id) {
    v <- state[[id]]
    if (!is.numeric(v) || length(v) != n_s) {
      .sg_stop("`state$", id, "` must be a numeric vector of length ", n_s, " (states per node).")
    }
    as.numeric(v)
  })
  unlist(parts, use.names = FALSE)
}

# "equilibrium" or "threshold" in every node of a landscape, each from its
# one-node model without movement (D-032). `temperature`: NULL (each node's
# constant forcing temperature), one value or one per node. Messages of the
# nodes are collected into one.
.sg_init_equilibrium_nodes <- function(model, temperature, method) {
  ids <- .sg_node_ids(model)
  temps <- .sg_node_temperatures(model, temperature, paste0("method = \"", method, "\""))
  x <- numeric(model$n_states)
  notes <- character(0)
  rows <- vector("list", length(ids))
  for (i in seq_along(ids)) {
    nm <- .sg_node_model(model, i)
    res <- withCallingHandlers(
      .sg_init_equilibrium(nm, temps[i],
        what = paste0("method = \"", method, "\""),
        root = if (method == "threshold") "lower" else "upper"
      ),
      message = function(m) {
        notes <<- c(notes, paste0("node ", ids[i], ": ", trimws(conditionMessage(m))))
        invokeRestart("muffleMessage")
      }
    )
    x[.sg_node_positions(model, i)] <- res$x
    inf <- res$info
    rows[[i]] <- data.frame(
      node = ids[i], temperature = inf$temperature, R0 = inf$R0, psi = inf$psi,
      eggs_per_day = inf$eggs_per_day, adult_females = inf$adult_females,
      stringsAsFactors = FALSE
    )
  }
  if (length(notes) > 0L) .sg_message(paste(notes, collapse = "\n"))
  info <- do.call(rbind, rows)
  rownames(info) <- NULL
  list(x = x, info = info)
}

# Temperature of each node for equilibria: given (one value or one per node)
# or each node's constant forcing temperature.
.sg_node_temperatures <- function(model, temperature, what) {
  ids <- .sg_node_ids(model)
  n <- length(ids)
  if (!is.null(temperature)) {
    ok <- is.numeric(temperature) && length(temperature) %in% c(1L, n) &&
      all(is.finite(temperature))
    if (!ok) {
      .sg_stop("`temperature` must be one finite value or one per node (", n, ") in degC.")
    }
    return(rep_len(as.numeric(temperature), n))
  }
  temps <- vapply(seq_len(n), function(i) {
    t <- .sg_forcing_constant_temp(.sg_node_model(model, i)$forcing, model$params)
    if (is.null(t)) NA_real_ else t
  }, 0)
  if (anyNA(temps)) {
    .sg_stop(
      what, " needs a `temperature` (degC): the forcing of node(s) ",
      paste(ids[is.na(temps)], collapse = ", "), " does not have one constant temperature."
    )
  }
  temps
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
  h <- model$hosts
  if (!is.null(h)) {
    ids <- .sg_node_ids(model)
    for (i in seq_along(ids)) {
      idx <- unlist(lapply(.sg_host_blocks, function(b) sg_state_index(model, b, node = i)))
      total <- sum(state[idx])
      if (abs(total - h$n[i]) > 1e-8 * max(1, h$n[i])) {
        .sg_warn(sprintf(
          paste(
            "The host states%s sum to %s, not host_density * area = %s;",
            "the number of hosts stays %s."
          ),
          if (length(ids) > 1L) paste0(" of node ", ids[i]) else "",
          format(total), format(h$n[i]), format(total)
        ))
      }
    }
  }
  as.numeric(state)
}
