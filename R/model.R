#' Assemble a model
#'
#' Combines parameters, an inheritance cube, weather forcing and releases into
#' one model object with everything the solvers need: the segmented adult
#' chains of every genotype, the layout of the state vector, counts per node
#' and the release events mapped to state positions.
#'
#' @param params An `sg_params` object ([sg_params()]).
#' @param cube An `sg_cube` object ([sg_cube()]).
#' @param forcing An `sg_forcing` object ([sg_forcing()]) built with the same
#'   temperature, development and mortality parameters as `params`. Leave it
#'   `NULL` when a `landscape` is given (each node's forcing comes from it).
#' @param releases Optional `sg_releases` object ([sg_releases()]).
#' @param area Node area in km^2 (> 0; default 1, so counts equal densities
#'   numerically; D-009). Wound density `wound_density` and the mate-finding
#'   half-saturation `male_density_half_mating` (per km^2) and releases given per km^2 are
#'   converted to counts per node with it. Not used with a `landscape`, which
#'   gives each node's area.
#' @param adult_segment_substages Optional list with elements `female` and/or `male`:
#'   the number of sub-stages in each segment of that sex's adult chain,
#'   overriding the proportional allocation (D-018). Each must sum to `adult_substages`.
#'   Female segments are delimited by `female_receptive_age_dd` and
#'   `oviposition_age_dd`, male segments by `male_mature_age_dd` (zero thresholds add no
#'   segment); `print()` of a model shows the
#'   segments.
#' @param adult_females Optional number of adult females (unmated plus mated)
#'   in the node at equilibrium. If given, the node area is chosen so that the
#'   wild-type equilibrium at the forcing's constant temperature has this many
#'   adult females (see Details). Give `area` or `adult_females`, not both.
#'   Needs a forcing with one constant temperature, e.g. from
#'   [sg_constant_weather()], and a cube with a single wild type that breeds
#'   true. Not available with a `landscape`.
#' @param landscape Optional `sg_landscape` object ([sg_landscape()]): several
#'   nodes, each with its own area, forcing, wound or host density and role
#'   (MODEL_SPEC section 8a, D-032). Without it the model has one node.
#' @param movement Optional `sg_movement_layer` object ([sg_movement_layer()])
#'   or a list with one such layer: how adults move between the landscape's
#'   nodes. Needs a `landscape`. (Several layers follow in a later version.)
#'
#' @details
#' **Adult chains** (model specification S6, D-018). Each adult chain is split
#' at its maturity thresholds into segments with exact mean durations. With
#' the defaults (`adult_substages` = 10, `oviposition_age_dd` = 61,
#' `female_receptive_age_dd` = `male_mature_age_dd` = 0) the female chain
#' has 3 sub-stages before first oviposition and 7 after; the male chain is a
#' single segment. A genotype with lifespan fraction l_g keeps the same
#' thresholds in degree-days; only the segment in which its lifespan l_g `adult_lifespan_dd`
#' ends changes length. A genotype whose lifespan ends before a threshold
#' never reaches it (reported by `print()`).
#'
#' **State layout.** All states are counts per node, stored in one vector of
#' blocks E, L, P, M, U, F (in this order). Within each block the index order
#' from slowest to fastest is node, host class, genotype, mate genotype
#' (F only), sub-stage. Use [sg_state_index()] and [sg_state_names()] rather
#' than computing positions by hand.
#'
#' **Several nodes** (D-032). With a landscape of N nodes the state vector is
#' the concatenation of N per-node vectors, each with the layout above (node
#' slowest), and state names end in `"@"` and the node id, e.g.
#' `"M[AA,1]@north"`. Each node's wound count (or host count), mate-finding
#' half-saturation and per-km^2 releases use that node's area and densities.
#' Releases need a `node` ([sg_releases()]).
#'
#' **Releases** enter the first sub-stage of their stage: eggs E, pupae P,
#' adult males M and adult females U (unmated).
#'
#' **Population size from `adult_females`** (D-028). The parameters set the
#' population density (per km^2); the node area sets how many individuals
#' that is. With `adult_females`, the area is the target divided by the
#' equilibrium density of adult females at the forcing's temperature (from
#' [sg_init()] with `method = "equilibrium"`). Densities in `params`
#' (`wound_density`, `male_density_half_mating`) stay as given, so the
#' number of wounds in the node (`wound_density` times the area, or with the
#' host-wound model the number of hosts, `host_density` times the area)
#' scales with the region. Every rate of the model depends on counts only
#' through densities, so the equilibrium is exactly proportional to the
#' area, also with larval crowding or an Allee effect. `print()` shows the
#' area chosen.
#'
#' **Fingerprint.** `model$meta$hash` is an MD5 hash of all model inputs
#' (parameters, cube, forcing, releases, area and adult chains), so two models
#' with the same hash give the same simulation.
#'
#' @return An object of class `sg_model`: a list with `params`, `cube`,
#'   `forcing` (with a landscape: the forcing of the first node; all nodes
#'   share its dates), `area` (one value per node), `dims` (numbers of
#'   genotypes, sub-stages, nodes and host classes), `adult` (female and male
#'   chains), `layout` (one row per state block of a node: name, first index
#'   within the node, length), `n_states` (all nodes), `counts` (wound
#'   number `H`, half-saturation `h_M`, one value per node, and the extinction
#'   threshold; internal names;
#'   `H` is `NA` with the host-wound model), `hosts` (`NULL` with
#'   `wound_model = "static"`; otherwise the host class, host count `n` of
#'   each node, infestation sub-stages `n_I` and per-host rates per day `wounding`,
#'   `healing`, `advance`, `treatment` and `attractiveness`, D-030),
#'   `releases` (events with day, state index and count; with several nodes
#'   also the node) and `meta` (package version and hash). With a landscape
#'   also `n_states_node` (states per node), `landscape` and `movement` (the
#'   merged movement groups, `dx_ref`, `t_ref` and `kappa_max`; `NULL`
#'   without movement).
#' @seealso [sg_state_index()], [sg_state_names()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   cube <- sg_cube(MGDrivE::cubeMendelian())
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1961-01-01", end = "1962-12-31")
#'   rel <- sg_releases("aa",
#'     stage = "adult", sex = "male",
#'     times = seq(60, by = 7, length.out = 10), number = 500
#'   )
#'   mod <- sg_model(p, cube, wx, rel, area = 1)
#'   mod
#'   head(sg_state_names(mod))
#'
#'   # three nodes in a row, the first a source, with adult movement (D-032)
#'   land <- sg_landscape(
#'     data.frame(
#'       id = c("south", "middle", "north"), area = 100,
#'       role = c("source", "study", "study")
#'     ),
#'     wx
#'   )
#'   chain <- sg_movement_layer(
#'     list(male = 0.05, virgin = 0.1, mated = 0.17),
#'     rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0))
#'   )
#'   mod3 <- sg_model(p, cube, landscape = land, movement = chain)
#'   mod3
#'   head(sg_state_names(mod3))
#' }
sg_model <- function(params, cube, forcing = NULL, releases = NULL, area = 1,
                     adult_segment_substages = NULL, adult_females = NULL,
                     landscape = NULL, movement = NULL) {
  if (!inherits(params, "sg_params")) {
    .sg_stop("`params` must be an sg_params object; create one with sg_params().")
  }
  if (!inherits(cube, "sg_cube")) {
    .sg_stop("`cube` must be an sg_cube object; create one with sg_cube().")
  }
  if (!is.null(landscape)) {
    if (!inherits(landscape, "sg_landscape")) {
      .sg_stop("`landscape` must be an sg_landscape object; create one with sg_landscape().")
    }
    if (!is.null(forcing)) {
      .sg_stop("Give `forcing` or `landscape`, not both: a landscape holds each node's forcing.")
    }
    if (!missing(area)) {
      .sg_stop("Give `area` or `landscape`, not both: a landscape holds each node's area.")
    }
    if (!is.null(adult_females)) {
      .sg_stop("`adult_females` is available only for one-node models without a landscape.")
    }
    forcing <- landscape$forcing[[landscape$nodes$forcing_series[1L]]]
    area <- landscape$nodes$area
  } else if (!is.null(movement)) {
    .sg_stop(
      "`movement` needs a `landscape` (sg_landscape()) whose nodes the flies move between."
    )
  }
  if (!inherits(forcing, "sg_forcing")) {
    .sg_stop("`forcing` must be an sg_forcing object; create one with sg_forcing().")
  }
  if (!is.null(adult_females)) {
    if (!missing(area)) .sg_stop("Give `area` or `adult_females`, not both.")
    ok <- is.numeric(adult_females) && length(adult_females) == 1L &&
      is.finite(adult_females) && adult_females > 0
    if (!ok) .sg_stop("`adult_females` must be a single positive number.")
    # build at 1 km^2, find the equilibrium density, then rebuild at the area
    # that holds `adult_females` (exact: the equilibrium scales with the area)
    unit <- sg_model(params, cube, forcing, releases,
      area = 1, adult_segment_substages = adult_segment_substages
    )
    sized <- .sg_area_for_females(unit, adult_females)
    model <- sg_model(params, cube, forcing, releases,
      area = sized$area, adult_segment_substages = adult_segment_substages
    )
    model$meta$adult_females <- sized
    return(model)
  }
  area_ok <- is.numeric(area) && length(area) == 1L && is.finite(area) && area > 0
  if (is.null(landscape) && !area_ok) {
    .sg_stop("`area` must be a single positive number (node area in km^2).")
  }
  .sg_check_forcing_params(forcing, params)
  movement <- .sg_model_movement(movement, landscape, params)
  adult_segment_substages <- .sg_check_adult_substages_arg(adult_segment_substages)

  p <- params
  n_g <- cube$n_genotypes
  female_thr <- if (p$mating_mode == "instant") {
    c(X_ovi = p$oviposition_age_dd)
  } else {
    c(X_rec = p$female_receptive_age_dd, X_ovi = p$oviposition_age_dd)
  }
  adult <- list(
    female = .sg_adult_chain(female_thr, p$adult_lifespan_dd, p$adult_substages, cube$ell,
      substages = adult_segment_substages$female, sex = "female"
    ),
    male = .sg_adult_chain(c(X_mat = p$male_mature_age_dd),
      p$adult_lifespan_dd, p$adult_substages, cube$ell,
      substages = adult_segment_substages$male, sex = "male"
    )
  )

  hosts <- p$wound_model == "hosts"
  n_node <- if (is.null(landscape)) 1L else landscape$n_nodes
  dims <- list(
    n_genotypes = n_g,
    n_E = p$egg_substages, n_L = p$larva_substages, n_P = p$pupa_substages,
    n_A = p$adult_substages,
    n_I = if (hosts) p$infestation_substages else 0L,
    n_node = n_node, n_host = 1L
  )
  # blocks of one node; the state vector holds n_node such vectors (node slowest)
  layout <- .sg_layout(utils::modifyList(dims, list(n_node = 1L)), hosts)
  if (n_node > 1L) layout$dims <- sub("^node 1 x ", "per node: ", layout$dims)
  n_states <- n_node * sum(layout$length)
  dens <- .sg_node_densities(landscape, p)

  model <- list(
    params = params,
    cube = cube,
    forcing = forcing,
    area = area,
    dims = dims,
    adult = adult,
    layout = layout,
    n_states = n_states,
    counts = list(
      H = if (hosts) rep(NA_real_, n_node) else dens$wound * area,
      h_M = p$male_density_half_mating * area,
      ext_threshold = p$extinction_threshold
    ),
    hosts = if (hosts) .sg_host_setup(p, area, dens$host) else NULL,
    releases = NULL,
    meta = list(package_version = as.character(utils::packageVersion("ScrewDrivE")), hash = NA)
  )
  # fields that one-node models without a landscape do not carry, so their
  # structure and hash stay as before (D-032)
  if (!is.null(landscape)) {
    model$n_states_node <- sum(layout$length)
    model$landscape <- landscape
    model$movement <- movement
  }
  class(model) <- "sg_model"
  model$releases <- .sg_release_events(model, releases)
  model$meta$hash <- .sg_model_hash(model)
  model
}

#' Positions and names of model states
#'
#' `sg_state_index()` returns the positions in the state vector of the states
#' of one block, optionally restricted to genotypes, mate genotypes and
#' sub-stages. `sg_state_names()` returns a name for every state, e.g.
#' `"E[AA,1]"` or `"F[AA|aa,3]"` (female genotype | mate genotype, sub-stage).
#'
#' @param model An `sg_model` object.
#' @param block State block: `"E"` (eggs), `"L"` (larvae), `"P"` (pupae),
#'   `"M"` (adult males), `"U"` (unmated females) or `"F"` (mated females);
#'   with `wound_model = "hosts"` also `"HS"` (hosts without an attractive
#'   wound), `"HW"` (wounded, not infested) and `"HI"` (infested, by
#'   sub-stage). Host states are named e.g. `"HW[hosts]"`, `"HI[hosts,2]"`.
#' @param genotype,mate,substage Optional genotype names, mate genotype names
#'   (block `"F"` only) and sub-stage numbers. `NULL` selects all. Host blocks
#'   have no genotype.
#' @param node Optional node ids or numbers (models with a landscape,
#'   [sg_landscape()]). `NULL` selects all nodes.
#' @return `sg_state_index()`: integer positions (1-based), ordered with the
#'   sub-stage varying fastest and the node slowest. `sg_state_names()`:
#'   character vector of length `model$n_states`; with more than one node each
#'   name ends in `"@"` and the node id, e.g. `"M[AA,1]@north"`.
#' @seealso [sg_model()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   mod <- sg_model(
#'     p, sg_cube(MGDrivE::cubeMendelian()),
#'     sg_forcing(sg_example_weather()[1:30, ], p)
#'   )
#'   sg_state_index(mod, "F", genotype = "AA", mate = "aa", substage = 1)
#'   sg_state_names(mod)[sg_state_index(mod, "M", genotype = "Aa")]
#' }
sg_state_index <- function(model, block = c("E", "L", "P", "M", "U", "F", "HS", "HW", "HI"),
                           genotype = NULL, mate = NULL, substage = NULL, node = NULL) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object.")
  }
  block <- .sg_match_choice(block, c("E", "L", "P", "M", "U", "F", "HS", "HW", "HI"), "block")
  geno <- model$cube$genotypes
  n_sub <- .sg_block_substages(model$dims, block)
  # first position of each selected node, minus 1 (node slowest, D-032)
  base <- (.sg_node_numbers(model, node) - 1L) * .sg_n_states_node(model)
  with_nodes <- function(within) as.integer(rep(base, each = length(within)) + within)
  if (block %in% .sg_host_blocks) {
    if (is.null(model$hosts)) {
      .sg_stop("Block ", block, " exists only with wound_model = \"hosts\".")
    }
    if (!is.null(genotype) || !is.null(mate)) {
      .sg_stop("Host block ", block, " has no genotype or mate.")
    }
    si <- if (is.null(substage)) seq_len(n_sub) else substage
    if (!is.numeric(si) || any(si < 1) || any(si > n_sub) || any(si != round(si))) {
      .sg_stop("`substage` must be whole numbers between 1 and ", n_sub, " for block ", block, ".")
    }
    return(with_nodes(model$layout$start[model$layout$block == block] + si - 1L))
  }
  pick_geno <- function(x, what) {
    if (is.null(x)) {
      return(seq_along(geno))
    }
    bad <- setdiff(x, geno)
    if (length(bad) > 0L) {
      .sg_stop("Unknown ", what, ": ", paste(bad, collapse = ", "), ".", .sg_suggest(bad[1L], geno))
    }
    match(x, geno)
  }
  gi <- pick_geno(genotype, "genotype")
  if (!is.null(mate) && block != "F") {
    .sg_stop("`mate` applies only to block \"F\" (mated females).")
  }
  mi <- if (block == "F") pick_geno(mate, "mate genotype") else 1L
  si <- if (is.null(substage)) seq_len(n_sub) else substage
  if (!is.numeric(si) || any(si < 1) || any(si > n_sub) || any(si != round(si))) {
    .sg_stop("`substage` must be whole numbers between 1 and ", n_sub, " for block ", block, ".")
  }
  start <- model$layout$start[model$layout$block == block]
  n_m <- if (block == "F") length(geno) else 1L
  grid <- expand.grid(s = si, m = mi, g = gi)
  with_nodes(start + ((grid$g - 1L) * n_m + (grid$m - 1L)) * n_sub + (grid$s - 1L))
}

#' @rdname sg_state_index
#' @export
sg_state_names <- function(model) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object.")
  }
  geno <- model$cube$genotypes
  out <- character(.sg_n_states_node(model)) # names of one node's states
  for (b in model$layout$block) {
    n_sub <- .sg_block_substages(model$dims, b)
    idx <- sg_state_index(model, b, node = 1L)
    if (b %in% .sg_host_blocks) {
      cls <- model$hosts$classes
      out[idx] <- if (b == "HI") {
        paste0("HI[", cls, ",", seq_len(n_sub), "]")
      } else {
        paste0(b, "[", cls, "]")
      }
    } else if (b == "F") {
      grid <- expand.grid(s = seq_len(n_sub), m = geno, g = geno, stringsAsFactors = FALSE)
      out[idx] <- paste0("F[", grid$g, "|", grid$m, ",", grid$s, "]")
    } else {
      grid <- expand.grid(s = seq_len(n_sub), g = geno, stringsAsFactors = FALSE)
      out[idx] <- paste0(b, "[", grid$g, ",", grid$s, "]")
    }
  }
  ids <- .sg_node_ids(model)
  if (length(ids) == 1L) {
    return(out)
  }
  paste0(rep(out, times = length(ids)), "@", rep(ids, each = length(out)))
}

#' @export
print.sg_model <- function(x, ...) {
  d <- x$dims
  f <- x$forcing
  cat(
    "<sg_model> ", format(x$n_states, big.mark = ","), " states; ", d$n_genotypes,
    " genotype", if (d$n_genotypes != 1L) "s", "; ", f$n_days, " days (",
    format(f$dates[1L]), " to ", format(f$dates[f$n_days]), ")\n",
    sep = ""
  )
  cat(sprintf(
    "  sub-stages: eggs %d, larvae %d, pupae %d, adults %d per sex\n",
    d$n_E, d$n_L, d$n_P, d$n_A
  ))
  for (sx in c("female", "male")) {
    ch <- x$adult[[sx]]
    b <- ch$boundaries
    fmt <- function(v) vapply(v, function(z) format(z, digits = 6L), "")
    segs <- paste0(
      "[", fmt(b[-length(b)]), ", ", fmt(b[-1L]), ") dd: ", ch$segment_counts,
      collapse = "; "
    )
    cat("  ", sx, " chain: ", segs, "\n", sep = "")
  }
  never <- .sg_never_reached(x)
  if (length(never) > 0L) {
    cat("  never reached:\n")
    cat(paste0("    ", never), sep = "\n")
  }
  p <- x$params
  h <- x$hosts
  if (d$n_node == 1L) {
    wounds <- if (is.null(h)) {
      sprintf("wounds H = %s per node", format(x$counts$H))
    } else {
      sprintf("hosts N_H = %s per node", format(h$n))
    }
    cat(sprintf(
      "  node area %s km^2; %s; mating \"%s\"; search competitors \"%s\"\n",
      format(x$area), wounds, p$mating_mode, p$wound_competitors
    ))
  } else {
    .sg_print_nodes(x)
    cat(sprintf("  mating \"%s\"; search competitors \"%s\"\n", p$mating_mode, p$wound_competitors))
  }
  .sg_print_movement(x)
  if (!is.null(h)) {
    cat(sprintf(
      paste0(
        "  host-wound model: %.3g of hosts wounded without flies, wounds heal in %g d; ",
        "infestations %g d in %d sub-stages; infested wounds x%g attractive; ",
        "treatment %g per day\n"
      ),
      p$wounded_fraction, p$wound_healing_days, p$infestation_days, h$n_I,
      h$attractiveness, h$treatment
    ))
  }
  sized <- x$meta$adult_females
  if (!is.null(sized)) {
    cat(sprintf(
      "  area chosen for %s adult females at the %g degC equilibrium (%.4g per km^2)\n",
      format(sized$adult_females), sized$temperature, sized$density
    ))
  }
  rel <- x$releases
  if (nrow(rel) > 0L) {
    cat(sprintf(
      "  releases: %d events on %d days (days %d to %d), %s individuals\n",
      nrow(rel), length(unique(rel$day)), min(rel$day), max(rel$day),
      format(sum(rel$count), big.mark = ",", scientific = FALSE)
    ))
  } else {
    cat("  releases: none\n")
  }
  cat("  hash: ", x$meta$hash, "\n", sep = "")
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# print() lines for a model with several nodes.
.sg_print_nodes <- function(x) {
  nd <- x$landscape$nodes
  rng <- function(v) {
    r <- range(v)
    if (r[1L] == r[2L]) format(r[1L]) else paste0(format(r[1L]), "-", format(r[2L]))
  }
  n_src <- sum(nd$role == "source")
  cat(sprintf(
    "  %d nodes (%d study, %d source), total area %s km^2; %d forcing series\n",
    nrow(nd), nrow(nd) - n_src, n_src, format(sum(nd$area), big.mark = ","),
    length(x$landscape$forcing)
  ))
  cat(if (is.null(x$hosts)) {
    sprintf("  wounds H per node: %s\n", rng(x$counts$H))
  } else {
    sprintf("  hosts N_H per node: %s\n", rng(x$hosts$n))
  })
}

# print() lines for the movement groups of a model (MODEL_SPEC 8a).
.sg_print_movement <- function(x) {
  mv <- x$movement
  if (is.null(mv)) {
    return(invisible(NULL))
  }
  for (g in mv$groups) {
    lost <- vapply(seq_len(ncol(g$loss)), function(i) any(g$loss[, i] > 0), TRUE)
    cat(sprintf(
      paste0(
        "  movement (%s clock): %d edges, %d node%s with open edges; ",
        "leave rates per day at %g degC: %s\n"
      ),
      g$clock, length(g$col), sum(lost), if (sum(lost) != 1L) "s" else "", mv$t_ref,
      paste0(rownames(g$rate), " ", vapply(seq_len(3L), function(s) {
        r <- range(g$rate[s, ])
        if (r[1L] == r[2L]) {
          format(r[1L], digits = 4L)
        } else {
          paste0(format(r[1L], digits = 4L), "-", format(r[2L], digits = 4L))
        }
      }, ""), collapse = ", ")
    ))
  }
  if (any(vapply(mv$groups, function(g) g$clock == "physiological", TRUE))) {
    cat(sprintf(
      "  physiological clock: kappa = min(%s, dx / %.4g)\n",
      format(mv$kappa_max), mv$dx_ref
    ))
  }
  invisible(NULL)
}

.sg_check_forcing_params <- function(forcing, params) {
  nms <- names(forcing$param_values)
  differ <- nms[!vapply(nms, function(n) identical(forcing$param_values[[n]], params[[n]]), TRUE)]
  if (length(differ) > 0L) {
    detail <- vapply(differ, function(n) {
      paste0(
        n, " (forcing ", .sg_fmt(forcing$param_values[[n]]),
        ", params ", .sg_fmt(params[[n]]), ")"
      )
    }, "")
    .sg_stop(
      "The forcing was built with different parameter values: ", paste(detail, collapse = "; "),
      ". Rebuild it with sg_forcing(weather, params)."
    )
  }
  invisible(TRUE)
}

.sg_check_adult_substages_arg <- function(x) {
  if (is.null(x)) {
    return(list())
  }
  if (!is.list(x) || is.null(names(x)) || !all(names(x) %in% c("female", "male"))) {
    .sg_stop("`adult_segment_substages` must be a list with elements `female` and/or `male`.")
  }
  x
}

# Host-wound blocks (D-030), after the fly blocks in host mode.
.sg_host_blocks <- c("HS", "HW", "HI")

.sg_block_substages <- function(dims, block) {
  if (block %in% .sg_host_blocks) {
    return(if (block == "HI") as.integer(dims$n_I) else 1L)
  }
  # EXPR is named so that the label E is not partially matched to it
  switch(EXPR = block,
    E = dims$n_E,
    L = dims$n_L,
    P = dims$n_P,
    dims$n_A
  )
}

# One row per block: name, first position (1-based), length, dimensions.
# With host = TRUE the host blocks HS, HW, HI follow the fly blocks (D-030).
.sg_layout <- function(dims, hosts = FALSE) {
  blocks <- c("E", "L", "P", "M", "U", "F", if (hosts) .sg_host_blocks)
  n_g <- dims$n_genotypes
  len <- vapply(blocks, function(b) {
    if (b %in% .sg_host_blocks) {
      return(as.integer(dims$n_node * dims$n_host * .sg_block_substages(dims, b)))
    }
    n_m <- if (b == "F") n_g else 1L
    as.integer(dims$n_node * dims$n_host * n_g * n_m * .sg_block_substages(dims, b))
  }, 1L)
  start <- cumsum(c(1L, len[-length(len)]))
  dims_txt <- vapply(blocks, function(b) {
    n_sub <- .sg_block_substages(dims, b)
    if (b %in% .sg_host_blocks) {
      sprintf("node %d x host %d x sub-stage %d", dims$n_node, dims$n_host, n_sub)
    } else if (b == "F") {
      sprintf(
        "node %d x host %d x genotype %d x mate %d x sub-stage %d",
        dims$n_node, dims$n_host, n_g, n_g, n_sub
      )
    } else {
      sprintf(
        "node %d x host %d x genotype %d x sub-stage %d",
        dims$n_node, dims$n_host, n_g, n_sub
      )
    }
  }, "")
  data.frame(
    block = blocks, start = as.integer(start), length = unname(len),
    dims = unname(dims_txt), stringsAsFactors = FALSE
  )
}

# Host-wound model (D-030, design note section 3.2): host count of each node
# (host density per km^2 times area) and per-host rates per day. Wounding is
# set so that, without flies, a fraction wounded_fraction of hosts carries an
# attractive wound: omega = gamma_W p0 / (1 - p0).
.sg_host_setup <- function(p, area, host_density = p$host_density) {
  healing <- 1 / p$wound_healing_days
  list(
    classes = "hosts",
    n = host_density * area,
    n_I = p$infestation_substages,
    wounding = healing * p$wounded_fraction / (1 - p$wounded_fraction),
    healing = healing,
    advance = p$infestation_substages / p$infestation_days,
    treatment = p$treatment_rate,
    attractiveness = p$infested_wound_attractiveness
  )
}

# Wound and host densities per km^2 of each node: the landscape's values where
# given, otherwise those of the parameters (D-032).
.sg_node_densities <- function(landscape, p) {
  if (is.null(landscape)) {
    return(list(wound = p$wound_density, host = p$host_density))
  }
  nd <- landscape$nodes
  hosts <- p$wound_model == "hosts"
  if (hosts && any(!is.na(nd$wound_density))) {
    .sg_warn(
      "The landscape's wound_density has no effect with wound_model = \"hosts\"; ",
      "wounds come from host_density and wounded_fraction."
    )
  }
  if (!hosts && any(!is.na(nd$host_density))) {
    .sg_warn("The landscape's host_density has no effect unless wound_model = \"hosts\".")
  }
  list(
    wound = ifelse(is.na(nd$wound_density), p$wound_density, nd$wound_density),
    host = ifelse(is.na(nd$host_density), p$host_density, nd$host_density)
  )
}

# Movement of a model (MODEL_SPEC 8a, N1-N3): the layers merged into groups,
# the reference temperature, dx at that temperature (S1, with its floor) and
# the cap of the physiological clock. NULL without movement.
.sg_model_movement <- function(movement, landscape, p) {
  if (is.null(movement)) {
    return(NULL)
  }
  layers <- if (inherits(movement, "sg_movement_layer")) list(movement) else movement
  ok <- is.list(layers) && length(layers) >= 1L &&
    all(vapply(layers, inherits, TRUE, "sg_movement_layer"))
  if (!ok) {
    .sg_stop(
      "`movement` must be an sg_movement_layer object (or a list of them); ",
      "create one with sg_movement_layer()."
    )
  }
  if (length(layers) > 1L) {
    .sg_stop(
      "This version accepts one movement layer per model; several layers (e.g. local ",
      "and long-distance movement) follow in a later version (design note phase 4)."
    )
  }
  list(
    groups = .sg_movement_groups(layers, landscape$nodes$id),
    t_ref = p$movement_t_ref,
    dx_ref = sg_degree_days(p$movement_t_ref, p),
    kappa_max = p$movement_kappa_max
  )
}

# Stop for functions that do not (yet) handle several nodes or movement.
.sg_one_node_only <- function(model, what, when = "in a later version") {
  if (model$dims$n_node > 1L || !is.null(model$movement)) {
    .sg_stop(
      what, " is not available for models with several nodes or movement yet (", when, ")."
    )
  }
  invisible(TRUE)
}

# Node i of a model as a one-node model (D-032): the node's forcing series,
# area, wounds H, mate finding h_M and host count, without landscape,
# movement or releases. Fields are replaced, not recomputed, so the numbers
# are exactly those the multi-node model uses. A model without a landscape
# is returned as it is.
.sg_node_model <- function(model, i) {
  land <- model$landscape
  if (is.null(land)) {
    return(model)
  }
  m <- model
  m$forcing <- land$forcing[[land$nodes$forcing_series[i]]]
  m$area <- model$area[i]
  m$dims$n_node <- 1L
  m$layout$dims <- sub("^per node: ", "node 1 x ", m$layout$dims)
  m$n_states <- model$n_states_node
  m$counts$H <- model$counts$H[i]
  m$counts$h_M <- model$counts$h_M[i]
  if (!is.null(m$hosts)) m$hosts$n <- model$hosts$n[i]
  m$releases <- model$releases[0L, c("day", "index", "count", "genotype", "stage", "sex"),
    drop = FALSE
  ]
  m$n_states_node <- NULL
  m$landscape <- NULL
  m$movement <- NULL
  m$meta$hash <- .sg_model_hash(m)
  m
}

# Positions of node i's states in the full state vector.
.sg_node_positions <- function(model, i) {
  n_s <- .sg_n_states_node(model)
  (i - 1L) * n_s + seq_len(n_s)
}

# Node ids ("1" for a model without a landscape) and states per node.
.sg_node_ids <- function(model) {
  if (is.null(model$landscape)) "1" else model$landscape$nodes$id
}
.sg_n_states_node <- function(model) {
  if (is.null(model$n_states_node)) model$n_states else model$n_states_node
}

# Node numbers (1-based) from ids or numbers; NULL means all nodes.
.sg_node_numbers <- function(model, node) {
  ids <- .sg_node_ids(model)
  if (is.null(node)) {
    return(seq_along(ids))
  }
  if (is.numeric(node)) {
    bad <- node < 1 | node > length(ids) | node != round(node) | is.na(node)
    if (any(bad)) {
      .sg_stop("`node` numbers must be whole numbers between 1 and ", length(ids), ".")
    }
    return(as.integer(node))
  }
  idx <- match(as.character(node), ids)
  if (anyNA(idx)) {
    bad <- as.character(node)[is.na(idx)][1L]
    .sg_stop(
      "Unknown node \"", bad, "\".", .sg_suggest(bad, ids),
      " Nodes: ", .sg_wrap_names(ids), "."
    )
  }
  idx
}

# Release events with the state position they are added to (first sub-stage).
# With several nodes the table also has the node id (D-032).
.sg_release_events <- function(model, releases) {
  ids <- .sg_node_ids(model)
  multi <- length(ids) > 1L
  ev <- .sg_resolve_releases(
    releases, model$forcing$dates, model$area, model$cube$genotypes, ids
  )
  cols <- c("day", "index", "count", "genotype", "stage", "sex", if (multi) "node")
  if (nrow(ev) == 0L) {
    ev$index <- integer(0)
    if (multi) ev$node <- character(0)
    return(ev[, cols])
  }
  block <- ifelse(ev$stage == "egg", "E", ifelse(ev$stage == "pupa", "P",
    ifelse(ev$sex == "male", "M", "U")
  ))
  node <- if (multi) ev$node else rep(1L, nrow(ev))
  ev$index <- vapply(seq_len(nrow(ev)), function(i) {
    sg_state_index(model, block[i], genotype = ev$genotype[i], substage = 1L, node = node[i])
  }, 1L)
  if (multi) ev$node <- ids[ev$node]
  key <- paste(ev$day, ev$index)
  agg <- ev[!duplicated(key), setdiff(cols, "count")]
  agg$count <- as.numeric(tapply(ev$count, factor(key, levels = unique(key)), sum))
  agg <- agg[order(agg$day, agg$index), cols]
  rownames(agg) <- NULL
  agg
}

# Genotypes that never reach a threshold of their chain.
.sg_never_reached <- function(model) {
  out <- character(0)
  for (sx in c("female", "male")) {
    ch <- model$adult[[sx]]
    for (thr in names(ch$flags)) {
      first <- which(ch$flags[[thr]])[1L]
      gs <- names(ch$last_class)[ch$last_class < first]
      if (length(gs) > 0L) {
        what <- switch(thr,
          X_ovi = "first oviposition (never lay eggs)",
          X_rec = "female receptivity (never mate)",
          X_mat = "male maturity (never available as mates)"
        )
        out <- c(out, paste0(paste(gs, collapse = ", "), ": ", what))
      }
    }
  }
  out
}

# Parameters that act only through movement (MODEL_SPEC N3, D-032).
.sg_movement_param_names <- c("movement_t_ref", "movement_kappa_max")

# MD5 hash of the model inputs, computed from a canonical text representation
# (independent of the R version, unlike serialize()). Movement parameters
# enter only when the model has movement, so models without movement keep the
# hashes they had before movement existed (D-032).
.sg_model_hash <- function(model) {
  f <- model$forcing
  params <- unclass(model$params)
  if (is.null(model$movement)) params[.sg_movement_param_names] <- NULL # keeps attribute preset
  inputs <- list(
    params = params,
    cube = model$cube[c(
      "genotypes", "wild_type", "release_type", "eta", "phi", "xiF", "xiM", "s",
      "ell", "egg_tensor"
    )],
    forcing = list(
      dates = as.numeric(f$dates), dx = f$dx, dx_egg = f$dx_egg, mu_T = f$mu_T,
      mu_E = f$mu_E, phi_ov = f$phi_ov, moisture = f$moisture
    ),
    area = model$area,
    adult = lapply(model$adult, function(ch) ch[c("segment_counts", "width", "last_class")]),
    releases = model$releases
  )
  # landscape and movement enter only when present, so one-node hashes are
  # unchanged (D-032)
  if (!is.null(model$landscape)) {
    inputs$landscape <- list(
      nodes = model$landscape$nodes,
      forcing = lapply(model$landscape$forcing, function(f) {
        list(
          dx = f$dx, dx_egg = f$dx_egg, mu_T = f$mu_T, mu_E = f$mu_E,
          phi_ov = f$phi_ov, moisture = f$moisture
        )
      })
    )
  }
  if (!is.null(model$movement)) inputs$movement <- model$movement
  txt <- paste(deparse(inputs, control = c("keepNA", "keepInteger", "digits17")),
    collapse = "\n"
  )
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(txt, tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}
