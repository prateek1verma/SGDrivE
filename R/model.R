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
#'   temperature, development and mortality parameters as `params`.
#' @param releases Optional `sg_releases` object ([sg_releases()]).
#' @param area Node area in km^2 (> 0; default 1, so counts equal densities
#'   numerically; D-009). Wound density `wound_density` and the mate-finding
#'   half-saturation `male_density_half_mating` (per km^2) and releases given per km^2 are
#'   converted to counts per node with it.
#' @param adult_segment_substages Optional list with elements `female` and/or `male`:
#'   the number of sub-stages in each segment of that sex's adult chain,
#'   overriding the proportional allocation (D-018). Each must sum to `adult_substages`.
#'   Female segments are delimited by `female_receptive_age_dd` and
#'   `oviposition_age_dd`, male segments by `male_mature_age_dd` (zero thresholds add no
#'   segment); `print()` of a model shows the
#'   segments.
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
#' **Releases** enter the first sub-stage of their stage: eggs E, pupae P,
#' adult males M and adult females U (unmated).
#'
#' **Fingerprint.** `model$meta$hash` is an MD5 hash of all model inputs
#' (parameters, cube, forcing, releases, area and adult chains), so two models
#' with the same hash give the same simulation.
#'
#' @return An object of class `sg_model`: a list with `params`, `cube`,
#'   `forcing`, `area`, `dims` (numbers of genotypes, sub-stages, nodes and
#'   host classes), `adult` (female and male chains), `layout` (one row per
#'   state block: name, first index, length), `n_states`, `counts` (wound
#'   number `H`, half-saturation `h_M`, extinction threshold, per node; internal names),
#'   `releases` (events with day, state index and count) and `meta`
#'   (package version and hash).
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
#' }
sg_model <- function(params, cube, forcing, releases = NULL, area = 1,
                     adult_segment_substages = NULL) {
  if (!inherits(params, "sg_params")) {
    .sg_stop("`params` must be an sg_params object; create one with sg_params().")
  }
  if (!inherits(cube, "sg_cube")) {
    .sg_stop("`cube` must be an sg_cube object; create one with sg_cube().")
  }
  if (!inherits(forcing, "sg_forcing")) {
    .sg_stop("`forcing` must be an sg_forcing object; create one with sg_forcing().")
  }
  area_ok <- is.numeric(area) && length(area) == 1L && is.finite(area) && area > 0
  if (!area_ok) {
    .sg_stop("`area` must be a single positive number (node area in km^2).")
  }
  .sg_check_forcing_params(forcing, params)
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

  dims <- list(
    n_genotypes = n_g,
    n_E = p$egg_substages, n_L = p$larva_substages, n_P = p$pupa_substages,
    n_A = p$adult_substages,
    n_node = 1L, n_host = 1L
  )
  layout <- .sg_layout(dims)
  n_states <- sum(layout$length)

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
      H = p$wound_density * area, h_M = p$male_density_half_mating * area,
      ext_threshold = p$extinction_threshold
    ),
    releases = NULL,
    meta = list(package_version = as.character(utils::packageVersion("SGDrivE")), hash = NA)
  )
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
#'   `"M"` (adult males), `"U"` (unmated females) or `"F"` (mated females).
#' @param genotype,mate,substage Optional genotype names, mate genotype names
#'   (block `"F"` only) and sub-stage numbers. `NULL` selects all.
#' @return `sg_state_index()`: integer positions (1-based), ordered with the
#'   sub-stage varying fastest. `sg_state_names()`: character vector of length
#'   `model$n_states`.
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
sg_state_index <- function(model, block = c("E", "L", "P", "M", "U", "F"),
                           genotype = NULL, mate = NULL, substage = NULL) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object.")
  }
  block <- .sg_match_choice(block, c("E", "L", "P", "M", "U", "F"), "block")
  geno <- model$cube$genotypes
  n_sub <- .sg_block_substages(model$dims, block)
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
  as.integer(start + ((grid$g - 1L) * n_m + (grid$m - 1L)) * n_sub + (grid$s - 1L))
}

#' @rdname sg_state_index
#' @export
sg_state_names <- function(model) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object.")
  }
  geno <- model$cube$genotypes
  out <- character(model$n_states)
  for (b in model$layout$block) {
    n_sub <- .sg_block_substages(model$dims, b)
    idx <- sg_state_index(model, b)
    if (b == "F") {
      grid <- expand.grid(s = seq_len(n_sub), m = geno, g = geno, stringsAsFactors = FALSE)
      out[idx] <- paste0("F[", grid$g, "|", grid$m, ",", grid$s, "]")
    } else {
      grid <- expand.grid(s = seq_len(n_sub), g = geno, stringsAsFactors = FALSE)
      out[idx] <- paste0(b, "[", grid$g, ",", grid$s, "]")
    }
  }
  out
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
  cat(sprintf(
    "  node area %s km^2; wounds H = %s per node; mating \"%s\"; search competitors \"%s\"\n",
    format(x$area), format(x$counts$H), p$mating_mode, p$wound_competitors
  ))
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

.sg_block_substages <- function(dims, block) {
  # EXPR is named so that the label E is not partially matched to it
  switch(EXPR = block,
    E = dims$n_E,
    L = dims$n_L,
    P = dims$n_P,
    dims$n_A
  )
}

# One row per block: name, first position (1-based), length, dimensions.
.sg_layout <- function(dims) {
  blocks <- c("E", "L", "P", "M", "U", "F")
  n_g <- dims$n_genotypes
  len <- vapply(blocks, function(b) {
    n_m <- if (b == "F") n_g else 1L
    as.integer(dims$n_node * dims$n_host * n_g * n_m * .sg_block_substages(dims, b))
  }, 1L)
  start <- cumsum(c(1L, len[-length(len)]))
  dims_txt <- vapply(blocks, function(b) {
    n_sub <- .sg_block_substages(dims, b)
    if (b == "F") {
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

# Release events with the state position they are added to (first sub-stage).
.sg_release_events <- function(model, releases) {
  ev <- .sg_resolve_releases(releases, model$forcing$dates, model$area, model$cube$genotypes)
  if (nrow(ev) == 0L) {
    ev$index <- integer(0)
    return(ev[, c("day", "index", "count", "genotype", "stage", "sex")])
  }
  block <- ifelse(ev$stage == "egg", "E", ifelse(ev$stage == "pupa", "P",
    ifelse(ev$sex == "male", "M", "U")
  ))
  ev$index <- vapply(seq_len(nrow(ev)), function(i) {
    sg_state_index(model, block[i], genotype = ev$genotype[i], substage = 1L)
  }, 1L)
  key <- paste(ev$day, ev$index)
  agg <- ev[!duplicated(key), c("day", "index", "genotype", "stage", "sex")]
  agg$count <- as.numeric(tapply(ev$count, factor(key, levels = unique(key)), sum))
  agg <- agg[order(agg$day, agg$index), c("day", "index", "count", "genotype", "stage", "sex")]
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

# MD5 hash of the model inputs, computed from a canonical text representation
# (independent of the R version, unlike serialize()).
.sg_model_hash <- function(model) {
  f <- model$forcing
  inputs <- list(
    params = unclass(model$params),
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
  txt <- paste(deparse(inputs, control = c("keepNA", "keepInteger", "digits17")),
    collapse = "\n"
  )
  tmp <- tempfile()
  on.exit(unlink(tmp), add = TRUE)
  writeLines(txt, tmp, useBytes = TRUE)
  unname(tools::md5sum(tmp))
}
