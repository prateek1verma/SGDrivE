# Output tables for deterministic (sg_sim) and stochastic (sg_sim_cb)
# simulations (D-026): long tables, totals by stage, sex and genotype, allele
# frequencies, female classes and the case proxy.

#' Totals by stage, sex and genotype
#'
#' Adds up a simulation's counts by any combination of life stage, sex,
#' genotype and (for mated females) mate genotype.
#'
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @param by Grouping variables, any of `"stage"`, `"sex"`, `"genotype"` and
#'   `"mate"`; `character(0)` gives the grand total.
#' @param stages Stages to include, any of `"eggs"`, `"larvae"`, `"pupae"`,
#'   `"males"`, `"unmated_females"`, `"mated_females"`, or the shortcuts
#'   `"adults"` (the last three), `"juveniles"` (the first three) and `"all"`.
#' @param per `"node"` (counts per node) or `"km2"` (divided by the node area).
#'
#' @details Stages are the model's blocks E, L, P, M, U and F; juveniles have
#'   sex `"unsexed"`. Grouping by `"mate"` needs the full state, so it is not
#'   available for stochastic runs with `keep = "totals"`; the mate of
#'   individuals other than mated females is `NA`.
#' @return A data frame with columns `day`, `date`, `replicate` (stochastic
#'   runs only), the `by` variables and `value`.
#' @seealso [sg_allele_freq()], [sg_females()], [sg_cases()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   head(sg_totals(sim))
#'   tail(sg_totals(sim, by = "sex", stages = "adults", per = "km2"))
#' }
sg_totals <- function(sim, by = c("stage", "genotype"), stages = "all",
                      per = c("node", "km2")) {
  per <- .sg_match_choice(per, c("node", "km2"), "per")
  by <- .sg_check_by(by)
  st <- .sg_stage_array(sim, need_mate = "mate" %in% by)
  keep <- st$cols$stage %in% .sg_expand_stages(stages)
  cols <- st$cols[keep, , drop = FALSE]
  a <- st$array[, keep, , drop = FALSE]
  key <- if (length(by)) do.call(paste, c(cols[by], sep = "\r")) else rep("all", nrow(cols))
  groups <- unique(key)
  agg <- outer(key, groups, `==`) * 1
  values <- .sg_apply_cols(a, agg)
  if (per == "km2") values <- values / .sg_sim_info(sim)$area
  levels_df <- cols[match(groups, key), by, drop = FALSE]
  .sg_long_table(sim, values, levels_df)
}

#' @rdname sg_totals
#' @param x A simulation from [simulate_ode()] or [simulate_cb()].
#' @param row.names,optional Not used.
#' @param ... Not used.
#' @export
# nolint start: object_name_linter. (row.names is the generic's argument)
as.data.frame.sg_sim <- function(x, row.names = NULL, optional = FALSE, ..., per = "node") {
  sg_totals(x, by = c("stage", "sex", "genotype", "mate"), per = per)
}

#' @rdname sg_totals
#' @export
as.data.frame.sg_sim_cb <- function(x, row.names = NULL, optional = FALSE, ...,
                                    per = "node") {
  by <- c("stage", "sex", "genotype", if (!is.null(x$state)) "mate")
  sg_totals(x, by = by, per = per)
}
# nolint end

#' Allele frequencies
#'
#' Frequencies of alleles among the individuals of chosen stages.
#'
#' @inheritParams sg_totals
#' @param stages Stages whose individuals are counted (see [sg_totals()]);
#'   default adults.
#' @param alleles Optional matrix of allele copies per genotype (rows named by
#'   genotype, columns by allele). By default each genotype name is split
#'   into single-letter alleles, as in MGDrivE cubes (`"Aa"` has one `A` and
#'   one `a`; `"XY"` one `X` and one `Y`).
#' @return A data frame with columns `day`, `date`, `replicate` (stochastic
#'   runs only), `allele`, `count` (allele copies) and `freq` (share of all
#'   copies; `NA` when no individuals are present).
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   sim <- simulate_ode(mod, init = sg_init(mod, "hardy_weinberg", alleles = c(A = 0.8, a = 0.2)))
#'   tail(sg_allele_freq(sim))
#' }
sg_allele_freq <- function(sim, stages = "adults", alleles = NULL) {
  geno <- .sg_sim_info(sim)$genotypes
  amat <- if (is.null(alleles)) {
    .sg_default_alleles(geno)
  } else {
    .sg_check_allele_matrix(alleles, geno)[geno, , drop = FALSE]
  }
  st <- .sg_stage_array(sim)
  keep <- st$cols$stage %in% .sg_expand_stages(stages)
  # each column (stage, genotype, mate) contributes its genotype's allele copies
  copies <- amat[st$cols$genotype[keep], , drop = FALSE]
  counts <- .sg_apply_cols(st$array[, keep, , drop = FALSE], copies)
  total <- apply(counts, c(1L, 3L), sum)
  freq <- sweep(counts, c(1L, 3L), total, "/")
  freq[!is.finite(freq)] <- NA_real_
  lev <- data.frame(allele = colnames(amat), stringsAsFactors = FALSE)
  out <- .sg_long_table(sim, counts, lev)
  names(out)[names(out) == "value"] <- "count"
  out$freq <- .sg_long_table(sim, freq, lev)$value
  out
}

#' Female classes
#'
#' Adult female numbers by class: the quantities that matter for genetic
#' control.
#'
#' @inheritParams sg_totals
#' @return A data frame with columns `day`, `date`, `replicate` (stochastic
#'   runs only) and
#'   \describe{
#'     \item{`females`}{All adult females (unmated and mated).}
#'     \item{`unmated`, `mated`}{Unmated and mated females.}
#'     \item{`wild_type`}{Females of the wild-type genotype(s).}
#'     \item{`egg_laying`}{Mated females old enough to lay eggs (`oviposition_age_dd`).}
#'     \item{`fertile_equivalent`}{Mated females weighted by the viability of
#'       their cross (viable eggs per egg relative to wild type x wild type):
#'       a female mated to a sterile male counts 0, one with half the viable
#'       eggs counts 0.5.}
#'   }
#'   `egg_laying` and `fertile_equivalent` need the full state, so they are
#'   `NA` for stochastic runs with `keep = "totals"`.
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   tail(sg_females(sim))
#' }
sg_females <- function(sim, per = c("node", "km2")) {
  per <- .sg_match_choice(per, c("node", "km2"), "per")
  info <- .sg_sim_info(sim)
  a <- .sg_sim_array(sim)
  nm <- dimnames(a$array)[[2L]]
  block <- substr(nm, 1L, 1L)
  geno <- sub("^.\\[([^],|]+).*$", "\\1", nm)
  # weighted column sums as one vector, time fastest then replicate (the
  # order of .sg_time_columns())
  col_sum <- function(w) {
    as.vector(.sg_apply_cols(a$array, matrix(as.numeric(w), ncol = 1L))[, 1L, ])
  }
  unmated <- col_sum(block == "U")
  mated <- col_sum(block == "F")
  out <- data.frame(
    females = unmated + mated, unmated = unmated, mated = mated,
    wild_type = col_sum(block %in% c("U", "F") & geno %in% info$wild_type)
  )
  if (a$level == "state") {
    sub <- as.integer(sub("^.*,([0-9]+)\\]$", "\\1", nm))
    is_f <- block == "F"
    out$egg_laying <- col_sum(is_f & info$laying[pmax(sub, 1L)])
    mate <- sub("^F\\[[^|]+\\|([^,]+),.*$", "\\1", nm)
    pair <- paste0(geno, "|", mate)
    wild_pairs <- outer(info$wild_type, info$wild_type, paste, sep = "|")
    v_ref <- max(info$viability[wild_pairs], na.rm = TRUE)
    weight <- ifelse(is_f, info$viability[pair] / v_ref, 0)
    weight[is.na(weight)] <- 0
    out$fertile_equivalent <- col_sum(weight)
  } else {
    out$egg_laying <- NA_real_
    out$fertile_equivalent <- NA_real_
  }
  if (per == "km2") out[] <- lapply(out, function(v) v / info$area)
  cbind(.sg_time_columns(sim), out)
}

#' Case proxy: wounds newly infested with viable eggs
#'
#' Myiasis cases by period: oviposition events (egg batches) that contain at
#' least one viable egg (D-026). Eggs of females mated to sterile males cannot
#' hatch, so their oviposition is not a case.
#'
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @param period `"day"`, `"week"` (seven-day blocks from the first day),
#'   `"month"` or `"year"` (calendar).
#' @param per `"km2"` (per km^2, default) or `"node"`.
#' @param host_density Optional density of hosts (animals per km^2); adds
#'   cases per host.
#' @details Cases come from the `fertile_batches` counter: a batch of C eggs
#'   (the clutch size, R phi_ov by default) laid by a female of genotype g
#'   mated to a male of genotype m holds at least one viable egg with
#'   probability 1 - exp(-C v), where v is the viable fraction of the cross's
#'   eggs (sum over offspring of s ih tau). In the ODE this gives the daily
#'   rate; in stochastic runs batches are drawn and those with viable eggs
#'   counted (whole numbers). Cases are counted at oviposition; viable eggs
#'   that later die of temperature are not removed.
#' @return A data frame with columns `period` (first day), `days`,
#'   `replicate` (stochastic runs only), `cases` (per km^2 or node), `events`
#'   (all oviposition events, also those with no viable egg),
#'   `wound_fraction` (cases per available wound and day: the mean daily
#'   fraction of wounds newly infested with viable eggs) and, with
#'   `host_density`, `cases_per_host`.
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   sg_cases(sim, period = "month", host_density = 20)
#' }
sg_cases <- function(sim, period = c("day", "week", "month", "year"), per = c("km2", "node"),
                     host_density = NULL) {
  period <- .sg_match_choice(period, c("day", "week", "month", "year"), "period")
  per <- .sg_match_choice(per, c("km2", "node"), "per")
  if (!is.null(host_density)) {
    ok <- is.numeric(host_density) && length(host_density) == 1L && is.finite(host_density) &&
      host_density > 0
    if (!ok) .sg_stop("`host_density` must be a single positive number (hosts per km^2).")
  }
  info <- .sg_sim_info(sim)
  k <- sim$counters
  if (!"fertile_batches" %in% dimnames(k)[[2L]]) {
    .sg_stop("This simulation was made by an older version of SGDrivE; run it again.")
  }
  column <- function(name) {
    v <- if (length(dim(k)) == 3L) k[, name, ] else k[, name]
    matrix(v, nrow = dim(k)[1L])
  }
  dates <- sim$date[-1L]
  key <- switch(EXPR = period,
    day = format(dates),
    week = format(dates[1L] + 7L * ((seq_along(dates) - 1L) %/% 7L)),
    month = format(dates, "%Y-%m-01"),
    year = format(dates, "%Y-01-01")
  )
  groups <- unique(key)
  agg <- outer(key, groups, `==`) * 1
  cases <- t(agg) %*% column("fertile_batches") # [period, replicate], per node
  events <- t(agg) %*% column("egg_batches")
  days <- colSums(agg)
  wounds <- info$H * info$area # wounds per node
  scale <- if (per == "km2") 1 / info$area else 1
  out <- lapply(seq_len(ncol(cases)), function(r) {
    d <- data.frame(
      period = as.Date(dates[match(groups, key)]), days = as.integer(days),
      cases = cases[, r] * scale, events = events[, r] * scale,
      wound_fraction = if (wounds > 0) cases[, r] / (wounds * days) else NA_real_
    )
    if (!is.null(host_density)) d$cases_per_host <- cases[, r] / info$area / host_density
    if (inherits(sim, "sg_sim_cb")) d <- cbind(d[1:2], replicate = r, d[-(1:2)])
    d
  })
  out <- do.call(rbind, out)
  rownames(out) <- NULL
  out
}

# ---- internal ------------------------------------------------------------------

# Output information stored in a simulation's metadata (D-026).
.sg_output_info <- function(model) {
  cb <- model$cube
  laying <- model$adult$female$flags$X_ovi
  if (is.null(laying)) laying <- rep(TRUE, model$dims$n_A)
  list(
    area = model$area, H = model$params$wound_density, genotypes = cb$genotypes,
    wild_type = cb$wild_type, laying = as.logical(laying), viability = cb$egg_total
  )
}

.sg_sim_info <- function(sim) {
  if (!inherits(sim, c("sg_sim", "sg_sim_cb"))) {
    .sg_stop("`sim` must be a simulation from simulate_ode() or simulate_cb().")
  }
  info <- sim$meta$info
  if (is.null(info)) {
    .sg_stop("This simulation was made by an older version of SGDrivE; run it again.")
  }
  info
}

# The simulation's counts as an array [time, column, replicate] and whether
# the columns are states or block x genotype totals.
.sg_sim_array <- function(sim) {
  .sg_sim_info(sim)
  if (inherits(sim, "sg_sim")) {
    a <- sim$state
    dim(a) <- c(dim(a), 1L)
    dimnames(a) <- list(NULL, colnames(sim$state), NULL)
    return(list(array = a, level = "state"))
  }
  if (!is.null(sim$state)) {
    list(array = sim$state, level = "state")
  } else {
    list(array = sim$totals, level = "totals")
  }
}

# Apply a column-aggregation matrix to every time x replicate slice.
.sg_apply_cols <- function(a, agg) {
  d <- dim(a)
  out <- vapply(seq_len(d[3L]), function(r) {
    m <- matrix(a[, , r], d[1L], d[2L])
    m %*% agg
  }, matrix(0, d[1L], ncol(agg)))
  array(out, c(d[1L], ncol(agg), d[3L]))
}

# Counts aggregated over sub-stages: array [time, (stage, genotype, mate),
# replicate] with a descriptor of the columns.
.sg_stage_array <- function(sim, need_mate = FALSE) {
  a <- .sg_sim_array(sim)
  nm <- dimnames(a$array)[[2L]]
  if (a$level == "totals" && need_mate) {
    .sg_stop("Mate genotypes need the full state; run simulate_cb() with keep = \"state\".")
  }
  block <- substr(nm, 1L, 1L)
  # names are "E[g,1]", "F[g|m,1]" (states) or "E[g]" (totals); a "]" first
  # in a bracket expression is literal
  geno <- sub("^.\\[([^],|]+).*$", "\\1", nm)
  mate <- ifelse(block == "F" & grepl("|", nm, fixed = TRUE),
    sub("^F\\[[^|]+\\|([^],]+).*$", "\\1", nm), NA_character_
  )
  key <- paste(block, geno, mate)
  groups <- unique(key)
  agg <- outer(key, groups, `==`) * 1
  first <- match(groups, key)
  stage_of <- c(
    E = "eggs", L = "larvae", P = "pupae", M = "males", U = "unmated_females",
    F = "mated_females"
  )
  sex_of <- c(E = "unsexed", L = "unsexed", P = "unsexed", M = "male", U = "female", F = "female")
  cols <- data.frame(
    stage = unname(stage_of[block[first]]), sex = unname(sex_of[block[first]]),
    genotype = geno[first], mate = mate[first], stringsAsFactors = FALSE
  )
  list(array = .sg_apply_cols(a$array, agg), cols = cols)
}

.sg_stage_names <- c("eggs", "larvae", "pupae", "males", "unmated_females", "mated_females")

.sg_expand_stages <- function(stages) {
  short <- list(
    all = .sg_stage_names, adults = .sg_stage_names[4:6], juveniles = .sg_stage_names[1:3]
  )
  out <- unique(unlist(lapply(stages, function(s) if (s %in% names(short)) short[[s]] else s)))
  bad <- setdiff(out, .sg_stage_names)
  if (!is.character(stages) || length(bad) > 0L) {
    .sg_stop(
      "Unknown stage(s): ", paste(bad, collapse = ", "), ". Use ",
      paste(c(.sg_stage_names, names(short)), collapse = ", "), "."
    )
  }
  out
}

.sg_check_by <- function(by) {
  allowed <- c("stage", "sex", "genotype", "mate")
  if (!is.character(by) || any(!by %in% allowed)) {
    .sg_stop("`by` must contain only ", paste0("\"", allowed, "\"", collapse = ", "), ".")
  }
  unique(by)
}

# Allele copies per genotype from single-letter alleles in genotype names.
.sg_default_alleles <- function(geno) {
  split <- strsplit(geno, "")
  alleles <- unique(unlist(split))
  m <- t(vapply(
    split, function(l) as.numeric(tabulate(match(l, alleles), length(alleles))),
    numeric(length(alleles))
  ))
  dimnames(m) <- list(geno, alleles)
  m
}

.sg_check_allele_matrix <- function(alleles, geno) {
  ok <- is.matrix(alleles) && is.numeric(alleles) && !is.null(rownames(alleles)) &&
    !is.null(colnames(alleles)) && all(geno %in% rownames(alleles)) && all(alleles >= 0)
  if (!ok) {
    .sg_stop(
      "`alleles` must be a non-negative numeric matrix with a row for every genotype (",
      .sg_wrap_names(geno), ") and one column per allele."
    )
  }
  alleles
}

# day, date and (stochastic runs) replicate for every time x replicate.
.sg_time_columns <- function(sim) {
  n_rep <- if (inherits(sim, "sg_sim_cb")) sim$meta$n_rep else 1L
  out <- data.frame(day = rep(sim$day, n_rep), date = rep(sim$date, n_rep))
  if (inherits(sim, "sg_sim_cb")) out$replicate <- rep(seq_len(n_rep), each = length(sim$day))
  out
}

# Long table from values [time, group, replicate] and the group levels.
.sg_long_table <- function(sim, values, levels_df) {
  d <- dim(values)
  tc <- .sg_time_columns(sim)
  out <- do.call(rbind, lapply(seq_len(d[2L]), function(j) {
    cbind(tc, levels_df[rep(j, nrow(tc)), , drop = FALSE], value = as.vector(values[, j, ]))
  }))
  rownames(out) <- NULL
  ord_cols <- c(if ("replicate" %in% names(out)) "replicate", "day")
  out[do.call(order, out[ord_cols]), , drop = FALSE]
}
