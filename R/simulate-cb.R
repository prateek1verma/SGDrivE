#' Stochastic simulation (chain-binomial)
#'
#' Simulates replicates of the model with whole numbers of individuals
#' (model specification S15-S19), using the compiled C++ solver.
#'
#' @param model An `sg_model` object ([sg_model()]).
#' @param init Initial state: an `sg_state` object from [sg_init()] built for
#'   this model (default: `sg_init(model)`), or a numeric vector of length
#'   `model$n_states`. Non-integer values are rounded without bias in each
#'   replicate (`floor(x)` plus a Bernoulli draw on the fraction), as are
#'   non-integer release sizes.
#' @param n_rep Number of replicates.
#' @param dt Time step in days; must divide the day into a whole number of
#'   steps. Default 0.1 day (D-022).
#' @param output_every Record the state every this many days (and always on
#'   the last day). Counters are always daily.
#' @param keep `"state"` (every compartment) or `"totals"` (totals by block
#'   and genotype, much smaller for long runs or many replicates).
#' @param seed A whole number: the seed of the replicates' random-number
#'   streams. `NULL` (default) draws one from R's current random-number
#'   generator, so `set.seed()` before the call also makes runs
#'   reproducible.
#' @param parallel If `TRUE`, replicates run in parallel with
#'   [future.apply::future_lapply()] under the plan set with
#'   `future::plan()` (for example `future::plan("multisession")`). Results
#'   are the same as with `FALSE`.
#' @param max_memory Largest allowed size of the results in bytes (default
#'   1e9). Use `keep = "totals"` or a larger `output_every` for big runs.
#'
#' @details
#' **Steps** (D-022). Each day is divided into steps of length `dt`; every
#' step of a day uses that day's forcing, and releases are added at the start
#' of their day. Within a step, individuals move along their sub-stage chains
#' with the exact probabilities of the chain over the step, die, mate, lay
#' eggs and emerge, with binomial, multinomial and Poisson draws. Arrivals
#' during a step (eggs laid, individuals moving on to the next stage,
#' emerged adults) are placed at their expected position. Eggs are laid in
#' batches by default (`oviposition` and `clutch_size` in [sg_params()],
#' D-024). The mean of a step is the deterministic map, which is within
#' about 0.5% of the ODE at `dt = 0.1`.
#'
#' **Random numbers.** Replicate `r` uses its own L'Ecuyer-CMRG stream
#' (stream `r` derived from `seed`), so it is identical whether it is run
#' alone, among other replicates, sequentially or in parallel. R's own
#' generators are used. The random-number generator of the session is left
#' as it was (apart from drawing the seed when `seed = NULL`).
#'
#' **Extinction.** A replicate whose population is zero at the end of a day
#' with no releases left stops early; its later days are zero.
#' `extinct_day` gives the first day ending with no individuals (`NA` if
#' none), also when later releases bring the population back.
#'
#' **Population size.** Counts are whole individuals per node, so the node
#' area matters: with the defaults (1 km^2, H = 5 wounds per km^2) a McAllen
#' population falls below one individual in most winters. Use a larger
#' `area` in [sg_model()] for populations that persist.
#'
#' @return An object of class `sg_sim_cb`, a list with:
#'   \describe{
#'     \item{`day`, `date`}{Output days (0 = initial state) and dates.}
#'     \item{`state` or `totals`}{Array of counts per node with dimensions
#'       output day x state x replicate, or output day x (block and genotype)
#'       x replicate with `keep = "totals"` (mated females summed over
#'       mates).}
#'     \item{`counters`}{Array of daily flows, day x counter x replicate.}
#'     \item{`releases`}{Array of individuals released, day x block x
#'       replicate.}
#'     \item{`extinct_day`}{First day ending with no individuals, per
#'       replicate.}
#'     \item{`meta`}{Settings, the seed and each replicate's stream, run
#'       time, model hash, package version, initial state (before rounding).}
#'   }
#' @seealso [simulate_ode()], [sg_init()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-03-01", end = "1971-06-30")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx, area = 1000)
#'   sim <- simulate_cb(mod, sg_init(mod, "burn_in", years = 1, weather = sg_example_weather()),
#'     n_rep = 5, seed = 1
#'   )
#'   sim
#' }
simulate_cb <- function(model, init = NULL, n_rep = 1L, dt = 0.1, output_every = 1L,
                        keep = c("state", "totals"), seed = NULL, parallel = FALSE,
                        max_memory = 1e9) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  whole <- function(x, lower) {
    is.numeric(x) && length(x) == 1L && is.finite(x) && x >= lower && x == round(x)
  }
  if (!whole(n_rep, 1)) .sg_stop("`n_rep` must be a whole number >= 1.")
  if (!whole(output_every, 1)) .sg_stop("`output_every` must be a whole number >= 1 (days).")
  keep <- .sg_match_choice(keep, c("state", "totals"), "keep")
  if (!is.null(seed) && !whole(seed, -.Machine$integer.max)) {
    .sg_stop("`seed` must be NULL or a single whole number.")
  }
  if (!is.logical(parallel) || length(parallel) != 1L || is.na(parallel)) {
    .sg_stop("`parallel` must be TRUE or FALSE.")
  }
  if (!is.numeric(max_memory) || length(max_memory) != 1L || !(max_memory > 0)) {
    .sg_stop("`max_memory` must be a positive number of bytes.")
  }
  if (parallel && !requireNamespace("future.apply", quietly = TRUE)) {
    .sg_stop(
      "parallel = TRUE needs the 'future.apply' package: install.packages(\"future.apply\")."
    )
  }
  n_rep <- as.integer(n_rep)
  n_days <- model$forcing$n_days
  x0 <- .sg_check_init(model, init)
  .sg_cb_inputs(model, x0, dt, stochastic = FALSE, n_days = n_days) # checks dt
  out_days <- as.integer(unique(c(seq.int(output_every, n_days, by = output_every), n_days)))
  record <- seq_len(n_days) %in% out_days
  view <- .sg_model_view(model)
  cnames <- attr(view, "counter_names")
  agg <- if (keep == "totals") .sg_totals_matrix(model) else NULL
  width <- if (keep == "state") model$n_states else ncol(agg)

  bytes <- 8 * n_rep * ((length(out_days) + 1) * width + n_days * (length(cnames) + 4))
  if (bytes > max_memory) {
    .sg_stop(sprintf(
      "The results would take about %.3g GB (limit `max_memory` = %.3g GB). %s",
      bytes / 1e9, max_memory / 1e9,
      "Use keep = \"totals\", a larger output_every, fewer replicates, or raise max_memory."
    ))
  }

  rng <- .sg_cb_streams(n_rep, seed)
  one <- function(r) {
    .sg_with_stream(rng$streams[[r]], {
      inp <- .sg_cb_inputs(model, x0, dt, stochastic = TRUE, n_days = n_days)
      res <- sg_cb_run_cpp(view, inp$x, as.numeric(inp$rel$count), dt, TRUE, n_days, record)
    })
    states <- rbind(inp$x, res$states)
    list(
      state = if (keep == "state") states else states %*% agg,
      counters = res$counters, releases = inp$released,
      extinct_day = if (res$extinct_day > 0L) res$extinct_day else NA_integer_
    )
  }
  t0 <- proc.time()[["elapsed"]]
  reps <- if (parallel) {
    future.apply::future_lapply(seq_len(n_rep), one, future.seed = TRUE)
  } else {
    lapply(seq_len(n_rep), one)
  }
  run_time <- proc.time()[["elapsed"]] - t0

  rep_names <- as.character(seq_len(n_rep))
  stack <- function(name, cols) {
    a <- simplify2array(lapply(reps, `[[`, name), higher = TRUE)
    dim(a) <- c(nrow(reps[[1L]][[name]]), length(cols), n_rep)
    dimnames(a) <- list(NULL, cols, rep_names)
    a
  }
  cols <- if (keep == "state") sg_state_names(model) else colnames(agg)
  dates <- model$forcing$dates
  out <- list(day = c(0L, out_days), date = c(dates[1L] - 1, dates[out_days]))
  out[[keep]] <- stack("state", cols)
  out$counters <- stack("counters", cnames)
  out$releases <- stack("releases", c("E", "P", "M", "U"))
  out$extinct_day <- stats::setNames(vapply(reps, `[[`, 0L, "extinct_day"), rep_names)
  p <- model$params
  out$meta <- list(
    backend = "cpp", dt = dt, n_rep = n_rep, keep = keep,
    output_every = as.integer(output_every), seed = rng$seed,
    streams = do.call(cbind, rng$streams), parallel = parallel, run_time = run_time,
    oviposition = p$oviposition, clutch_size = p$clutch_size, mating = p$mating,
    model_hash = model$meta$hash,
    package_version = as.character(utils::packageVersion("SGDrivE")), init = x0
  )
  structure(out, class = "sg_sim_cb")
}

#' @export
print.sg_sim_cb <- function(x, ...) {
  m <- x$meta
  n_days <- dim(x$counters)[1L]
  cat(sprintf(
    "<sg_sim_cb> %d replicate%s x %d days (%s to %s), dt = %g day, eggs in %s\n",
    m$n_rep, if (m$n_rep == 1L) "" else "s", n_days, format(x$date[2L]),
    format(x$date[length(x$date)]), m$dt, m$oviposition
  ))
  tot <- .sg_cb_block_totals(x)
  last <- tot[dim(tot)[1L], , , drop = FALSE]
  dim(last) <- dim(tot)[2:3]
  fmt <- function(v) format(v, digits = 4L, big.mark = ",")
  cat("  final totals, mean [min, max] over replicates:\n")
  for (b in seq_len(nrow(last))) {
    cat(sprintf(
      "    %s %s [%s, %s]\n", dimnames(tot)[[2L]][b], fmt(mean(last[b, ])),
      fmt(min(last[b, ])), fmt(max(last[b, ]))
    ))
  }
  ext <- x$extinct_day
  cat(sprintf(
    "  extinct: %d of %d replicates%s\n", sum(!is.na(ext)), length(ext),
    if (any(!is.na(ext))) {
      sprintf(" (median day %g)", stats::median(ext, na.rm = TRUE))
    } else {
      ""
    }
  ))
  egg_cols <- grepl("^eggs_laid", dimnames(x$counters)[[2L]])
  eggs <- apply(x$counters[, egg_cols, , drop = FALSE], 3L, sum)
  cat(sprintf(
    "  eggs laid %s and released %s per replicate (mean)\n",
    fmt(mean(eggs)), fmt(mean(apply(x$releases, 3L, sum)))
  ))
  cat(sprintf("  seed %s; run time %.2f s\n", format(m$seed), m$run_time))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# Matrix summing states into totals by block and genotype (mated females
# over mates): columns "E[g]", ..., "F[g]".
.sg_totals_matrix <- function(model) {
  nm <- sg_state_names(model)
  block <- substr(nm, 1L, 1L)
  geno <- sub("^.\\[([^,|]+).*$", "\\1", nm)
  cols <- as.vector(outer(
    model$cube$genotypes, c("E", "L", "P", "M", "U", "F"),
    function(g, b) paste0(b, "[", g, "]")
  ))
  agg <- matrix(0, length(nm), length(cols), dimnames = list(NULL, cols))
  agg[cbind(seq_along(nm), match(paste0(block, "[", geno, "]"), cols))] <- 1
  agg
}

# Totals by block [output day, block, replicate] of an sg_sim_cb object.
.sg_cb_block_totals <- function(sim) {
  a <- if (!is.null(sim$state)) sim$state else sim$totals
  block <- substr(dimnames(a)[[2L]], 1L, 1L)
  blocks <- c("E", "L", "P", "M", "U", "F")
  out <- vapply(blocks, function(b) {
    apply(a[, block == b, , drop = FALSE], c(1L, 3L), sum)
  }, matrix(0, dim(a)[1L], dim(a)[3L]))
  out <- aperm(array(out, c(dim(a)[1L], dim(a)[3L], length(blocks))), c(1L, 3L, 2L))
  dimnames(out) <- list(NULL, blocks, dimnames(a)[[3L]])
  out
}

# Seeds of the replicates' L'Ecuyer-CMRG streams (D-022). Without `seed`,
# one is drawn from the session's generator (which advances it, so repeated
# calls differ); the session's generator kind and state are then restored.
.sg_cb_streams <- function(n_rep, seed) {
  if (is.null(seed)) seed <- sample.int(.Machine$integer.max, 1L)
  seed <- as.integer(seed)
  streams <- .sg_with_stream(NULL, {
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    s <- get(".Random.seed", envir = globalenv())
    out <- vector("list", n_rep)
    for (r in seq_len(n_rep)) {
      out[[r]] <- s
      s <- parallel::nextRNGStream(s)
    }
    out
  })
  list(seed = seed, streams = streams)
}

# Evaluate `expr` with the L'Ecuyer-CMRG stream `stream` as R's generator
# state (or, with NULL, leave it to `expr`), restoring the session's
# generator kind and state afterwards.
.sg_with_stream <- function(stream, expr) {
  genv <- globalenv()
  had_seed <- exists(".Random.seed", envir = genv, inherits = FALSE)
  old_seed <- if (had_seed) get(".Random.seed", envir = genv) else NULL
  old_kind <- RNGkind()
  on.exit({
    RNGkind(old_kind[1L], old_kind[2L], old_kind[3L])
    if (had_seed) {
      assign(".Random.seed", old_seed, envir = genv) # nolint: object_name_linter.
    } else if (exists(".Random.seed", envir = genv, inherits = FALSE)) {
      rm(".Random.seed", envir = genv)
    }
  })
  if (!is.null(stream)) {
    RNGkind("L'Ecuyer-CMRG")
    assign(".Random.seed", stream, envir = genv) # nolint: object_name_linter.
  }
  expr
}
