#' Deterministic simulation (ODE)
#'
#' Solves the model's ordinary differential equations over the forcing period
#' (model specification S10-S14).
#'
#' @param model An `sg_model` object ([sg_model()]).
#' @param init Initial state: an `sg_state` object from [sg_init()] built for
#'   this model (default: `sg_init(model)`), or a numeric vector of length
#'   `model$n_states`.
#' @param backend `"R"`: the pure-R reference implementation, solved with
#'   'deSolve' (milestone M2). A compiled C++ backend follows in milestone M3.
#' @param method Integration method passed to [deSolve::ode()]; default
#'   `"lsoda"`, which switches automatically between a non-stiff (Adams) and a
#'   stiff (BDF) method. `"radau"`, `"bdf"`, `"ode45"` and other deSolve
#'   methods can be used to check that results do not depend on the solver.
#' @param rtol,atol Relative and absolute error tolerances of the integrator.
#' @param output_every Record the state every this many days (and always on
#'   the last day). Counters and diagnostics are always daily.
#'
#' @details
#' **Day-by-day integration** (D-019). Each day is integrated separately from
#' the start of the day (time d - 1) to its end (time d) with that day's
#' forcing held constant, so the change of weather at midnight is an exact
#' breakpoint. Releases (D-005) are added to the state at the start of their
#' day, before it is integrated: a release is an instantaneous jump, so the
#' solver only ever sees smooth problems and releases add no stiffness.
#'
#' **Counters.** Alongside the states, the flows of the model are integrated
#' over each day (eggs laid by genotype, egg batches, larval and pupal
#' entries, emergences, matings by female x mate genotype, deaths by block and
#' by old age), giving daily totals used for flow accounting and outputs.
#'
#' **Diagnostics.** For each day: integration steps, derivative evaluations,
#' whether `lsoda` finished on its stiff method, the largest rate (transfer
#' plus mortality) of any compartment, and the smallest state value. Values
#' slightly below zero (round-off near extinction) are reported, not removed.
#'
#' @return An object of class `sg_sim`, a list with:
#'   \describe{
#'     \item{`day`, `date`}{Output days (0 = initial state) and dates.}
#'     \item{`state`}{Matrix of states (rows: output days; columns:
#'       [sg_state_names()]), counts per node at the end of each day.}
#'     \item{`counters`}{Matrix of daily flow totals (rows: days 1 to n).}
#'     \item{`releases`}{Matrix of individuals released per day and block.}
#'     \item{`diagnostics`}{Data frame of daily solver diagnostics.}
#'     \item{`meta`}{Backend, method, tolerances, run time, model hash,
#'       package and deSolve versions, initial state.}
#'   }
#' @seealso [sg_model()], [sg_init()]
#' @export
#' @examples
#' if (requireNamespace("deSolve", quietly = TRUE) &&
#'   requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-04-30")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   sim <- simulate_ode(mod)
#'   sim
#' }
simulate_ode <- function(model, init = NULL, backend = "R", method = "lsoda",
                         rtol = 1e-8, atol = 1e-10, output_every = 1L) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  backend <- .sg_match_choice(backend, c("R", "cpp"), "backend")
  if (backend == "cpp") {
    .sg_stop("The C++ backend is not available yet (milestone M3); use backend = \"R\".")
  }
  if (!requireNamespace("deSolve", quietly = TRUE)) {
    .sg_stop("The R backend needs the 'deSolve' package: install.packages(\"deSolve\").")
  }
  x0 <- .sg_check_init(model, init)
  tol_ok <- is.numeric(rtol) && is.numeric(atol) && length(rtol) == 1L &&
    length(atol) == 1L && rtol > 0 && atol > 0
  if (!tol_ok) {
    .sg_stop("`rtol` and `atol` must be single positive numbers.")
  }
  every_ok <- is.numeric(output_every) && length(output_every) == 1L &&
    output_every >= 1 && output_every == round(output_every)
  if (!every_ok) {
    .sg_stop("`output_every` must be a whole number >= 1 (days).")
  }

  setup <- .sg_reference_setup(model)
  n <- model$n_states
  n_days <- model$forcing$n_days
  cnames <- setup$counter_names
  n_c <- length(cnames)
  out_days <- as.integer(unique(c(seq.int(output_every, n_days, by = output_every), n_days)))

  state <- matrix(NA_real_, length(out_days) + 1L, n,
    dimnames = list(NULL, sg_state_names(model))
  )
  state[1L, ] <- x0
  counters <- matrix(0, n_days, n_c, dimnames = list(NULL, cnames))
  rel_blocks <- c("E", "P", "M", "U")
  released <- matrix(0, n_days, length(rel_blocks), dimnames = list(NULL, rel_blocks))
  diag <- data.frame(
    day = seq_len(n_days), steps = NA_integer_, evaluations = NA_integer_,
    stiff = NA, max_rate = NA_real_, min_state = NA_real_
  )
  rel <- model$releases
  rel_block <- substr(sg_state_names(model)[rel$index], 1L, 1L)

  func <- function(t, y, parms) {
    r <- .sg_rates_r(parms$day, y[seq_len(n)], setup)
    list(c(r$deriv, r$counters))
  }

  x <- x0
  row <- 1L
  t0 <- proc.time()[["elapsed"]]
  for (d in seq_len(n_days)) {
    # releases: instantaneous jumps at the start of the day (D-005, D-019)
    on_day <- which(rel$day == d)
    if (length(on_day) > 0L) {
      for (j in on_day) {
        x[rel$index[j]] <- x[rel$index[j]] + rel$count[j]
        released[d, rel_block[j]] <- released[d, rel_block[j]] + rel$count[j]
      }
    }
    diag$max_rate[d] <- .sg_max_rate(d, x, setup)
    sol <- deSolve::ode(
      y = c(x, numeric(n_c)), times = c(d - 1, d), func = func,
      parms = list(day = d), method = method, rtol = rtol, atol = atol
    )
    if (nrow(sol) < 2L || anyNA(sol[2L, -1L])) {
      .sg_stop("The integration failed on day ", d, " (", format(model$forcing$dates[d]), ").")
    }
    y <- sol[2L, -1L]
    x <- unname(y[seq_len(n)])
    counters[d, ] <- y[n + seq_len(n_c)]
    ist <- attr(sol, "istate")
    diag$steps[d] <- ist[2L]
    diag$evaluations[d] <- ist[3L]
    if (method == "lsoda") diag$stiff[d] <- isTRUE(ist[15L] == 2L)
    diag$min_state[d] <- min(x)
    if (d %in% out_days) {
      row <- row + 1L
      state[row, ] <- x
    }
  }
  run_time <- proc.time()[["elapsed"]] - t0

  dates <- model$forcing$dates
  structure(
    list(
      day = c(0L, out_days),
      date = c(dates[1L] - 1, dates[out_days]),
      state = state,
      counters = counters,
      releases = released,
      diagnostics = diag,
      meta = list(
        backend = backend, method = method, rtol = rtol, atol = atol,
        output_every = as.integer(output_every), run_time = run_time,
        model_hash = model$meta$hash,
        package_version = as.character(utils::packageVersion("SGDrivE")),
        desolve_version = as.character(utils::packageVersion("deSolve")),
        init = x0
      )
    ),
    class = "sg_sim"
  )
}

#' @export
print.sg_sim <- function(x, ...) {
  n_days <- nrow(x$counters)
  m <- x$meta
  cat(
    "<sg_sim> ", n_days, " days (", format(x$date[2L]), " to ",
    format(x$date[length(x$date)]), "), backend \"", m$backend, "\", method \"",
    m$method, "\"\n",
    sep = ""
  )
  tot <- .sg_sim_block_totals(x)
  last <- tot[nrow(tot), ]
  cat(
    "  final totals: ",
    paste0(names(last), " ", format(last, digits = 4L), collapse = ", "), "\n",
    sep = ""
  )
  cat(sprintf(
    "  eggs laid %s; adults emerged %s; released %s\n",
    format(sum(x$counters[, grepl("^eggs_laid", colnames(x$counters))]), digits = 5L),
    format(sum(x$counters[, c("emerged_male", "emerged_female")]), digits = 5L),
    format(sum(x$releases), digits = 5L)
  ))
  d <- x$diagnostics
  cat(sprintf(
    "  solver: %.0f steps/day (max %d), %s, largest rate %.3g per day, min state %.3g\n",
    mean(d$steps), max(d$steps),
    if (all(is.na(d$stiff))) {
      "stiffness not reported for this method"
    } else {
      sprintf("stiff method on %d of %d days", sum(d$stiff), nrow(d))
    },
    max(d$max_rate), min(d$min_state)
  ))
  cat(sprintf("  run time %.2f s\n", m$run_time))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

.sg_check_init <- function(model, init) {
  if (is.null(init)) {
    return(as.numeric(sg_init(model)))
  }
  if (inherits(init, "sg_state")) {
    if (!identical(attr(init, "model_hash"), model$meta$hash)) {
      .sg_stop(
        "`init` was built for a different model (its model hash differs). ",
        "Rebuild it with sg_init(model)."
      )
    }
    return(as.numeric(init))
  }
  as.numeric(sg_init(model, "state", state = init))
}

# Largest total rate (transfer + mortality) of any compartment that is in use.
.sg_max_rate <- function(day, x, setup) {
  a <- .sg_rates_r(day, x, setup, aux = TRUE)$aux
  q_f <- a$q_f[is.finite(setup$w_f)]
  q_m <- a$q_m[is.finite(setup$w_m)]
  max(
    a$q[["E"]] + a$mu[["E"]], a$q[["L"]] + a$mu[["L"]], a$q[["P"]] + a$mu[["P"]],
    if (length(q_f)) max(q_f) + a$mu[["F"]] else 0,
    if (length(q_m)) max(q_m) + a$mu[["M"]] else 0
  )
}

# Totals by block at each output day (rows) from the state matrix.
.sg_sim_block_totals <- function(sim) {
  block <- substr(colnames(sim$state), 1L, 1L)
  tot <- vapply(c("E", "L", "P", "M", "U", "F"), function(b) {
    rowSums(sim$state[, block == b, drop = FALSE])
  }, numeric(nrow(sim$state)))
  if (!is.matrix(tot)) tot <- matrix(tot, nrow = 1L, dimnames = list(NULL, names(tot)))
  tot
}
