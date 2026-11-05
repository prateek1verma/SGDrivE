#' Deterministic simulation (ODE)
#'
#' Solves the model's ordinary differential equations over the forcing period
#' (model specification S10-S14).
#'
#' @param model An `sg_model` object ([sg_model()]).
#' @param init Initial state: an `sg_state` object from [sg_init()] built for
#'   this model (default: `sg_init(model)`), or a numeric vector of length
#'   `model$n_states`.
#' @param backend `"cpp"` (default): the compiled C++ rate kernel and day
#'   loop, 12 to 70 times faster than the R backend. `"R"`: the pure-R
#'   reference implementation, slow but easy to read, used to test the C++
#'   code. The two agree to within the integration tolerance.
#' @param method Integration method. `"dopri5"` (C++ backend only) is the
#'   package's own Dormand-Prince 5(4) integrator, run entirely in C++. Any
#'   method of [deSolve::ode()] can be used with either backend: `"lsoda"`,
#'   which switches automatically between a non-stiff (Adams) and a stiff (BDF)
#'   method, `"radau"`, `"bdf"`, `"ode45"` and others, e.g. to check that
#'   results do not depend on the solver. Default: `"dopri5"` for the C++
#'   backend, `"lsoda"` for the R backend.
#' @param rtol,atol Relative and absolute error tolerances of the integrator.
#' @param output_every Record the state every this many days (and always on
#'   the last day). Counters and diagnostics are always daily.
#' @param stiff_fallback Only for `method = "dopri5"`: if `TRUE`, a
#'   day on which the explicit integrator detects stiffness (or fails) is
#'   integrated again with `"lsoda"`, using the C++ rates. This needs the
#'   'deSolve' package, but only when such a day occurs.
#'   If `FALSE`, stiff days are integrated with Dormand-Prince regardless
#'   (accurate but slow), and a failure is an error. The default `NULL` means
#'   `TRUE` for one node and `FALSE` for several nodes, where `TRUE` is an
#'   error: `lsoda` builds a dense Jacobian of all states (D-032).
#' @param max_memory Largest allowed size in bytes of the stored results
#'   (states on output days and daily counters; default 1e9) and, with more
#'   than one node, of the dense Jacobian of the deSolve methods `"lsoda"`,
#'   `"lsode"`, `"bdf"` and `"radau"`. Larger runs stop before they start.
#'
#' @details
#' **Day-by-day integration** (D-019). Each day is integrated separately from
#' the start of the day (time d - 1) to its end (time d) with that day's
#' forcing held constant, so the change of weather at midnight is an exact
#' breakpoint. Releases (D-005) are added to the state at the start of their
#' day, before it is integrated: a release is an instantaneous jump, so the
#' solver only ever sees smooth problems and releases add no stiffness.
#'
#' **Stiffness** (D-021). The model is rarely stiff (larval crowding and very
#' short sub-stages can make it so), which favours an explicit integrator.
#' Dormand-Prince estimates the stiffness of the problem at every step
#' (Hairer's test); when it is detected, that one day is handed to `"lsoda"`
#' and the C++ day loop resumes on the next day.
#'
#' **Counters.** Alongside the states, the flows of the model are integrated
#' over each day (eggs laid by genotype, egg batches, larval and pupal
#' entries, emergences, matings by female x mate genotype, deaths by block and
#' by old age), giving daily totals used for flow accounting and outputs.
#'
#' **Diagnostics.** For each day: integration steps, rejected steps
#' (Dormand-Prince only), derivative evaluations, whether stiffness was
#' detected (`lsoda`: finished on its stiff method; Dormand-Prince: Hairer's
#' test; `NA` for other methods), whether the day was integrated by the
#' `lsoda` fallback (steps and evaluations are then those of `lsoda`), the
#' largest rate of any compartment (transfer plus mortality, plus mating for
#' unmated females; D-026; with several nodes plus the node's total leave
#' rate), and the smallest
#' state value. Values slightly below zero (round-off near extinction) are
#' reported, not removed. With several nodes and `"dopri5"`, the column
#' `limiting_node` gives the node whose states most often had the largest
#' error in the day's accepted steps, i.e. the node that limits the step
#' size (all nodes share one step size).
#'
#' **Several nodes** (D-032). A model with a landscape is integrated as one
#' system; movement terms are part of the rates (MODEL_SPEC 8a). States and
#' counters are named with the node, e.g. `"M[AA,1]@north"`,
#' `"emigrants_F@north"`.
#'
#' @return An object of class `sg_sim`, a list with:
#'   \describe{
#'     \item{`day`, `date`}{Output days (0 = initial state) and dates.}
#'     \item{`state`}{Matrix of states (rows: output days; columns:
#'       [sg_state_names()]), counts per node at the end of each day.}
#'     \item{`counters`}{Matrix of daily flow totals (rows: days 1 to n).}
#'     \item{`releases`}{Matrix of individuals released per day and block.}
#'     \item{`diagnostics`}{Data frame of daily solver diagnostics.}
#'     \item{`meta`}{Backend, method, tolerances, stiff fallback setting and
#'       number of fallback days, run time, model hash, package and deSolve
#'       versions, initial state.}
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
#'   # the pure-R reference backend gives the same result, more slowly
#'   sim_r <- simulate_ode(mod, backend = "R")
#'   max(abs(sim_r$state - sim$state)) / max(sim$state)
#' }
simulate_ode <- function(model, init = NULL, backend = "cpp", method = NULL,
                         rtol = 1e-8, atol = 1e-10, output_every = 1L,
                         stiff_fallback = NULL, max_memory = 1e9) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  backend <- .sg_match_choice(backend, c("R", "cpp"), "backend")
  if (is.null(method)) method <- if (backend == "cpp") "dopri5" else "lsoda"
  if (!is.character(method) || length(method) != 1L || is.na(method)) {
    .sg_stop("`method` must be a single string, e.g. \"dopri5\" or \"lsoda\".")
  }
  if (method == "dopri5" && backend == "R") {
    .sg_stop(
      "Method \"dopri5\" is the C++ backend's own integrator: use backend = \"cpp\", ",
      "or a deSolve method such as \"ode45\" with backend = \"R\"."
    )
  }
  n_node <- model$dims$n_node
  if (is.null(stiff_fallback)) stiff_fallback <- n_node == 1L
  if (!is.logical(stiff_fallback) || length(stiff_fallback) != 1L || is.na(stiff_fallback)) {
    .sg_stop("`stiff_fallback` must be TRUE or FALSE (or NULL for the default).")
  }
  if (stiff_fallback && n_node > 1L && method == "dopri5") {
    .sg_stop(
      "stiff_fallback = TRUE is not available with several nodes: its lsoda method ",
      "builds a dense Jacobian of all ", format(model$n_states, big.mark = ","),
      " states. Use stiff_fallback = FALSE (the default for several nodes)."
    )
  }
  if (!is.numeric(max_memory) || length(max_memory) != 1L || !(max_memory > 0)) {
    .sg_stop("`max_memory` must be a single positive number (bytes).")
  }
  dense <- c("lsoda", "lsode", "bdf", "radau")
  jac_bytes <- 8 * as.numeric(model$n_states)^2
  if (n_node > 1L && method %in% dense && jac_bytes > max_memory) {
    .sg_stop(
      "Method \"", method, "\" needs a dense Jacobian of ", format(jac_bytes / 1e6, digits = 3L),
      " MB for ", format(model$n_states, big.mark = ","), " states, more than max_memory (",
      format(max_memory / 1e6, digits = 3L), " MB). Use method = \"dopri5\" (backend = \"cpp\")",
      " or a non-stiff deSolve method such as \"ode45\", or raise max_memory."
    )
  }
  if (method != "dopri5") .sg_need_desolve(paste0("Method \"", method, "\""))
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

  n <- model$n_states
  n_days <- model$forcing$n_days
  out_days <- as.integer(unique(c(seq.int(output_every, n_days, by = output_every), n_days)))
  n_counters <- length(.sg_model_counter_names(model))
  out_bytes <- 8 * (as.numeric(length(out_days) + 1L) * n + as.numeric(n_days) * n_counters)
  if (out_bytes > max_memory) {
    .sg_stop(
      "The results would take ", format(out_bytes / 1e6, digits = 3L), " MB (states on ",
      length(out_days) + 1L, " output days and daily counters), more than max_memory (",
      format(max_memory / 1e6, digits = 3L), " MB). Use a larger output_every or raise max_memory."
    )
  }

  # rate functions of the chosen backend: derivatives and counter rates
  # c(dydt, counters) for y = c(state, counters), and the largest rate
  if (backend == "R") {
    setup <- .sg_reference_setup(model)
    cnames <- setup$counter_names
    ptr <- NULL
    rates <- function(day, y) {
      r <- .sg_rates_r(day, y[seq_len(n)], setup)
      c(r$deriv, r$counters)
    }
    max_rate <- function(day, x) .sg_max_rate(day, x, setup)
  } else {
    view <- .sg_model_view(model)
    cnames <- attr(view, "counter_names")
    ptr <- sg_kernel_new_cpp(view)
    rates <- function(day, y) sg_kernel_rates_cpp(ptr, day, y)
    max_rate <- function(day, x) sg_kernel_max_rate_cpp(ptr, day, x)
  }
  func <- function(t, y, parms) list(rates(parms$day, y))

  res <- list(
    state = matrix(NA_real_, length(out_days) + 1L, n,
      dimnames = list(NULL, sg_state_names(model))
    ),
    counters = matrix(0, n_days, length(cnames), dimnames = list(NULL, cnames)),
    diag = data.frame(
      day = seq_len(n_days), steps = NA_integer_, rejected = NA_integer_,
      evaluations = NA_integer_, stiff = NA, fallback = FALSE, max_rate = NA_real_,
      min_state = NA_real_
    )
  )
  if (n_node > 1L) res$diag$limiting_node <- NA_character_
  res$state[1L, ] <- x0
  ctx <- list(
    model = model, n = n, n_c = length(cnames), n_days = n_days, out_days = out_days, func = func,
    max_rate = max_rate, rtol = rtol, atol = atol, node_ids = .sg_node_ids(model)
  )

  t0 <- proc.time()[["elapsed"]]
  res <- if (method == "dopri5") {
    .sg_run_dopri5(res, ctx, x0, ptr, stiff_fallback)
  } else {
    .sg_run_desolve(res, ctx, x0, method)
  }
  run_time <- proc.time()[["elapsed"]] - t0

  dates <- model$forcing$dates
  structure(
    list(
      day = c(0L, out_days),
      date = c(dates[1L] - 1, dates[out_days]),
      state = res$state,
      counters = res$counters,
      releases = .sg_release_matrix(model),
      diagnostics = res$diag,
      meta = list(
        backend = backend, method = method, rtol = rtol, atol = atol,
        output_every = as.integer(output_every),
        stiff_fallback = method == "dopri5" && stiff_fallback,
        fallback_days = sum(res$diag$fallback), run_time = run_time,
        model_hash = model$meta$hash, info = .sg_output_info(model),
        package_version = as.character(utils::packageVersion("ScrewDrivE")),
        desolve_version = if (.sg_has_desolve()) {
          as.character(utils::packageVersion("deSolve"))
        } else {
          NA_character_
        },
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
  if (!is.null(m$info$hosts)) {
    hs <- sg_hosts(x)
    cat(sprintf(
      "  hosts: %s infested of %s at the end (prevalence %.3g); %s infestations\n",
      format(hs$infested[nrow(hs)], digits = 4L), format(hs$hosts[nrow(hs)], digits = 6L),
      hs$prevalence[nrow(hs)], format(.sg_counter_total(x, "infestations"), digits = 5L)
    ))
  }
  cat(sprintf(
    "  eggs laid %s; adults emerged %s; released %s\n",
    format(.sg_counter_total(x, "eggs_laid\\[[^]]*\\]"), digits = 5L),
    format(.sg_counter_total(x, "emerged_(male|female)"), digits = 5L),
    format(sum(x$releases), digits = 5L)
  ))
  n_node <- length(unique(sub("^[^@]*@?", "", colnames(x$state))))
  if (n_node > 1L) {
    nm <- colnames(x$state)
    fly <- !grepl("^H[SWI]\\[", nm)
    last <- x$state[nrow(x$state), fly]
    node_of <- sub("^[^@]*@", "", nm[fly])
    node_tot <- tapply(last, factor(node_of, levels = unique(node_of)), sum)
    cat(sprintf(
      "  nodes: %d, %d with flies at the end (total >= 1)%s\n", n_node, sum(node_tot >= 1),
      if (any(grepl("^emigrants_", colnames(x$counters)))) {
        sprintf(
          "; movers %s, lost from the landscape %s",
          format(.sg_counter_total(x, "emigrants_[MUF]"), digits = 5L),
          format(.sg_counter_total(x, "lost_[MUF]"), digits = 5L)
        )
      } else {
        ""
      }
    ))
  }
  d <- x$diagnostics
  cat(sprintf(
    "  solver: %.0f steps/day (max %d), %s, largest rate %.3g per day, min state %.3g\n",
    mean(d$steps), max(d$steps),
    if (m$method == "dopri5") {
      sprintf(
        "stiffness detected on %d of %d days, %d integrated by the lsoda fallback",
        sum(d$stiff), nrow(d), sum(d$fallback)
      )
    } else if (all(is.na(d$stiff))) {
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

# Total over days (and nodes) of the counters whose name, without the node
# suffix, matches the regular expression `pattern`.
.sg_counter_total <- function(sim, pattern) {
  cn <- colnames(sim$counters)
  sum(sim$counters[, grepl(paste0("^", pattern, "(@.*)?$"), cn), drop = FALSE])
}

# Day loop with a deSolve method (either backend's rates), D-019.
.sg_run_desolve <- function(res, ctx, x, method) {
  rel <- ctx$model$releases
  row <- 1L
  for (d in seq_len(ctx$n_days)) {
    # releases: instantaneous jumps at the start of the day (D-005, D-019)
    for (j in which(rel$day == d)) x[rel$index[j]] <- x[rel$index[j]] + rel$count[j]
    res$diag$max_rate[d] <- ctx$max_rate(d, x)
    step <- .sg_desolve_day(ctx, x, d, method)
    x <- step$x
    res <- .sg_store_day(res, d, step)
    if (d %in% ctx$out_days) {
      row <- row + 1L
      res$state[row, ] <- x
    }
  }
  res
}

# Day loop in C++ (Dormand-Prince); days on which it stops (stiff or failed)
# are integrated with lsoda using the C++ rates, then the C++ loop resumes.
.sg_run_dopri5 <- function(res, ctx, x, ptr, stiff_fallback) {
  record <- seq_len(ctx$n_days) %in% ctx$out_days
  status_text <- c("too many steps", "step size too small", "stiff", "non-finite state")
  d <- 1L
  while (d <= ctx$n_days) {
    r <- sg_ode_run_cpp(
      ptr, x, d, ctx$n_days, ctx$rtol, ctx$atol, record, stiff_fallback,
      .sg_dopri5_max_steps
    )
    done <- seq.int(d, length.out = length(r$steps))
    if (length(done) > 0L) {
      res$counters[done, ] <- r$counters
      res$diag$steps[done] <- r$steps
      res$diag$rejected[done] <- r$rejected
      res$diag$evaluations[done] <- r$evaluations
      res$diag$stiff[done] <- r$stiff
      res$diag$max_rate[done] <- r$max_rate
      res$diag$min_state[done] <- r$min_state
      if (length(r$limiting_node) > 0L) {
        res$diag$limiting_node[done] <- ifelse(r$limiting_node > 0L,
          ctx$node_ids[pmax(1L, r$limiting_node)], NA_character_
        )
      }
      res$state[match(r$state_days, ctx$out_days) + 1L, ] <- r$states
    }
    x <- r$x
    if (r$stopped_day == 0L) break
    d <- r$stopped_day
    if (r$stop_status == 4L) {
      .sg_stop(
        "The solution became non-finite on day ", d, " (", format(ctx$model$forcing$dates[d]),
        "); check the parameters."
      )
    }
    if (!stiff_fallback) {
      .sg_stop(
        "The Dormand-Prince integration failed on day ", d, " (",
        format(ctx$model$forcing$dates[d]), "): ", status_text[r$stop_status],
        ". Use stiff_fallback = TRUE or a stiff deSolve method such as \"lsoda\"."
      )
    }
    .sg_need_desolve(paste0(
      "Day ", d, " (", format(ctx$model$forcing$dates[d]), ") is stiff; its fallback"
    ), " (or set stiff_fallback = FALSE)")
    # x is the state at the start of day d, after its releases
    res$diag$max_rate[d] <- ctx$max_rate(d, x)
    step <- .sg_desolve_day(ctx, x, d, "lsoda")
    step$stiff <- r$stop_status == 3L
    step$rejected <- NA_integer_
    x <- step$x
    res <- .sg_store_day(res, d, step)
    res$diag$fallback[d] <- TRUE
    if (record[d]) res$state[match(d, ctx$out_days) + 1L, ] <- x
    d <- d + 1L
  }
  res
}

.sg_need_desolve <- function(what, alternative = "") {
  if (!.sg_has_desolve()) {
    .sg_stop(what, " needs the 'deSolve' package: install.packages(\"deSolve\")", alternative, ".")
  }
}
.sg_has_desolve <- function() requireNamespace("deSolve", quietly = TRUE)

# Maximum number of Dormand-Prince steps per day.
.sg_dopri5_max_steps <- 1e5

# Integrate one day (from d - 1 to d) with a deSolve method, counters reset.
.sg_desolve_day <- function(ctx, x, d, method) {
  sol <- deSolve::ode(
    y = c(x, numeric(ctx$n_c)), times = c(d - 1, d), func = ctx$func,
    parms = list(day = d), method = method, rtol = ctx$rtol, atol = ctx$atol
  )
  ist <- attr(sol, "istate")
  rst <- attr(sol, "rstate")
  # deSolve does not always flag a failure: it may return a negative code,
  # missing values, or (lsoda, when the step size underflows) a success code
  # with the starting values, having reached an earlier time (rstate[3])
  failed <- nrow(sol) < 2L || any(!is.finite(sol[2L, -1L])) ||
    (!is.null(ist) && !is.na(ist[1L]) && ist[1L] < 0) ||
    (method == "lsoda" && length(rst) >= 3L && is.finite(rst[3L]) && rst[3L] < d - 1e-9)
  if (failed) {
    .sg_stop(
      "The integration failed on day ", d, " (", format(ctx$model$forcing$dates[d]), "): ",
      "the solution may have become non-finite or grown without bound; check the parameters."
    )
  }
  y <- unname(sol[2L, -1L])
  list(
    x = y[seq_len(ctx$n)], counters = y[ctx$n + seq_len(ctx$n_c)],
    steps = ist[2L], evaluations = ist[3L],
    stiff = if (method == "lsoda") isTRUE(ist[15L] == 2L) else NA
  )
}

# Store one day's counters and diagnostics.
.sg_store_day <- function(res, d, step) {
  res$counters[d, ] <- step$counters
  res$diag$steps[d] <- step$steps
  res$diag$evaluations[d] <- step$evaluations
  res$diag$stiff[d] <- step$stiff
  if (!is.null(step$rejected)) res$diag$rejected[d] <- step$rejected
  res$diag$min_state[d] <- min(step$x)
  res
}

# Individuals released per day (rows) and block (columns E, P, M, U).
.sg_release_matrix <- function(model) {
  rel <- model$releases
  blocks <- c("E", "P", "M", "U")
  out <- matrix(0, model$forcing$n_days, length(blocks), dimnames = list(NULL, blocks))
  rel_block <- substr(sg_state_names(model)[rel$index], 1L, 1L)
  for (j in seq_along(rel$day)) {
    out[rel$day[j], rel_block[j]] <- out[rel$day[j], rel_block[j]] + rel$count[j]
  }
  out
}

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

# Largest total rate (transfer + mortality, and mating for unmated females;
# D-026) of any compartment that is in use. With several nodes: the largest
# over nodes of the node's largest rate plus its largest total leave rate
# (MODEL_SPEC 8a).
.sg_max_rate <- function(day, x, setup) {
  if (is.null(setup$node_setups)) {
    return(.sg_max_rate_node(day, x, setup))
  }
  n_s <- setup$n_states_node
  mv <- .sg_movement_r(day, x, setup)
  max(vapply(seq_len(setup$n_node), function(i) {
    r <- .sg_max_rate_node(day, x[(i - 1L) * n_s + seq_len(n_s)], setup$node_setups[[i]])
    r + if (is.null(mv)) 0 else max(mv$lambda[, i])
  }, 0))
}

.sg_max_rate_node <- function(day, x, setup) {
  a <- .sg_rates_node(day, x, setup, aux = TRUE)$aux
  p <- setup$p
  nu_i <- if (p$mating_mode == "instant") {
    rep(p$mating_rate * a$g_find, setup$n_A)
  } else {
    p$mating_rate * setup$receptive * a$g_find
  }
  mate_ok <- rowSums(a$pi) > 0
  # unmated females [substage, genotype]: ageing, death and mating
  u_rate <- a$q_f + a$mu[["F"]] + outer(nu_i, as.numeric(mate_ok))
  u_rate <- u_rate[is.finite(setup$w_f)]
  q_m <- a$q_m[is.finite(setup$w_m)]
  h <- setup$hosts
  host_rate <- if (is.null(h)) {
    0
  } else {
    max(h$wounding, h$healing + a$infestation_hazard, h$advance + h$treatment)
  }
  max(
    a$q[["E"]] + a$mu[["E"]], a$q[["L"]] + a$mu[["L"]], a$q[["P"]] + a$mu[["P"]],
    if (length(u_rate)) max(u_rate) else 0,
    if (length(q_m)) max(q_m) + a$mu[["M"]] else 0,
    host_rate
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
