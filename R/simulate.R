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
#' @param node_extinction If `TRUE`, the node extinction rule (MODEL_SPEC N5,
#'   D-032) is applied at the end of each day: every node whose flies (all
#'   stages, genotypes and sub-stages) total less than `extinction_threshold`
#'   ([sg_params()], default 1 individual) is set to zero; with the host-wound
#'   model its infested hosts return to the hosts without an attractive wound
#'   and are counted in `infestations_ended`. The removed flies are counted by
#'   node and stage (`removed_E` ... `removed_F`). Default `FALSE`, also for
#'   several nodes (D-032); see Details.
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
#' **Node extinction** (MODEL_SPEC N5). In the ODE, movement puts tiny amounts
#' of flies into every connected node, and a node dying out in winter decays
#' towards zero without reaching it; without an Allee effect such fractions
#' of a fly grow again in spring. For questions of spread, die-back and
#' recolonisation use `node_extinction = TRUE` or the chain-binomial solver.
#' The rule is applied after each day is integrated and before the next day's
#' releases, so recorded states show the zeros; the flow counters stay
#' balanced through the `removed_*` counters.
#'
#' @param record What the run stores: `"state"` (default: full states on
#'   output days and all daily counters), `"totals"` (totals by stage and
#'   genotype per node, and all counters), or an [sg_record()] object with
#'   chosen quantities and counter sums, computed in C++ during the run.
#'
#' @return An object of class `sg_sim`, a list with:
#'   \describe{
#'     \item{`day`, `date`}{Output days (0 = initial state) and dates.}
#'     \item{`state`}{Matrix of states (rows: output days; columns:
#'       [sg_state_names()]), counts per node at the end of each day; `NULL`
#'       when not recorded.}
#'     \item{`totals`, `recorded`}{With [sg_record()]: totals by stage and
#'       genotype, and the recorded quantities (rows: output days).}
#'     \item{`counters`}{Matrix of daily flow totals (rows: days 1 to n);
#'       `NULL` when only counter sums are recorded.}
#'     \item{`counter_sums`}{With `sg_counter()` terms: daily counter sums.}
#'     \item{`final_state`}{The state at the end of the run (named), to
#'       continue from.}
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
                         stiff_fallback = NULL, node_extinction = FALSE, max_memory = 1e9,
                         record = "state") {
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
  if (!is.logical(node_extinction) || length(node_extinction) != 1L || is.na(node_extinction)) {
    .sg_stop("`node_extinction` must be TRUE or FALSE.")
  }
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
  ext <- if (node_extinction) .sg_extinction_setup(model) else NULL
  all_counters <- c(.sg_model_counter_names(model), ext$names)
  rec <- .sg_recorder(record, model, all_counters)
  n_out <- as.numeric(length(out_days) + 1L)
  width_out <- (if (rec$state) n else 0) + length(rec$names_x)
  width_day <- (if (rec$counters) length(all_counters) else 0) + length(rec$names_c)
  out_bytes <- 8 * (n_out * width_out + as.numeric(n_days) * width_day + n)
  if (out_bytes > max_memory) {
    .sg_stop(
      "The results would take ", format(out_bytes / 1e6, digits = 3L), " MB (",
      if (rec$state) "states" else "recorded quantities", " on ", n_out,
      " output days and daily counters), more than max_memory (",
      format(max_memory / 1e6, digits = 3L), " MB). Use a larger output_every, record less ",
      "(sg_record()) or raise max_memory."
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
    diag = data.frame(
      day = seq_len(n_days), steps = NA_integer_, rejected = NA_integer_,
      evaluations = NA_integer_, stiff = NA, fallback = FALSE, max_rate = NA_real_,
      min_state = NA_real_
    ),
    rec = rec
  )
  if (rec$state) {
    res$state <- matrix(NA_real_, n_out, n, dimnames = list(NULL, sg_state_names(model)))
  }
  if (length(rec$names_x) > 0L) {
    res$recorded <- matrix(NA_real_, n_out, length(rec$names_x), dimnames = list(NULL, rec$names_x))
  }
  if (rec$counters) {
    res$counters <- matrix(0, n_days, length(all_counters), dimnames = list(NULL, all_counters))
  }
  if (length(rec$names_c) > 0L) {
    res$counter_sums <- matrix(0, n_days, length(rec$names_c), dimnames = list(NULL, rec$names_c))
  }
  if (n_node > 1L) res$diag$limiting_node <- NA_character_
  res <- .sg_record_row(res, 1L, x0)
  ctx <- list(
    model = model, n = n, n_c = length(cnames), n_days = n_days, out_days = out_days, func = func,
    max_rate = max_rate, rtol = rtol, atol = atol, node_ids = .sg_node_ids(model),
    ext = ext
  )

  t0 <- proc.time()[["elapsed"]]
  res <- if (method == "dopri5") {
    .sg_run_dopri5(res, ctx, x0, ptr, stiff_fallback)
  } else {
    .sg_run_desolve(res, ctx, x0, method)
  }
  run_time <- proc.time()[["elapsed"]] - t0

  dates <- model$forcing$dates
  totals <- recorded <- NULL
  if (!is.null(res$recorded)) {
    is_total <- seq_len(ncol(res$recorded)) <= rec$n_totals
    if (any(is_total)) totals <- res$recorded[, is_total, drop = FALSE]
    if (any(!is_total)) recorded <- res$recorded[, !is_total, drop = FALSE]
  }
  final_state <- stats::setNames(res$final_state, sg_state_names(model))
  structure(
    list(
      day = c(0L, out_days),
      date = c(dates[1L] - 1, dates[out_days]),
      state = res$state,
      totals = totals,
      recorded = recorded,
      counters = res$counters,
      counter_sums = res$counter_sums,
      final_state = final_state,
      releases = .sg_release_matrix(model),
      diagnostics = res$diag,
      meta = list(
        backend = backend, method = method, rtol = rtol, atol = atol,
        output_every = as.integer(output_every),
        stiff_fallback = method == "dopri5" && stiff_fallback,
        node_extinction = node_extinction,
        record = list(
          sets = rec$record$sets, quantities = setdiff(rec$names_x, colnames(totals)),
          counter_sums = rec$names_c
        ),
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
  n_days <- nrow(x$diagnostics)
  m <- x$meta
  cat(
    "<sg_sim> ", n_days, " days (", format(x$date[2L]), " to ",
    format(x$date[length(x$date)]), "), backend \"", m$backend, "\", method \"",
    m$method, "\"\n",
    sep = ""
  )
  # what the run kept (phase 1 WP9): states or totals, all counters
  stock <- if (!is.null(x$state)) x$state else x$totals
  has_counters <- !is.null(x$counters)
  if (!is.null(stock)) {
    tot <- .sg_sim_block_totals(x)
    last <- tot[nrow(tot), ]
    cat(
      "  final totals", if (.sg_sim_nodes(x)$multi) " (study nodes)", ": ",
      paste0(names(last), " ", format(last, digits = 4L), collapse = ", "), "\n",
      sep = ""
    )
  }
  if (!is.null(m$info$hosts) && !is.null(stock) && has_counters) {
    hs <- sg_hosts(x)
    cat(sprintf(
      "  hosts: %s infested of %s at the end (prevalence %.3g); %s infestations\n",
      format(hs$infested[nrow(hs)], digits = 4L), format(hs$hosts[nrow(hs)], digits = 6L),
      hs$prevalence[nrow(hs)], format(sum(sg_cases(x, "year", per = "node")$cases), digits = 5L)
    ))
  }
  if (has_counters) {
    cat(sprintf(
      "  eggs laid %s; adults emerged %s; released %s\n",
      format(.sg_counter_total(x, "eggs_laid\\[[^]]*\\]"), digits = 5L),
      format(.sg_counter_total(x, "emerged_(male|female)"), digits = 5L),
      format(sum(x$releases), digits = 5L)
    ))
  }
  n_node <- if (is.null(stock)) 1L else length(unique(sub("^[^@]*@?", "", colnames(stock))))
  if (n_node > 1L) {
    nm <- colnames(stock)
    fly <- !grepl("^H[SWI]\\[", nm)
    last <- stock[nrow(stock), fly]
    node_of <- sub("^[^@]*@", "", nm[fly])
    node_tot <- tapply(last, factor(node_of, levels = unique(node_of)), sum)
    cat(sprintf(
      "  nodes: %d, %d with flies at the end (total >= 1)%s\n", n_node, sum(node_tot >= 1),
      if (has_counters && any(grepl("^emigrants_", colnames(x$counters)))) {
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
  .sg_print_recorded("recorded (last output day)", x$recorded, function(v) v[length(v)])
  .sg_print_recorded("counter sums (whole run)", x$counter_sums, sum)
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

# print() lines for recorded quantities or counter sums: up to six values.
.sg_print_recorded <- function(label, m, fun) {
  if (is.null(m) || ncol(m) == 0L) {
    return(invisible(NULL))
  }
  v <- apply(m, 2L, fun)
  show <- utils::head(v, 6L)
  cat(
    "  ", label, ": ",
    paste0(names(show), " ", format(show, digits = 4L), collapse = ", "),
    if (length(v) > 6L) paste0(", ... (", length(v), " in all)"), "\n",
    sep = ""
  )
  invisible(NULL)
}

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
    step <- .sg_extinction_step(.sg_desolve_day(ctx, x, d, method), ctx)
    x <- step$x
    res <- .sg_store_day(res, d, step)
    if (d %in% ctx$out_days) {
      row <- row + 1L
      res <- .sg_record_row(res, row, x)
    }
  }
  res$final_state <- x
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
      .sg_dopri5_max_steps, !is.null(ctx$ext),
      if (is.null(ctx$ext)) 1 else ctx$ext$threshold,
      if (is.null(ctx$ext)) -1L else ctx$ext$ended_index,
      res$rec[c("state", "counters", "wx", "wc")]
    )
    done <- seq.int(d, length.out = length(r$steps))
    if (length(done) > 0L) {
      if (!is.null(res$counters)) res$counters[done, ] <- r$counters
      if (!is.null(res$counter_sums)) res$counter_sums[done, ] <- r$counter_sums
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
      rows <- match(r$state_days, ctx$out_days) + 1L
      if (!is.null(res$state)) res$state[rows, ] <- r$states
      if (!is.null(res$recorded)) res$recorded[rows, ] <- r$recorded
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
    step <- .sg_extinction_step(.sg_desolve_day(ctx, x, d, "lsoda"), ctx)
    step$stiff <- r$stop_status == 3L
    step$rejected <- NA_integer_
    x <- step$x
    res <- .sg_store_day(res, d, step)
    res$diag$fallback[d] <- TRUE
    if (record[d]) res <- .sg_record_row(res, match(d, ctx$out_days) + 1L, x)
    d <- d + 1L
  }
  res$final_state <- x
  res
}

# Record state x on output row `row` (R paths): the full state and/or the
# recorded weighted sums W_x x (phase 1 WP9).
.sg_record_row <- function(res, row, x) {
  if (!is.null(res$state)) res$state[row, ] <- x
  if (!is.null(res$recorded)) res$recorded[row, ] <- .sg_csr_multiply(res$rec$wx, x)
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

# Node extinction rule (MODEL_SPEC N5): positions of each node's fly stages
# and host blocks, the counter names of the removed flies, the threshold and
# the position of infestations_ended within a node's counters (0-based for
# C++, -1 without hosts).
.sg_extinction_setup <- function(model) {
  ids <- .sg_node_ids(model)
  stages <- c("E", "L", "P", "M", "U", "F")
  base <- paste0("removed_", stages)
  node_counters <- .sg_counter_names(model$cube$genotypes, hosts = !is.null(model$hosts))
  ended <- match("infestations_ended", node_counters)
  list(
    threshold = model$counts$ext_threshold,
    names = if (length(ids) == 1L) {
      base
    } else {
      paste0(rep(base, length(ids)), "@", rep(ids, each = length(base)))
    },
    stages = lapply(seq_along(ids), function(i) {
      lapply(stats::setNames(stages, stages), function(b) sg_state_index(model, b, node = i))
    }),
    hs = if (is.null(model$hosts)) NULL else sg_state_index(model, "HS"),
    hi = if (is.null(model$hosts)) {
      NULL
    } else {
      lapply(seq_along(ids), function(i) {
        sg_state_index(model, "HI", node = i)
      })
    },
    ended_col = if (is.na(ended)) NULL else (seq_along(ids) - 1L) * length(node_counters) + ended,
    ended_index = if (is.na(ended)) -1L else ended - 1L
  )
}

# Apply the node extinction rule to one day's result (R paths; the C++ day
# loop does the same in node_extinction()). Adds step$removed.
.sg_extinction_step <- function(step, ctx) {
  ext <- ctx$ext
  if (is.null(ext)) {
    return(step)
  }
  removed <- numeric(length(ext$names))
  for (i in seq_along(ext$stages)) {
    sums <- vapply(ext$stages[[i]], function(idx) sum(step$x[idx]), 0)
    if (!(sum(sums) < ext$threshold)) next
    for (idx in ext$stages[[i]]) step$x[idx] <- 0
    removed[(i - 1L) * 6L + seq_len(6L)] <- sums
    if (!is.null(ext$hi)) {
      infested <- sum(step$x[ext$hi[[i]]])
      step$x[ext$hi[[i]]] <- 0
      step$x[ext$hs[i]] <- step$x[ext$hs[i]] + infested
      step$counters[ext$ended_col[i]] <- step$counters[ext$ended_col[i]] + infested
    }
  }
  step$removed <- removed
  step
}

# Store one day's counters and diagnostics.
.sg_store_day <- function(res, d, step) {
  full <- c(step$counters, step$removed)
  if (!is.null(res$counters)) res$counters[d, ] <- full
  if (!is.null(res$counter_sums)) res$counter_sums[d, ] <- .sg_csr_multiply(res$rec$wc, full)
  res$diag$steps[d] <- step$steps
  res$diag$evaluations[d] <- step$evaluations
  res$diag$stiff[d] <- step$stiff
  if (!is.null(step$rejected)) res$diag$rejected[d] <- step$rejected
  res$diag$min_state[d] <- min(step$x)
  res
}

# Individuals released per day (rows) and block (columns E, P, M, U); with
# several nodes per block and node ("M@north", D-032).
.sg_release_matrix <- function(model) {
  rel <- model$releases
  cols <- .sg_release_columns(model)
  out <- matrix(0, model$forcing$n_days, length(cols), dimnames = list(NULL, cols))
  nm <- sg_state_names(model)[rel$index]
  rel_col <- substr(nm, 1L, 1L)
  if (length(cols) > 4L) rel_col <- paste0(rel_col, "@", .sg_name_node(nm))
  for (j in seq_along(rel$day)) {
    out[rel$day[j], rel_col[j]] <- out[rel$day[j], rel_col[j]] + rel$count[j]
  }
  out
}

# Columns of the release matrix: blocks E, P, M, U, with several nodes per
# block and node ("M@north").
.sg_release_columns <- function(model) {
  blocks <- c("E", "P", "M", "U")
  ids <- .sg_node_ids(model)
  if (length(ids) == 1L) {
    return(blocks)
  }
  paste0(rep(blocks, length(ids)), "@", rep(ids, each = length(blocks)))
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

# Totals by block at each output day (rows) from the state matrix; with
# several nodes over the study nodes (D-032).
.sg_sim_block_totals <- function(sim) {
  st <- sim$state
  if (inherits(sim, "sg_sim") && (is.null(st) || .sg_sim_nodes(sim)$multi)) {
    a <- .sg_sim_array(sim)
    st <- matrix(a$array[, , 1L], dim(a$array)[1L], dimnames = list(NULL, dimnames(a$array)[[2L]]))
  }
  block <- substr(colnames(st), 1L, 1L)
  tot <- vapply(c("E", "L", "P", "M", "U", "F"), function(b) {
    rowSums(st[, block == b, drop = FALSE])
  }, numeric(nrow(st)))
  if (!is.matrix(tot)) tot <- matrix(tot, nrow = 1L, dimnames = list(NULL, names(tot)))
  tot
}
