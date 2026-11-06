#' What a simulation stores
#'
#' `sg_record()` says which results [simulate_ode()] and [simulate_cb()]
#' keep (phases 1 and 2 of the multi-node work, D-032, design note section
#' 7.1): full states, totals,
#' chosen weighted sums of states (`sg_quantity()`, `sg_quantity_weights()`),
#' all daily counters, or chosen sums of counters (`sg_counter()`). The sums
#' are computed in C++ during the run, so a long run of many nodes stores
#' only what is needed.
#'
#' @param ... Any of the built-in sets `"state"` (full states on output days),
#'   `"totals"` (totals by stage and genotype, and host blocks, for each node),
#'   `"counters"` (all daily counters) and `"cases"` (the counter sums
#'   [sg_cases()] needs, for every node), and objects from `sg_quantity()`,
#'   `sg_quantity_weights()` and `sg_counter()`.
#' @param every Period of the counter sums (`sg_counter()` terms and the
#'   `"cases"` set): `"day"` (default), `"week"` (7-day blocks from the first
#'   simulated day, as `sg_cases(period = "week")`), `"month"` or `"year"`
#'   (calendar months and years). The sums are added up by period during the
#'   run, so only one row per period is stored. Full counters
#'   (`"counters"`) stay daily.
#'
#' @details
#' * `simulate_ode(record = "totals")` (the default since phase 3, decision
#'   13) keeps totals and all counters; `record = "state"` keeps full states
#'   and all counters, as before. `simulate_cb(record = )` works the same
#'   way for each replicate, with an extra last dimension for the replicates
#'   (`keep = "totals"` is a shortcut for `record = "totals"`).
#' * Quantities are stored in `sim$recorded` (one row per output day, day 0
#'   included), totals in `sim$totals`, states in `sim$state`.
#' * With `sg_counter()` terms or the `"cases"` set only those counter sums
#'   are kept (in `sim$counter_sums`, one row per period of `every`, with the
#'   periods' first days and lengths in `sim$counter_periods`) unless
#'   `"counters"` is also given; without them all counters are kept in
#'   `sim$counters`.
#' * The `"cases"` set records, for every node, the counters `sg_cases()`
#'   reads (`fertile_batches` and `egg_batches`; with the host-wound model
#'   `infestations`, `egg_batches` and `superinfestations`), named like the
#'   counters, so `sg_cases()` works without full counters: weekly county
#'   cases over decades take a few megabytes.
#' * `sim$final_state` (the state at the end of the run; for `simulate_cb()`
#'   a matrix state x replicate) is always kept, so a run can be continued
#'   from it.
#' * Summaries such as [sg_totals()] need `"state"` or `"totals"`; [sg_cases()]
#'   and [sg_flow_balance()] need `"counters"`. They say so when it is missing.
#'
#' @return An object of class `sg_record`.
#' @seealso [simulate_ode()], [simulate_cb()]
#' @export
#' @examples
#' rec <- sg_record(
#'   "totals",
#'   sg_quantity("adults", stage = "adults", by = "node"),
#'   sg_quantity("wild_mated", stage = "F", genotype = "AA"),
#'   sg_counter("cases", counter = "fertile_batches")
#' )
#' rec
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx)
#'   sim <- simulate_ode(mod, record = rec)
#'   tail(sim$recorded)
#'   tail(sim$counter_sums)
#'
#'   # three nodes with movement: adults and cases per node
#'   land <- sg_landscape(data.frame(id = c("a", "b", "c"), area = 100), wx)
#'   chain <- sg_movement_layer(
#'     list(male = 0.05, virgin = 0.1, mated = 0.17),
#'     rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0))
#'   )
#'   mod3 <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), landscape = land, movement = chain)
#'   sim3 <- simulate_ode(mod3,
#'     init = sg_init(mod3, nodes = "a"),
#'     record = sg_record(
#'       sg_quantity("adults", stage = "adults", by = "node"),
#'       sg_counter("cases", counter = "fertile_batches", by = "node")
#'     )
#'   )
#'   tail(sim3$recorded)
#' }
sg_record <- function(..., every = c("day", "week", "month", "year")) {
  every <- .sg_match_choice(every, c("day", "week", "month", "year"), "every")
  parts <- list(...)
  sets <- character(0)
  terms <- list()
  for (x in parts) {
    if (is.character(x)) {
      bad <- setdiff(x, .sg_record_sets)
      if (length(bad) > 0L) {
        .sg_stop(
          "Unknown record set \"", bad[1L], "\"; use ",
          paste0("\"", .sg_record_sets, "\"", collapse = ", "), ".",
          .sg_suggest(bad[1L], .sg_record_sets)
        )
      }
      sets <- union(sets, x)
    } else if (inherits(x, c("sg_quantity", "sg_quantity_weights", "sg_counter", "sg_alleles"))) {
      terms[[length(terms) + 1L]] <- x
    } else {
      .sg_stop(
        "sg_record() takes \"state\", \"totals\", \"counters\", \"cases\" and objects from ",
        "sg_quantity(), sg_quantity_weights(), sg_alleles() and sg_counter()."
      )
    }
  }
  if (length(sets) == 0L && length(terms) == 0L) {
    .sg_stop("sg_record() needs at least one set or quantity.")
  }
  structure(list(sets = sets, terms = terms, every = every), class = "sg_record")
}

#' @export
print.sg_record <- function(x, ...) {
  cat(
    "<sg_record>",
    paste0(
      if (length(x$sets)) paste0("sets: ", paste(x$sets, collapse = ", ")),
      if (!identical(x$every, "day") && !is.null(x$every)) paste0("; counter sums by ", x$every)
    ), "\n"
  )
  for (t in x$terms) {
    kind <- switch(class(t)[1L],
      sg_quantity = "quantity",
      sg_quantity_weights = "weights",
      sg_alleles = "allele copies",
      sg_counter = "counter sum"
    )
    cat("  ", kind, " \"", t$name, "\"",
      if (length(t$by)) paste0(" by ", paste(t$by, collapse = ", ")), "\n",
      sep = ""
    )
  }
  invisible(x)
}

#' @rdname sg_record
#' @param name Name of the quantity or counter sum (no `"["`, `"]"` or `"@"`).
#' @param stage States to add up: blocks `"E"`, `"L"`, `"P"`, `"M"`, `"U"`,
#'   `"F"`, or `"adults"` (M, U, F), `"juveniles"` (E, L, P), `"all"`
#'   (default when `host` is not given).
#' @param sex Optional `"male"` (M) or `"female"` (U, F).
#' @param genotype,mate Optional genotype names; `mate` for mated females (F).
#' @param substage Optional sub-stage numbers.
#' @param node Optional node ids or numbers (models with a landscape).
#'   Default: the study nodes, as in the summaries ([sg_totals()]).
#' @param host Host blocks `"HS"`, `"HW"`, `"HI"` (host-wound model) instead of
#'   fly stages.
#' @param by Split into one quantity per group: any of `"node"`, `"genotype"`,
#'   `"stage"` and `"mate"` (`sg_counter()`: `"node"`). Names follow the state
#'   names, e.g. `"adults[AA]@north"`.
#' @export
sg_quantity <- function(name, stage = NULL, sex = NULL, genotype = NULL, mate = NULL,
                        substage = NULL, node = NULL, host = NULL, by = character(0)) {
  .sg_check_record_name(name)
  if (!is.null(host) && !is.null(stage)) .sg_stop("Give `stage` or `host`, not both.")
  if (!is.null(host)) {
    bad <- setdiff(host, .sg_host_blocks)
    if (length(bad) > 0L || !is.character(host)) {
      .sg_stop("`host` must be any of \"HS\", \"HW\", \"HI\".")
    }
  } else {
    stage <- .sg_record_blocks(if (is.null(stage)) "all" else stage)
  }
  if (!is.null(sex)) {
    sex <- .sg_match_choice(sex, c("male", "female"), "sex")
    if (!is.null(host)) .sg_stop("`sex` applies to fly stages, not to hosts.")
    stage <- intersect(stage, if (sex == "male") "M" else c("U", "F"))
    if (length(stage) == 0L) .sg_stop("No stage of the chosen `sex` in `stage`.")
  }
  by <- .sg_check_record_by(by, c("node", "genotype", "stage", "mate"))
  if ("mate" %in% by && !identical(stage, "F")) {
    .sg_stop("`by = \"mate\"` needs stage = \"F\" (only mated females have a mate).")
  }
  if (!is.null(host) && any(c("genotype", "mate") %in% by)) {
    .sg_stop("Host blocks have no genotype or mate to group by.")
  }
  if (!is.null(mate) && !identical(stage, "F")) .sg_stop("`mate` needs stage = \"F\".")
  structure(
    list(
      name = name, stage = stage, host = host, genotype = genotype, mate = mate,
      substage = substage, node = node, by = by
    ),
    class = "sg_quantity"
  )
}

#' @rdname sg_record
#' @param W For `sg_quantity_weights()`: a numeric matrix with one row per
#'   quantity and one column per state ([sg_state_names()]), or a data frame
#'   with columns `quantity`, `state` (state names or positions) and
#'   `weight`.
#' @param names Names of the quantities of a matrix `W` (default its row
#'   names).
#' @export
sg_quantity_weights <- function(W, names = NULL) { # nolint: object_name_linter.
  if (is.data.frame(W)) {
    miss <- setdiff(c("quantity", "state", "weight"), names(W))
    if (length(miss) > 0L) {
      .sg_stop("A weight table needs columns `quantity`, `state` and `weight`.")
    }
    if (!is.numeric(W$weight) || any(!is.finite(W$weight))) {
      .sg_stop("`W$weight` must be finite numbers.")
    }
    q <- unique(as.character(W$quantity))
    for (nm in q) .sg_check_record_name(nm)
    return(structure(list(name = q, table = W, by = character(0)), class = "sg_quantity_weights"))
  }
  if (!is.matrix(W) || !is.numeric(W) || any(!is.finite(W))) {
    .sg_stop("`W` must be a finite numeric matrix (quantities x states) or a weight table.")
  }
  if (is.null(names)) names <- rownames(W)
  if (is.null(names) || length(names) != nrow(W) || anyDuplicated(names)) {
    .sg_stop("Give one distinct name per row of `W` (`names` or row names).")
  }
  for (nm in names) .sg_check_record_name(nm)
  structure(list(name = names, matrix = unname(W), by = character(0)),
    class = "sg_quantity_weights"
  )
}

#' @rdname sg_record
#' @param alleles For `sg_alleles()`: optional matrix of allele copies per
#'   genotype (rows named by genotype, columns by allele), as in
#'   [sg_allele_freq()]; by default each genotype name is split into
#'   single-letter alleles.
#' @details `sg_alleles(name, ...)` records the allele copies carried by the
#'   individuals of `stage` (default adults; mated females count their own
#'   genotype), one quantity per allele named `name[allele]` (with
#'   `by = "node"` also `@node`). [sg_allele_freq()] reads them when the run
#'   kept neither states nor totals.
#' @export
sg_alleles <- function(name, alleles = NULL, stage = "adults", sex = NULL, node = NULL,
                       by = character(0)) {
  .sg_check_record_name(name)
  stage <- .sg_record_blocks(stage)
  if (!is.null(sex)) {
    sex <- .sg_match_choice(sex, c("male", "female"), "sex")
    stage <- intersect(stage, if (sex == "male") "M" else c("U", "F"))
    if (length(stage) == 0L) .sg_stop("No stage of the chosen `sex` in `stage`.")
  }
  if (!is.null(alleles) && (!is.matrix(alleles) || !is.numeric(alleles))) {
    .sg_stop("`alleles` must be a numeric matrix (rows: genotypes, columns: alleles).")
  }
  structure(
    list(
      name = name, alleles = alleles, stage = stage, node = node,
      by = .sg_check_record_by(by, "node")
    ),
    class = "sg_alleles"
  )
}

#' @rdname sg_record
#' @param counter Base names of the daily counters to add up, e.g.
#'   `"fertile_batches"` or `c("deaths_M", "senescence_M")` (see the columns
#'   of `sim$counters`, without the `"@"` node suffix).
#' @export
sg_counter <- function(name, counter, node = NULL, by = character(0)) {
  .sg_check_record_name(name)
  if (!is.character(counter) || length(counter) < 1L || anyNA(counter)) {
    .sg_stop("`counter` must be counter names, e.g. \"fertile_batches\".")
  }
  structure(
    list(name = name, counter = counter, node = node, by = .sg_check_record_by(by, "node")),
    class = "sg_counter"
  )
}

# ---- internal ------------------------------------------------------------------

.sg_record_sets <- c("state", "totals", "counters", "cases")

.sg_check_record_name <- function(name) {
  ok <- is.character(name) && length(name) == 1L && !is.na(name) && nzchar(name) &&
    !grepl("[][@]", name)
  if (!ok) .sg_stop("A record name must be one non-empty string without \"[\", \"]\" or \"@\".")
  invisible(TRUE)
}

.sg_check_record_by <- function(by, allowed) {
  if (!is.character(by) || any(!by %in% allowed)) {
    .sg_stop("`by` must contain only ", paste0("\"", allowed, "\"", collapse = ", "), ".")
  }
  unique(by)
}

.sg_record_blocks <- function(stage) {
  short <- list(
    all = c("E", "L", "P", "M", "U", "F"), adults = c("M", "U", "F"),
    juveniles = c("E", "L", "P")
  )
  out <- unique(unlist(lapply(stage, function(s) if (s %in% names(short)) short[[s]] else s)))
  bad <- setdiff(out, short$all)
  if (!is.character(stage) || length(bad) > 0L) {
    .sg_stop(
      "Unknown stage(s) in `stage`: ", paste(bad, collapse = ", "),
      ". Use E, L, P, M, U, F, \"adults\", \"juveniles\" or \"all\"."
    )
  }
  out
}

# The recorder of a run (phase 1 WP9): from `record` ("state", "totals", or
# an sg_record) and the model, whether full states and counters are kept and
# the weighted sums W_x (totals first, then quantities) and W_c in CSR form
# (0-based for C++), with their names. `counter_names` are all counters of
# the run (with removed_* under the node extinction rule).
.sg_recorder <- function(record, model, counter_names) {
  if (is.character(record)) {
    if (length(record) != 1L || !record %in% c("state", "totals")) {
      .sg_stop("`record` must be \"state\", \"totals\" or an sg_record() object.")
    }
    record <- sg_record(record, "counters")
  }
  if (!inherits(record, "sg_record")) {
    .sg_stop("`record` must be \"state\", \"totals\" or an sg_record() object.")
  }
  sums <- Filter(function(t) inherits(t, "sg_counter"), record$terms)
  if ("cases" %in% record$sets) { # the counters sg_cases() reads, every node
    cases <- if (is.null(model$hosts)) {
      c("fertile_batches", "egg_batches")
    } else {
      c("infestations", "egg_batches", "superinfestations")
    }
    all_nodes <- seq_along(.sg_node_ids(model))
    sums <- c(sums, lapply(cases, function(cn) {
      sg_counter(cn, counter = cn, node = all_nodes, by = "node")
    }))
  }
  quantities <- Filter(function(t) !inherits(t, "sg_counter"), record$terms)
  rows_x <- list()
  if ("totals" %in% record$sets) {
    agg <- .sg_totals_matrix(model)
    rows_x <- lapply(seq_len(ncol(agg)), function(j) {
      list(name = colnames(agg)[j], index = which(agg[, j] != 0), weight = 1)
    })
    n_totals <- ncol(agg)
  } else {
    n_totals <- 0L
  }
  allele_info <- list()
  for (q in quantities) {
    if (inherits(q, "sg_alleles")) {
      a <- .sg_allele_rows(q, model)
      rows_x <- c(rows_x, a$rows)
      allele_info[[length(allele_info) + 1L]] <- a$info
    } else {
      rows_x <- c(rows_x, .sg_quantity_rows(q, model))
    }
  }
  rows_c <- unlist(lapply(sums, .sg_counter_rows, model = model, counter_names = counter_names),
    recursive = FALSE
  )
  names_x <- vapply(rows_x, `[[`, "", "name")
  names_c <- vapply(rows_c, `[[`, "", "name")
  dup <- c(names_x[duplicated(names_x)], names_c[duplicated(names_c)])
  if (length(dup) > 0L) .sg_stop("Recorded names must be distinct; repeated: ", dup[1L], ".")
  every <- if (is.null(record$every)) "day" else record$every
  per <- .sg_counter_periods(model$forcing$dates, every)
  list(
    state = "state" %in% record$sets,
    counters = "counters" %in% record$sets || length(sums) == 0L,
    n_totals = n_totals,
    names_x = names_x, names_c = names_c,
    wx = .sg_rows_csr(rows_x), wc = .sg_rows_csr(rows_c),
    every = every, c_period = as.integer(per$index - 1L),
    counter_periods = per$periods,
    alleles = allele_info,
    record = record
  )
}

# Nodes of a quantity or counter sum: those given, or the study nodes.
.sg_record_nodes <- function(model, node) {
  if (!is.null(node) || is.null(model$landscape)) {
    return(.sg_node_numbers(model, node))
  }
  which(model$landscape$nodes$role == "study")
}

# Periods of counter sums over the days of a run with these dates: the
# 1-based period of each day and a data frame of each period's first day
# (`period`) and number of days (`days`). Weeks are 7-day blocks from the
# first day (as sg_cases(period = "week")); months and years are calendar
# periods.
.sg_counter_periods <- function(dates, every) {
  key <- switch(EXPR = every,
    day = seq_along(dates),
    week = (seq_along(dates) - 1L) %/% 7L,
    month = format(dates, "%Y-%m"),
    year = format(dates, "%Y")
  )
  index <- match(key, unique(key))
  first <- !duplicated(index)
  list(
    index = index,
    periods = data.frame(period = dates[first], days = as.integer(tabulate(index)))
  )
}

# CSR form (0-based columns) of rows given as list(name, index, weight).
.sg_rows_csr <- function(rows) {
  len <- vapply(rows, function(r) length(r$index), 1L)
  list(
    row_ptr = as.integer(c(0L, cumsum(len))),
    col = as.integer(unlist(lapply(rows, `[[`, "index"))) - 1L,
    w = as.numeric(unlist(lapply(rows, function(r) rep_len(r$weight, length(r$index)))))
  )
}

# W_x . x with the CSR form of .sg_rows_csr() (R paths of simulate_ode()).
.sg_csr_multiply <- function(m, v) {
  n <- length(m$row_ptr) - 1L
  if (n == 0L) {
    return(numeric(0))
  }
  row <- rep.int(seq_len(n), diff(m$row_ptr))
  out <- numeric(n)
  if (length(m$col) > 0L) {
    s <- tapply(m$w * v[m$col + 1L], factor(row, levels = seq_len(n)), sum)
    out <- as.numeric(s)
    out[is.na(out)] <- 0
  }
  out
}

# Rows of W_x for one sg_alleles term: for each allele (and node with
# by = "node") the states of the chosen stages and nodes weighted by their
# genotype's copies of the allele (mated females by their own genotype).
# Also the term's description for sg_allele_freq(): name, blocks, alleles,
# node ids and whether it is split by node.
.sg_allele_rows <- function(a, model) {
  geno <- model$cube$genotypes
  amat <- if (is.null(a$alleles)) {
    .sg_default_alleles(geno)
  } else {
    .sg_check_allele_matrix(a$alleles, geno)[geno, , drop = FALSE]
  }
  ids <- .sg_node_ids(model)
  nodes <- .sg_record_nodes(model, a$node)
  by_node <- "node" %in% a$by && length(ids) > 1L
  groups <- if (by_node) as.list(nodes) else list(nodes)
  rows <- list()
  for (grp in groups) {
    idx <- integer(0)
    gen <- character(0)
    for (k in grp) {
      for (b in a$stage) {
        for (g in geno) {
          i <- sg_state_index(model, b, genotype = g, node = k)
          idx <- c(idx, i)
          gen <- c(gen, rep(g, length(i)))
        }
      }
    }
    for (al in colnames(amat)) {
      w <- amat[gen, al]
      keep <- w > 0
      rows[[length(rows) + 1L]] <- list(
        name = paste0(a$name, "[", al, "]", if (by_node) paste0("@", ids[grp]) else ""),
        index = idx[keep], weight = unname(w[keep])
      )
    }
  }
  list(
    rows = rows,
    info = list(
      name = a$name, blocks = a$stage, alleles = colnames(amat), nodes = ids[nodes],
      by_node = by_node
    )
  )
}

# Rows of W_x for one sg_quantity or sg_quantity_weights.
.sg_quantity_rows <- function(q, model) {
  if (inherits(q, "sg_quantity_weights")) {
    return(.sg_weight_rows(q, model))
  }
  ids <- .sg_node_ids(model)
  nodes <- .sg_record_nodes(model, q$node)
  geno <- model$cube$genotypes
  pick <- function(x, what) {
    if (is.null(x)) {
      return(geno)
    }
    bad <- setdiff(x, geno)
    if (length(bad) > 0L) {
      .sg_stop(
        "Unknown ", what, " \"", bad[1L], "\" in sg_quantity(\"", q$name, "\").",
        .sg_suggest(bad[1L], geno)
      )
    }
    x
  }
  genos <- pick(q$genotype, "genotype")
  mates <- if (identical(q$stage, "F")) pick(q$mate, "mate genotype") else NA_character_
  blocks <- if (is.null(q$host)) q$stage else q$host
  if (!is.null(q$host) && is.null(model$hosts)) {
    .sg_stop("sg_quantity(\"", q$name, "\") uses host blocks, but the model has no hosts.")
  }
  # every selected state with its group labels
  cells <- list()
  for (b in blocks) {
    for (k in nodes) {
      if (b %in% .sg_host_blocks) {
        idx <- sg_state_index(model, b, substage = q$substage, node = k)
        cells[[length(cells) + 1L]] <- data.frame(
          index = idx, stage = b, genotype = NA, mate = NA, node = ids[k]
        )
        next
      }
      for (g in genos) {
        for (m in mates) {
          idx <- sg_state_index(model, b,
            genotype = g, mate = if (b == "F" && !is.na(m)) m else NULL,
            substage = q$substage, node = k
          )
          cells[[length(cells) + 1L]] <- data.frame(
            index = idx, stage = b, genotype = g, mate = if (b == "F") m else NA, node = ids[k]
          )
        }
      }
    }
  }
  cells <- do.call(rbind, cells)
  if (identical(mates, NA_character_) || !identical(q$stage, "F")) cells$mate <- NA
  cells <- cells[!duplicated(cells$index), , drop = FALSE]
  inside <- setdiff(q$by, "node")
  key <- if (length(inside)) do.call(paste, c(cells[inside], sep = ",")) else rep("", nrow(cells))
  nd <- if ("node" %in% q$by && length(ids) > 1L) paste0("@", cells$node) else rep("", nrow(cells))
  label <- paste0(q$name, ifelse(nzchar(key), paste0("[", key, "]"), ""), nd)
  lapply(unique(label), function(l) {
    list(name = l, index = sort(cells$index[label == l]), weight = 1)
  })
}

.sg_weight_rows <- function(q, model) {
  if (!is.null(q$matrix)) {
    if (ncol(q$matrix) != model$n_states) {
      .sg_stop(
        "sg_quantity_weights(): `W` has ", ncol(q$matrix), " columns but the model has ",
        model$n_states, " states."
      )
    }
    return(lapply(seq_len(nrow(q$matrix)), function(r) {
      idx <- which(q$matrix[r, ] != 0)
      list(name = q$name[r], index = idx, weight = q$matrix[r, idx])
    }))
  }
  tab <- q$table
  st <- if (is.numeric(tab$state)) {
    as.integer(tab$state)
  } else {
    match(as.character(tab$state), sg_state_names(model))
  }
  if (anyNA(st) || any(st < 1L | st > model$n_states)) {
    .sg_stop(
      "sg_quantity_weights(): unknown state(s) in the weight table, e.g. \"",
      tab$state[is.na(st) | st < 1L][1L], "\"."
    )
  }
  lapply(q$name, function(nm) {
    sel <- as.character(tab$quantity) == nm
    w <- tapply(tab$weight[sel], st[sel], sum)
    list(name = nm, index = as.integer(names(w)), weight = as.numeric(w))
  })
}

# Rows of W_c for one sg_counter: the named counters, of the chosen nodes,
# summed (by node with by = "node").
.sg_counter_rows <- function(t, model, counter_names) {
  ids <- .sg_node_ids(model)
  nodes <- ids[.sg_record_nodes(model, t$node)]
  base <- sub("@.*$", "", counter_names)
  node_of <- if (length(ids) > 1L) .sg_name_node(counter_names) else rep(ids, length(counter_names))
  bad <- setdiff(t$counter, base)
  if (length(bad) > 0L) {
    .sg_stop(
      "Unknown counter \"", bad[1L], "\" in sg_counter(\"", t$name, "\").",
      .sg_suggest(bad[1L], unique(base)), " Counters: ", .sg_wrap_names(unique(base)), "."
    )
  }
  sel <- base %in% t$counter & node_of %in% nodes
  if ("node" %in% t$by && length(ids) > 1L) {
    lapply(nodes, function(id) {
      list(name = paste0(t$name, "@", id), index = which(sel & node_of == id), weight = 1)
    })
  } else {
    list(list(name = t$name, index = which(sel), weight = 1))
  }
}
