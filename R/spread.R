# Spread metrics of a simulation on a landscape (multi-node work phase 3,
# WP6; design note dev/MOVEMENT_DESIGN.md section 7.2, D-032).

#' Spread of flies through a landscape
#'
#' When each node was first reached and how far the flies have spread: per
#' node the first day with flies, the last day with flies and the number of
#' output days with flies; over time the number of occupied nodes and the
#' position of the front.
#'
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @param threshold A node counts as occupied on an output day when its
#'   flies (eggs, larvae, pupae and adults of all genotypes) total at least
#'   `threshold`. Default the model's `extinction_threshold` (1 individual),
#'   as the node extinction rule.
#' @param nodes Nodes to include (ids or numbers; default the study nodes).
#' @param position Position of each node along the spread, e.g. km along a
#'   transect: one value per node of the landscape (in node order, or named
#'   by node id). Default the nodes' `x` coordinates ([sg_landscape()]); the
#'   front is `NA` when there are none.
#' @param quantity Optional name of a recorded quantity with one value per
#'   node (`sg_quantity(quantity, stage = "all", by = "node")` in
#'   [sg_record()]), used instead of the states or totals.
#'
#' @details Days are output days (`output_every` in the simulation), with
#'   day 0 the initial state: nodes occupied at the start have
#'   `first_arrival` 0. The front on a day is the largest `position` of an
#'   occupied node. For stochastic runs every replicate is reported.
#'
#'   Fractions of a fly in the ODE (without the node extinction rule) count
#'   once they reach `threshold`.
#'
#' @return An object of class `sg_spread`, a list with
#'   \describe{
#'     \item{`nodes`}{Data frame: `replicate` (stochastic runs), `node`,
#'       `first_arrival` (first output day occupied, `NA` if never),
#'       `last_occupied` (last output day occupied) and `days_occupied`
#'       (number of output days occupied, day 0 excluded).}
#'     \item{`time`}{Data frame: `day`, `date`, `replicate` (stochastic
#'       runs), `occupied` (number of occupied nodes) and `front`.}
#'     \item{`threshold`, `position`}{The threshold and the positions used.}
#'   }
#' @seealso [sg_totals()], [plot.sg_sim()] (`what = "nodes"` and `"front"`)
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_constant_weather(27, days = 120), p)
#'   n <- 6
#'   edges <- data.frame(from = 1:(n - 1), to = 2:n, border = 10, dist = 10)
#'   land <- sg_landscape(data.frame(id = 1:n, area = 100, x = (1:n - 1) * 10), wx)
#'   layer <- sg_movement_diffusion(rep(100, n), edges, D = 5)
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), landscape = land, movement = layer)
#'   x0 <- sg_init(mod, values = c(E = 50, L = 50, P = 50, M = 20, U = 5, F = 20), nodes = 1)
#'   sim <- simulate_cb(mod, x0, n_rep = 3, seed = 1)
#'   sp <- sg_spread(sim)
#'   sp
#'   head(sp$nodes)
#' }
sg_spread <- function(sim, threshold = NULL, nodes = NULL, position = NULL, quantity = NULL) {
  info <- .sg_sim_info(sim)
  if (is.null(threshold)) {
    threshold <- if (is.null(info$extinction_threshold)) 1 else info$extinction_threshold
  }
  ok <- is.numeric(threshold) && length(threshold) == 1L && is.finite(threshold) && threshold > 0
  if (!ok) .sg_stop("`threshold` must be one positive number of individuals.")
  sel <- .sg_sim_nodes(sim, nodes)
  flies <- .sg_node_flies(sim, sel, nodes, quantity) # [time, node, replicate]
  pos <- .sg_spread_position(info, sel, position)
  occ <- flies >= threshold
  n_t <- dim(occ)[1L]
  n_rep <- dim(occ)[3L]
  stochastic <- inherits(sim, "sg_sim_cb")
  day <- sim$day
  node_rows <- lapply(seq_len(n_rep), function(r) {
    o <- matrix(occ[, , r], n_t)
    d <- data.frame(
      node = sel$ids,
      first_arrival = apply(o, 2L, function(v) if (any(v)) day[which(v)[1L]] else NA_integer_),
      last_occupied = apply(o, 2L, function(v) if (any(v)) day[max(which(v))] else NA_integer_),
      days_occupied = as.integer(colSums(o[-1L, , drop = FALSE])),
      stringsAsFactors = FALSE
    )
    if (stochastic) d <- cbind(replicate = r, d)
    d
  })
  time_rows <- lapply(seq_len(n_rep), function(r) {
    o <- matrix(occ[, , r], n_t)
    front <- apply(o, 1L, function(v) if (any(v) && !all(is.na(pos))) max(pos[v]) else NA_real_)
    d <- data.frame(day = day, date = sim$date, occupied = as.integer(rowSums(o)), front = front)
    if (stochastic) d <- cbind(d[1:2], replicate = r, d[-(1:2)])
    d
  })
  out <- list(
    nodes = do.call(rbind, node_rows), time = do.call(rbind, time_rows),
    threshold = threshold, position = stats::setNames(pos, sel$ids)
  )
  rownames(out$nodes) <- NULL
  rownames(out$time) <- NULL
  structure(out, class = "sg_spread")
}

#' @export
print.sg_spread <- function(x, ...) {
  nd <- x$nodes
  n_node <- length(unique(nd$node))
  n_rep <- if (is.null(nd$replicate)) 1L else length(unique(nd$replicate))
  rep_of <- if (n_rep > 1L) nd$replicate else rep(1L, nrow(nd))
  reached <- tapply(!is.na(nd$first_arrival), rep_of, sum)
  cat(sprintf(
    "<sg_spread> %d node%s%s; occupied = at least %g flies\n", n_node,
    if (n_node == 1L) "" else "s",
    if (n_rep > 1L) sprintf(", %d replicates", n_rep) else "", x$threshold
  ))
  cat(sprintf(
    "  nodes ever occupied: %s of %d%s\n",
    if (n_rep > 1L) {
      sprintf("%.3g (mean; %d-%d)", mean(reached), min(reached), max(reached))
    } else {
      reached
    },
    n_node, if (n_rep > 1L) " per replicate" else ""
  ))
  tt <- x$time
  last <- tt[tt$day == max(tt$day), , drop = FALSE]
  cat(sprintf("  occupied at the end: %s\n", paste(range(last$occupied), collapse = "-")))
  if (!all(is.na(tt$front))) {
    cat(sprintf(
      "  front at the end: %s (position of the farthest occupied node)\n",
      paste(format(range(last$front, na.rm = TRUE), digits = 4L), collapse = " to ")
    ))
  }
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# Flies per node [time, node, replicate] of the chosen nodes: from a recorded
# per-node quantity, or the fly columns (not hosts) of the states or totals.
.sg_node_flies <- function(sim, sel, nodes, quantity) {
  if (!is.null(quantity)) {
    rec <- sim$recorded
    cols <- if (sel$multi) paste0(quantity, "@", sel$ids) else quantity
    if (is.null(rec) || !all(cols %in% dimnames(rec)[[2L]])) {
      .sg_stop(
        "No recorded quantity \"", quantity, "\" for every chosen node; record it with ",
        "sg_quantity(\"", quantity, "\", stage = \"all\", by = \"node\")."
      )
    }
    if (length(dim(rec)) == 2L) {
      rec <- array(rec, c(dim(rec), 1L), dimnames = list(NULL, colnames(rec), NULL))
    }
    return(rec[, cols, , drop = FALSE])
  }
  if (is.null(sim$state) && is.null(sim$totals)) {
    .sg_stop(
      "sg_spread() needs the states or totals of the run, or a recorded per-node quantity ",
      "(`quantity`, e.g. sg_quantity(\"flies\", stage = \"all\", by = \"node\"))."
    )
  }
  a <- .sg_sim_array(sim, nodes)
  fly <- !.sg_is_host_name(dimnames(a$array)[[2L]])
  per_node <- vapply(sel$ids, function(id) {
    apply(a$array[, fly & a$node == id, , drop = FALSE], c(1L, 3L), sum)
  }, matrix(0, dim(a$array)[1L], dim(a$array)[3L]))
  per_node <- array(per_node, c(dim(a$array)[1L], dim(a$array)[3L], length(sel$ids)))
  aperm(per_node, c(1L, 3L, 2L))
}

# Positions of the chosen nodes: `position` (per node of the landscape, in
# node order or named by id) or the nodes' x coordinates; NA without either.
.sg_spread_position <- function(info, sel, position) {
  all_ids <- info$nodes$id
  if (is.null(position)) {
    x <- info$nodes$x
    if (is.null(x)) x <- rep(NA_real_, length(all_ids))
    return(as.numeric(x[match(sel$ids, all_ids)]))
  }
  if (!is.numeric(position) || anyNA(position)) {
    .sg_stop("`position` must be numbers, one per node of the landscape.")
  }
  if (!is.null(names(position))) {
    p <- position[sel$ids]
    if (anyNA(p)) {
      .sg_stop("`position` has no value for node(s) ", .sg_wrap_names(sel$ids[is.na(p)]), ".")
    }
    return(unname(as.numeric(p)))
  }
  if (length(position) != length(all_ids)) {
    .sg_stop("`position` must have one value per node of the landscape (", length(all_ids), ").")
  }
  as.numeric(position[match(sel$ids, all_ids)])
}
