#' Movement layer
#'
#' Defines how adult flies move between the nodes of a landscape
#' (MODEL_SPEC section 8a, N1-N3, D-032). A layer gives, for each class of
#' adults (males, unmated females, mated females) and node, the rate at which
#' flies leave the node, and for each node the probabilities of the
#' destinations. Movers keep their genotype, mate genotype and age. Pass the
#' layer to [sg_model()] with `movement = `.
#'
#' @param rate Leave rates per day at the reference temperature
#'   `movement_t_ref` ([sg_params()], default 27 degC): a named list (or named
#'   numeric vector) with elements `male`, `virgin` (unmated females) and
#'   `mated` (mated females), each one value for all nodes or one per node in
#'   node order; >= 0. A daily leave fraction f measured in the field
#'   converts to a rate as -log(1 - f).
#' @param dest Destination probabilities pi(i, j) that a fly leaving node i
#'   goes to node j, shared by the three classes: either a numeric square
#'   matrix (rows: from, columns: to; with the node ids as dimnames, or in
#'   node order) or a data frame of edges with columns `from`, `to` (node ids)
#'   and `prob`. Entries are in \[0, 1\], the diagonal is 0 and each row sums
#'   to at most 1. The remainder 1 - sum_j pi(i, j) of the flies leaving node
#'   i leaves the landscape (an open edge) and is counted as lost; rows that
#'   sum to 1 are closed edges.
#' @param clock How the rates change with temperature: `"physiological"`
#'   (default): kappa = min(`movement_kappa_max`, dx(T) / dx(`movement_t_ref`))
#'   times the given rates, so flies move less in the cold and not at all at
#'   or below the development threshold; `"calendar"`: the given rates every
#'   day (for wind-borne movement and comparisons with MGDrivE2).
#'
#' @details No study relates screwworm movement to temperature; the cap
#'   `movement_kappa_max` (1 by default, `Inf` for uncapped degree-day
#'   scaling) and the clock are best treated as a sensitivity analysis
#'   (D-032, design note dev/MOVEMENT_DESIGN.md section 5.3). Version 0.2
#'   accepts one layer per model; several layers (e.g. local and long-distance
#'   movement) follow.
#'
#' @return An object of class `sg_movement_layer`: a list with `rate`
#'   (list of three numeric vectors), `edges` (data frame `from`, `to`,
#'   `prob`, ordered by `from` then `to`; ids as text, or as numbers when
#'   `dest` had no names), `n_nodes` (from a matrix `dest`, else `NA`) and
#'   `clock`.
#' @seealso [sg_landscape()], [sg_model()]
#' @export
#' @examples
#' # three nodes in a row; flies leave to their neighbours, closed ends
#' dest <- rbind(
#'   c(0, 1, 0),
#'   c(0.5, 0, 0.5),
#'   c(0, 1, 0)
#' )
#' layer <- sg_movement_layer(
#'   rate = list(male = 0.05, virgin = 0.1, mated = -log(1 - 0.16)),
#'   dest = dest
#' )
#' layer
#'
#' # the same as an edge table with node ids
#' edges <- data.frame(
#'   from = c("a", "b", "b", "c"), to = c("b", "a", "c", "b"),
#'   prob = c(1, 0.5, 0.5, 1)
#' )
#' sg_movement_layer(list(male = 0.05, virgin = 0.1, mated = 0.17), edges, clock = "calendar")
sg_movement_layer <- function(rate, dest, clock = c("physiological", "calendar")) {
  clock <- .sg_match_choice(clock, c("physiological", "calendar"), "clock")
  d <- .sg_dest_edges(dest)
  rate <- .sg_check_rates(rate, d$n_nodes)
  structure(
    list(rate = rate, edges = d$edges, n_nodes = d$n_nodes, clock = clock),
    class = "sg_movement_layer"
  )
}

#' @export
print.sg_movement_layer <- function(x, ...) {
  e <- x$edges
  out_sum <- tapply(e$prob, factor(e$from, levels = unique(e$from)), sum)
  open <- sum(out_sum < 1 - 1e-12)
  fmt_rate <- function(v) {
    if (length(v) == 1L) {
      return(format(v, digits = 4L))
    }
    paste0(format(min(v), digits = 4L), "-", format(max(v), digits = 4L))
  }
  cat(sprintf(
    "<sg_movement_layer> clock \"%s\"; %d edges from %d nodes (%d with open edges)\n",
    x$clock, nrow(e), length(out_sum), open
  ))
  cat(sprintf(
    "  leave rates per day at movement_t_ref: male %s, virgin %s, mated %s\n",
    fmt_rate(x$rate$male), fmt_rate(x$rate$virgin), fmt_rate(x$rate$mated)
  ))
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

.sg_movement_classes <- c("male", "virgin", "mated")

# Leave rates: list(male, virgin, mated), each length 1 or n (n = NA: not
# known yet, checked against the landscape in .sg_movement_groups()).
.sg_check_rates <- function(rate, n) {
  if (is.numeric(rate) && !is.null(names(rate))) rate <- as.list(rate)
  if (!is.list(rate) || is.null(names(rate))) {
    .sg_stop("`rate` must be a named list with elements `male`, `virgin` and `mated`.")
  }
  missing <- setdiff(.sg_movement_classes, names(rate))
  extra <- setdiff(names(rate), .sg_movement_classes)
  if (length(missing) > 0L || length(extra) > 0L) {
    .sg_stop(
      "`rate` must have exactly the elements `male`, `virgin` and `mated`",
      if (length(missing) > 0L) paste0("; missing: ", paste(missing, collapse = ", ")) else "",
      if (length(extra) > 0L) {
        paste0(
          "; unknown: ", paste(extra, collapse = ", "), ".",
          .sg_suggest(extra[1L], .sg_movement_classes)
        )
      } else {
        "."
      }
    )
  }
  out <- lapply(.sg_movement_classes, function(cl) {
    v <- rate[[cl]]
    if (!is.numeric(v) || length(v) < 1L || any(!is.finite(v)) || any(v < 0)) {
      .sg_stop("`rate$", cl, "` must be finite numbers >= 0 (leave rates per day).")
    }
    if (!is.na(n) && !length(v) %in% c(1L, n)) {
      .sg_stop(
        "`rate$", cl, "` must have one value or one per node (", n, "); got ", length(v), "."
      )
    }
    as.numeric(unname(v))
  })
  stats::setNames(out, .sg_movement_classes)
}

# Destinations as a sorted edge table (from, to, prob > 0) and the number of
# nodes when `dest` is a matrix.
.sg_dest_edges <- function(dest) {
  tol <- 1e-12
  if (is.matrix(dest)) {
    if (!is.numeric(dest) || nrow(dest) != ncol(dest) || nrow(dest) < 1L) {
      .sg_stop("`dest` must be a numeric square matrix (rows: from, columns: to).")
    }
    if (any(!is.finite(dest)) || any(dest < 0) || any(dest > 1)) {
      .sg_stop("Entries of `dest` must be finite probabilities in [0, 1].")
    }
    rn <- rownames(dest)
    cn <- colnames(dest)
    if (!identical(rn, cn)) {
      .sg_stop("Row and column names of `dest` must be the same node ids in the same order.")
    }
    if (any(diag(dest) != 0)) {
      .sg_stop("The diagonal of `dest` must be 0 (staying is not a move).")
    }
    nz <- which(dest > 0, arr.ind = TRUE)
    nz <- nz[order(nz[, 1L], nz[, 2L]), , drop = FALSE]
    ids <- if (is.null(rn)) seq_len(nrow(dest)) else rn
    edges <- data.frame(
      from = ids[nz[, 1L]], to = ids[nz[, 2L]], prob = dest[nz],
      stringsAsFactors = FALSE
    )
    n <- nrow(dest)
  } else if (is.data.frame(dest)) {
    miss <- setdiff(c("from", "to", "prob"), names(dest))
    if (length(miss) > 0L) {
      .sg_stop(
        "An edge table `dest` needs columns `from`, `to` and `prob`; missing: ",
        paste(miss, collapse = ", "), "."
      )
    }
    p <- dest$prob
    if (!is.numeric(p) || any(!is.finite(p)) || any(p < 0) || any(p > 1)) {
      .sg_stop("`dest$prob` must be finite probabilities in [0, 1].")
    }
    num_ids <- is.numeric(dest$from) && is.numeric(dest$to)
    from <- if (num_ids) as.integer(dest$from) else as.character(dest$from)
    to <- if (num_ids) as.integer(dest$to) else as.character(dest$to)
    if (anyNA(from) || anyNA(to)) .sg_stop("`dest$from` and `dest$to` must not be NA.")
    not_whole <- any(from != dest$from) || any(to != dest$to) || any(from < 1L) || any(to < 1L)
    if (num_ids && not_whole) {
      .sg_stop("Numeric node numbers in `dest` must be whole numbers >= 1.")
    }
    if (any(from == to)) .sg_stop("`dest` must not contain edges from a node to itself.")
    key <- paste(from, to, sep = "\r")
    if (anyDuplicated(key)) {
      dup <- sub("\r", " -> ", key[duplicated(key)][1L])
      .sg_stop("`dest` lists an edge more than once: ", dup, ".")
    }
    keep <- p > 0
    edges <- data.frame(
      from = from[keep], to = to[keep], prob = as.numeric(p[keep]),
      stringsAsFactors = FALSE
    )
    edges <- edges[order(edges$from, edges$to), , drop = FALSE]
    rownames(edges) <- NULL
    n <- NA_integer_
  } else {
    .sg_stop("`dest` must be a numeric matrix or a data frame with columns `from`, `to`, `prob`.")
  }
  if (nrow(edges) > 0L) {
    rs <- tapply(edges$prob, factor(edges$from, levels = unique(edges$from)), sum)
    if (any(rs > 1 + tol)) {
      bad <- names(rs)[rs > 1 + tol][1L]
      .sg_stop(
        "Destination probabilities from node ", bad, " sum to ", format(max(rs), digits = 6L),
        "; each row of `dest` must sum to at most 1 (the remainder leaves the landscape)."
      )
    }
  }
  list(edges = edges, n_nodes = n)
}

# Movement groups (MODEL_SPEC 8a "Layers and groups"): layers that share a
# clock (and, from phase 4, a schedule) are merged, with
#   delta^G[s, i] = sum_l delta^l[s, i]
#   pi^G[s, i, j] = sum_l delta^l[s, i] pi^l[i, j] / delta^G[s, i]
# A group with one layer keeps that layer's probabilities exactly. Each group:
# clock; rate [3, n] (rows male, virgin, mated); edges in CSR form (row_ptr,
# 1-based, length n + 1; col, 1-based); prob [3, n_edges]; loss [3, n], the
# share of leavers that leaves the landscape.
.sg_movement_groups <- function(layers, node_ids) {
  n <- length(node_ids)
  resolve <- function(v, what) {
    if (is.numeric(v)) {
      if (any(v > n)) {
        .sg_stop("Movement layer refers to node ", max(v), " but the landscape has ", n, " nodes.")
      }
      return(as.integer(v))
    }
    idx <- match(v, node_ids)
    if (anyNA(idx)) {
      bad <- v[is.na(idx)][1L]
      .sg_stop("Movement layer refers to unknown node \"", bad, "\".", .sg_suggest(bad, node_ids))
    }
    idx
  }
  prep <- lapply(layers, function(l) {
    if (!is.na(l$n_nodes) && l$n_nodes != n) {
      .sg_stop(
        "The movement layer's `dest` has ", l$n_nodes, " nodes but the landscape has ", n, "."
      )
    }
    rate <- vapply(.sg_movement_classes, function(cl) {
      v <- l$rate[[cl]]
      if (!length(v) %in% c(1L, n)) {
        .sg_stop(
          "`rate$", cl, "` must have one value or one per node (", n, "); got ", length(v), "."
        )
      }
      rep_len(v, n)
    }, numeric(n))
    from <- resolve(l$edges$from)
    to <- resolve(l$edges$to)
    list(clock = l$clock, rate = t(rate), key = (from - 1) * n + to, prob = l$edges$prob)
  })
  clocks <- vapply(prep, `[[`, "", "clock")
  lapply(unique(clocks), function(ck) {
    members <- prep[clocks == ck]
    rate <- Reduce(`+`, lapply(members, `[[`, "rate"))
    key <- sort(unique(unlist(lapply(members, `[[`, "key"))))
    from <- as.integer((key - 1) %/% n + 1)
    to <- as.integer((key - 1) %% n + 1)
    prob <- matrix(0, 3L, length(key), dimnames = list(.sg_movement_classes, NULL))
    loss <- matrix(0, 3L, n, dimnames = list(.sg_movement_classes, node_ids))
    for (s in 1:3) {
      moving <- rate[s, from] > 0
      if (length(members) == 1L) {
        m <- members[[1L]]
        pr <- m$prob[match(key, m$key)]
      } else {
        num <- numeric(length(key))
        for (m in members) {
          hit <- match(m$key, key)
          num[hit] <- num[hit] + m$rate[s, from[hit]] * m$prob
        }
        pr <- num / rate[s, from]
      }
      pr[!moving] <- 0
      prob[s, ] <- pr
      out <- numeric(n)
      out[sort(unique(from))] <- tapply(pr, from, sum)
      ls <- 1 - out
      ls[ls < 1e-12 | rate[s, ] == 0] <- 0 # rows that sum to 1 up to round-off are closed
      loss[s, ] <- ls
    }
    dimnames(rate) <- list(.sg_movement_classes, node_ids)
    list(
      clock = ck, rate = rate,
      row_ptr = as.integer(c(0L, cumsum(tabulate(from, n))) + 1L),
      col = to, prob = prob, loss = loss
    )
  })
}
