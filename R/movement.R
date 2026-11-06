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

#' Local diffusion between nodes
#'
#' Builds a movement layer for local, diffusion-like movement on any network
#' of nodes (a transect, a grid, counties) from a diffusion coefficient and
#' the geometry of the nodes, by the finite-volume formula (design note
#' dev/MOVEMENT_DESIGN.md section 6.1, D-032):
#' \deqn{\textrm{rate}_{ij} = D L_{ij} / (d_{ij} A_i)}{rate[i, j] = D L[i, j] / (d[i, j] A[i])}
#' per day at the reference temperature `movement_t_ref`, with \eqn{A_i}{A[i]}
#' the area of node i (km^2), \eqn{L_{ij}}{L[i, j]} the length of the border
#' it shares with node j (km) and \eqn{d_{ij}}{d[i, j]} the distance between
#' their centroids (km). The leave rate of node i is the sum of its rates and
#' the destination probabilities are the rates divided by it.
#'
#' @param area Node areas in km^2, in node order (or named by node id).
#' @param border Shared border lengths in km: a symmetric square matrix
#'   (zero for nodes that do not touch, zero diagonal; with node ids as
#'   dimnames or in node order), or an edge table with columns `from`, `to`
#'   (node ids or numbers), `border` and `dist`, each pair of neighbours
#'   listed once or in both directions with equal values.
#' @param dist Distances between node centroids in km: a square matrix like
#'   `border` (only entries of touching nodes are used), or `NULL` when
#'   `border` is an edge table.
#' @param D Diffusion coefficient in km^2 per day at `movement_t_ref`: one
#'   value for all adults, or a named list (or vector) with `male`, `virgin`
#'   and `mated`.
#' @param outer_border Optional length of each node's border with the
#'   outside of the landscape (km, in node order; 0 for inner nodes). Flies
#'   crossing it leave the landscape (open edges, counted as lost). `NULL`
#'   (default): closed edges.
#' @param outer_dist Distance from each node's centroid to a mirrored node
#'   outside (km). Default `area / outer_border`, twice the
#'   centroid-to-border distance of a rectangular cell whose outer side is
#'   `outer_border` long.
#' @param clock Movement clock, as in [sg_movement_layer()].
#'
#' @details
#' The package computes nothing from shapes: compute areas, shared border
#' lengths and centroid distances with your own tools (for example
#' `sf::st_area()`, `sf::st_touches()` with the lengths of the shared
#' boundaries, and `sf::st_distance()` between centroids).
#'
#' On a transect of cells with spacing h (and any width), each neighbour
#' gets \eqn{D/h^2}{D / h^2}: an interior cell leaves at \eqn{2D/h^2}{2 D / h^2},
#' a closed end at \eqn{D/h^2}{D / h^2} and an open end (`outer_border` = the
#' cell width) at \eqn{2D/h^2}{2 D / h^2}, half of it lost. On a regular grid of
#' square cells an interior cell leaves at \eqn{4D/h^2}{4 D / h^2}. The
#' destination probabilities do not depend on D, so the three classes share
#' them, as a layer requires.
#'
#' Rates depend on the grid: the same D gives \eqn{D/h^2}{D / h^2}, so a
#' 50 km grid has 2,500 times smaller rates than a 1 km grid. Use
#' `check_movement()` (phase 3) to see the lifetime dispersal a layer
#' implies.
#'
#' @return An `sg_movement_layer` ([sg_movement_layer()]) with `n_nodes`
#'   set, which can be inspected and changed like any other layer.
#' @seealso [sg_movement_layer()], [sg_landscape()], [sg_model()]
#' @export
#' @examples
#' # a transect of 5 cells, 50 km long and 50 km wide, open at the north end
#' h <- 50
#' n <- 5
#' border <- matrix(0, n, n)
#' border[cbind(1:(n - 1), 2:n)] <- h
#' border <- border + t(border)
#' dist <- matrix(h, n, n)
#' layer <- sg_movement_diffusion(
#'   area = rep(h * h, n), border = border, dist = dist,
#'   D = list(male = 0.25, virgin = 0.4, mated = 0.4),
#'   outer_border = c(0, 0, 0, 0, h)
#' )
#' layer
#' layer$rate$mated # 0.4 / 50^2 at the south end, twice that inside
#'
#' # the same transect as an edge table
#' edges <- data.frame(from = 1:4, to = 2:5, border = h, dist = h)
#' sg_movement_diffusion(rep(h * h, n), edges, D = 0.4, outer_border = c(0, 0, 0, 0, h))
sg_movement_diffusion <- function(area, border, dist = NULL,
                                  D, # nolint: object_name_linter.
                                  outer_border = NULL, outer_dist = NULL,
                                  clock = c("physiological", "calendar")) {
  clock <- .sg_match_choice(clock, c("physiological", "calendar"), "clock")
  if (!is.numeric(area) || length(area) < 1L || any(!is.finite(area)) || any(area <= 0)) {
    .sg_stop("`area` must be positive node areas in km^2.")
  }
  n <- length(area)
  geo <- .sg_diffusion_geometry(border, dist, n, names(area))
  ids <- geo$ids
  coef <- .sg_check_diffusion_coef(D)
  # outer border: open edges
  w_out <- numeric(n)
  if (!is.null(outer_border)) {
    ok <- is.numeric(outer_border) && length(outer_border) == n && all(is.finite(outer_border)) &&
      all(outer_border >= 0)
    if (!ok) .sg_stop("`outer_border` must be ", n, " lengths >= 0 in km (one per node).")
    if (is.null(outer_dist)) {
      outer_dist <- ifelse(outer_border > 0, area / pmax(outer_border, 1e-300), Inf)
    }
    ok <- is.numeric(outer_dist) && length(outer_dist) == n && !anyNA(outer_dist) &&
      all(outer_dist[outer_border > 0] > 0 & is.finite(outer_dist[outer_border > 0]))
    if (!ok) {
      .sg_stop("`outer_dist` must be ", n, " positive distances in km (one per node).")
    }
    open <- outer_border > 0
    w_out[open] <- outer_border[open] / (outer_dist[open] * area[open])
  } else if (!is.null(outer_dist)) {
    .sg_stop("`outer_dist` needs `outer_border` (the border with the outside).")
  }
  # per unit D: rate i -> j = L_ij / (d_ij A_i)
  e <- geo$edges
  e$w <- e$border / (e$dist * area[e$i])
  w_in <- vapply(seq_len(n), function(i) sum(e$w[e$i == i]), 0)
  w_tot <- w_in + w_out
  prob <- e$w / w_tot[e$i]
  dest <- data.frame(from = ids[e$i], to = ids[e$j], prob = prob, stringsAsFactors = FALSE)
  rate <- lapply(coef, function(d) d * w_tot)
  layer <- sg_movement_layer(rate, dest, clock = clock)
  layer$n_nodes <- n
  layer
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

# Diffusion coefficients of sg_movement_diffusion(): list(male, virgin, mated)
# of single positive numbers.
.sg_check_diffusion_coef <- function(d) {
  if (is.numeric(d) && length(d) == 1L && is.null(names(d))) {
    d <- list(male = d, virgin = d, mated = d)
  }
  if (is.numeric(d) && !is.null(names(d))) d <- as.list(d)
  ok <- is.list(d) && setequal(names(d), .sg_movement_classes) && length(d) == 3L &&
    all(vapply(d, function(v) is.numeric(v) && length(v) == 1L && is.finite(v) && v >= 0, TRUE))
  if (!ok) {
    .sg_stop(
      "`D` must be one diffusion coefficient >= 0 (km^2 per day) or a named list with ",
      "`male`, `virgin` and `mated`."
    )
  }
  d[.sg_movement_classes]
}

# Geometry of sg_movement_diffusion(): node ids and the directed edges (i, j
# as positions, border, dist) between touching nodes, from matrices or an
# edge table.
.sg_diffusion_geometry <- function(border, dist, n, area_names) {
  sym <- function(m) isTRUE(all.equal(m, t(m), tolerance = 1e-12, check.attributes = FALSE))
  if (is.matrix(border)) {
    if (!is.numeric(border) || !identical(dim(border), c(n, n))) {
      .sg_stop("`border` must be a ", n, " x ", n, " numeric matrix (one row per node of `area`).")
    }
    if (any(!is.finite(border)) || any(border < 0) || any(diag(border) != 0) || !sym(border)) {
      .sg_stop("`border` must hold finite lengths >= 0, be symmetric and have a zero diagonal.")
    }
    if (!is.matrix(dist) || !is.numeric(dist) || !identical(dim(dist), c(n, n))) {
      .sg_stop("`dist` must be a ", n, " x ", n, " numeric matrix of centroid distances.")
    }
    ids <- rownames(border)
    if (!identical(ids, colnames(border))) {
      .sg_stop("Row and column names of `border` must be the same node ids.")
    }
    if (is.null(ids)) ids <- area_names
    if (!is.null(area_names) && !identical(ids, area_names)) {
      .sg_stop("The names of `area` and the dimnames of `border` must be the same node ids.")
    }
    nz <- which(border > 0, arr.ind = TRUE)
    nz <- nz[order(nz[, 1L], nz[, 2L]), , drop = FALSE]
    d <- dist[nz]
    if (any(!is.finite(d)) || any(d <= 0) || !sym(dist * (border > 0))) {
      .sg_stop("`dist` must be positive, finite and symmetric for every pair of touching nodes.")
    }
    edges <- data.frame(i = nz[, 1L], j = nz[, 2L], border = border[nz], dist = d)
  } else if (is.data.frame(border)) {
    if (!is.null(dist)) .sg_stop("With an edge table `border`, leave `dist` NULL (use its column).")
    miss <- setdiff(c("from", "to", "border", "dist"), names(border))
    if (length(miss) > 0L) {
      .sg_stop(
        "An edge table `border` needs columns `from`, `to`, `border` and `dist`; missing: ",
        paste(miss, collapse = ", "), "."
      )
    }
    ids <- area_names
    pos <- function(v) {
      if (is.numeric(v)) {
        if (any(v < 1 | v > n | v != round(v), na.rm = TRUE) || anyNA(v)) {
          .sg_stop("Node numbers in `border` must be whole numbers between 1 and ", n, ".")
        }
        return(as.integer(v))
      }
      if (is.null(ids)) {
        .sg_stop("Node ids in `border` need `area` named by node id (or use node numbers).")
      }
      p <- match(as.character(v), ids)
      if (anyNA(p)) {
        bad <- as.character(v)[is.na(p)][1L]
        .sg_stop("Unknown node \"", bad, "\" in `border`.", .sg_suggest(bad, ids))
      }
      p
    }
    i <- pos(border$from)
    j <- pos(border$to)
    b <- border$border
    d <- border$dist
    ok <- is.numeric(b) && is.numeric(d) && all(is.finite(b)) && all(b >= 0) &&
      all(is.finite(d)) && all(d > 0)
    if (!ok) .sg_stop("`border$border` must be lengths >= 0 and `border$dist` positive (km).")
    if (any(i == j)) .sg_stop("`border` must not pair a node with itself.")
    both <- data.frame(i = c(i, j), j = c(j, i), border = c(b, b), dist = c(d, d))
    key <- paste(both$i, both$j)
    first <- !duplicated(key)
    # a pair given in both directions must agree
    agree <- vapply(split(seq_len(nrow(both)), key), function(r) {
      length(unique(both$border[r])) == 1L && length(unique(both$dist[r])) == 1L
    }, TRUE)
    if (!all(agree)) {
      .sg_stop("A pair of nodes appears in `border` with different `border` or `dist` values.")
    }
    edges <- both[first & both$border > 0, , drop = FALSE]
    edges <- edges[order(edges$i, edges$j), , drop = FALSE]
  } else {
    .sg_stop(
      "`border` must be a numeric matrix or a data frame with columns `from`, `to`, ",
      "`border`, `dist`."
    )
  }
  rownames(edges) <- NULL
  if (is.null(ids)) ids <- seq_len(n)
  list(ids = ids, edges = edges)
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
