# check_movement(): leave rates translated back into biology (multi-node
# work phase 3, WP3; design note dev/MOVEMENT_DESIGN.md section 6.2, D-032).

#' Check a model's movement
#'
#' Translates a model's leave rates back into biology, so that rates on the
#' wrong scale for the grid show up (design note section 6.2, D-032): for
#' each node and class of adults at the reference temperature
#' `movement_t_ref`, the daily and lifetime probability of leaving the node,
#' the expected number of moves per lifetime and, when the landscape has
#' coordinates, the expected path length per lifetime and the diffusion
#' coefficient the rates imply.
#'
#' @param model An `sg_model` object with movement ([sg_model()] with
#'   `landscape` and `movement`).
#' @param temperatures Constant temperatures (degC) for a table of lifetime
#'   movement against temperature with the model's clock and
#'   `movement_kappa_max` (design note 5.3); `NULL` for none.
#' @param genotype Genotype whose adult lifespan is used (default the first
#'   wild-type genotype).
#'
#' @details
#' With the total leave rate \eqn{\lambda}{lambda} of a class in a node
#' (summed over the movement groups, each scaled by its clock kappa) and the
#' adult lifespan \eqn{\tau}{tau} (death by the adult hazard or old age at the
#' end of the adult chain, exact from the chain's sub-stage rates; moving
#' does not change it):
#' * `daily_leave` = \eqn{1 - e^{-\lambda}}{1 - exp(-lambda)};
#' * `lifetime_leave` = \eqn{1 - E[e^{-\lambda\tau}]}{1 - E[exp(-lambda tau)]},
#'   the probability of leaving at least once (for unmated females as if
#'   they never mated);
#' * `moves_per_life` = \eqn{\lambda E[\tau]}{lambda E[tau]};
#' * `lost_share`: share of the movers that leaves the landscape;
#' * with coordinates `x` (and `y`) in [sg_landscape()]: `path_per_life` =
#'   moves per lifetime times the mean distance of a move, and `D` =
#'   \eqn{\lambda \sum_j p_{ij} d_{ij}^2 / (2k)}{lambda sum_j p[i, j] d[i, j]^2
#'   / (2 k)}, with k = 1 when all nodes lie on a line and 2 otherwise
#'   (only moves within the landscape count).
#'
#' The temperature table gives, for each class and temperature, the clock's
#' kappa, the adult lifespan and the moves per lifetime and lifetime leave
#' probability in the node with the largest leave rate of that class.
#'
#' A lifetime leave probability above 0.99 is reported as a warning: flies
#' almost surely leave their node, which usually means that a rate meant for
#' a finer grid is used on a coarser one (or vice versa). [check_model()]
#' includes this check and also warns when movement within one step of the
#' chain-binomial solver is large.
#'
#' @return An object of class `sg_movement_check` (printed as a report), a
#'   list with `nodes` (one row per node and class), `life` (adult lifespan
#'   by class at `movement_t_ref`), `temperatures` (or `NULL`), `t_ref`,
#'   `kappa_max`, `clocks`, `dimension` and `warnings`.
#' @seealso [sg_movement_layer()], [sg_movement_diffusion()], [check_model()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_constant_weather(27, days = 30), p)
#'   # 5 cells of 10 km in a row, D = 0.4 km^2 per day for all adults
#'   h <- 10
#'   edges <- data.frame(from = 1:4, to = 2:5, border = h, dist = h)
#'   layer <- sg_movement_diffusion(rep(h^2, 5), edges, D = 0.4)
#'   land <- sg_landscape(data.frame(id = 1:5, area = h^2, x = (0:4) * h), wx)
#'   mod <- sg_model(p, sg_cube(MGDrivE::cubeMendelian()), landscape = land, movement = layer)
#'   check_movement(mod)
#' }
check_movement <- function(model, temperatures = c(18, 22, 27, 30, 35), genotype = NULL) {
  if (!inherits(model, "sg_model")) {
    .sg_stop("`model` must be an sg_model object; create one with sg_model().")
  }
  mv <- model$movement
  if (is.null(mv)) {
    .sg_stop("The model has no movement; give `movement` to sg_model() (see sg_movement_layer()).")
  }
  if (!is.null(temperatures)) .sg_check_temperature_input(temperatures, "temperatures")
  nd <- model$landscape$nodes
  ids <- nd$id
  n <- length(ids)
  p <- model$params
  t_ref <- mv$t_ref
  classes <- .sg_movement_classes
  life <- lapply(stats::setNames(classes, classes), function(cl) {
    .sg_adult_life(model, cl, t_ref, genotype)
  })
  life_days <- vapply(life, .sg_expected_life, 0)

  # distances between nodes (km) when the landscape has coordinates
  has_xy <- !is.null(nd$x) && !anyNA(nd$x)
  xy <- if (has_xy) cbind(nd$x, if (!is.null(nd$y) && !anyNA(nd$y)) nd$y else 0) else NULL
  dimension <- if (has_xy) .sg_landscape_dimension(xy) else NA_integer_
  dist <- if (has_xy) as.matrix(stats::dist(xy)) else NULL

  sh <- .sg_movement_shares(mv, n, function(g) .sg_clock_kappa(g$clock, t_ref, p, mv))
  rows <- lapply(seq_along(classes), function(k) {
    lam <- sh$lambda[k, ]
    path <- rep(NA_real_, n)
    d_impl <- rep(NA_real_, n)
    if (has_xy) {
      for (i in seq_len(n)) {
        pr <- sh$prob[[i]]
        if (length(pr$to) == 0L) {
          path[i] <- 0
          d_impl[i] <- 0
          next
        }
        dd <- dist[i, pr$to]
        path[i] <- sum(pr$p[k, ] * dd)
        d_impl[i] <- lam[i] * sum(pr$p[k, ] * dd^2) / (2 * dimension)
      }
    }
    data.frame(
      node = ids, role = nd$role, class = classes[k], rate = lam,
      daily_leave = -expm1(-lam),
      lifetime_leave = .sg_lifetime_leave(lam, life[[k]]),
      moves_per_life = lam * life_days[[k]],
      lost_share = ifelse(lam > 0, sh$lost[k, ], NA_real_),
      path_per_life = lam * life_days[[k]] * path,
      D = d_impl,
      stringsAsFactors = FALSE
    )
  })
  nodes <- do.call(rbind, rows)
  rownames(nodes) <- NULL

  temps <- NULL
  if (!is.null(temperatures)) {
    temps <- do.call(rbind, lapply(seq_along(classes), function(k) {
      i_max <- which.max(sh$lambda[k, ])
      do.call(rbind, lapply(temperatures, function(tc) {
        sht <- .sg_movement_shares(mv, n, function(g) .sg_clock_kappa(g$clock, tc, p, mv))
        lt <- .sg_adult_life(model, classes[k], tc, genotype)
        lam <- sht$lambda[k, i_max]
        el <- .sg_expected_life(lt)
        data.frame(
          temperature = tc, class = classes[k], node = ids[i_max],
          kappa = if (sh$lambda[k, i_max] > 0) lam / sh$lambda[k, i_max] else NA_real_,
          life = el, moves_per_life = lam * el,
          lifetime_leave = .sg_lifetime_leave(lam, lt), stringsAsFactors = FALSE
        )
      }))
    }))
    rownames(temps) <- NULL
  }

  warn <- character(0)
  high <- nodes[nodes$lifetime_leave > 0.99, , drop = FALSE]
  if (nrow(high) > 0L) {
    warn <- c(warn, paste0(
      "Lifetime leave probability above 0.99 (", paste(unique(high$class), collapse = ", "),
      " in node(s) ", .sg_wrap_names(unique(high$node)), "): flies almost surely leave their ",
      "node; check that the rates suit the grid (a rate for a finer grid used on a coarser one?)."
    ))
  }
  structure(
    list(
      nodes = nodes,
      life = data.frame(class = classes, life = unname(life_days), stringsAsFactors = FALSE),
      temperatures = temps, t_ref = t_ref, kappa_max = mv$kappa_max,
      clocks = vapply(mv$groups, `[[`, "", "clock"), dimension = dimension, warnings = warn
    ),
    class = "sg_movement_check"
  )
}

#' @export
print.sg_movement_check <- function(x, ...) {
  cat(sprintf(
    "<sg_movement_check> %d nodes; clock %s, kappa_max %s; rates at %g degC\n",
    length(unique(x$nodes$node)), paste(unique(x$clocks), collapse = " + "),
    format(x$kappa_max), x$t_ref
  ))
  cat(sprintf(
    "  adult life at %g degC: male %.3g d, female %.3g d\n", x$t_ref,
    x$life$life[x$life$class == "male"], x$life$life[x$life$class == "mated"]
  ))
  rng <- function(v) {
    v <- v[!is.na(v)]
    if (length(v) == 0L) {
      return("-")
    }
    r <- range(v)
    if (isTRUE(all.equal(r[1L], r[2L]))) {
      format(r[1L], digits = 3L)
    } else {
      paste0(format(r[1L], digits = 3L), "-", format(r[2L], digits = 3L))
    }
  }
  cat("  by class, range over nodes:\n")
  for (cl in unique(x$nodes$class)) {
    d <- x$nodes[x$nodes$class == cl, ]
    cat(sprintf(
      "    %-6s leave per day %s, in a lifetime %s; moves per life %s%s\n", cl,
      rng(d$daily_leave), rng(d$lifetime_leave), rng(d$moves_per_life),
      if (all(is.na(d$D))) {
        ""
      } else {
        sprintf("; path %s km; D %s km^2/d", rng(d$path_per_life), rng(d$D))
      }
    ))
  }
  if (!is.null(x$temperatures)) {
    cat("  moves per lifetime against temperature (node with the largest rate):\n")
    tt <- stats::reshape(x$temperatures[, c("temperature", "class", "moves_per_life")],
      idvar = "temperature", timevar = "class", direction = "wide"
    )
    names(tt) <- sub("^moves_per_life\\.", "", names(tt))
    print(tt, row.names = FALSE, digits = 3L)
  }
  if (length(x$warnings)) {
    cat("  warnings:\n")
    for (w in x$warnings) cat("   - ", w, "\n", sep = "")
  }
  invisible(x)
}

# ---- internal ------------------------------------------------------------------

# Kappa of a movement group's clock at a constant temperature (MODEL_SPEC N3).
.sg_clock_kappa <- function(clock, temp, p, mv) {
  if (clock == "calendar") {
    return(1)
  }
  min(mv$kappa_max, sg_degree_days(temp, p) / mv$dx_ref)
}

# Total leave rates and movement shares over the groups (N4) with kappa per
# group from `kappa_of(group)` (a number, or one value per node): lambda
# [3, n]; for each node its destinations (`to`, positions) with
# probabilities p [3, deg]; lost [3, n], the share of movers lost.
.sg_movement_shares <- function(mv, n, kappa_of) {
  lam_g <- lapply(mv$groups, function(g) sweep(g$rate, 2L, rep_len(kappa_of(g), n), `*`))
  lam <- Reduce(`+`, lam_g)
  prob <- vector("list", n)
  lost <- matrix(0, 3L, n)
  for (i in seq_len(n)) {
    edges <- lapply(mv$groups, function(g) {
      seq.int(g$row_ptr[i], length.out = g$row_ptr[i + 1L] - g$row_ptr[i])
    })
    to <- sort(unique(unlist(Map(function(g, e) g$col[e], mv$groups, edges))))
    pr <- matrix(0, 3L, length(to))
    for (k in seq_along(mv$groups)) {
      g <- mv$groups[[k]]
      for (cl in 1:3) {
        if (lam[cl, i] == 0) next
        w <- lam_g[[k]][cl, i] / lam[cl, i]
        for (e in edges[[k]]) {
          d <- match(g$col[e], to)
          pr[cl, d] <- pr[cl, d] + w * g$prob[cl, e]
        }
        lost[cl, i] <- lost[cl, i] + w * g$loss[cl, i]
      }
    }
    prob[[i]] <- list(to = to, p = pr)
  }
  list(lambda = lam, prob = prob, lost = lost)
}

# 1 when the coordinates lie on a line, else 2.
.sg_landscape_dimension <- function(xy) {
  if (nrow(xy) < 3L) {
    return(1L)
  }
  centred <- sweep(xy, 2L, colMeans(xy))
  sv <- svd(centred)$d
  if (length(sv) < 2L || sv[2L] <= 1e-9 * max(sv[1L], 1)) 1L else 2L
}
