#' Landscape of several nodes
#'
#' Defines the nodes of a multi-node model (MODEL_SPEC section 8a, D-032):
#' each node has an area, its own weather (forcing), optional coordinates,
#' optional wound or host densities and a role. Pass the landscape to
#' [sg_model()] with `landscape = `; adult flies move between nodes as given
#' by a movement layer ([sg_movement_layer()]).
#'
#' @param nodes A data frame with one row per node and columns:
#'   \describe{
#'     \item{`id`}{Node names: unique, non-empty, without `"@"` (used in
#'       state names such as `"M[AA,1]@n3"`).}
#'     \item{`area`}{Node area in km^2 (> 0).}
#'     \item{`x`, `y`}{Optional coordinates in km (used by later diagnostics
#'       and constructors; not by the model equations).}
#'     \item{`role`}{Optional: `"study"` (default) or `"source"`. Source
#'       nodes are simulated like study nodes and feed them through movement
#'       (e.g. an overwintering area outside the study region); summaries
#'       report study nodes by default.}
#'     \item{`wound_density`, `host_density`}{Optional per-node densities per
#'       km^2; `NA` (or a missing column) means the value in [sg_params()].
#'       `wound_density` is used with `wound_model = "static"`,
#'       `host_density` with `wound_model = "hosts"`.}
#'   }
#' @param forcing One `sg_forcing` object ([sg_forcing()]) for all nodes, or a
#'   list with one per node (named by node id, or in node order). All must
#'   cover the same dates and be built with the same parameters.
#'
#' @details Identical forcing objects are stored once, with an index per node,
#'   so nodes that share a weather station share one forcing series.
#'
#' @return An object of class `sg_landscape`: a list with `nodes` (the node
#'   table with all columns filled in and a column `forcing_series`),
#'   `forcing` (the distinct forcing series) and `n_nodes`.
#' @seealso [sg_movement_layer()], [sg_model()]
#' @export
#' @examples
#' p <- sg_params()
#' warm <- sg_forcing(sg_constant_weather(27, days = 30), p)
#' cool <- sg_forcing(sg_constant_weather(20, days = 30), p)
#' land <- sg_landscape(
#'   data.frame(id = c("south", "middle", "north"), area = c(100, 100, 150)),
#'   forcing = list(south = warm, middle = warm, north = cool)
#' )
#' land
#' as.data.frame(land)
sg_landscape <- function(nodes, forcing) {
  nodes <- .sg_check_nodes(nodes)
  n <- nrow(nodes)
  series <- .sg_forcing_series(.sg_node_forcings(forcing, nodes$id))
  nodes$forcing_series <- series$index
  structure(
    list(nodes = nodes, forcing = series$forcing, n_nodes = n),
    class = "sg_landscape"
  )
}

#' @export
print.sg_landscape <- function(x, ...) {
  nd <- x$nodes
  f <- x$forcing[[1L]]
  n_src <- sum(nd$role == "source")
  cat(sprintf(
    "<sg_landscape> %d node%s (%d study, %d source), total area %s km^2\n",
    x$n_nodes, if (x$n_nodes != 1L) "s" else "", x$n_nodes - n_src, n_src,
    format(sum(nd$area), big.mark = ",")
  ))
  cat(sprintf(
    "  forcing: %d distinct series, %d days (%s to %s)\n",
    length(x$forcing), f$n_days, format(f$dates[1L]), format(f$dates[f$n_days])
  ))
  show <- utils::head(nd, 10L)
  print(show, row.names = FALSE)
  if (x$n_nodes > 10L) cat("  ... ", x$n_nodes - 10L, " more nodes\n", sep = "")
  invisible(x)
}

#' @export
as.data.frame.sg_landscape <- function(x, ...) {
  x$nodes
}

# ---- internal ------------------------------------------------------------------

# Node table with every column filled in: id, area, x, y, role,
# wound_density, host_density (NA = the value in the parameters).
.sg_check_nodes <- function(nodes) {
  if (!is.data.frame(nodes)) {
    .sg_stop("`nodes` must be a data frame with columns `id` and `area` (one row per node).")
  }
  if (nrow(nodes) < 1L) .sg_stop("`nodes` must have at least one row.")
  missing <- setdiff(c("id", "area"), names(nodes))
  if (length(missing) > 0L) {
    .sg_stop("`nodes` lacks column(s) ", paste0("`", missing, "`", collapse = ", "), ".")
  }
  known <- c("id", "area", "x", "y", "role", "wound_density", "host_density")
  extra <- setdiff(names(nodes), known)
  if (length(extra) > 0L) {
    .sg_stop(
      "Unknown column(s) in `nodes`: ", paste0("`", extra, "`", collapse = ", "),
      ". Allowed: ", paste0("`", known, "`", collapse = ", "), ".",
      .sg_suggest(extra[1L], known)
    )
  }
  n <- nrow(nodes)
  id <- as.character(nodes$id)
  if (anyNA(id) || any(!nzchar(id))) .sg_stop("Node ids in `nodes$id` must be non-empty.")
  dup <- unique(id[duplicated(id)])
  if (length(dup) > 0L) {
    .sg_stop("Node ids must be unique; repeated: ", paste(dup, collapse = ", "), ".")
  }
  if (any(grepl("@", id, fixed = TRUE))) {
    .sg_stop("Node ids must not contain \"@\" (it separates the node in state names).")
  }
  num_col <- function(nm, default, lower, strict, allow_na) {
    v <- if (nm %in% names(nodes)) nodes[[nm]] else rep(default, n)
    if (!is.numeric(v) && !all(is.na(v))) .sg_stop("`nodes$", nm, "` must be numeric.")
    v <- as.numeric(v)
    bad <- if (allow_na) !is.na(v) & !is.finite(v) else !is.finite(v)
    if (!is.na(lower)) {
      below <- if (strict) v <= lower else v < lower
      bad <- bad | (!is.na(v) & below)
    }
    if (any(bad)) {
      .sg_stop(
        "`nodes$", nm, "` must be ", if (allow_na) "NA or ",
        if (is.na(lower)) "finite" else if (strict) "> 0" else ">= 0",
        "; node ", id[which(bad)[1L]], " has ", .sg_fmt(v[which(bad)[1L]]), "."
      )
    }
    v
  }
  area <- num_col("area", NA, 0, TRUE, FALSE)
  x <- num_col("x", NA_real_, NA, FALSE, TRUE)
  y <- num_col("y", NA_real_, NA, FALSE, TRUE)
  wd <- num_col("wound_density", NA_real_, 0, FALSE, TRUE)
  hd <- num_col("host_density", NA_real_, 0, FALSE, TRUE)
  role <- if ("role" %in% names(nodes)) as.character(nodes$role) else rep("study", n)
  bad_role <- is.na(role) | !role %in% c("study", "source")
  if (any(bad_role)) {
    .sg_stop(
      "`nodes$role` must be \"study\" or \"source\"; node ", id[which(bad_role)[1L]],
      " has ", .sg_fmt(role[which(bad_role)[1L]]), "."
    )
  }
  if (!any(role == "study")) .sg_stop("A landscape needs at least one study node.")
  data.frame(
    id = id, area = area, x = x, y = y, role = role,
    wound_density = wd, host_density = hd, stringsAsFactors = FALSE
  )
}

# One sg_forcing per node, in node order, from a single forcing or a list.
.sg_node_forcings <- function(forcing, ids) {
  if (inherits(forcing, "sg_forcing")) {
    return(rep(list(forcing), length(ids)))
  }
  if (!is.list(forcing) || !all(vapply(forcing, inherits, TRUE, "sg_forcing"))) {
    .sg_stop(
      "`forcing` must be one sg_forcing object or a list of them (one per node); ",
      "create them with sg_forcing()."
    )
  }
  if (length(forcing) != length(ids)) {
    .sg_stop(
      "`forcing` has ", length(forcing), " elements but the landscape has ",
      length(ids), " nodes; give one forcing per node or one for all."
    )
  }
  nm <- names(forcing)
  if (!is.null(nm) && any(nzchar(nm))) {
    unknown <- setdiff(nm, ids)
    if (length(unknown) > 0L || anyDuplicated(nm)) {
      .sg_stop(
        "Names of `forcing` must be the node ids, each once",
        if (length(unknown) > 0L) paste0("; unknown: ", paste(unknown, collapse = ", ")) else "",
        "."
      )
    }
    forcing <- forcing[ids]
  }
  unname(forcing)
}

# Distinct forcing series and the series index (1-based) of each node. All
# forcings must share dates and forcing parameters.
.sg_forcing_series <- function(forcings) {
  first <- forcings[[1L]]
  distinct <- list()
  index <- integer(length(forcings))
  for (i in seq_along(forcings)) {
    f <- forcings[[i]]
    if (!identical(f$dates, first$dates)) {
      .sg_stop(
        "All forcings of a landscape must cover the same dates; node ", i,
        " runs from ", format(f$dates[1L]), " to ", format(f$dates[f$n_days]),
        ", node 1 from ", format(first$dates[1L]), " to ", format(first$dates[first$n_days]), "."
      )
    }
    if (!identical(f$param_values, first$param_values)) {
      .sg_stop(
        "All forcings of a landscape must be built with the same parameters; node ", i,
        " differs from node 1. Rebuild them with sg_forcing(weather, params)."
      )
    }
    hit <- 0L
    for (k in seq_along(distinct)) {
      if (identical(unclass(distinct[[k]]), unclass(f))) {
        hit <- k
        break
      }
    }
    if (hit == 0L) {
      distinct[[length(distinct) + 1L]] <- f
      hit <- length(distinct)
    }
    index[i] <- hit
  }
  list(forcing = distinct, index = index)
}
