# Plots of simulations (D-026): ggplot2, a suggested package. Stochastic
# replicates are drawn as separate lines; no summary bands.

#' Plot a simulation
#'
#' Time series of a deterministic or stochastic simulation as a ggplot2
#' plot: totals by stage, totals by genotype, allele frequencies, female
#' classes or cases. Stochastic replicates are drawn as separate thin lines.
#'
#' @param x A simulation from [simulate_ode()] or [simulate_cb()].
#' @param what What to plot: `"stages"` (totals by stage), `"genotypes"`
#'   (totals by genotype of the chosen `stages`), `"alleles"` (allele
#'   frequencies, [sg_allele_freq()]), `"females"` (female classes,
#'   [sg_females()]), `"cases"` (cases per period, [sg_cases()]) or `"hosts"`
#'   (hosts without a wound, wounded and infested, [sg_hosts()]; host-wound
#'   model only).
#' @param stages Stages included for `"stages"`, `"genotypes"` and
#'   `"alleles"` (see [sg_totals()]); default all stages for `"stages"`,
#'   adults otherwise.
#' @param per `"km2"` (default) or `"node"`.
#' @param period For `"cases"`: `"day"`, `"week"` (default), `"month"` or
#'   `"year"`.
#' @param host_density For `"cases"`: hosts per km^2; plots cases per host
#'   (with the host-wound model cases per host are plotted without it).
#' @param replicates For stochastic runs: which replicates to draw (default
#'   all).
#' @param palette Colour palette, all suitable for colour-blind readers:
#'   `"okabe_ito"` (default; Okabe & Ito 2008, up to 8 series), `"brewer"`
#'   (ColorBrewer Dark2, up to 8 series) or `"viridis"` (any number). With 9
#'   to 12 series both `"okabe_ito"` and `"brewer"` use ColorBrewer Paired,
#'   with more than 12 viridis.
#' @param nodes For models with a landscape: the nodes to sum over (default
#'   the study nodes); see [sg_totals()].
#' @param ... Not used.
#' @return A ggplot object (add layers, scales or themes with `+`). The
#'   plotted data are returned by [sg_plot_data()].
#' @seealso [sg_plot_data()]
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE) &&
#'   requireNamespace("ggplot2", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-01-01", end = "1971-12-31")
#'   mod <- sg_model(
#'     p, sg_cube(MGDrivE::cubeMendelian()), wx,
#'     sg_releases("aa", "adult", "male", times = seq(60, 200, 7), number = 30)
#'   )
#'   sim <- simulate_ode(mod)
#'   plot(sim)
#'   plot(sim, "alleles")
#'   plot(sim, "cases", period = "month")
#' }
plot.sg_sim <- function(x, what = c("stages", "genotypes", "alleles", "females", "cases", "hosts"),
                        stages = NULL, per = c("km2", "node"), period = "week",
                        host_density = NULL, replicates = NULL,
                        palette = c("okabe_ito", "brewer", "viridis"), nodes = NULL, ...) {
  .sg_plot(x, what, stages, per, period, host_density, replicates, palette, nodes)
}

#' @rdname plot.sg_sim
#' @export
plot.sg_sim_cb <- function(x,
                           what = c("stages", "genotypes", "alleles", "females", "cases", "hosts"),
                           stages = NULL, per = c("km2", "node"), period = "week",
                           host_density = NULL, replicates = NULL,
                           palette = c("okabe_ito", "brewer", "viridis"), nodes = NULL, ...) {
  .sg_plot(x, what, stages, per, period, host_density, replicates, palette, nodes)
}

#' Data behind a simulation plot
#'
#' The long table that [plot.sg_sim()] draws, for custom plots.
#'
#' @inheritParams plot.sg_sim
#' @param sim A simulation from [simulate_ode()] or [simulate_cb()].
#' @return A data frame with columns `date`, `replicate` (stochastic runs
#'   only), `series` (stage, genotype, allele or female class, with spaces
#'   for underscores; `"cases"` or `"cases per host"` for cases) and `value`,
#'   with attributes `ylab` and `title` used by the plot.
#' @export
#' @examples
#' if (requireNamespace("MGDrivE", quietly = TRUE)) {
#'   p <- sg_params()
#'   wx <- sg_forcing(sg_example_weather(), p, start = "1971-05-01", end = "1971-05-31")
#'   sim <- simulate_ode(sg_model(p, sg_cube(MGDrivE::cubeMendelian()), wx))
#'   head(sg_plot_data(sim, "genotypes"))
#' }
sg_plot_data <- function(sim,
                         what = c("stages", "genotypes", "alleles", "females", "cases", "hosts"),
                         stages = NULL, per = c("km2", "node"), period = "week",
                         host_density = NULL, replicates = NULL, nodes = NULL) {
  what <- .sg_match_choice(
    what, c("stages", "genotypes", "alleles", "females", "cases", "hosts"), "what"
  )
  per <- .sg_match_choice(per, c("km2", "node"), "per")
  .sg_sim_info(sim)
  unit <- if (per == "km2") "per km^2" else "per node"
  if (is.null(stages)) stages <- if (what == "stages") "all" else "adults"
  d <- switch(EXPR = what,
    stages = {
      t <- sg_totals(sim, by = "stage", stages = stages, per = per, nodes = nodes)
      out <- .sg_series(t, "stage", "value")
      lev <- intersect(.sg_stage_names, unique(out$series))
      out$series <- factor(gsub("_", " ", out$series), levels = gsub("_", " ", lev))
      structure(out, ylab = paste("individuals", unit), title = "Totals by stage")
    },
    genotypes = {
      t <- sg_totals(sim, by = "genotype", stages = stages, per = per, nodes = nodes)
      structure(.sg_series(t, "genotype", "value"),
        ylab = paste("individuals", unit),
        title = paste("Totals by genotype:", paste(stages, collapse = ", "))
      )
    },
    alleles = {
      t <- sg_allele_freq(sim, stages = stages, nodes = nodes)
      structure(.sg_series(t, "allele", "freq"),
        ylab = "allele frequency",
        title = paste("Allele frequencies:", paste(stages, collapse = ", "))
      )
    },
    females = {
      f <- sg_females(sim, per = per, nodes = nodes)
      classes <- c("females", "unmated", "mated", "wild_type", "egg_laying", "fertile_equivalent")
      classes <- classes[vapply(classes, function(cl) !all(is.na(f[[cl]])), TRUE)]
      ids <- intersect(c("date", "replicate"), names(f))
      out <- do.call(rbind, lapply(classes, function(cl) {
        cbind(f[ids], series = cl, value = f[[cl]])
      }))
      out$series <- factor(gsub("_", " ", out$series), levels = gsub("_", " ", classes))
      structure(out, ylab = paste("females", unit), title = "Female classes")
    },
    cases = {
      cs <- sg_cases(sim, period = period, per = per, host_density = host_density, nodes = nodes)
      col <- if (is.null(host_density)) "cases" else "cases_per_host"
      ids <- intersect(c("period", "replicate"), names(cs))
      out <- cbind(cs[ids], series = gsub("_", " ", col), value = cs[[col]])
      names(out)[1L] <- "date"
      structure(out,
        ylab = if (is.null(host_density)) {
          paste("cases", unit, "per", period)
        } else {
          paste("cases per host per", period)
        },
        title = if (is.null(.sg_sim_info(sim)$hosts)) {
          "Wounds newly infested with viable eggs"
        } else {
          "New infestations of hosts"
        }
      )
    },
    hosts = {
      hs <- sg_hosts(sim, per = per, nodes = nodes)
      classes <- c("unwounded", "wounded", "infested")
      ids <- intersect(c("date", "replicate"), names(hs))
      out <- do.call(rbind, lapply(classes, function(cl) {
        cbind(hs[ids], series = cl, value = hs[[cl]])
      }))
      out$series <- factor(out$series, levels = classes)
      structure(out, ylab = paste("hosts", unit), title = "Hosts")
    }
  )
  if (!is.null(replicates) && "replicate" %in% names(d)) {
    keep_attr <- attributes(d)[c("ylab", "title")]
    d <- d[d$replicate %in% replicates, , drop = FALSE]
    attributes(d)[c("ylab", "title")] <- keep_attr
  }
  rownames(d) <- NULL
  d
}

# ---- internal ------------------------------------------------------------------

.sg_has_ggplot2 <- function() requireNamespace("ggplot2", quietly = TRUE)

.sg_plot <- function(sim, what, stages, per, period, host_density, replicates,
                     palette = "okabe_ito", nodes = NULL) {
  palette <- .sg_match_choice(palette, c("okabe_ito", "brewer", "viridis"), "palette")
  if (!.sg_has_ggplot2()) {
    .sg_stop("Plots need the 'ggplot2' package: install.packages(\"ggplot2\").")
  }
  d <- sg_plot_data(sim, what, stages, per, period, host_density, replicates, nodes)
  stochastic <- "replicate" %in% names(d)
  if (stochastic) d$line <- interaction(d$series, d$replicate, drop = TRUE)
  aes_of <- function(...) do.call(ggplot2::aes, lapply(list(...), as.name))
  map <- if (stochastic) {
    aes_of(x = "date", y = "value", colour = "series", group = "line")
  } else {
    aes_of(x = "date", y = "value", colour = "series")
  }
  n_rep <- if (stochastic) length(unique(d$replicate)) else 1L
  one_series <- length(unique(d$series)) == 1L
  ggplot2::ggplot(d, map) +
    ggplot2::geom_line(alpha = if (n_rep > 1L) max(0.15, 1 / sqrt(n_rep)) else 1, linewidth = 0.5) +
    ggplot2::labs(
      x = NULL, y = attr(d, "ylab"), colour = NULL,
      title = attr(d, "title"),
      subtitle = if (stochastic) sprintf("%d replicate%s", n_rep, if (n_rep == 1L) "" else "s")
    ) +
    .sg_colour_scale(palette, length(unique(d$series))) +
    ggplot2::theme_minimal() +
    ggplot2::guides(
      colour = if (one_series) "none" else ggplot2::guide_legend(override.aes = list(alpha = 1))
    )
}

# Okabe & Ito (2008) colour-blind-safe palette, yellow last (low contrast on
# white).
.sg_okabe_ito <- c(
  "#000000", "#E69F00", "#56B4E9", "#009E73", "#0072B2", "#D55E00", "#CC79A7", "#F0E442"
)

# Discrete colour scale for n series from a colour-blind-friendly palette:
# Okabe-Ito or ColorBrewer Dark2 up to 8 series, ColorBrewer Paired for 9 to
# 12 (viridis tells many categories apart poorly), viridis beyond.
.sg_colour_scale <- function(palette, n) {
  if (palette != "viridis" && n <= 8L) {
    if (palette == "okabe_ito") {
      return(ggplot2::scale_colour_manual(values = unname(.sg_okabe_ito[seq_len(n)])))
    }
    return(ggplot2::scale_colour_brewer(palette = "Dark2"))
  }
  if (palette != "viridis" && n <= 12L) {
    return(ggplot2::scale_colour_brewer(palette = "Paired"))
  }
  ggplot2::scale_colour_viridis_d(end = 0.9)
}

# Rename a long table's grouping column to `series` and its value column to
# `value`, keeping date and replicate.
.sg_series <- function(t, col, value_col) {
  ids <- intersect(c("date", "replicate"), names(t))
  out <- t[c(ids, col, value_col)]
  names(out) <- c(ids, "series", "value")
  out
}
