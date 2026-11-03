# Scenario 17: sensitivity analysis of pgSIT control (BRT and figures)
#
# Reads the Latin hypercube that scenario 16 ran and fits boosted regression
# trees to two responses:
#
#   suppression  the lowest population reached: log10 of the smallest
#                number of adult wild females over the run, floored at
#                min_floor because below one female the deterministic
#                value is no longer meaningful (lower means better control);
#   cases        myiasis cases averted over the run (%).
#
# The relative influence of each parameter is then drawn with the
# biological parameters, which field data would have to pin down, separated
# from the engineering ones, which a release programme controls. A second
# figure checks the proxy the whole analysis rests on: that the lowest
# deterministic density predicts how often the stochastic model eliminates
# the wild type.
#
# Needs gbm, and uses dismo::gbm.step when it is installed (the method of
# Elith, Leathwick & Hastie 2008) and gbm::gbm with cross-validation
# otherwise. Install with install.packages("gbm").
#
# Run time: a minute or two, most of it fitting the trees.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Settings and data --------------------------------------------------------
sensitivity_dir <- file.path(getwd(), "scenario16_sensitivity") # from scenario 16
out_dir <- file.path(getwd(), "scenario17_figures") # figures are saved here
min_floor <- 1e-2 # lowest wild-female number the response distinguishes
learning_rate <- 0.005 # BRT shrinkage, lowered automatically if needed
tree_complexity <- 3 # interaction depth
bag_fraction <- 0.6
max_trees <- 6000
min_trees <- 1000 # fits below this are refitted more slowly (Elith et al. 2008)
seed <- 20260811
# UC Berkeley palette: Berkeley Blue for biology, Gold Dark for engineering
group_colours <- c(Biology = "#002676", Engineering = "#FC9313")

data <- readRDS(file.path(sensitivity_dir, "sensitivity_data.rds"))
parameters <- data$parameters
d <- data$sensitivity[data$sensitivity$ok, ]
d$suppression <- log10(pmax(d$min_wt, min_floor))
d$cases <- d$cases_averted
nrow(d) # samples that had a wild equilibrium

font <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  installed <- unique(systemfonts::system_fonts()$family)
  font <- c(intersect(c("Calibri", "Arial"), installed), "sans")[1]
}
save_figure <- function(plot, name, width, height) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(out_dir, paste0(name, ".png")), plot, width = width, height = height,
         units = "in", dpi = 300, device = png_device, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot, width = width, height = height,
         units = "in", device = grDevices::cairo_pdf)
}

# 2. Fit the boosted regression trees -----------------------------------------
# gbm.step chooses the number of trees by cross-validation itself; with plain
# gbm the number comes from gbm.perf() on the cross-validation error.
if (!requireNamespace("gbm", quietly = TRUE)) {
  stop("package 'gbm' is needed: install.packages(\"gbm\")")
}
# One fit at a given learning rate.
fit_once <- function(frame, rate) {
  set.seed(seed)
  if (requireNamespace("dismo", quietly = TRUE)) {
    fit <- dismo::gbm.step(
      data = frame, gbm.x = parameters$name, gbm.y = "y", family = "gaussian",
      tree.complexity = tree_complexity, learning.rate = rate,
      bag.fraction = bag_fraction, max.trees = max_trees, silent = TRUE,
      plot.main = FALSE
    )
    if (is.null(fit)) return(NULL)
    list(fit = fit, influence = summary(fit, plotit = FALSE), n_trees = fit$n.trees,
         explained = 1 - fit$cv.statistics$deviance.mean / fit$self.statistics$mean.null,
         rate = rate)
  } else {
    fit <- gbm::gbm(
      y ~ ., data = frame, distribution = "gaussian", n.trees = max_trees,
      interaction.depth = tree_complexity, shrinkage = rate,
      bag.fraction = bag_fraction, cv.folds = 5, n.minobsinnode = 10,
      n.cores = 1 # gbm splits the cross-validation over cores by default,
      # which would make the fit, and so the figure, differ from run to run
    )
    best <- which.min(fit$cv.error)
    null_dev <- mean((frame$y - mean(frame$y))^2)
    list(fit = fit, influence = summary(fit, n.trees = best, plotit = FALSE),
         n_trees = best, explained = 1 - fit$cv.error[best] / null_dev, rate = rate)
  }
}

# The learning rate is chosen so that the cross-validated optimum falls
# between min_trees and max_trees (Elith, Leathwick & Hastie 2008): a fit
# that stops below min_trees has learned too fast, and one that runs into
# max_trees never reached its optimum at all, so the rate is lowered or
# raised and the fit repeated.
fit_brt <- function(response) {
  frame <- d[c(parameters$name, response)]
  names(frame)[ncol(frame)] <- "y"
  rate <- learning_rate
  best <- NULL
  for (attempt in 1:5) {
    out <- fit_once(frame, rate)
    if (is.null(out)) {
      rate <- rate / 4 # the fit failed; try again more slowly
      next
    }
    if (out$n_trees >= min_trees && out$n_trees < max_trees) return(out) # converged
    # keep the closest fit so far, in case none of the attempts converges
    if (is.null(best) || abs(out$n_trees - min_trees) < abs(best$n_trees - min_trees)) {
      best <- out
    }
    rate <- if (out$n_trees < min_trees) rate / 4 else rate * 3
  }
  if (is.null(best)) stop("the boosted regression tree did not converge for ", response)
  warning(response, ": the optimum was not bracketed (", best$n_trees,
          " trees at rate ", signif(best$rate, 3), "); check the fit")
  best
}

responses <- c(suppression = "Lowest population reached",
               cases = "Myiasis cases averted (%)")
fits <- lapply(names(responses), fit_brt)
names(fits) <- names(responses)
vapply(fits, function(f) c(trees = f$n_trees, explained = f$explained), numeric(2))

# 3. Relative influence -------------------------------------------------------
influence <- do.call(rbind, lapply(names(fits), function(r) {
  inf <- fits[[r]]$influence
  data.frame(response = r, name = as.character(inf$var), influence = inf$rel.inf)
}))
influence <- merge(influence, parameters[c("name", "group", "label")], by = "name")
order_by <- stats::aggregate(influence ~ label, data = influence, FUN = mean)
influence$label <- factor(influence$label,
                          levels = order_by$label[order(order_by$influence)])
influence$response <- factor(responses[influence$response], levels = responses)
influence[order(-influence$influence), c("response", "label", "group", "influence")][1:10, ]

# One panel per response, sharing the parameter names on the left. The
# figure is wide and short (figure_width x figure_height below), so the
# type is sized to leave every one of the 15 rows legible.
figure_width <- 14.72
figure_height <- 7

fig_influence <- ggplot(influence, aes(influence, label, fill = group)) +
  geom_col(width = 0.74) +
  geom_text(aes(label = sprintf("%.0f", influence)), hjust = -0.28, size = 5,
            family = font, fontface = "bold", colour = "grey25") +
  facet_wrap(~response, nrow = 1) +
  scale_fill_manual(
    values = group_colours, name = NULL,
    labels = c(Biology = "Biology: must be measured",
               Engineering = "Engineering: can be improved")
  ) +
  scale_x_continuous(name = "Relative influence (%)",
                     expand = expansion(mult = c(0, 0.14))) +
  labs(y = NULL, title = "What decides whether pgSIT releases suppress the population",
       subtitle = sprintf(
         "Boosted regression trees over %d Latin hypercube samples of %d parameters",
         nrow(d), nrow(parameters)
       )) +
  theme_classic(base_size = 19, base_family = font) +
  theme(
    axis.title.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.text = element_text(colour = "black"),
    axis.text.y = element_text(size = 17, hjust = 1),
    axis.line = element_line(linewidth = 0.8),
    axis.ticks = element_line(linewidth = 0.8),
    axis.ticks.length = unit(0.18, "cm"),
    # faint vertical rules make the bar lengths easier to compare
    panel.grid.major.x = element_line(colour = "grey91", linewidth = 0.4),
    strip.background = element_rect(fill = "grey93", colour = NA),
    strip.text = element_text(face = "bold", size = 19, margin = margin(7, 6, 7, 6)),
    panel.spacing = unit(1.1, "cm"),
    legend.position = "top",
    legend.justification = "left",
    legend.text = element_text(size = 18),
    legend.margin = margin(0, 0, 4, 0),
    legend.key.size = unit(0.75, "cm"),
    plot.title.position = "plot", # align the title to the figure, not the panel
    plot.title = element_text(face = "bold", size = 24, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 17, colour = "grey30", margin = margin(b = 2)),
    plot.margin = margin(12, 16, 10, 12)
  )
fig_influence
save_figure(fig_influence, "sensitivity_influence",
            width = figure_width, height = figure_height)

# 4. The proxy check ----------------------------------------------------------
# Elimination in the stochastic model against the lowest density the
# deterministic model reached, with a logistic fit and McFadden's pseudo-R2.
v <- data$validation
v <- v[!is.na(v$p_elimination), ]
v$suppression <- log10(pmax(v$min_wt, min_floor))
n_rep <- data$settings$n_rep_validate
v$successes <- round(v$p_elimination * n_rep)
link <- stats::glm(cbind(successes, n_rep - successes) ~ suppression,
                   data = v, family = stats::binomial())
pseudo_r2 <- 1 - link$deviance / link$null.deviance
half <- (stats::qlogis(0.5) - stats::coef(link)[1]) / stats::coef(link)[2]
c(pseudo_r2 = unname(pseudo_r2), half_at = unname(half))

curve <- data.frame(suppression = seq(min(v$suppression), max(v$suppression),
                                      length.out = 200))
curve$p <- stats::predict(link, curve, type = "response")

fig_validation <- ggplot(v, aes(suppression, p_elimination)) +
  geom_line(data = curve, aes(suppression, p), colour = "grey30", linewidth = 1.2) +
  geom_point(size = 4, colour = group_colours[["Biology"]], alpha = 0.75) +
  annotate("text", x = min(v$suppression), y = 0.08, hjust = 0, size = 6.2, family = font,
           label = sprintf("pseudo-R² = %.2f", pseudo_r2)) +
  scale_y_continuous(labels = function(x) sprintf("%.0f%%", 100 * x)) +
  labs(
    x = expression(bold(paste("Lowest population reached, ", log[10], " adult wild females"))),
    y = "Replicates eliminated",
    title = "Deterministic minimum predicts elimination",
    subtitle = sprintf("%d parameter sets, %d replicates each", nrow(v), n_rep)
  ) +
  theme_classic(base_size = 19, base_family = font) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.line = element_line(linewidth = 0.8),
    plot.title.position = "plot",
    plot.title = element_text(face = "bold", size = 20, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 16, colour = "grey30", margin = margin(b = 10)),
    plot.margin = margin(12, 16, 10, 12)
  )
fig_validation
save_figure(fig_validation, "sensitivity_validation", width = 7.5, height = 6)

# 5. Numbers for the caption --------------------------------------------------
totals <- stats::aggregate(influence ~ response + group, data = influence, FUN = sum)
totals$influence <- round(totals$influence, 1)
totals
