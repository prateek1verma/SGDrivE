# Scenario 15: the release-effort figure (poster panel)
#
# Reads the thresholds that scenario 14 saved and draws one 14 x 8 inch
# panel in which each strategy has two bars:
#
#   left axis   how many released males per wild male it needs for a 90%
#               chance of elimination;
#   right axis  the myiasis cases averted over two years at that effort,
#               as a percentage of a no-release population.
#
# Both carry a 95% confidence interval: the delta-method interval of the
# fitted threshold, and the bootstrap interval of the median cases averted.
# Each whisker is drawn twice, a white line under a dark one, so that it
# stays visible against the dark bar as well as against the light one.
#
# The threshold over the plausible range of mating competitiveness eta is
# printed at the end for the caption; it is a different kind of uncertainty
# and would only crowd the figure.
#
# The panel has two y axes, so the axis titles and the values are coloured
# to match what they describe. Colours are from the UC Berkeley palette
# (brand.berkeley.edu/visual-identity/colors).
#
# Run time: a few seconds.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Settings and data --------------------------------------------------------
threshold_dir <- file.path(getwd(), "scenario14_threshold") # written by scenario 14
out_dir <- file.path(getwd(), "scenario15_figures") # figures are saved here
# UC Berkeley palette: Berkeley Blue and Gold Dark. Blue against gold stays
# distinguishable with any common form of colour blindness, because the two
# differ in lightness as well as in hue.
colour_threshold <- "#002676"
colour_averted <- "#FC9313"
# Short strategy names: the full labels do not fit a bar chart axis.
# "Sorted SIT" is the sex-sorted release of scenario 10.
short_labels <- c(classic = "Classic\nSIT", sorted = "Sorted\nSIT",
                  fsridl = "fsRIDL", pgsit = "pgSIT")

data <- readRDS(file.path(threshold_dir, "threshold_summary.rds"))
thresholds <- data$thresholds
outcomes <- data$outcomes
target_p <- data$settings$target_p
n_rep <- data$settings$n_rep_final

# The release schedule and the follow-up, read from the run so that the
# subtitle always describes the data in hand. Periods under a year are
# given in months.
period_text <- function(days) {
  if (days >= 365) {
    years <- days / 365
    sprintf("%s year%s", format(round(years, 1)), if (years >= 1.95) "s" else "")
  } else {
    sprintf("%s months", format(round(days / 30.44)))
  }
}
release_text <- period_text(data$settings$n_releases * data$settings$release_interval)
follow_text <- period_text(data$settings$days)

central <- thresholds[thresholds$role == "central", ]
panel <- merge(central[c("strategy", "eta", "threshold", "lower", "upper")],
               outcomes[c("strategy", "median_averted", "averted_lower",
                          "averted_upper", "p_elimination")], by = "strategy")
panel <- panel[order(panel$threshold), ] # cheapest strategy first
panel$x <- seq_len(nrow(panel))
keys <- as.character(panel$strategy)
panel$label <- unname(ifelse(keys %in% names(short_labels), short_labels[keys], keys))
panel[c("strategy", "eta", "threshold", "median_averted")]

# 2. Figure style -------------------------------------------------------------
font <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  installed <- unique(systemfonts::system_fonts()$family)
  font <- c(intersect(c("Calibri", "Arial"), installed), "sans")[1]
}

# 3. The two axes -------------------------------------------------------------
# Cases averted are drawn on the threshold scale and read off a secondary
# axis: averted_scale() maps them onto it, averted_unscale() back for the
# labels. The secondary axis covers the range of the data, padded, rather
# than 0 to 100, so that the differences between the strategies are visible.
y_max <- 1.16 * max(panel$upper)
averted_range <- range(c(panel$averted_lower, panel$averted_upper))
pad <- 0.45 * diff(averted_range)
averted_range <- c(averted_range[1] - pad, min(100, averted_range[2] + pad))
averted_scale <- function(a) (a - averted_range[1]) / diff(averted_range) * y_max
averted_unscale <- function(y) y / y_max * diff(averted_range) + averted_range[1]
bar <- 0.34 # width of each bar
offset <- 0.19 # the two bars sit either side of the strategy's position
cap <- 0.13 # width of the caps on the confidence intervals

# A confidence interval as a dark capped whisker on a white one, which keeps
# it legible whatever the bar underneath it is filled with.
interval <- function(mapping) {
  list(
    geom_errorbar(mapping, width = cap, linewidth = 2.6, colour = "white"),
    geom_errorbar(mapping, width = cap, linewidth = 1, colour = "grey10")
  )
}

fig <- ggplot(panel) +
  # release effort, left axis
  geom_col(aes(x = x - offset, y = threshold), fill = colour_threshold, width = bar) +
  interval(aes(x = x - offset, ymin = lower, ymax = upper)) +
  geom_text(aes(x = x - offset, y = upper, label = sprintf("%.1f", threshold)),
            vjust = -0.6, size = 8, family = font, fontface = "bold",
            colour = colour_threshold) +
  # cases averted, right axis
  geom_col(aes(x = x + offset, y = averted_scale(median_averted)),
           fill = colour_averted, width = bar) +
  interval(aes(x = x + offset, ymin = averted_scale(averted_lower),
               ymax = averted_scale(averted_upper))) +
  geom_text(aes(x = x + offset, y = averted_scale(averted_upper),
                label = sprintf("%.1f%%", median_averted)),
            vjust = -0.6, size = 8, family = font, fontface = "bold",
            colour = colour_averted) +
  scale_x_continuous(breaks = panel$x, labels = panel$label,
                     expand = expansion(add = 0.5)) +
  scale_y_continuous(
    name = "Released males per wild male",
    expand = expansion(mult = c(0, 0.06)),
    sec.axis = sec_axis(transform = averted_unscale, name = "Myiasis cases averted (%)")
  ) +
  coord_cartesian(ylim = c(0, y_max)) +
  labs(
    x = NULL,
    title = sprintf("Release effort for %.0f%% elimination, and the cases it averts",
                    100 * target_p),
    subtitle = sprintf(
      "Weekly releases for %s, population followed for %s  ·  %d replicates per strategy",
      release_text, follow_text, n_rep
    )
  ) +
  theme_classic(base_size = 27, base_family = font) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(colour = "black"),
    axis.text.x = element_text(face = "bold", margin = margin(t = 8)),
    axis.title.y.left = element_text(colour = colour_threshold),
    axis.text.y.left = element_text(colour = colour_threshold),
    axis.title.y.right = element_text(colour = colour_averted, angle = -90),
    axis.text.y.right = element_text(colour = colour_averted),
    axis.line = element_line(linewidth = 0.9),
    axis.ticks = element_line(linewidth = 0.9),
    axis.ticks.length = unit(0.22, "cm"),
    plot.title = element_text(face = "bold", size = 30, margin = margin(b = 6)),
    plot.subtitle = element_text(size = 21, colour = "grey30", margin = margin(b = 16)),
    plot.margin = margin(16, 22, 12, 16)
  )
fig

# 4. Save ---------------------------------------------------------------------
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
ggsave(file.path(out_dir, "release_threshold.png"), fig, width = 14, height = 8,
       units = "in", dpi = 300, device = png_device, bg = "white")
ggsave(file.path(out_dir, "release_threshold.pdf"), fig, width = 14, height = 8,
       units = "in", device = grDevices::cairo_pdf)

# 5. Numbers for the caption --------------------------------------------------
# The eta range is the threshold at the low and high plausible mating
# competitiveness of each strategy (scenario 14), left out of the figure.
span <- stats::aggregate(threshold ~ strategy, data = thresholds, FUN = range)
i <- match(panel$strategy, span$strategy)
data.frame(
  strategy = panel$strategy, eta = panel$eta,
  threshold = round(panel$threshold, 2),
  times_cheapest = round(panel$threshold / min(panel$threshold), 1),
  eta_range = sprintf("%.1f-%.1f", span$threshold[, 1][i], span$threshold[, 2][i]),
  averted = round(panel$median_averted, 1),
  elimination = panel$p_elimination, n_rep = n_rep
)
