# Scenario 13: four genetic-control strategies compared, with poster figures
#
# Stochastic population dynamics of a wild screwworm population at a
# constant 25 degC under four release strategies, each with the settings of
# its own scenario:
#
# 1. classic SIT (scenario 07): irradiated sterile males and females;
# 2. sex-sorted SIT (scenario 10): irradiated sterile males with 0.1% female
#    contamination;
# 3. fsRIDL (scenario 12): fertile males carrying a female-lethal transgene,
#    not irradiated;
# 4. pgSIT (scenario 11): sterile pgSIT males, not irradiated.
#
# All start from the same wild equilibrium (wild_adults adults, half of them
# female) and release every week for 26 weeks. Replicates run in chunks, in
# parallel, and each chunk is saved to out_dir/runs as it finishes, so memory
# stays small and an interrupted run resumes where it stopped. The script
# saves the trajectories (strategy_comparison_data.rds) and two figures in
# out_dir: one panel per strategy (28 x 8 inches, for a poster) and all
# strategies in one plot. Colours are from the UC Berkeley palette
# (brand.berkeley.edu/visual-identity/colors).
#
# Run time: about a minute for 100 replicates per strategy on a laptop.

library(ScrewDrivE)
library(ggplot2) # plots

# 1. Settings -----------------------------------------------------------------
temperature <- 25 # constant temperature (degC)
days <- 730 # simulated days (two years)
wild_adults <- 20000 # adults at equilibrium (mgdrive-nws AdPopEQ)
n_rep <- 100 # stochastic replicates per strategy
chunk_size <- 100 # replicates per chunk saved to disk (memory per chunk ~30 MB)
n_workers <- max(1, parallel::detectCores() - 1) # parallel workers (1 = sequential)
seed <- 20260811 # random seed (as in mgdrive-nws)
n_releases <- 26 # weekly releases (half a year)
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)
out_dir <- file.path(getwd(), "scenario13_figures") # figures and data are saved here

# Adult field mortality (per day) beyond the temperature mortality, as in
# scenarios 07, 10, 11 and 12
p <- sg_params(female_background_mortality = 0.1, male_background_mortality = 0.12)

# Release strategies. eta: mating competitiveness of released males;
# male_ratio: released males per wild adult per release (before the
# sequestration correction); female_percent: released females, % of the male
# release; mating_capacity: females a released male can mate with (Inf: no
# correction); lifespan: fraction of the wild adult lifespan kept by the
# released genotype.
strategies <- list(
  classic = list(
    label = "Classic SIT (M + F, Rad)", cube = "sit", eta = 0.4, male_ratio = 1,
    female_percent = 100, mating_capacity = 2, lifespan = 1
  ),
  sorted = list(
    label = "Sex-sorted SIT (0.1% F, Rad)", cube = "sit", eta = 0.5,
    male_ratio = 1, female_percent = 0.1, mating_capacity = Inf, lifespan = 1
  ),
  fsridl = list(
    label = "fsRIDL (Male-only, No Rad)", cube = "fsridl", eta = 0.9, male_ratio = 1,
    female_percent = 0, mating_capacity = Inf, lifespan = 1, female_lethality = 1
  ),
  pgsit = list(
    label = "pgSIT (Male-only, No Rad)", cube = "pgsit", eta = 0.95, male_ratio = 1,
    female_percent = 0, mating_capacity = Inf, lifespan = 1, penetrance = 1
  )
)
# UC Berkeley palette: Gold Dark, Green Medium, Rose Medium, Berkeley Blue
colours <- c(classic = "#FC9313", sorted = "#018943", fsridl = "#E7115E", pgsit = "#002676")

# 2. Inheritance cubes --------------------------------------------------------
# Two genotypes for SIT and pgSIT (W wild, S released), three for fsRIDL (WW,
# WH, HH); see scenarios 07, 11 and 12 for the reasoning.
two_genotype_cube <- function(eta_s, lifespan, s_s, offspring_survival) {
  g <- c("W", "S")
  dims <- list(g, g, g)
  ih <- array(0, c(2, 2, 2), dimnames = dims)
  ih[, , "W"] <- 1 # offspring are recorded as wild type
  tau <- array(1, c(2, 2, 2), dimnames = dims)
  tau[, "S", ] <- offspring_survival # offspring sired by released males
  eta <- matrix(1, 2, 2, dimnames = list(g, g))
  eta[, "S"] <- eta_s
  sg_cube(list(
    genotypesID = g, ih = ih, tau = tau, eta = eta, phi = c(W = 0.5, S = 0.5),
    xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), s = c(W = 1, S = s_s),
    omega = c(W = 1, S = 1), wildType = "W", releaseType = "S"
  ), lifespan_fraction = c(W = 1, S = lifespan))
}

fsridl_cube <- function(eta_h, lifespan, female_lethality) {
  g <- c("WW", "WH", "HH")
  dims <- list(g, g, g)
  gametes <- rbind(WW = c(W = 1, H = 0), WH = c(W = 0.5, H = 0.5), HH = c(W = 0, H = 1))
  ih <- array(0, c(3, 3, 3), dimnames = dims)
  for (f in g) {
    for (m in g) {
      a <- gametes[f, ]
      b <- gametes[m, ]
      ih[f, m, ] <- c(a[["W"]] * b[["W"]], a[["W"]] * b[["H"]] + a[["H"]] * b[["W"]],
                      a[["H"]] * b[["H"]])
    }
  }
  carriers <- c("WH", "HH")
  tau <- array(1, c(3, 3, 3), dimnames = dims)
  tau[, , carriers] <- 1 - 0.5 * female_lethality # carrier daughters die as embryos
  phi <- c(WW = 0.5, WH = 0, HH = 0)
  phi[carriers] <- 0.5 * (1 - female_lethality) / (1 - 0.5 * female_lethality)
  eta <- matrix(1, 3, 3, dimnames = list(g, g))
  eta[, carriers] <- eta_h
  ones <- c(WW = 1, WH = 1, HH = 1)
  sg_cube(list(
    genotypesID = g, ih = ih, tau = tau, eta = eta, phi = phi, xiF = ones,
    xiM = ones, s = ones, omega = ones, wildType = "WW", releaseType = "HH"
  ), lifespan_fraction = c(WW = 1, WH = lifespan, HH = lifespan))
}

make_cube <- function(st) {
  switch(st$cube,
    sit = two_genotype_cube(st$eta, st$lifespan, s_s = 0, offspring_survival = 0),
    pgsit = two_genotype_cube(st$eta, st$lifespan, s_s = 1, offspring_survival = 1 - st$penetrance),
    fsridl = fsridl_cube(st$eta, st$lifespan, st$female_lethality)
  )
}

# 3. Simulate every strategy --------------------------------------------------
# For each strategy: the wild population sized at equilibrium, the releases
# (with the sequestration correction where mating_capacity is finite), then
# an ODE run and n_rep chain-binomial replicates from the wild equilibrium.
# Replicates keep only totals (keep = "totals", enough for wild-type females)
# and run in chunks of chunk_size, in parallel; each chunk has its own seed
# (seed + chunk), as replicate streams derive from the seed. Chunks already
# in out_dir/runs are reused: delete that folder after changing settings.
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
chunk_dir <- file.path(out_dir, "runs")
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)

setup_strategy <- function(key) {
  st <- strategies[[key]]
  cube <- make_cube(st)
  wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
  female_ratio <- st$male_ratio * st$female_percent / 100
  male_ratio <- max(0, st$male_ratio - female_ratio / st$mating_capacity)
  rel <- sg_releases(cube$release_type, stage = "adult", sex = "male",
                     times = release_days, number = male_ratio * wild_adults)
  if (female_ratio > 0) {
    rel <- c(rel, sg_releases(cube$release_type, stage = "adult", sex = "female",
                              times = release_days, number = female_ratio * wild_adults))
  }
  model <- sg_model(p, cube, forcing, releases = rel, area = wild$area)
  list(model = model, init = sg_init(model, "equilibrium"))
}

parallel_run <- n_workers > 1 && requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)
if (parallel_run) future::plan("multisession", workers = n_workers)

ode_list <- list()
n_chunks <- ceiling(n_rep / chunk_size)
for (key in names(strategies)) {
  s <- setup_strategy(key)
  o <- sg_females(simulate_ode(s$model, s$init))
  ode_list[[key]] <- data.frame(strategy = key, day = o$day, females = o$wild_type)
  for (k in seq_len(n_chunks)) {
    file <- file.path(chunk_dir, sprintf("%s_chunk%03d.rds", key, k))
    if (file.exists(file)) next
    n_k <- min(chunk_size, n_rep - (k - 1) * chunk_size)
    sim <- simulate_cb(s$model, s$init, n_rep = n_k, seed = seed + k,
                       keep = "totals", parallel = parallel_run)
    f <- sg_females(sim)
    saveRDS(data.frame(strategy = key, replicate = f$replicate + (k - 1) * chunk_size,
                       day = f$day, females = f$wild_type), file)
    rm(sim, f)
    gc()
    message(key, ": chunk ", k, " of ", n_chunks, " saved")
  }
}
if (parallel_run) future::plan("sequential")

# Trajectories of all chunks of this run (n_chunks per strategy)
chunk_files <- file.path(chunk_dir, sprintf("%s_chunk%03d.rds",
                                            rep(names(strategies), each = n_chunks),
                                            rep(seq_len(n_chunks), length(strategies))))
traj <- do.call(rbind, lapply(chunk_files, readRDS))
ode <- do.call(rbind, ode_list)

# Save the data needed to redraw the figures (reload with readRDS())
saveRDS(
  list(
    traj = traj, ode = ode, strategies = strategies, colours = colours,
    settings = list(
      temperature = temperature, days = days, wild_adults = wild_adults, n_rep = n_rep,
      chunk_size = chunk_size, seed = seed, n_releases = n_releases,
      release_interval = release_interval, first_release = first_release, params = p
    )
  ),
  file.path(out_dir, "strategy_comparison_data.rds")
)

labels <- vapply(strategies, `[[`, "", "label")
traj$strategy <- factor(traj$strategy, levels = names(strategies))
ode$strategy <- factor(ode$strategy, levels = names(strategies))

# Replicates with no wild-type adult female left at the end
final <- traj[traj$day == max(traj$day), ]
eliminated <- tapply(final$females < 1, final$strategy, sum)
eliminated

# 4. Figure style -------------------------------------------------------------
# Calibri or Arial when installed (checked with systemfonts), else the
# default sans-serif font. PNGs are written with ragg and PDFs with cairo so
# that system fonts are embedded.
font <- "sans"
if (requireNamespace("systemfonts", quietly = TRUE)) {
  installed <- unique(systemfonts::system_fonts()$family)
  font <- c(intersect(c("Calibri", "Arial"), installed), "sans")[1]
}
release_times <- release_days - 1 # a release on day d is added at time d - 1
line_alpha <- min(0.35, 15 / n_rep) # lighter lines when there are many replicates
y_max <- 1.25 * wild_adults / 2
k_labels <- function(x) ifelse(x == 0, "0", paste0(x / 1000, "k"))
save_figure <- function(plot, name, width, height) {
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  png_device <- if (requireNamespace("ragg", quietly = TRUE)) ragg::agg_png else "png"
  ggsave(file.path(out_dir, paste0(name, ".png")), plot,
         width = width, height = height, units = "in", dpi = 300, device = png_device, bg = "white")
  ggsave(file.path(out_dir, paste0(name, ".pdf")), plot,
         width = width, height = height, units = "in", device = grDevices::cairo_pdf)
  file.path(out_dir, paste0(name, c(".png", ".pdf")))
}
poster_theme <- function(base_size) {
  theme_classic(base_size = base_size, base_family = font) +
    theme(
      axis.title = element_text(face = "bold"),
      axis.text = element_text(colour = "black"),
      axis.line = element_line(linewidth = 0.9),
      axis.ticks = element_line(linewidth = 0.9),
      axis.ticks.length = unit(0.25, "cm"),
      plot.margin = margin(20, 30, 15, 15)
    )
}

# 5. Figure 1: one panel per strategy (28 x 8 inches) --------------------------
panel_labels <- paste0(seq_along(labels), ". ", labels)
names(panel_labels) <- names(labels)
counts <- data.frame(
  strategy = factor(names(eliminated), levels = names(strategies)),
  text = paste0("Eliminated: ", eliminated, "/", n_rep)
)
fig_panels <- ggplot(traj, aes(day, females, group = replicate, colour = strategy)) +
  geom_vline(xintercept = release_times, linetype = "dashed", colour = "grey60",
             linewidth = 0.35) +
  geom_hline(yintercept = wild_adults / 2, linetype = "dashed", colour = "grey40",
             linewidth = 0.6) +
  geom_line(alpha = line_alpha, linewidth = 0.5) +
  geom_text(data = counts, aes(x = days, y = 0.96 * y_max, label = text),
            inherit.aes = FALSE, hjust = 1, size = 9, family = font) +
  facet_wrap(~strategy, nrow = 1, labeller = as_labeller(panel_labels)) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_continuous(breaks = c(0, 180, 365, 545, 730), expand = expansion(mult = 0.02)) +
  scale_y_continuous(labels = k_labels, breaks = seq(0, y_max, by = 2500),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, days), ylim = c(0, y_max)) +
  labs(x = "Time (days)", y = "Adult wild-type females") +
  poster_theme(base_size = 30) +
  theme(
    strip.background = element_rect(fill = "grey92", colour = "black", linewidth = 0.9),
    strip.text = element_text(face = "bold", size = 27, margin = margin(10, 6, 10, 6)),
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.9),
    panel.spacing = unit(1.2, "cm")
  )
fig_panels
save_figure(fig_panels, "strategy_panels", width = 28, height = 8)

# 6. Figure 2: all strategies in one plot -------------------------------------
# The legend sits in the lower right, which the trajectories leave empty once
# the populations have recovered (or died out at zero).
legend_labels <- paste0(labels, "  (eliminated ", eliminated, "/", n_rep, ")")
fig_all <- ggplot(traj, aes(day, females, group = interaction(strategy, replicate),
                            colour = strategy)) +
  geom_vline(xintercept = release_times, linetype = "dashed", colour = "grey60",
             linewidth = 0.35) +
  geom_hline(yintercept = wild_adults / 2, linetype = "dashed", colour = "grey40",
             linewidth = 0.6) +
  geom_line(alpha = line_alpha, linewidth = 0.5) +
  scale_colour_manual(values = colours, labels = legend_labels, name = "Release strategy") +
  scale_x_continuous(breaks = c(0, 56, 112, 182, 300, 365, 500, 600, 730), # 182: end of releases
                     expand = expansion(mult = 0.01)) +
  scale_y_continuous(labels = k_labels, breaks = seq(0, y_max, by = 2500),
                     expand = expansion(mult = c(0, 0.02))) +
  coord_cartesian(xlim = c(0, days), ylim = c(0, y_max)) +
  labs(x = "Time (days)", y = "Adult wild-type females") +
  guides(colour = guide_legend(byrow = TRUE, override.aes = list(alpha = 1, linewidth = 5))) +
  poster_theme(base_size = 40) +
  theme(
    panel.border = element_rect(fill = NA, colour = "black", linewidth = 0.9),
    axis.line = element_blank(),
    legend.position = "inside",
    legend.position.inside = c(0.985, 0.035),
    legend.justification = c(1, 0),
    legend.title = element_text(face = "bold", size = 30, margin = margin(b = 14)),
    legend.text = element_text(size = 26),
    legend.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    legend.margin = margin(16, 20, 16, 16),
    legend.key.width = unit(1.8, "cm"),
    legend.key.height = unit(0.9, "cm"),
    legend.key.spacing.y = unit(0.5, "cm") # white space between legend entries
  )
fig_all
save_figure(fig_all, "strategy_all", width = 18, height = 10)
