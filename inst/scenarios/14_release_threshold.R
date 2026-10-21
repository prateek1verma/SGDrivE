# Scenario 14: the release effort each strategy needs (simulation)
#
# How many sterile or transgenic males per wild male must be released every
# week for six months to eliminate the wild population with probability
# target_p? For the four strategies of scenario 13 (classic SIT, sex-sorted
# SIT, fsRIDL, pgSIT) this script finds that threshold and, at it, the time
# to elimination and the myiasis cases averted.
#
# Unlike a full grid, the mating competitiveness eta is not swept: it is a
# property of the technology, so each strategy keeps its own value (step 1)
# and the threshold is found again at a low and a high plausible value, which
# gives the range drawn on the figure.
#
# Steps for each strategy and eta:
#
#   1. a coarse geometric ladder of release ratios brackets the threshold;
#   2. a fine ladder inside the bracket is added, a logistic regression of
#      elimination on log(release ratio) is fitted to the points near it, and
#      the threshold is the ratio at which the fitted probability is target_p
#      (95% confidence interval by the delta method);
# For fsRIDL the eta range is narrow, because the sons that carry the
# transgene inherit the released males' competitiveness, so raising it helps
# the released males and their carrier sons alike.
#
#   3. at the threshold of the central eta, n_rep_final replicates give the
#      median days to elimination and the cases averted, both with
#      bootstrap confidence intervals, against a no-release baseline.
#
# Release ratios are per wild MALE at equilibrium (the overflooding ratio of
# the SIT literature), not per wild adult as in scenarios 07 and 10 to 13.
#
# Each simulated cell is cached in out_dir/cells, so memory stays small and
# an interrupted run resumes. Delete that folder after changing a setting.
#
# Run time: about three minutes with the defaults below.

library(SGDrivE)

# 1. Settings -----------------------------------------------------------------
temperature <- 25 # constant temperature (degC)
days <- 365*5 # simulated days (two years)
wild_adults <- 20000 # adults at equilibrium (mgdrive-nws AdPopEQ)
n_releases <- 52*2 # weekly releases (half a year)
release_interval <- 7 # days between releases
first_release <- 1 # day of the first release (day 1 = first simulated day)
target_p <- 0.9 # elimination probability the threshold is defined by
seed <- 20260811 # random seed (as in mgdrive-nws)
n_rep_coarse <- 20 # replicates per cell, bracketing ladder
n_rep_fine <- 50 # replicates per cell, ladder near the threshold
n_rep_final <- 150 # replicates at the threshold, and for the baseline
n_boot <- 2000 # bootstrap samples for the medians
n_workers <- max(1, parallel::detectCores() - 2) # parallel workers (1 = sequential)
out_dir <- file.path(getwd(), "scenario14_threshold") # results are saved here

# Released males per wild male: the ladder that brackets the threshold
coarse_ladder <- 0.5 * 2^(0:8) # 0.5, 1, 2, ... 128
n_fine <- 7 # points of the ladder inside the bracket

# Strategies as in scenario 13. eta: mating competitiveness of released males
# relative to wild males, as central, low and high plausible values;
# irradiated males compete poorly, pgSIT and fsRIDL males are not irradiated.
# female_percent: released females, % of the male release; mating_capacity:
# females a released male can mate with (Inf: no sequestration correction);
# lifespan: fraction of the wild adult lifespan kept by the released genotype.
strategies <- list(
  classic = list(
    label = "Classic SIT (M + F, Rad)", cube = "sit", eta = c(0.4),
    female_percent = 100, mating_capacity = 2, lifespan = 1
  ),
  sorted = list(
    label = "Sex-sorted SIT (0.1% F, Rad)", cube = "sit", eta = c(0.5),
    female_percent = 0.1, mating_capacity = Inf, lifespan = 1
  ),
  fsridl = list(
    label = "fsRIDL (Male-only, No Rad)", cube = "fsridl", eta = c(0.90),
    female_percent = 0, mating_capacity = Inf, lifespan = 1, female_lethality = 1
  ),
  pgsit = list(
    label = "pgSIT (Male-only, No Rad)", cube = "pgsit", eta = c(0.95),
    female_percent = 0, mating_capacity = Inf, lifespan = 1, penetrance = 1
  )
)

# Adult field mortality (per day) beyond the temperature mortality, as in
# scenarios 07 and 10 to 13
p <- sg_params(female_background_mortality = 0.1, male_background_mortality = 0.12)

# 2. Inheritance cubes (as in scenario 13) ------------------------------------
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

make_cube <- function(st, eta) {
  switch(st$cube,
    sit = two_genotype_cube(eta, st$lifespan, s_s = 0, offspring_survival = 0),
    pgsit = two_genotype_cube(eta, st$lifespan, s_s = 1,
                              offspring_survival = 1 - st$penetrance),
    fsridl = fsridl_cube(eta, st$lifespan, st$female_lethality)
  )
}

# 3. One simulated cell -------------------------------------------------------
# The node area comes from the wild equilibrium, which does not involve the
# released genotype, so it is the same for every strategy and eta. Wild males
# are half of the adults at equilibrium, and the release ratio is per wild
# male. Elimination means no individual of the wild-type genotype is left in
# any life stage (simulate_cb()'s wild_extinct_day).
forcing <- sg_forcing(sg_constant_weather(temperature, days = days), p)
release_days <- seq(first_release, by = release_interval, length.out = n_releases)
wild_males <- wild_adults / 2
area <- sg_model(p, make_cube(strategies$pgsit, 1), forcing,
                 adult_females = wild_adults / 2)$area

cell_dir <- file.path(out_dir, "cells")
dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)

strategy_model <- function(key, eta, ratio) {
  st <- strategies[[key]]
  cube <- make_cube(st, eta)
  female_ratio <- ratio * st$female_percent / 100
  male_ratio <- max(0, ratio - female_ratio / st$mating_capacity)
  rel <- NULL
  if (male_ratio > 0) {
    rel <- sg_releases(cube$release_type, stage = "adult", sex = "male",
                       times = release_days, number = male_ratio * wild_males)
  }
  if (female_ratio > 0) {
    rel <- c(rel, sg_releases(cube$release_type, stage = "adult", sex = "female",
                              times = release_days, number = female_ratio * wild_males))
  }
  sg_model(p, cube, forcing, releases = rel, area = area)
}

# Cases per replicate: myiasis cases (oviposition events with a viable egg)
# over the whole run, in the node.
replicate_cases <- function(sim) {
  cases <- sg_cases(sim, period = "year", per = "node")
  as.numeric(tapply(cases$cases, cases$replicate, sum))
}

run_cell <- function(key, eta, ratio, n_rep, cell_seed, parallel) {
  file <- file.path(cell_dir, sprintf("%s_e%03.0f_r%09.4f_n%04d.rds",
                                      key, 100 * eta, ratio, n_rep))
  if (file.exists(file)) return(readRDS(file))
  model <- strategy_model(key, eta, ratio)
  sim <- simulate_cb(model, sg_init(model, "equilibrium"), n_rep = n_rep,
                     seed = cell_seed, keep = "totals", parallel = parallel)
  out <- list(
    strategy = key, eta = eta, ratio = ratio, n_rep = n_rep,
    extinct_day = as.numeric(sim$wild_extinct_day), # NA where the wild type survives
    cases = replicate_cases(sim)
  )
  saveRDS(out, file)
  message(sprintf("  %s, eta %.2f, %.3f males per wild male: eliminated %d of %d",
                  key, eta, ratio, sum(!is.na(out$extinct_day)), n_rep))
  out
}

# 4. Threshold for one strategy and eta ---------------------------------------
# The threshold is where the fitted logistic curve crosses target_p. Its
# confidence interval comes from the delta method on log(ratio); the fit uses
# the ladder points around the bracket, where the curve is informative.
fit_threshold <- function(cells, bracket) {
  d <- data.frame(
    ratio = vapply(cells, `[[`, 0, "ratio"),
    n = vapply(cells, `[[`, 0, "n_rep"),
    k = vapply(cells, function(x) sum(!is.na(x$extinct_day)), 0)
  )
  d <- d[d$ratio >= bracket[1] / 3 & d$ratio <= bracket[2] * 3, ]
  d <- stats::aggregate(cbind(k, n) ~ ratio, data = d, FUN = sum)
  fit <- try(stats::glm(cbind(k, n - k) ~ log(ratio), data = d,
                        family = stats::binomial()), silent = TRUE)
  if (inherits(fit, "try-error") || stats::coef(fit)[2] <= 0) {
    return(list(threshold = NA_real_, lower = NA_real_, upper = NA_real_, fit = NULL))
  }
  b <- unname(stats::coef(fit))
  target <- stats::qlogis(target_p)
  log_threshold <- (target - b[1]) / b[2]
  grad <- c(-1 / b[2], -(target - b[1]) / b[2]^2)
  se <- sqrt(drop(t(grad) %*% stats::vcov(fit) %*% grad))
  list(
    threshold = exp(log_threshold), lower = exp(log_threshold - 1.96 * se),
    upper = exp(log_threshold + 1.96 * se), fit = fit
  )
}

find_threshold <- function(key, eta, parallel) {
  cells <- list()
  take <- function(ratio, n_rep, tag) {
    # own random streams per strategy, eta, ratio and ladder
    cell_seed <- seed + sum(utf8ToInt(paste0(key, tag))) + round(1000 * ratio) +
      round(1e4 * eta)
    run_cell(key, eta, ratio, n_rep, cell_seed, parallel)
  }
  # coarse ladder: stop as soon as a ratio reaches the target
  p_hat <- rep(NA_real_, length(coarse_ladder))
  for (i in seq_along(coarse_ladder)) {
    cell <- take(coarse_ladder[i], n_rep_coarse, "coarse")
    cells[[length(cells) + 1L]] <- cell
    p_hat[i] <- mean(!is.na(cell$extinct_day))
    if (p_hat[i] >= target_p) break
  }
  hit <- which(p_hat >= target_p)
  if (!length(hit)) {
    return(list(strategy = key, eta = eta, threshold = NA_real_, lower = NA_real_,
                upper = NA_real_, above_max = TRUE, cells = cells))
  }
  i <- hit[1]
  bracket <- c(if (i == 1) coarse_ladder[1] / 2 else coarse_ladder[i - 1], coarse_ladder[i])
  # fine ladder inside the bracket, then the logistic fit
  fine <- exp(seq(log(bracket[1]), log(bracket[2]), length.out = n_fine))
  for (r in fine) cells[[length(cells) + 1L]] <- take(r, n_rep_fine, "fine")
  fit <- fit_threshold(cells, bracket)
  if (is.na(fit$threshold)) fit$threshold <- bracket[2] # fallback: the first ratio that passed
  c(list(strategy = key, eta = eta, above_max = FALSE, cells = cells), fit[1:3])
}

# 5. Run: thresholds at three eta, then replicates at the central threshold ---
parallel_run <- n_workers > 1 && requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)
if (parallel_run) future::plan("multisession", workers = n_workers)

boot_ci <- function(x, probs = c(0.025, 0.975)) {
  x <- x[is.finite(x)]
  if (length(x) < 2) return(c(NA_real_, NA_real_))
  set.seed(seed)
  bs <- replicate(n_boot, stats::median(sample(x, replace = TRUE)))
  unname(stats::quantile(bs, probs))
}

# No-release baseline: the cases a wild population causes in two years
baseline <- run_cell("pgsit", 1, 0, n_rep_final, seed + 1, parallel_run)
baseline_cases <- mean(baseline$cases)
baseline_cases

thresholds <- NULL
outcomes <- NULL
for (key in names(strategies)) {
  message("strategy: ", key)
  etas <- strategies[[key]]$eta # central, low, high
  found <- lapply(etas, function(e) find_threshold(key, e, parallel_run))
  thresholds <- rbind(thresholds, data.frame(
    strategy = key, role = c("central", "low", "high"), eta = etas,
    threshold = vapply(found, function(x) x$threshold, 0),
    lower = vapply(found, function(x) x$lower, 0),
    upper = vapply(found, function(x) x$upper, 0),
    above_max = vapply(found, function(x) isTRUE(x$above_max), TRUE)
  ))
  central <- found[[1]]
  if (is.na(central$threshold)) next
  # replicates at the threshold: time to elimination and cases averted
  at <- run_cell(key, etas[1], central$threshold, n_rep_final, seed + 2, parallel_run)
  averted <- 100 * (1 - at$cases / baseline_cases)
  days_ci <- boot_ci(at$extinct_day)
  averted_ci <- boot_ci(averted)
  outcomes <- rbind(outcomes, data.frame(
    strategy = key, eta = etas[1], threshold = central$threshold,
    n_rep = n_rep_final, p_elimination = mean(!is.na(at$extinct_day)),
    median_days = stats::median(at$extinct_day, na.rm = TRUE),
    days_lower = days_ci[1], days_upper = days_ci[2],
    median_averted = stats::median(averted),
    averted_lower = averted_ci[1], averted_upper = averted_ci[2],
    median_cases = stats::median(at$cases), baseline_cases = baseline_cases
  ))
}
if (parallel_run) future::plan("sequential")

# 6. Save ---------------------------------------------------------------------
thresholds$strategy <- factor(thresholds$strategy, levels = names(strategies))
outcomes$strategy <- factor(outcomes$strategy, levels = names(strategies))
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
saveRDS(
  list(
    thresholds = thresholds, outcomes = outcomes,
    labels = vapply(strategies, `[[`, "", "label"), strategies = strategies,
    settings = list(
      temperature = temperature, days = days, wild_adults = wild_adults,
      wild_males = wild_males, area = area, target_p = target_p, seed = seed,
      n_releases = n_releases, release_interval = release_interval,
      first_release = first_release, n_rep_coarse = n_rep_coarse,
      n_rep_fine = n_rep_fine, n_rep_final = n_rep_final,
      baseline_cases = baseline_cases, params = p
    )
  ),
  file.path(out_dir, "threshold_summary.rds")
)
write.csv(thresholds, file.path(out_dir, "thresholds.csv"), row.names = FALSE)
write.csv(outcomes, file.path(out_dir, "outcomes.csv"), row.names = FALSE)

thresholds
outcomes
