# Scenario 18: stochastic sensitivity analysis of pgSIT control (simulations)
#
# The same question as scenario 16 - which parameters decide whether a pgSIT
# release campaign works - but every sample of the Latin hypercube is run
# with the chain-binomial solver instead of the deterministic one, so that
# the probability of eliminating the wild type is measured directly rather
# than inferred from how low the deterministic population fell.
#
# For each of n_samples parameter sets:
#
#   * n_rep stochastic replicates with the releases, giving the share of
#     replicates in which no wild-type individual was left (the elimination
#     probability), the days it took, the lowest number of adult wild
#     females reached and the myiasis cases;
#   * one deterministic run without releases as the baseline for cases
#     averted. The baseline has no releases and a large population, so its
#     stochastic mean is well approximated by the ODE and running it
#     stochastically would double the cost for nothing.
#
# Scenario 19 fits the boosted regression trees and draws the figures.
# The design (15 parameters, their ranges, the re-sizing of every sample to
# the same 10,000 adult wild females, releases per wild male) is that of
# scenario 16; see its comments for why wound search matters only because
# the mate-finding scale is swept as well.
#
# Samples are cached in out_dir/chunks, so an interrupted run resumes.
# Delete that folder after changing a setting.
#
# Run time: about twenty minutes with the defaults below on eight workers,
# a few hours if run sequentially.

library(SGDrivE)

# 1. Settings and the sampled parameters --------------------------------------
n_samples <- 5000 # Latin hypercube samples
n_rep <- 50 # stochastic replicates per sample
chunk_size <- 250 # samples per cached chunk
days <- 1825 # simulated days (five years)
wild_adults <- 20000 # adults at equilibrium, so 10,000 wild females
seed <- 20260811
n_workers <- max(1, parallel::detectCores() - 2) # parallel workers (1 = sequential)
out_dir <- file.path(getwd(), "scenario18_stochastic") # results are saved here

# Each row: the parameter, its range, whether it is sampled on the log scale,
# and the group it belongs to in the figure. The ranges are what the analysis
# rests on: a parameter given a wide range will look influential, so each one
# should be defensible and is quoted in the figure caption.
parameters <- rbind(
  # biology: what field data would have to pin down
  data.frame(name = "temperature", low = 18, high = 35, log = FALSE,
             group = "Biology", label = "Temperature"),
  data.frame(name = "wound_search_rate", low = 5e-4, high = 0.03, log = TRUE,
             group = "Biology", label = "Wound search rate"),
  data.frame(name = "male_density_half_mating", low = 0, high = 10, log = FALSE,
             group = "Biology", label = "Mate-finding (Allee) scale"),
  data.frame(name = "female_background_mortality", low = 0.02, high = 0.25, log = FALSE,
             group = "Biology", label = "Female field mortality"),
  data.frame(name = "male_background_mortality", low = 0.02, high = 0.25, log = FALSE,
             group = "Biology", label = "Male field mortality"),
  data.frame(name = "larva_survival_daily", low = 0.85, high = 0.98, log = FALSE,
             group = "Biology", label = "Larval survival"),
  data.frame(name = "eggs_per_day_max", low = 60, high = 120, log = FALSE,
             group = "Biology", label = "Eggs per female per day"),
  data.frame(name = "larval_crowding_coef", low = 0, high = 0.08, log = FALSE,
             group = "Biology", label = "Larval crowding"),
  data.frame(name = "oviposition_age_dd", low = 45, high = 75, log = FALSE,
             group = "Biology", label = "Pre-oviposition period"),
  # engineering: what a release programme controls
  data.frame(name = "release_ratio", low = 0.5, high = 50, log = TRUE,
             group = "Engineering", label = "Released males per wild male"),
  data.frame(name = "eta", low = 0.3, high = 1, log = FALSE,
             group = "Engineering", label = "Mating competitiveness"),
  # sampled as the leak, 1 - penetrance: elimination needs the leak below
  # 1/R0, so a linear scale over 0.9 to 1 would put most samples in a region
  # where no release size can work
  data.frame(name = "leakage", low = 1e-4, high = 0.1, log = TRUE,
             group = "Engineering", label = "Sterility leak (1 - penetrance)"),
  data.frame(name = "lifespan_released", low = 0.4, high = 1, log = FALSE,
             group = "Engineering", label = "Released-male lifespan"),
  data.frame(name = "release_interval", low = 7, high = 21, log = FALSE,
             group = "Engineering", label = "Release interval"),
  data.frame(name = "campaign_days", low = 182, high = 730, log = FALSE,
             group = "Engineering", label = "Campaign length")
)

# 2. The Latin hypercube ------------------------------------------------------
# One stratified draw per parameter, in base R: each column is a random
# permutation of n strata, jittered inside its stratum.
latin_hypercube <- function(n, k) {
  m <- matrix(0, n, k)
  for (j in seq_len(k)) m[, j] <- (sample.int(n) - stats::runif(n)) / n
  m
}
set.seed(seed)
unit <- latin_hypercube(n_samples, nrow(parameters))
design <- as.data.frame(mapply(function(u, low, high, log_scale) {
  if (log_scale) exp(log(low) + u * (log(high) - log(low))) else low + u * (high - low)
}, split(unit, col(unit)), parameters$low, parameters$high, parameters$log))
names(design) <- parameters$name
design$release_interval <- round(design$release_interval)
design$sample <- seq_len(n_samples)

# 3. One sample ---------------------------------------------------------------
# The pgSIT cube of scenario 11: W wild, S the released sterile male; the
# offspring of S survive with probability 1 - penetrance.
pgsit_cube <- function(eta_s, lifespan, penetrance) {
  g <- c("W", "S")
  dims <- list(g, g, g)
  ih <- array(0, c(2, 2, 2), dimnames = dims)
  ih[, , "W"] <- 1
  tau <- array(1, c(2, 2, 2), dimnames = dims)
  tau[, "S", ] <- 1 - penetrance
  eta <- matrix(1, 2, 2, dimnames = list(g, g))
  eta[, "S"] <- eta_s
  sg_cube(list(
    genotypesID = g, ih = ih, tau = tau, eta = eta, phi = c(W = 0.5, S = 0.5),
    xiF = c(W = 1, S = 1), xiM = c(W = 1, S = 1), s = c(W = 1, S = 1),
    omega = c(W = 1, S = 1), wildType = "W", releaseType = "S"
  ), lifespan_fraction = c(W = 1, S = lifespan))
}

sample_params <- function(row) {
  sg_params(
    wound_search_rate = row$wound_search_rate,
    male_density_half_mating = row$male_density_half_mating,
    female_background_mortality = row$female_background_mortality,
    male_background_mortality = row$male_background_mortality,
    larva_survival_daily = row$larva_survival_daily,
    eggs_per_day_max = row$eggs_per_day_max,
    oviposition_age_dd = row$oviposition_age_dd,
    larval_crowding_mode = "per_wound", # a coefficient of 0 is "none"
    larval_crowding_coef = row$larval_crowding_coef
  )
}

sample_models <- function(row) {
  p <- sample_params(row)
  forcing <- sg_forcing(sg_constant_weather(row$temperature, days = days), p)
  cube <- pgsit_cube(row$eta, row$lifespan_released, 1 - row$leakage)
  wild <- sg_model(p, cube, forcing, adult_females = wild_adults / 2)
  times <- seq(1, by = row$release_interval,
               length.out = max(1, floor(row$campaign_days / row$release_interval)))
  rel <- sg_releases("S", stage = "adult", sex = "male", times = times,
                     number = row$release_ratio * wild_adults / 2)
  list(wild = wild, treated = sg_model(p, cube, forcing, releases = rel, area = wild$area),
       n_releases = length(times))
}

# Myiasis cases per replicate over the whole run, in the node.
replicate_cases <- function(sim) {
  cases <- sg_cases(sim, period = "year", per = "node")
  as.numeric(tapply(cases$cases, cases$replicate, sum))
}

# Stochastic responses. A sample is dropped when the wild population has no
# positive equilibrium (R0 <= 1), which the ranges allow. simulate_cb() is
# left sequential here because the samples themselves run in parallel.
run_sample <- function(row) {
  out <- data.frame(
    sample = row$sample, ok = FALSE, r0 = NA_real_, area = NA_real_,
    n_releases = NA_integer_, n_rep = n_rep, n_eliminated = NA_integer_,
    p_elimination = NA_real_, median_days = NA_real_, median_min_wt = NA_real_,
    median_cases = NA_real_, baseline_cases = NA_real_, cases_averted = NA_real_
  )
  res <- try(
    {
      m <- sample_models(row)
      base_sim <- simulate_ode(m$wild, sg_init(m$wild, "equilibrium"))
      sim <- simulate_cb(m$treated, sg_init(m$treated, "equilibrium"), n_rep = n_rep,
                         seed = seed + row$sample, keep = "totals")
      females <- sg_females(sim)
      list(m = m, base_sim = base_sim, sim = sim, females = females)
    },
    silent = TRUE
  )
  if (inherits(res, "try-error")) return(out)
  extinct_day <- as.numeric(res$sim$wild_extinct_day) # NA where the wild type survives
  minima <- as.numeric(tapply(res$females$wild_type, res$females$replicate, min))
  cases <- replicate_cases(res$sim)
  baseline <- sum(sg_cases(res$base_sim, period = "year", per = "node")$cases)
  out$ok <- TRUE
  out$r0 <- sg_r0(res$m$wild, row$temperature)
  out$area <- res$m$wild$area
  out$n_releases <- res$m$n_releases
  out$n_eliminated <- sum(!is.na(extinct_day))
  out$p_elimination <- mean(!is.na(extinct_day))
  out$median_days <- if (all(is.na(extinct_day))) {
    NA_real_
  } else {
    stats::median(extinct_day, na.rm = TRUE)
  }
  out$median_min_wt <- stats::median(minima)
  out$median_cases <- stats::median(cases)
  out$baseline_cases <- baseline
  out$cases_averted <- 100 * (1 - stats::median(cases) / baseline)
  out
}

# 4. Run the hypercube, one cached chunk at a time ----------------------------
chunk_dir <- file.path(out_dir, "chunks")
dir.create(chunk_dir, showWarnings = FALSE, recursive = TRUE)
parallel_run <- n_workers > 1 && requireNamespace("future", quietly = TRUE) &&
  requireNamespace("future.apply", quietly = TRUE)
if (parallel_run) future::plan("multisession", workers = n_workers)
apply_over <- function(x, f) {
  if (parallel_run) future.apply::future_lapply(x, f, future.seed = TRUE) else lapply(x, f)
}

chunks <- split(seq_len(n_samples), ceiling(seq_len(n_samples) / chunk_size))
started <- Sys.time()
for (k in seq_along(chunks)) {
  file <- file.path(chunk_dir, sprintf("stochastic_%03d.rds", k))
  if (file.exists(file)) next
  rows <- lapply(chunks[[k]], function(i) design[i, ])
  saveRDS(do.call(rbind, apply_over(rows, run_sample)), file)
  elapsed <- as.numeric(difftime(Sys.time(), started, units = "mins"))
  message(sprintf("chunk %d of %d (%.1f min elapsed)", k, length(chunks), elapsed))
}
if (parallel_run) future::plan("sequential")

results <- do.call(rbind, lapply(
  file.path(chunk_dir, sprintf("stochastic_%03d.rds", seq_along(chunks))), readRDS
))
sensitivity <- merge(design, results, by = "sample")
sum(!sensitivity$ok) # samples without a positive wild equilibrium

# 5. Save ---------------------------------------------------------------------
saveRDS(
  list(
    sensitivity = sensitivity, parameters = parameters,
    settings = list(
      n_samples = n_samples, n_rep = n_rep, days = days, wild_adults = wild_adults,
      seed = seed, strategy = "pgSIT", wound_model = "static", solver = "chain-binomial"
    )
  ),
  file.path(out_dir, "stochastic_sensitivity_data.rds")
)
write.csv(sensitivity, file.path(out_dir, "stochastic_sensitivity_samples.csv"),
          row.names = FALSE)

summary(sensitivity$p_elimination)
summary(sensitivity$cases_averted)
