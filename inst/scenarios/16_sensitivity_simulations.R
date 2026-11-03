# Scenario 16: sensitivity analysis of pgSIT control (simulations)
#
# Which parameters decide whether a pgSIT release campaign suppresses a
# screwworm population: the biology we would have to measure, or the
# engineering we can improve? This script draws a Latin hypercube over both
# groups (step 1), runs the deterministic model once with releases and once
# without for every sample, and records how far the wild population fell and
# how many myiasis cases were averted. Scenario 17 fits the boosted
# regression trees and draws the figure.
#
# Deterministic runs are used because they are fast (about 0.2 s a sample)
# and because the lowest wild-female density they reach is what governs
# stochastic extinction. To show that, step 5 re-runs a subsample with the
# chain-binomial solver and records how often the wild type was eliminated,
# so that scenario 17 can relate the two.
#
# Two things to know about the design:
#
# * every sample is re-sized with `adult_females`, so each starts from the
#   same 10,000 adult wild females and the responses measure how hard the
#   population is to control, not how large it is;
# * releases are set per wild male. Without an Allee effect the model is
#   then scale-free, and `wound_search_rate` only sets the area; it matters
#   here because `male_density_half_mating` (mate finding, a density) is
#   swept as well.
#
# Samples are cached in out_dir/chunks, so an interrupted run resumes.
# Delete that folder after changing a setting.
#
# Run time: about ten minutes with the defaults below, less in parallel.

library(ScrewDrivE)

# 1. Settings and the sampled parameters --------------------------------------
n_samples <- 5000 # Latin hypercube samples
chunk_size <- 500 # samples per cached chunk
n_validate <- 120 # samples re-run with the stochastic solver (step 5)
n_rep_validate <- 50 # stochastic replicates per validation sample
days <- 1825 # simulated days (five years)
wild_adults <- 20000 # adults at equilibrium, so 10,000 wild females
seed <- 20260811
n_workers <- max(1, parallel::detectCores() - 2) # parallel workers (1 = sequential)
out_dir <- file.path(getwd(), "scenario16_sensitivity") # results are saved here

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
  data.frame(name = "release_ratio", low = 0.5, high = 20, log = TRUE,
             group = "Engineering", label = "Released males per wild male"),
  data.frame(name = "eta", low = 0.2, high = 1, log = FALSE,
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
summary(design[parameters$name[1:4]])

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

# Everything one sample needs: the wild population it starts from and the
# same model with the pgSIT releases added.
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

total_cases <- function(sim) sum(sg_cases(sim, period = "year", per = "node")$cases)

# Deterministic responses. A sample is dropped when the wild population has
# no positive equilibrium (R0 <= 1), which the ranges allow.
run_sample <- function(row) {
  out <- data.frame(
    sample = row$sample, ok = FALSE, r0 = NA_real_, area = NA_real_,
    n_releases = NA_integer_, start_wt = NA_real_, min_wt = NA_real_,
    final_wt = NA_real_, cases = NA_real_, baseline_cases = NA_real_,
    cases_averted = NA_real_
  )
  res <- try(
    {
      m <- sample_models(row)
      base_sim <- simulate_ode(m$wild, sg_init(m$wild, "equilibrium"))
      sim <- simulate_ode(m$treated, sg_init(m$treated, "equilibrium"))
      females <- sg_females(sim)
      list(m = m, base_sim = base_sim, sim = sim, females = females)
    },
    silent = TRUE
  )
  if (inherits(res, "try-error")) return(out)
  out$ok <- TRUE
  out$r0 <- sg_r0(res$m$wild, row$temperature)
  out$area <- res$m$wild$area
  out$n_releases <- res$m$n_releases
  out$start_wt <- res$females$wild_type[1]
  out$min_wt <- min(res$females$wild_type)
  out$final_wt <- res$females$wild_type[nrow(res$females)]
  out$cases <- total_cases(res$sim)
  out$baseline_cases <- total_cases(res$base_sim)
  out$cases_averted <- 100 * (1 - out$cases / out$baseline_cases)
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
for (k in seq_along(chunks)) {
  file <- file.path(chunk_dir, sprintf("deterministic_%03d.rds", k))
  if (file.exists(file)) next
  rows <- lapply(chunks[[k]], function(i) design[i, ])
  saveRDS(do.call(rbind, apply_over(rows, run_sample)), file)
  message("deterministic chunk ", k, " of ", length(chunks))
}
results <- do.call(rbind, lapply(
  file.path(chunk_dir, sprintf("deterministic_%03d.rds", seq_along(chunks))), readRDS
))
sensitivity <- merge(design, results, by = "sample")
sum(!sensitivity$ok) # samples without a positive wild equilibrium

# 5. Stochastic check of the deterministic proxy ------------------------------
# The subsample is spread over the range of the deterministic minimum, so
# that scenario 17 can fit elimination probability against it across the
# whole range rather than only where most samples happen to fall.
ok <- sensitivity[sensitivity$ok, ]
set.seed(seed + 1)
strata <- cut(log10(pmax(ok$min_wt, 1e-6)), breaks = 12, labels = FALSE)
per_stratum <- ceiling(n_validate / length(unique(strata)))
picked <- unlist(lapply(split(ok$sample, strata), function(s) {
  s[sample.int(length(s), min(per_stratum, length(s)))]
}), use.names = FALSE)
spare <- setdiff(ok$sample, picked) # top up if the thin strata ran short
if (length(picked) < n_validate && length(spare)) {
  picked <- c(picked, spare[sample.int(length(spare),
                                       min(n_validate - length(picked), length(spare)))])
}
picked <- sort(unique(picked))[seq_len(min(n_validate, length(unique(picked))))]

run_stochastic <- function(row) {
  res <- try(
    {
      m <- sample_models(row)
      sim <- simulate_cb(m$treated, sg_init(m$treated, "equilibrium"),
                         n_rep = n_rep_validate, seed = seed + row$sample,
                         keep = "totals")
      mean(!is.na(as.numeric(sim$wild_extinct_day)))
    },
    silent = TRUE
  )
  data.frame(sample = row$sample,
             p_elimination = if (inherits(res, "try-error")) NA_real_ else res)
}

validation_file <- file.path(chunk_dir, "stochastic.rds")
if (!file.exists(validation_file)) {
  rows <- lapply(picked, function(i) design[design$sample == i, ])
  saveRDS(do.call(rbind, apply_over(rows, run_stochastic)), validation_file)
}
validation <- merge(readRDS(validation_file),
                    sensitivity[c("sample", "min_wt", "cases_averted")], by = "sample")
if (parallel_run) future::plan("sequential")

# 6. Save ---------------------------------------------------------------------
saveRDS(
  list(
    sensitivity = sensitivity, validation = validation, parameters = parameters,
    settings = list(
      n_samples = n_samples, n_validate = n_validate, n_rep_validate = n_rep_validate,
      days = days, wild_adults = wild_adults, seed = seed, strategy = "pgSIT",
      wound_model = "static"
    )
  ),
  file.path(out_dir, "sensitivity_data.rds")
)
write.csv(sensitivity, file.path(out_dir, "sensitivity_samples.csv"), row.names = FALSE)

summary(sensitivity$cases_averted)
summary(log10(pmax(sensitivity$min_wt, 1e-2)))
