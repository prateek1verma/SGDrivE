# Texas 1962, weeks 1-26: model fits with an explicit observation model
#
# Observation model taken from the case-inversion exercise
# (MGDrivE1-screwworm/inference/case_inversion/nws_cases.stan), kept here in the
# script; the SGDrivE package is not changed:
#
#   cases_w ~ NegBin(mu_w, phi),  mu_w = rho * sum over week w of (g * batches)
#
# where `batches` are the daily egg-laying events simulated by SGDrivE
# (counter "egg_batches"), g is a discretised gamma delay from egg laying to
# report (mean 6 d, sd 3 d, up to 20 d), and rho is reported cases per egg batch
# (under-reporting, reportable hosts, egg masses per wound).
#
# Parameters are fitted by maximum likelihood (no MCMC): wound density, rho,
# phi, and depending on the scenario the weekly immigration rate (adults of both
# sexes from 1 February) or a development-rate scale. Development rates are at
# their literature value unless stated.
#
# Run from the repository root:  Rscript data-raw/texas_nws/analysis/obs_model_fits_1962.R
# Run time: about 15 minutes.

suppressMessages(devtools::load_all(".", quiet = TRUE))
suppressMessages(library(ggplot2))

data_dir <- "data-raw/texas_nws"
out_dir <- file.path(data_dir, "analysis", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Settings. Cases identify only rho x wound density (as in the case-inversion
# model), so rho needs a prior. Beta(2, 198) has mean 0.01 (about 1 reported case
# per 100 egg batches, the order of the 1976 Texas confirmed-case / infested-animal
# ratio) and a 95% range of about 0.001-0.03. NULL fits rho by likelihood alone.
# The immigration bound scales with rho: at rho near 0.01 about 100 times more
# flies are needed than at rho = 1, so the cap is 100,000 adults a week.
rho_prior <- c(2, 198)
imm_upper <- 1e5
run_tag <- if (is.null(rho_prior)) "" else "_rho_prior"

# 1. Observation model (from the case-inversion exercise) ---------------------
delay_kernel <- function(mean = 6, sd = 3, max_days = 20) {
  g <- diff(pgamma(0:max_days, shape = (mean / sd)^2, rate = mean / sd^2))
  g / sum(g) # g[j] = P(delay = j - 1 days)
}
g_delay <- delay_kernel()

case_mean <- function(batches, rho, n_weeks = 26) {
  lam <- vapply(seq_along(batches), function(t) {
    j <- seq_len(min(length(g_delay), t))
    sum(g_delay[j] * batches[t - j + 1])
  }, numeric(1))
  rho * vapply(seq_len(n_weeks), function(w) sum(lam[7 * (w - 1) + 1:7]), numeric(1))
}

nb_loglik <- function(cases, mu, phi) {
  sum(dnbinom(cases, mu = mu + 1e-9, size = phi, log = TRUE))
}

# 2. Process model (SGDrivE, unchanged) ---------------------------------------
wx_all <- read.csv(file.path(data_dir, "weather_1961_1962", "station_daily_1961_1962_filled.csv"),
  stringsAsFactors = FALSE)
station_weather <- function(station) {
  wx <- subset(wx_all, station_id == station)
  wx <- data.frame(date = as.Date(wx$date), tmax = wx$tmax, tmin = wx$tmin, rain = wx$prcp)
  wx[order(wx$date), ]
}
cube <- sg_cube(list(
  genotypesID = "W", wildType = "W",
  ih = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  tau = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  eta = matrix(1, 1, 1, dimnames = list("W", "W")),
  phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1), omega = c(W = 1)
))

counties <- list(
  Uvalde = list(station = "USC00419265", area = 4040, cases = "uvalde_weekly_1962.csv"),
  Goliad = list(station = "USC00413618", area = 2220, cases = "goliad_weekly_1962.csv")
)

# Daily egg batches over 1 January to 1 July 1962 (182 days = weeks 1-26)
simulate_batches <- function(county, wd, dev = 1, n_imm = 0, egg_temp = "ambient", start = 32) {
  cfg <- counties[[county]]
  wx <- station_weather(cfg$station)
  run <- wx[wx$date >= as.Date("1962-01-01") & wx$date <= as.Date("1962-07-01"), ]
  burn <- wx[wx$date < as.Date("1962-01-01"), ]
  p <- sg_params(preset = "gp2014", egg_temp_source = egg_temp)
  base_slope <- p$dev_rate_slope
  p$wound_density <- wd
  p$dev_rate_slope <- base_slope * dev
  rel <- if (n_imm > 0) {
    sg_releases("W", stage = "adult", sex = "both", times = seq(start, 182, by = 7),
      number = n_imm, per = "node", split = 0.5)
  }
  m <- sg_model(p, cube, sg_forcing(run, p), releases = rel, area = cfg$area)
  sim <- simulate_ode(m, suppressWarnings(sg_init(m, "burn_in", years = 1L, weather = burn)))
  sim$counters[, "egg_batches"]
}

# 3. Scenarios ----------------------------------------------------------------
# par: log wound density, logit rho, log phi, then the scenario's extra parameter
scenarios <- list(
  list(id = "U1", county = "Uvalde", label = "Uvalde: literature development, no immigration", extra = "none"),
  list(id = "U2", county = "Uvalde", label = "Uvalde: literature development + immigration from 1 Feb", extra = "immigration"),
  list(id = "U3", county = "Uvalde", label = "Uvalde: immigration + eggs at host temperature", extra = "immigration", egg = "host"),
  list(id = "U4", county = "Uvalde", label = "Uvalde: development-rate scale fitted (reference)", extra = "dev"),
  list(id = "G1", county = "Goliad", label = "Goliad: literature development, no immigration", extra = "none"),
  list(id = "G2", county = "Goliad", label = "Goliad: literature development + immigration from 1 Feb", extra = "immigration")
)

unpack <- function(par, sc) {
  out <- list(wd = exp(par[1]), rho = plogis(par[2]), phi = exp(par[3]), n_imm = 0, dev = 1)
  if (sc$extra == "immigration") out$n_imm <- exp(par[4])
  if (sc$extra == "dev") out$dev <- exp(par[4])
  out
}

fit_scenario <- function(sc) {
  obs <- subset(read.csv(file.path(data_dir, counties[[sc$county]]$cases)), week <= 26)$cases
  egg <- if (is.null(sc$egg)) "ambient" else sc$egg
  data_negll <- function(par) {
    th <- unpack(par, sc)
    b <- simulate_batches(sc$county, th$wd, th$dev, th$n_imm, egg)
    if (length(b) != 182 || anyNA(b)) return(1e10)
    -nb_loglik(obs, case_mean(b, th$rho), th$phi)
  }
  # objective: data negative log-likelihood, plus the rho prior if set (MAP)
  negll <- function(par) {
    v <- data_negll(par)
    if (!is.null(rho_prior)) v <- v - dbeta(plogis(par[2]), rho_prior[1], rho_prior[2], log = TRUE)
    v
  }
  lower <- c(log(1e-3), qlogis(1e-5), log(0.1))
  upper <- c(log(1e4), qlogis(0.999), log(1e4))
  starts <- if (is.null(rho_prior)) {
    list(c(log(1), qlogis(0.5), log(5)), c(log(10), qlogis(0.02), log(5)))
  } else {
    list(c(log(10), qlogis(0.01), log(5)), c(log(100), qlogis(0.005), log(5)))
  }
  if (sc$extra == "immigration") {
    lower <- c(lower, log(0.1)); upper <- c(upper, log(imm_upper))
    imm_start <- if (is.null(rho_prior)) c(5, 5) else c(100, 2000)
    starts <- Map(function(s, n) c(s, log(n)), starts, imm_start)
  }
  if (sc$extra == "dev") {
    lower <- c(lower, log(0.2)); upper <- c(upper, log(20))
    starts <- lapply(starts, function(s) c(s, log(2)))
  }
  fits <- lapply(starts, function(s) optim(s, negll, method = "L-BFGS-B", lower = lower, upper = upper))
  best <- fits[[which.min(vapply(fits, `[[`, numeric(1), "value"))]]
  th <- unpack(best$par, sc)
  b <- simulate_batches(sc$county, th$wd, th$dev, th$n_imm, egg)
  mu <- case_mean(b, th$rho)
  k <- length(best$par)
  dnll <- data_negll(best$par) # AIC uses the data likelihood only
  list(sc = sc, theta = th, negll = dnll, neg_log_post = best$value, aic = 2 * dnll + 2 * k, k = k,
    at_bound = any(abs(best$par - lower) < 1e-3 | abs(best$par - upper) < 1e-3),
    weekly = data.frame(
      scenario = sc$label, county = sc$county, week = 1:26, observed = obs, expected = mu,
      lo = qnbinom(0.05, mu = mu + 1e-9, size = th$phi),
      hi = qnbinom(0.95, mu = mu + 1e-9, size = th$phi),
      true_batches = vapply(1:26, function(w) sum(b[7 * (w - 1) + 1:7]), numeric(1))
    ))
}

results <- lapply(scenarios, function(sc) {
  cat("fitting", sc$id, "\n")
  fit_scenario(sc)
})

# 4. Summary table ------------------------------------------------------------
summary_tab <- do.call(rbind, lapply(results, function(r) data.frame(
  id = r$sc$id, scenario = r$sc$label, k = r$k,
  neg_loglik = round(r$negll, 2), neg_log_posterior = round(r$neg_log_post, 2), AIC = round(r$aic, 1),
  wound_density = signif(r$theta$wd, 3), rho = signif(r$theta$rho, 3),
  phi = signif(r$theta$phi, 3), immigrants_week = signif(r$theta$n_imm, 3),
  dev_scale = signif(r$theta$dev, 3),
  jan_jun_expected = round(sum(r$weekly$expected)), jan_jun_observed = sum(r$weekly$observed),
  egg_batches_jan_jun = round(sum(r$weekly$true_batches)),
  any_parameter_at_bound = r$at_bound
)))
print(summary_tab, row.names = FALSE)
write.csv(summary_tab, file.path(out_dir, paste0("obs_model_fit_summary", run_tag, ".csv")), row.names = FALSE)

# 5. Figures ------------------------------------------------------------------
blue <- "#2a78d6"; ink <- "#3d3d3a"
weekly <- do.call(rbind, lapply(results, function(r) {
  w <- r$weekly
  w$panel <- sprintf("%s\nAIC %.1f; rho = %.3g; wounds %.3g/km^2%s%s", r$sc$label, r$aic,
    r$theta$rho, r$theta$wd,
    if (r$theta$n_imm > 0) sprintf("; immigrants %.3g/wk", r$theta$n_imm) else "",
    if (r$theta$dev != 1) sprintf("; dev x%.2f", r$theta$dev) else "")
  w
}))
write.csv(weekly, file.path(out_dir, paste0("obs_model_fit_weekly", run_tag, ".csv")), row.names = FALSE)
theme_fit <- theme_minimal(base_size = 11) +
  theme(
    panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(colour = "#e6e5e0", linewidth = 0.3),
    axis.text = element_text(colour = ink), axis.title = element_text(colour = ink),
    strip.text = element_text(colour = ink, face = "bold", hjust = 0, size = 9.5),
    plot.title = element_text(colour = ink, face = "bold"),
    plot.subtitle = element_text(colour = "#6b6a63"),
    legend.position = "top", legend.justification = "left", legend.title = element_blank(),
    plot.background = element_rect(fill = "#fcfcfb", colour = NA)
  )
month_axis <- scale_x_continuous(breaks = c(1, 5, 9, 13, 18, 22, 26),
  labels = c("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul 1"))

plot_county <- function(cty, file, height) {
  d <- subset(weekly, county == cty)
  d$panel <- factor(d$panel, levels = unique(d$panel))
  g <- ggplot(d, aes(week)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, alpha = "90% prediction interval"), fill = blue) +
    geom_col(aes(y = observed, fill = "Observed (USDA weekly cases)"), width = 0.75) +
    geom_line(aes(y = expected, colour = "Expected reported cases"), linewidth = 0.9) +
    geom_point(aes(y = expected, colour = "Expected reported cases"), size = 1.2) +
    scale_fill_manual(values = c("Observed (USDA weekly cases)" = "#cfcdc4")) +
    scale_colour_manual(values = c("Expected reported cases" = blue)) +
    scale_alpha_manual(values = c("90% prediction interval" = 0.18)) +
    facet_wrap(~panel, ncol = 1, scales = "free_y") + month_axis +
    labs(x = NULL, y = "Screwworm cases per week",
      title = sprintf("%s County, January–June 1962: fits with the observation model", cty),
      subtitle = paste0(paste0("Cases ~ NegBin(rho x delayed egg batches, phi); delay from laying to report: gamma, mean 6 d, sd 3 d.", if (!is.null(rho_prior)) sprintf("\nPrior rho ~ Beta(%g, %g).", rho_prior[1], rho_prior[2]) else "", "\n"),
        "Lower AIC is better. Each panel has its own y-axis.")) +
    theme_fit
  ggsave(file.path(out_dir, file), g, width = 9, height = height, dpi = 160)
}
plot_county("Uvalde", paste0("obs_fit_uvalde", run_tag, ".png"), 12)
plot_county("Goliad", paste0("obs_fit_goliad", run_tag, ".png"), 6.8)
cat("saved to", out_dir, "\n")
