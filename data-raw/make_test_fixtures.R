# Build the small reference fixtures used by the WP3 (weather and forcing)
# tests in tests/testthat/fixtures/.
#
# Run from the package root:  Rscript data-raw/make_test_fixtures.R
#
# Inputs:
# * the unchanged legacy Gutierrez & Ponti R translation in
#   data-raw/reference_impl/scripts/ (temperature functions);
# * data-raw/reference/gp2014_mcallen_discrete.rds (legacy McAllen run);
# * the raw NOAA GHCN-Daily station files USC00415701.csv and
#   USW00012959.csv, from the directory in the environment variable
#   SGDRIVE_RAW_WEATHER (default: the screwworm project folder on Box). They
#   are not stored in this repository (13 MB); download them from
#   https://www.ncei.noaa.gov/data/global-historical-climatology-network-daily/access/

stopifnot(file.exists("DESCRIPTION"))
fix_dir <- file.path("tests", "testthat", "fixtures")
dir.create(fix_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. legacy temperature functions on a grid --------------------------------
legacy <- new.env()
sys.source("data-raw/reference_impl/scripts/gutierrez_ponti_2014_pbdm.R", envir = legacy)
gp <- legacy$gp_params()
temp <- seq(-15, 50, by = 0.25)
grid <- data.frame(
  temp = temp,
  degree_days = legacy$gp_degree_days(temp, gp),
  temp_mortality = legacy$gp_temp_mortality(temp, gp),
  temp_survival = legacy$gp_temp_survival(temp, gp),
  egg_viability = legacy$gp_egg_viability(temp, gp),
  ovip_scalar = legacy$gp_ff_temperature(temp, gp$TL, gp$TU_ovip)
)
utils::write.csv(grid, file.path(fix_dir, "gp2014_temperature_grid.csv"), row.names = FALSE)

# ---- 2. legacy daily forcing terms on McAllen weather -------------------------
ref <- readRDS(file.path("data-raw", "reference", "gp2014_mcallen_discrete.rds"))
forcing <- ref[, c("date", "tmean", "aver5", "dd", "temp_mort", "temp_surv", "egg_viab", "phi_temp")]
saveRDS(forcing, file.path(fix_dir, "gp2014_mcallen_forcing.rds"), compress = "xz")

# ---- 3. raw NOAA extracts around filled periods -------------------------------
raw_dir <- Sys.getenv(
  "SGDRIVE_RAW_WEATHER",
  file.path(
    "/Users/prateekverma/Library/CloudStorage/Box-Box/Research Projects/Screwworm",
    "Simulations/MGDrivE1-screwworm/data/weather"
  )
)
# Windows containing interpolated days (1960-09-24) and backup-filled days
# (April 1975) of the bundled McAllen series.
windows <- list(
  c("1960-08-15", "1960-10-31"),
  c("1975-03-15", "1975-05-31")
)
keep_cols <- c("STATION", "DATE", "NAME", "PRCP", "TMAX", "TMIN")
for (st in c("USC00415701", "USW00012959")) {
  raw <- utils::read.csv(file.path(raw_dir, paste0(st, ".csv")),
    colClasses = "character", check.names = FALSE
  )
  d <- as.Date(raw$DATE)
  sel <- Reduce(`|`, lapply(windows, function(w) d >= as.Date(w[1]) & d <= as.Date(w[2])))
  sub <- raw[sel, intersect(keep_cols, names(raw))]
  utils::write.csv(sub, file.path(fix_dir, paste0("ghcn_", st, "_subset.csv")), row.names = FALSE)
  message(st, ": ", nrow(sub), " rows")
}
message("Fixtures written to ", fix_dir)
