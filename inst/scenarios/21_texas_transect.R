# Scenario 21: a south-north transect through Texas, 1961-1962
#
# The first landscape of the multi-node model (design note
# dev/MOVEMENT_DESIGN.md, decision 5): 18 cells of 50 x 50 km in a line along
# 98.2 degW from the Rio Grande valley (26.5 degN) to the Red River (34.1
# degN), each with the daily weather of the GHCN station assigned to the
# nearest county, and a source node south of them with Brownsville weather
# standing for the overwintering area in the lower Rio Grande valley and
# northern Mexico (design note 5.5). Adults move by local diffusion only, at
# the rate measured in the field (Matlock & Skoda 2009: about 16% of females
# and 10% of males leave a station array per day, D of about 0.4 km^2 per
# day). The 1961-62 winter was severe; the questions are:
#   1. do northern cells die out in winter and get recolonised in spring?
#   2. how fast does the front move north in 1962, against the 80-160 km per
#      generation observed in the 1960s (FAO)?
# Expected (design note section 2): local diffusion alone is far too slow;
# fast spread needs long-distance movement (wind, livestock), phase 4.
#
# The ODE runs with the node extinction rule (nodes below one fly are set
# to zero at the end of each day, design note 5.7) and the chain-binomial
# solver with 10 replicates. Scenario 22 varies the movement clock, D and
# the time step, and shows the ODE without the extinction rule.
#
# Data (not shipped with the package; see data-raw/texas_nws/): run from
# the repository root. Run time: about a minute.

library(ScrewDrivE)

# 0. Paths ---------------------------------------------------------------------
data_dir <- "data-raw/texas_nws/weather_1961_1962"
fig_dir <- "inst/scenarios/scenario21_figures"
dir.create(fig_dir, showWarnings = FALSE)

# 1. Cells and their weather -----------------------------------------------------
# Cell i is centred i * 50 km north of the source node, along 98.2 degW. Each
# cell takes the station of the county whose centroid is nearest its centre
# (great-circle distance); neighbouring cells may share a station.
h <- 50 # km
n_cells <- 18
counties <- read.csv(file.path(data_dir, "county_station_map_1961_1962.csv"))
great_circle <- function(lat1, lon1, lat2, lon2) {
  r <- pi / 180
  cos_angle <- sin(lat1 * r) * sin(lat2 * r) +
    cos(lat1 * r) * cos(lat2 * r) * cos((lon1 - lon2) * r)
  6371 * acos(pmin(1, cos_angle))
}
cells <- data.frame(
  id = sprintf("c%02d", seq_len(n_cells)),
  lat = 26.25 + (seq_len(n_cells) - 0.5) * h / 111.2, lon = -98.2,
  x = seq_len(n_cells) * h # km north of the source node
)
nearest <- vapply(seq_len(n_cells), function(i) {
  which.min(great_circle(cells$lat[i], cells$lon[i], counties$county_lat, counties$county_lon))
}, 1L)
cells$county <- counties$county[nearest]
cells$station <- counties$station_id[nearest]
source_station <- counties$station_id[counties$county == "Cameron"] # BROWNSVILLE
print(cells[c("id", "lat", "county", "station")], row.names = FALSE)

weather_all <- read.csv(file.path(data_dir, "station_daily_1961_1962_filled.csv"))
station_weather <- function(station) {
  w <- weather_all[weather_all$station_id == station, ]
  w <- data.frame(date = as.Date(w$date), tmax = w$tmax, tmin = w$tmin, rain = w$prcp)
  w[order(w$date), ]
}

# 2. Model -------------------------------------------------------------------------
# Default parameters, one wild genotype, static wounds (5 per km^2, so 12,500
# wounds per cell: realistic numbers, not a smoke test).
p <- sg_params()
cube <- sg_cube(list(
  genotypesID = "W", wildType = "W",
  ih = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  tau = array(1, c(1, 1, 1), dimnames = list("W", "W", "W")),
  eta = matrix(1, 1, 1, dimnames = list("W", "W")),
  phi = c(W = 0.5), xiF = c(W = 1), xiM = c(W = 1), s = c(W = 1), omega = c(W = 1)
))
stations <- c(source_station, cells$station)
forcing <- lapply(stations, function(st) sg_forcing(station_weather(st), p)) # 1961-1962
nodes <- data.frame(
  id = c("source", cells$id), area = h^2, x = c(0, cells$x),
  role = c("source", rep("study", n_cells))
)
land <- sg_landscape(nodes, forcing)

# Local diffusion between neighbouring cells (shared border 50 km, centroids
# 50 km apart), so each neighbour gets D / h^2 per day at 27 degC; the north
# end is open (flies crossing it leave the transect), the source's south side
# and the transect's sides are closed. The physiological clock with the
# default cap (movement_kappa_max = 1): less movement below 27 degC, none
# below the development threshold.
d_coef <- list(male = 0.25, virgin = 0.4, mated = 0.4) # km^2 per day
layer <- sg_movement_diffusion(
  stats::setNames(nodes$area, nodes$id),
  data.frame(from = nodes$id[-nrow(nodes)], to = nodes$id[-1L], border = h, dist = h),
  D = d_coef, outer_border = c(rep(0, n_cells), h)
)
model <- sg_model(p, cube, landscape = land, movement = layer)
model

# What these rates mean for a fly: the chance to leave its cell in a lifetime
# is well below 1%, so cells exchange very few flies.
mc <- check_movement(model, temperatures = NULL)
mc

# 3. Start and runs ------------------------------------------------------------
# Every node starts on 1 January 1961 with the same densities; 1961 then sets
# the populations that enter the 1961-62 winter (a burn-in year), and 1962 is
# the year reported.
x0 <- sg_init(model, values = c(E = 2, L = 2, P = 4, M = 1, U = 0.2, F = 1), per = "km2")
ode <- simulate_ode(model, x0, node_extinction = TRUE, rtol = 1e-6, atol = 1e-6)
cb <- simulate_cb(model, x0, n_rep = 10, seed = 21)

# 4. Die-back, recolonisation and the front ------------------------------------
# A cell is occupied with at least one fly (any stage). The front is the
# northernmost occupied cell, in km north of the source node. Generation time
# at 27 degC: egg, larva, pupa and the pre-oviposition period.
gen_days <- (p$egg_duration_dd + p$pupa_duration_dd + p$oviposition_age_dd) /
  sg_degree_days(27, p) + p$larva_duration_days
winter <- as.Date(c("1961-11-01", "1962-04-30"))

transect_summary <- function(tt) {
  # tt: the time table of sg_spread() for one run or replicate
  w <- tt[tt$date >= winter[1L] & tt$date <= winter[2L], ]
  low <- w[which.min(w$front), ]
  y62 <- tt[tt$date > low$date & tt$date <= as.Date("1962-12-31"), ]
  top <- max(y62$front)
  reach <- y62$date[which(y62$front >= top)[1L]]
  climb <- y62[y62$date <= reach, ]
  speed <- if (top > low$front && nrow(climb) > 2L) {
    unname(stats::coef(stats::lm(front ~ as.numeric(date), climb))[2L])
  } else {
    0
  }
  data.frame(
    winter_min_occupied = min(w$occupied), front_winter_min_km = low$front,
    date_winter_min = low$date, front_1962_08_01_km = tt$front[tt$date == as.Date("1962-08-01")],
    front_speed_km_per_day = speed, front_speed_km_per_generation = speed * gen_days
  )
}

all_nodes <- seq_len(nrow(nodes))
sp_ode <- sg_spread(ode, nodes = all_nodes)
sp_cb <- sg_spread(cb, nodes = all_nodes)
results <- rbind(
  cbind(run = "ODE, extinction rule", transect_summary(sp_ode$time)),
  do.call(rbind, lapply(seq_len(10), function(r) {
    tt <- sp_cb$time[sp_cb$time$replicate == r, ]
    cbind(run = paste("chain-binomial, replicate", r), transect_summary(tt))
  }))
)
print(results, row.names = FALSE, digits = 3)
cat(sprintf(
  paste0(
    "\nGeneration time at 27 degC: %.1f days. Front speed in 1962 (chain-binomial, mean): ",
    "%.2f km per day = %.0f km per generation, against 80-160 km per generation observed.\n"
  ),
  gen_days, mean(results$front_speed_km_per_day[-1L]),
  mean(results$front_speed_km_per_generation[-1L])
))
utils::write.csv(results, file.path(fig_dir, "transect_summary.csv"), row.names = FALSE)

# 5. Figures ---------------------------------------------------------------------
if (requireNamespace("ggplot2", quietly = TRUE)) {
  heat <- plot(ode, "nodes", nodes = all_nodes, stages = "adults") +
    ggplot2::labs(subtitle = "ODE with the node extinction rule; cells c01 (south) to c18 (north)")
  ggplot2::ggsave(file.path(fig_dir, "fig1_adults_by_cell.png"), heat,
    width = 9, height = 5, dpi = 150
  )
  fronts <- rbind(
    data.frame(date = sp_ode$time$date, front = sp_ode$time$front, run = "ODE", line = "ODE"),
    data.frame(
      date = sp_cb$time$date, front = sp_cb$time$front, run = "chain-binomial",
      line = paste("cb", sp_cb$time$replicate)
    )
  )
  front_plot <- ggplot2::ggplot(fronts, ggplot2::aes(date, front, group = line, colour = run)) +
    ggplot2::geom_step(alpha = 0.7) +
    ggplot2::scale_colour_manual(values = c(ODE = "#000000", `chain-binomial` = "#E69F00")) +
    ggplot2::labs(
      x = NULL, y = "km north of the source node", colour = NULL,
      title = "Northernmost cell with flies",
      subtitle = "Local diffusion only (D = 0.4 km^2/day); 10 chain-binomial replicates"
    ) +
    ggplot2::theme_minimal()
  ggplot2::ggsave(file.path(fig_dir, "fig2_front.png"), front_plot,
    width = 9, height = 5, dpi = 150
  )
}
