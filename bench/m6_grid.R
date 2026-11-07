# Shared setup of the multi-node benchmarks (bench/m6_*.R; phase 3 WP7):
# regular grids of 50 km cells with the SIT model (2 genotypes) and one year
# of McAllen weather (1971). Sourced from the package root after
# devtools::load_all().

source("tests/testthat/helper-cubes.R")

m6_params <- sg_params()
m6_year <- local({
  wx <- sg_example_weather()
  wx[format(wx$date, "%Y") == "1971", ]
})
m6_cube <- sg_cube(sg_test_cube_sit(residual = 0))
m6_h <- 50 # km

# n cells in a near-square grid (rows filled left to right), 50 km apart.
# weather: "distinct" (McAllen shifted from -1 to +1 degC across the nodes,
# one forcing series per node) or "shared" (one series for all); d_coef:
# diffusion coefficient (km^2 per day) for all adults; movement = FALSE for
# a landscape without movement; params: model parameters.
m6_grid_model <- function(n, weather = c("distinct", "shared"), d_coef = 0.4, movement = TRUE,
                          params = m6_params) {
  weather <- match.arg(weather)
  h <- m6_h
  k <- ceiling(sqrt(n))
  xy <- data.frame(x = ((seq_len(n) - 1L) %% k) * h, y = ((seq_len(n) - 1L) %/% k) * h)
  ids <- paste0("c", seq_len(n))
  dd <- as.matrix(stats::dist(xy))
  nb <- which(abs(dd - h) < 1e-9 & upper.tri(dd), arr.ind = TRUE)
  forcing <- if (weather == "shared" || n == 1L) {
    sg_forcing(m6_year, params)
  } else {
    lapply(seq_len(n), function(i) {
      w <- m6_year
      shift <- 2 * (i - 1) / (n - 1) - 1
      w[c("tmax", "tmin", "tmean")] <- w[c("tmax", "tmin", "tmean")] + shift
      sg_forcing(w, params)
    })
  }
  land <- sg_landscape(data.frame(id = ids, area = h^2, x = xy$x, y = xy$y), forcing)
  if (n == 1L || !movement) {
    return(sg_model(params, m6_cube, landscape = land))
  }
  edges <- data.frame(from = ids[nb[, 1L]], to = ids[nb[, 2L]], border = h, dist = h)
  layer <- sg_movement_diffusion(stats::setNames(rep(h^2, n), ids), edges, D = d_coef)
  sg_model(params, m6_cube, landscape = land, movement = layer)
}

# A start with the same densities in every node (per km^2), or in node 1 only.
m6_init <- function(m, where = c("all", "corner")) {
  where <- match.arg(where)
  v <- c(E = 5, L = 5, P = 5, M = 2, U = 0.5, F = 2)
  if (where == "all") sg_init(m, values = v, per = "km2") else sg_init(m, values = v, per = "km2", nodes = 1)
}
