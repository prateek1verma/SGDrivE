# Validation of the chain-binomial solver on landscapes (phase 2 WP6; D-032):
# the deterministic map converges to the multi-node ODE as dt -> 0, the
# replicate mean is close to the deterministic map, a sink node goes extinct
# and is recolonised, and a node behind a barrier stays empty.

cbv_cube <- sg_cube(sg_test_cube_mendelian())

cbv_model <- function(clock = "physiological", hosts = FALSE, days = 30) {
  p <- if (hosts) {
    sg_params(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  } else {
    sg_params()
  }
  start <- as.Date("1971-05-01")
  wx <- sg_forcing(sg_example_weather(), p, start = start, end = start + days - 1)
  mild <- sg_forcing(sg_constant_weather(22, days = days, start = start), p)
  nd <- data.frame(id = c("a", "b", "c"), area = c(40, 60, 80))
  lay <- sg_movement_layer(
    list(male = 0.3, virgin = 0.5, mated = 0.4),
    rbind(c(0, 1, 0), c(0.5, 0, 0.3), c(0, 1, 0)),
    clock = clock
  )
  sg_model(p, cbv_cube, landscape = sg_landscape(nd, list(wx, wx, mild)), movement = lay)
}

cbv_init <- function(m) {
  v <- c(E = 200, L = 200, P = 200, M = 100, U = 20, F = 100)
  as.numeric(sg_init(m, values = v, nodes = "a"))
}

test_that("the deterministic map converges to the multi-node ODE as dt -> 0", {
  for (clock in c("physiological", "calendar")) {
    for (hosts in c(FALSE, TRUE)) {
      m <- cbv_model(clock, hosts)
      x0 <- cbv_init(m)
      ode <- simulate_ode(m, x0)
      nm <- sg_state_names(m)
      fly <- !.sg_is_host_name(nm)
      move <- grepl("^(emigrants|immigrants|lost)_", colnames(ode$counters))
      err <- function(dt) {
        d <- .sg_cb_run_cpp(m, init = x0, dt = dt, stochastic = FALSE)
        node_err <- max(vapply(c("@a", "@b", "@c"), function(id) {
          cols <- fly & endsWith(nm, id)
          ref <- rowSums(ode$state[, cols])
          max(abs(rowSums(d$state[, cols]) - ref)) / max(ref)
        }, 0))
        k_ode <- colSums(ode$counters[, move])
        move_err <- max(abs(colSums(d$counters[, move]) - k_ode) / pmax(k_ode, 1))
        c(node = node_err, move = move_err)
      }
      e1 <- err(0.1)
      e2 <- err(0.05)
      info <- paste(clock, hosts)
      # first order in dt: halving the step about halves the differences
      expect_lt(e2[["node"]], 0.6 * e1[["node"]], label = paste("node", info))
      expect_lt(e2[["move"]], 0.6 * e1[["move"]], label = paste("movement", info))
      expect_lt(e2[["node"]], 0.01, label = info)
      # movement counters are first order too; unmated females are the
      # slowest (they emerge after the movement of a step and stay unmated
      # for about a day), about 4.5% at dt = 0.05
      expect_lt(e2[["move"]], 0.06, label = info)
    }
  }
})

test_that("the replicate mean is close to the deterministic map", {
  m <- cbv_model(days = 12)
  x0 <- cbv_init(m)
  det <- .sg_cb_run_cpp(m, init = x0, dt = 0.1, stochastic = FALSE)
  sim <- simulate_cb(m, x0, n_rep = 150, seed = 3, record = sg_record(
    sg_quantity("adults", stage = "adults", by = "node", node = c("a", "b", "c"))
  ))
  nm <- sg_state_names(m)
  for (id in c("a", "b", "c")) {
    adult <- substr(nm, 1L, 1L) %in% c("M", "U", "F") & endsWith(nm, paste0("@", id))
    ref <- sum(det$state[13L, adult])
    v <- sim$recorded[13L, paste0("adults@", id), ]
    se <- stats::sd(v) / sqrt(length(v))
    # within sampling error and the small bias of a nonlinear model
    expect_lt(abs(mean(v) - ref), 4 * se + 0.01 * ref, label = id)
  }
})

# Node a: a large population at 25 degC; node b: a sink without wounds (no
# eggs), reached by rare emigrants from a; node c: behind a barrier (no
# edges, no movement), empty at the start.
cbv_metapop <- function() {
  p <- sg_params()
  warm <- sg_forcing(sg_constant_weather(25, days = 150, start = as.Date("1971-05-01")), p)
  nd <- data.frame(id = c("a", "b", "c"), area = c(20, 2, 5), wound_density = c(5, 0, 5))
  rate <- c(3e-4, 0.3, 0)
  lay <- sg_movement_layer(list(male = rate, virgin = rate, mated = rate),
    rbind(c(0, 1, 0), c(1, 0, 0), c(0, 0, 0)),
    clock = "calendar"
  )
  sg_model(p, cbv_cube, landscape = sg_landscape(nd, warm), movement = lay)
}

test_that("a sink node goes extinct and is recolonised; a barrier node stays empty", {
  m <- cbv_metapop()
  x0 <- suppressMessages(sg_init(m, "equilibrium"))
  x0[!endsWith(names(x0), "@a")] <- 0 # flies in node a only
  sim <- simulate_cb(m, x0, n_rep = 5, seed = 1, record = sg_record(
    "counters", sg_quantity("flies", stage = "all", by = "node", node = c("a", "b", "c"))
  ))
  # node a persists
  expect_gt(min(sim$recorded[, "flies@a", ]), 0)
  # node b: occupied on some days and empty on others in every replicate,
  # with extinctions followed by recolonisation
  b <- sim$recorded[, "flies@b", ]
  occ <- colMeans(b > 0)
  expect_true(all(occ > 0.05 & occ < 0.95))
  ext <- apply(b > 0, 2L, function(o) sum(o[-length(o)] & !o[-1L]))
  recol <- apply(b > 0, 2L, function(o) {
    first_ext <- which(o[-length(o)] & !o[-1L])[1L]
    !is.na(first_ext) && any(o[seq.int(first_ext + 1L, length(o))])
  })
  expect_true(all(ext >= 2L))
  expect_true(all(recol))
  # node c: no flies on any day, no immigrants
  expect_identical(max(sim$recorded[, "flies@c", ]), 0)
  expect_identical(sum(sim$counters[, paste0("immigrants_", c("M", "U", "F"), "@c"), ]), 0)
})
