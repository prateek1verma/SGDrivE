# Host-wound model (D-030): state layout, model assembly and initial states.

cube3 <- sg_cube(sg_test_cube_mendelian())
host_params <- function(...) {
  base <- list(wound_model = "hosts", host_density = 50, wounded_fraction = 0.1)
  do.call(sg_params, utils::modifyList(base, list(...)))
}
model_of <- function(p, area = 1) {
  sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 10), p), area = area)
}

test_that("static models keep the fly layout and have no host blocks", {
  m <- model_of(sg_params())
  expect_identical(m$layout$block, c("E", "L", "P", "M", "U", "F"))
  expect_null(m$hosts)
  expect_identical(m$dims$n_I, 0L)
  expect_error(sg_state_index(m, "HW"), "exists only with wound_model = \"hosts\"")
  expect_false(any(grepl("^H", sg_state_names(m))))
})

test_that("host models append HS, HW and HI to the layout", {
  m_static <- model_of(sg_params())
  m <- model_of(host_params(infestation_substages = 4))
  expect_identical(m$layout$block, c("E", "L", "P", "M", "U", "F", "HS", "HW", "HI"))
  expect_identical(m$layout$length[7:9], c(1L, 1L, 4L))
  expect_identical(m$n_states, m_static$n_states + 6L)
  # fly blocks are where they were
  expect_identical(m$layout[1:6, ], m_static$layout)
  nm <- sg_state_names(m)
  expect_identical(
    tail(nm, 6),
    c("HS[hosts]", "HW[hosts]", paste0("HI[hosts,", 1:4, "]"))
  )
  expect_identical(nm[seq_len(m_static$n_states)], sg_state_names(m_static))
  expect_identical(sg_state_index(m, "HS"), m_static$n_states + 1L)
  expect_identical(sg_state_index(m, "HI", substage = 2:3), m_static$n_states + 4:5)
  expect_error(sg_state_index(m, "HI", substage = 5), "between 1 and 4")
  expect_error(sg_state_index(m, "HW", genotype = "AA"), "no genotype")
})

test_that("host counts and rates follow the parameters", {
  m <- model_of(host_params(
    host_density = 40, wounded_fraction = 0.2, wound_healing_days = 5,
    infestation_days = 8, infestation_substages = 4,
    infested_wound_attractiveness = 3, treatment_rate = 0.1
  ), area = 2.5)
  h <- m$hosts
  expect_identical(h$classes, "hosts")
  expect_equal(h$n, 100) # 40 per km^2 on 2.5 km^2
  expect_equal(h$healing, 1 / 5)
  expect_equal(h$wounding, (1 / 5) * 0.2 / 0.8) # wounded fraction 0.2 without flies
  expect_equal(h$wounding / (h$wounding + h$healing), 0.2)
  expect_equal(h$advance, 4 / 8)
  expect_identical(h$treatment, 0.1)
  expect_identical(h$attractiveness, 3)
  expect_true(is.na(m$counts$H))
  out <- capture.output(print(m))
  expect_true(any(grepl("hosts N_H = 100 per node", out)))
  expect_true(any(grepl("host-wound model: 0.2 of hosts wounded", out)))
  # the hash sees the host parameters
  expect_false(identical(m$meta$hash, model_of(host_params())$meta$hash))
})

test_that("constant starts put hosts at their fly-free steady state", {
  m <- model_of(host_params(host_density = 30, wounded_fraction = 0.1))
  m_static <- model_of(sg_params())
  x <- sg_init(m)
  expect_equal(unname(x[sg_state_index(m, "HS")]), 27)
  expect_equal(unname(x[sg_state_index(m, "HW")]), 3)
  expect_true(all(x[sg_state_index(m, "HI")] == 0))
  # the fly part is the same as in the static model
  fly <- seq_len(m_static$n_states)
  expect_equal(unname(x[fly]), unname(as.numeric(sg_init(m_static))))
  printed <- capture.output(print(x))
  expect_true(any(grepl("hosts: 27 without a wound, 3 wounded, 0 infested", printed)))
  hw <- sg_init(m, "hardy_weinberg", alleles = c(A = 0.5, a = 0.5))
  expect_equal(sum(hw[SGDrivE:::.sg_host_index(m)]), 30)
})

test_that("a full state is checked against the number of hosts", {
  m <- model_of(host_params(host_density = 30, wounded_fraction = 0.1))
  x <- as.numeric(sg_init(m))
  expect_silent(sg_init(m, "state", state = x))
  x[sg_state_index(m, "HI", substage = 1)] <- 2
  expect_warning(sg_init(m, "state", state = x), "host states sum to 32")
  expect_error(sg_init(m, "state", state = x[-1]), "length")
})

test_that("check_model() refuses host models until WP7", {
  m <- model_of(host_params())
  expect_error(check_model(m), "not available yet")
})

test_that("the model view carries the host block to C++", {
  m <- model_of(host_params(infestation_substages = 3, treatment_rate = 0.2))
  v <- SGDrivE:::.sg_model_view(m)
  info <- SGDrivE:::sg_model_view_info_cpp(v)
  expect_true(info$host_model)
  expect_identical(info$n_I, 3L)
  expect_identical(info$host_offsets, as.integer(sg_state_index(m, "HS") - 1L) + 0:2)
  h <- m$hosts
  expect_equal(
    info$host_rates,
    c(h$n, h$wounding, h$healing, h$advance, h$treatment, h$attractiveness)
  )
  expect_identical(info$dims[6], m$n_states)
  # static models: no host block
  info0 <- SGDrivE:::sg_model_view_info_cpp(SGDrivE:::.sg_model_view(model_of(sg_params())))
  expect_false(info0$host_model)
  expect_identical(info0$host_offsets, rep(-1L, 3))
  # the C++ side checks the host layout
  bad <- v
  bad$host_offsets <- v$host_offsets + 1L
  expect_error(SGDrivE:::sg_model_view_info_cpp(bad), "host block")
})

# ---- reference ODE (WP3) ------------------------------------------------------

host_cols <- function(sim, block) grep(paste0("^", block, "\\["), colnames(sim$state))
host_total <- function(sim, block) rowSums(sim$state[, host_cols(sim, block), drop = FALSE])
run_r <- function(p, init = NULL, days = 120, area = 10, temp = 25) {
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(temp, days = days), p), area = area)
  x0 <- if (is.null(init)) sg_init(m) else init(m)
  list(model = m, sim = simulate_ode(m, init = x0, backend = "R"))
}

test_that("hosts are conserved and host counters account for every move", {
  r <- run_r(host_params(infested_wound_attractiveness = 2, treatment_rate = 0.05))
  s <- r$sim
  tot <- host_total(s, "HS") + host_total(s, "HW") + host_total(s, "HI")
  expect_lt(max(abs(tot - 500)), 1e-8)
  k <- colSums(s$counters)
  change <- function(b) {
    v <- host_total(s, b)
    v[length(v)] - v[1]
  }
  flow_s <- -k[["hosts_wounded"]] + k[["wounds_healed"]] + k[["infestations_ended"]] +
    k[["infestations_treated"]]
  flow_w <- k[["hosts_wounded"]] - k[["wounds_healed"]] - k[["infestations"]]
  flow_i <- k[["infestations"]] - k[["infestations_ended"]] - k[["infestations_treated"]]
  expect_equal(change("HS"), flow_s, tolerance = 1e-6)
  expect_equal(change("HW"), flow_w, tolerance = 1e-6)
  expect_equal(change("HI"), flow_i, tolerance = 1e-6)
  expect_gt(k[["infestations"]], 0)
  # every fertile batch is either a new infestation or a superinfestation
  expect_equal(
    k[["infestations"]] + k[["superinfestations"]], k[["fertile_batches"]],
    tolerance = 1e-8
  )
  # the fly flow balance still closes
  expect_lt(max(abs(sg_flow_balance(s)$residual)), 1e-6)
})

test_that("without flies hosts relax to the fly-free steady state", {
  p <- host_params(host_density = 20, wounded_fraction = 0.2, wound_healing_days = 4)
  r <- run_r(p, days = 200, area = 1, init = function(m) {
    x <- numeric(m$n_states)
    x[sg_state_index(m, "HI", substage = 1)] <- 20 # all hosts infested, no flies
    x
  })
  last <- nrow(r$sim$state)
  expect_equal(host_total(r$sim, "HW")[last], 4, tolerance = 1e-6)
  expect_equal(host_total(r$sim, "HS")[last], 16, tolerance = 1e-6)
  expect_lt(host_total(r$sim, "HI")[last], 1e-8)
  expect_identical(sum(r$sim$counters[, "infestations"]), 0)
})

test_that("infestations last infestation_days on average, with Erlang spread", {
  for (n_i in c(1L, 4L)) {
    p <- host_params(infestation_days = 6, infestation_substages = n_i)
    r <- run_r(p, days = 150, area = 1, init = function(m) {
      x <- numeric(m$n_states)
      x[sg_state_index(m, "HI", substage = 1)] <- m$hosts$n
      x
    })
    ended <- r$sim$counters[, "infestations_ended"]
    expect_equal(sum(ended), 50, tolerance = 1e-6)
    t_mid <- seq_along(ended) - 0.5
    mean_t <- sum(t_mid * ended) / 50
    var_t <- sum((t_mid - mean_t)^2 * ended) / 50 - 1 / 12 # minus daily binning
    expect_equal(mean_t, 6, tolerance = 0.01)
    expect_equal(var_t, 36 / n_i, tolerance = 0.03)
  }
})

test_that("treatment shortens infestations and kills larvae", {
  start_inf <- function(m) {
    x <- numeric(m$n_states)
    x[sg_state_index(m, "HI", substage = 1)] <- m$hosts$n
    x
  }
  kappa <- 0.2
  r <- run_r(host_params(treatment_rate = kappa, infestation_substages = 4),
    days = 150, area = 1, init = start_inf
  )
  # mean time infested D = (1 - (q / (q + kappa))^n) / kappa (design note 3.7)
  q <- 4 / 6
  d_expected <- (1 - (q / (q + kappa))^4) / kappa
  i_t <- host_total(r$sim, "HI")
  d_sim <- sum((i_t[-1] + i_t[-length(i_t)]) / 2) / 50 # trapezoid over days
  expect_equal(d_sim, d_expected, tolerance = 0.01)
  k <- colSums(r$sim$counters)
  expect_equal(k[["infestations_treated"]] + k[["infestations_ended"]], 50, tolerance = 1e-6)
  # with flies: treatment lowers prevalence and larvae
  base <- run_r(host_params())$sim
  treated <- run_r(host_params(treatment_rate = 0.3))$sim
  last <- nrow(base$state)
  expect_lt(host_total(treated, "HI")[last], host_total(base, "HI")[last])
  larvae <- function(s) sum(s$state[last, grep("^L\\[", colnames(s$state))])
  expect_lt(larvae(treated), larvae(base))
})

test_that("rates: attractive sites replace H, crowding per infested host", {
  p_h <- host_params(host_density = 50, wounded_fraction = 0.1)
  p_s <- sg_params(wound_density = 5) # = host_density * wounded_fraction
  m_h <- model_of(p_h, area = 2)
  m_s <- model_of(p_s, area = 2)
  x_s <- as.numeric(sg_init(m_s, values = c(E = 3, L = 3, P = 3, M = 2, F = 2)))
  x_h <- as.numeric(sg_init(m_h, values = c(E = 3, L = 3, P = 3, M = 2, F = 2)))
  fly <- seq_len(m_s$n_states)
  expect_identical(x_h[fly], x_s[fly])
  r_s <- SGDrivE:::.sg_rates_r(1, x_s, SGDrivE:::.sg_reference_setup(m_s), aux = TRUE)
  r_h <- SGDrivE:::.sg_rates_r(1, x_h, SGDrivE:::.sg_reference_setup(m_h), aux = TRUE)
  # no infested hosts: A_s = wounded hosts = H, so the flies see the same rates
  expect_equal(r_h$aux$sites, 10)
  expect_equal(r_h$aux$psi, r_s$aux$psi, tolerance = 1e-14)
  expect_equal(r_h$deriv[fly], r_s$deriv[fly], tolerance = 1e-14)
  expect_equal(r_h$counters[names(r_s$counters)], r_s$counters, tolerance = 1e-14)
  expect_false(any(SGDrivE:::.sg_host_counter_names %in% names(r_s$counters)))
  # crowding per infested host, plus treatment in the larval hazard
  p_c <- host_params(
    larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.02,
    treatment_rate = 0.1, infested_wound_attractiveness = 3
  )
  m_c <- model_of(p_c, area = 2)
  x_c <- x_h
  x_c[sg_state_index(m_c, "HI", substage = 2)] <- 4
  x_c[sg_state_index(m_c, "HS")] <- x_c[sg_state_index(m_c, "HS")] - 4
  a <- SGDrivE:::.sg_rates_r(1, x_c, SGDrivE:::.sg_reference_setup(m_c), aux = TRUE)$aux
  l_tot <- sum(x_c[sg_state_index(m_c, "L")])
  expect_equal(
    a$mu[["L"]], -log(0.95) + 0.1 + 0.02 * l_tot / (4 + 1e-6 * 100),
    tolerance = 1e-14
  )
  expect_equal(a$sites, 10 + 3 * 4)
})

test_that("at most N_H hosts are infested, however many flies", {
  # a wound receives at most wound_finding_max * wound_search_rate * dx batches
  # per day; with the GP2014 search rate (0.001 per dd) that is about 0.01 per
  # day, so a fast search rate is used to saturate the two hosts
  p <- host_params(
    host_density = 2, wounded_fraction = 0.9, infested_wound_attractiveness = 5,
    wound_search_rate = 1
  )
  r <- run_r(p, days = 60, area = 1, init = function(m) {
    sg_init(m, values = c(E = 500, L = 500, P = 500, M = 300, F = 300))
  })
  i_t <- host_total(r$sim, "HI")
  expect_lte(max(i_t), 2 + 1e-9)
  # hosts are infested as soon as they are wounded: infested share about
  # 6 / (6 + 7 / 9) = 0.89 of the two hosts at wounded_fraction 0.9
  expect_gt(max(i_t), 1.5)
  k <- colSums(r$sim$counters)
  expect_gt(k[["superinfestations"]], 10 * k[["infestations"]])
  expect_gt(min(r$sim$state), -1e-6)
})

# ---- C++ kernel (WP4) ---------------------------------------------------------

host_configs <- list(
  plain = host_params(),
  full = host_params(
    infested_wound_attractiveness = 2.5, treatment_rate = 0.05, infestation_substages = 4,
    larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01
  ),
  no_attraction = host_params(
    infested_wound_attractiveness = 0, infestation_substages = 1,
    male_density_half_mating = 2, wound_competitors = "all_adults"
  )
)

test_that("C++ rates match the R reference in host mode", {
  set.seed(30)
  for (nm in names(host_configs)) {
    m <- model_of(host_configs[[nm]], area = 3)
    setup <- SGDrivE:::.sg_reference_setup(m)
    view <- SGDrivE:::.sg_model_view(m)
    for (rep in 1:3) {
      x <- stats::runif(m$n_states, 0, 20)
      if (rep == 3) x[sg_state_index(m, "HW")] <- 0 # no uninfested wound
      r <- SGDrivE:::.sg_rates_r(5, x, setup, aux = TRUE)
      cc <- SGDrivE:::sg_rates_cpp(view, 5L, x)
      expect_equal(cc$deriv, r$deriv, tolerance = 1e-12, info = nm)
      expect_equal(cc$counters, unname(r$counters), tolerance = 1e-12, info = nm)
      agg <- SGDrivE:::sg_aggregates_cpp(view, 5L, x)
      expect_equal(agg$sites, r$aux$sites, tolerance = 1e-14, info = nm)
      expect_equal(agg$psi, r$aux$psi, tolerance = 1e-14, info = nm)
      expect_equal(agg$mu[2], r$aux$mu[["L"]], tolerance = 1e-14, info = nm)
      expect_equal(
        SGDrivE:::sg_max_rate_cpp(view, 5L, x), SGDrivE:::.sg_max_rate(5, x, setup),
        tolerance = 1e-12, info = nm
      )
    }
  }
})

test_that("C++ and R backends give the same host-mode simulation", {
  skip_if_not_installed("deSolve")
  for (nm in names(host_configs)) {
    p <- host_configs[[nm]]
    m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(25, days = 120), p), area = 10)
    x0 <- sg_init(m)
    s_r <- simulate_ode(m, x0, backend = "R", rtol = 1e-10, atol = 1e-12)
    s_c <- simulate_ode(m, x0, rtol = 1e-10, atol = 1e-12)
    expect_lt(max(abs(s_c$state - s_r$state)) / max(s_r$state), 1e-8)
    expect_lt(max(abs(s_c$counters - s_r$counters)) / max(abs(s_r$counters)), 1e-8)
    # deSolve driven by the C++ rates agrees too
    s_l <- simulate_ode(m, x0, method = "lsoda", rtol = 1e-10, atol = 1e-12)
    expect_lt(max(abs(s_l$state - s_r$state)) / max(s_r$state), 1e-6)
    tot <- rowSums(s_c$state[, grep("^H", colnames(s_c$state))])
    expect_lt(max(abs(tot - m$hosts$n)), 1e-8 * m$hosts$n)
  }
})

test_that("the C++ view refuses crowding per infested host without hosts", {
  v <- SGDrivE:::.sg_model_view(model_of(sg_params()))
  v$larval_per_infested_host <- 1L
  expect_error(SGDrivE:::sg_model_view_info_cpp(v), "needs the host model")
})

# ---- chain-binomial step (WP5) ------------------------------------------------

cb_configs <- list(
  batches = host_params(host_density = 5, wounded_fraction = 0.3, wound_search_rate = 0.05),
  full_eggs = host_params(
    host_density = 5, wounded_fraction = 0.3, wound_search_rate = 0.05,
    infested_wound_attractiveness = 3, treatment_rate = 0.1, infestation_substages = 3,
    larval_crowding_mode = "per_infested_host", larval_crowding_coef = 0.01,
    oviposition_mode = "eggs"
  ),
  clutch_two_hosts = host_params(
    host_density = 2, wounded_fraction = 0.6, wound_search_rate = 1, clutch_size = 200,
    infestation_substages = 1
  )
)
cb_model <- function(p, days = 20, area = 4) {
  sg_model(p, cube3, sg_forcing(sg_constant_weather(25, days = days), p), area = area)
}
cb_init <- function(m) sg_init(m, values = c(E = 5, L = 5, P = 5, M = 3, U = 1, F = 3))

test_that("C++ host steps equal the R reference draw for draw", {
  ns <- asNamespace("SGDrivE")
  for (nm in names(cb_configs)) {
    m <- cb_model(cb_configs[[nm]], area = if (nm == "clutch_two_hosts") 1 else 4)
    set.seed(3)
    a <- ns$.sg_cb_run_r(m, cb_init(m))
    set.seed(3)
    b <- ns$.sg_cb_run_cpp(m, cb_init(m))
    expect_identical(b$state, a$state, info = nm)
    expect_identical(b$counters, a$counters, info = nm)
    hosts <- grep("^H", colnames(b$state))
    expect_true(all(rowSums(b$state[, hosts]) == sum(b$state[1, hosts])), info = nm)
    expect_true(all(b$state == round(b$state)), info = nm)
    k <- b$counters
    expect_identical(
      k[, "infestations"] + k[, "superinfestations"], k[, "fertile_batches"],
      info = nm
    )
    da <- ns$.sg_cb_run_r(m, cb_init(m), stochastic = FALSE)
    db <- ns$.sg_cb_run_cpp(m, cb_init(m), stochastic = FALSE)
    expect_lt(max(abs(da$state - db$state)) / max(da$state), 1e-12)
    expect_lt(max(abs(da$counters - db$counters)) / max(abs(da$counters)), 1e-12)
  }
})

test_that("a step never infests more hosts than are wounded", {
  p <- cb_configs$clutch_two_hosts
  m <- cb_model(p, area = 1)
  x <- as.numeric(sg_init(m, values = c(E = 50, L = 50, P = 50, M = 30, F = 300)))
  x[sg_state_index(m, "HS")] <- 0
  x[sg_state_index(m, "HW")] <- 2
  x[sg_state_index(m, "HI")] <- 0
  x <- round(x)
  set.seed(5)
  out <- SGDrivE:::sg_cb_step_rep_cpp(SGDrivE:::.sg_model_view(m), 1L, 0.1, x, 2000L)
  n <- m$n_states
  flows <- out[, n + seq_along(SGDrivE:::.sg_counter_names(cube3$genotypes, TRUE))]
  colnames(flows) <- SGDrivE:::.sg_counter_names(cube3$genotypes, TRUE)
  expect_true(all(flows[, "infestations"] <= 2))
  expect_gt(mean(flows[, "infestations"] == 2), 0.05) # both hosts hit in some steps
  expect_gt(sum(flows[, "superinfestations"]), 0)
  hosts <- grep("^H", sg_state_names(m))
  expect_true(all(rowSums(out[, hosts]) == 2))
  expect_true(all(out[, sg_state_index(m, "HI")] <= 2))
})

test_that("the occupancy draw has the exact mean of random hits", {
  set.seed(9)
  occ <- SGDrivE:::.sg_cb_draws(TRUE)$occupancy
  k <- vapply(1:20000, function(i) occ(5, 3), 0)
  expected <- 3 * (1 - (2 / 3)^5)
  expect_lt(abs(mean(k) - expected), 4 * stats::sd(k) / sqrt(20000))
  expect_true(all(k >= 1 & k <= 3))
  expect_identical(occ(0, 3), 0)
  expect_identical(occ(4, 0), 0)
  det <- SGDrivE:::.sg_cb_draws(FALSE)$occupancy
  expect_equal(det(0.01, 3), 0.01, tolerance = 1e-2) # few batches: about one host each
  expect_lt(det(1e6, 3), 3 + 1e-12)
})

test_that("the deterministic host map converges to the ODE as dt shrinks", {
  ns <- asNamespace("SGDrivE")
  p <- cb_configs$full_eggs
  m <- cb_model(p, days = 15, area = 40)
  x0 <- cb_init(m)
  ode <- simulate_ode(m, x0, rtol = 1e-10, atol = 1e-12)
  hosts <- grep("^H", sg_state_names(m))
  err <- vapply(c(0.1, 0.02), function(dt) {
    d <- ns$.sg_cb_run_cpp(m, x0, dt = dt, stochastic = FALSE)
    max(abs(d$state[16, hosts] - ode$state[16, hosts])) / max(ode$state[16, hosts])
  }, 0)
  expect_lt(err[2], err[1])
  expect_lt(err[2], 0.01)
})

test_that("simulate_cb() runs host models", {
  p <- host_params(host_density = 20, wounded_fraction = 0.2, wound_search_rate = 0.02)
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(25, days = 60), p), area = 25)
  x0 <- sg_init(m, values = c(E = 20, L = 20, P = 20, M = 10, F = 10))
  sim <- simulate_cb(m, x0, n_rep = 40, seed = 2)
  hosts <- grep("^H", dimnames(sim$state)[[2]])
  tot <- apply(sim$state[, hosts, ], c(1, 3), sum)
  expect_true(all(tot == rep(tot[1, ], each = nrow(tot))))
  # replicate mean of infested hosts follows the ODE
  ode <- simulate_ode(m, x0)
  inf <- grep("^HI", dimnames(sim$state)[[2]])
  cb_inf <- mean(apply(sim$state[31:61, inf, ], 3, sum)) / 31
  ode_inf <- sum(ode$state[31:61, grep("^HI", colnames(ode$state))]) / 31
  expect_equal(cb_inf, ode_inf, tolerance = 0.15)
  totals <- simulate_cb(m, x0, n_rep = 2, seed = 2, keep = "totals")
  expect_true(all(c("HS[hosts]", "HW[hosts]", "HI[hosts]") %in% dimnames(totals$totals)[[2]]))
  expect_equal(totals$totals[, "HI[hosts]", ], apply(simulate_cb(m, x0, n_rep = 2, seed = 2)$state[
    , inf, ,
    drop = FALSE
  ], c(1, 3), sum))
})

test_that("flies can die out while the hosts go on", {
  p <- host_params(host_density = 10, wounded_fraction = 0.2)
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(10, days = 120), p), area = 1)
  x0 <- sg_init(m, values = c(E = 1, L = 1, P = 1, M = 1, F = 1))
  sim <- simulate_cb(m, x0, n_rep = 3, seed = 4)
  expect_true(all(!is.na(sim$extinct_day)))
  hosts <- grep("^H", dimnames(sim$state)[[2]])
  last <- dim(sim$state)[1]
  expect_true(all(apply(sim$state[last, hosts, , drop = FALSE], 3, sum) > 0))
  fly <- seq_len(min(hosts) - 1)
  expect_true(all(sim$state[last, fly, ] == 0))
})

# ---- equilibria, starts, scaling and calibration (WP6) --------------------------
# A realistic region: 1000 km^2 with 30 animals per km^2 (30,000 hosts).

region_params <- function(...) {
  base <- list(
    wound_model = "hosts", host_density = 30, wounded_fraction = 0.05,
    wound_search_rate = 0.02
  )
  do.call(sg_params, utils::modifyList(base, list(...)))
}
region_model <- function(p, temp = 25, days = 30, area = 1000) {
  sg_model(p, cube3, sg_forcing(sg_constant_weather(temp, days = days), p), area = area)
}

test_that("the host-mode equilibrium is a fixed point of the ODE", {
  p <- region_params(
    infested_wound_attractiveness = 3, treatment_rate = 0.05, infestation_substages = 4,
    larval_crowding_mode = "per_infested_host", larval_crowding_coef = 1e-4
  )
  m <- region_model(p)
  e <- sg_init(m, "equilibrium")
  info <- attr(e, "equilibrium")
  expect_gt(info$adult_females, 1000)
  hosts <- SGDrivE:::.sg_host_index(m)
  expect_equal(sum(e[hosts]), 30000, tolerance = 1e-10)
  sim <- simulate_ode(m, e, rtol = 1e-10, atol = 1e-8)
  rel <- max(abs(sweep(sim$state, 2, as.numeric(e))) / pmax(as.numeric(e), 1))
  expect_lt(rel, 1e-6)
  eq <- sg_equilibria(m)
  expect_identical(nrow(eq), 2L)
  expect_identical(eq$stable, c(FALSE, TRUE))
  expect_equal(eq$adult_females[2], info$adult_females)
  expect_equal(eq$infested_hosts[2], sum(e[sg_state_index(m, "HI")]), tolerance = 1e-10)
  expect_equal(eq$prevalence[2], eq$infested_hosts[2] / 30000)
})

test_that("equilibrium prevalence does not depend on infested-wound attractiveness", {
  # egg balance fixes the strike rate per wound, so alpha scales the flies only
  eq1 <- sg_equilibria(region_model(region_params()))
  eq10 <- sg_equilibria(region_model(region_params(infested_wound_attractiveness = 10)))
  expect_equal(eq10$prevalence[2], eq1$prevalence[2], tolerance = 1e-8)
  expect_gt(eq10$adult_females[2], 2 * eq1$adult_females[2])
})

test_that("with an Allee effect the hosts' region has two stable equilibria", {
  p <- region_params(male_density_half_mating = 20)
  m <- region_model(p, days = 365)
  eq <- sg_equilibria(m)
  expect_identical(eq$stable, c(TRUE, FALSE, TRUE)) # extinction, threshold, upper
  thr <- sg_init(m, "threshold")
  up <- sg_init(m, "equilibrium")
  expect_equal(attr(thr, "threshold")$adult_females, eq$adult_females[2])
  expect_equal(attr(up, "equilibrium")$adult_females, eq$adult_females[3])
  fly <- setdiff(seq_len(m$n_states), SGDrivE:::.sg_host_index(m))
  start <- function(k) {
    x <- as.numeric(thr)
    x[fly] <- k * x[fly]
    x
  }
  females <- function(k) tail(sg_females(simulate_ode(m, start(k)))$females, 1)
  expect_lt(females(0.9), 0.01 * eq$adult_females[2])
  expect_equal(females(1.1), eq$adult_females[3], tolerance = 1e-3)
})

test_that("several stable positive equilibria give a message and the upper one", {
  m <- region_model(region_params())
  real <- SGDrivE:::.sg_equilibria_core(m, 25, "test")
  fake <- real
  fake$table <- real$table[c(1, 2, 2, 2), ]
  fake$table$eggs_per_day <- c(0, 1, 2, 3) * real$table$eggs_per_day[2]
  fake$table$adult_females <- c(0, 10, 20, 30)
  fake$table$stable <- c(FALSE, TRUE, FALSE, TRUE)
  fake$states <- real$states[c(1, 2, 2, 2)]
  local_mocked_bindings(.sg_equilibria_core = function(...) fake)
  expect_message(x <- sg_init(m, "equilibrium"), "2 stable equilibria at 25 degC")
  expect_equal(attr(x, "equilibrium")$adult_females, 30)
  expect_equal(suppressMessages(attr(sg_init(m, "threshold"), "threshold"))$adult_females, 20)
})

test_that("adult_females sizes a host region exactly", {
  p <- region_params(treatment_rate = 0.02)
  m <- sg_model(p, cube3, sg_forcing(sg_constant_weather(25, days = 30), p),
    adult_females = 50000
  )
  e <- sg_init(m, "equilibrium")
  expect_equal(attr(e, "equilibrium")$adult_females, 50000, tolerance = 1e-8)
  expect_equal(m$hosts$n, 30 * m$area)
})

test_that("R0 is unchanged by hosts and lowered by treatment", {
  m_static <- region_model(sg_params())
  m_hosts <- region_model(region_params())
  expect_equal(sg_r0(m_hosts, c(20, 25)), sg_r0(m_static, c(20, 25)))
  m_treat <- region_model(region_params(treatment_rate = 0.3))
  expect_true(all(sg_r0(m_treat, c(20, 25)) < sg_r0(m_hosts, c(20, 25))))
})

test_that("a burn-in runs the hosts with the flies", {
  m <- region_model(region_params(), days = 365)
  b <- sg_init(m, "burn_in")
  expect_lt(attr(b, "burn_in")$start_effect, 1e-3)
  expect_equal(sum(b[SGDrivE:::.sg_host_index(m)]), 30000, tolerance = 1e-8)
  e <- sg_init(m, "equilibrium")
  expect_lt(max(abs(b - e) / pmax(e, 1)), 1e-3) # a year at 25 degC reaches equilibrium
})

test_that("calibrate_population() calibrates each population-size parameter", {
  m_static <- region_model(sg_params(), days = 365)
  a <- calibrate_population(m_static, target = 50)
  b <- calibrate_wound_density(m_static, target = 50)
  expect_identical(a$parameter, "wound_density")
  expect_equal(a$wound_density, b$wound_density)
  m <- region_model(region_params(), days = 365)
  wf <- calibrate_population(m, target = 50)
  expect_identical(wf$parameter, "wounded_fraction")
  expect_identical(wf$method, "root search")
  expect_equal(wf$achieved, 50, tolerance = 1e-6)
  expect_true(wf$value > 0 && wf$value < 1)
  expect_identical(wf$params$host_density, 30)
  hd <- calibrate_population(m, target = 50, parameter = "host_density")
  expect_identical(hd$method, "proportional")
  expect_identical(hd$evaluations, 2L)
  expect_equal(hd$achieved, 50, tolerance = 1e-6)
  expect_equal(hd$model$hosts$n, hd$value * 1000)
  expect_true(any(grepl("host_density = .* animals per km\\^2", capture.output(print(hd)))))
  expect_error(calibrate_population(m, 50, parameter = "wound_density"), "is not used")
  expect_error(calibrate_population(m_static, 50, parameter = "host_density"), "is not used")
  expect_error(calibrate_wound_density(m, 50), "use calibrate_population")
  expect_error(calibrate_population(m, 1e9), "No wounded fraction")
})
