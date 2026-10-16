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

test_that("the stochastic solver and dynamic starts refuse host models until implemented", {
  m <- model_of(host_params())
  expect_error(simulate_cb(m, n_rep = 1, seed = 1), "not available yet for wound_model = \"hosts\"")
  expect_error(sg_init(m, "equilibrium"), "not available yet")
  expect_error(sg_init(m, "threshold"), "not available yet")
  expect_error(check_model(m), "not available yet")
  expect_error(sg_r0(m, 20), "not available yet")
  p <- host_params()
  expect_error(
    sg_model(p, cube3, sg_forcing(sg_constant_weather(20, days = 10), p), adult_females = 100),
    "not available yet"
  )
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
