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

test_that("solvers and dynamic starts refuse host models until implemented", {
  m <- model_of(host_params())
  expect_error(simulate_ode(m), "not available yet for wound_model = \"hosts\"")
  expect_error(simulate_cb(m, n_rep = 1, seed = 1), "not available yet")
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
