# ggplot2 plots of simulations (D-026) and the data behind them.

ns <- asNamespace("SGDrivE")
p <- sg_params()
wx <- sg_forcing(sg_example_weather(), p, "1971-05-01", "1971-06-29")
mod <- sg_model(p, sg_cube(sg_test_cube_mendelian()), wx,
  sg_releases("aa", "adult", "male", times = c(5, 20), number = 200),
  area = 50
)
init <- sg_init(mod, values = c(E = 20, L = 20, P = 20, M = 10, F = 10), per = "km2")
sim <- simulate_ode(mod, init = init)
cb <- simulate_cb(mod, init, n_rep = 3, seed = 1)
whats <- c("stages", "genotypes", "alleles", "females", "cases")

test_that("plot data: one row per date, series (and replicate)", {
  d <- sg_plot_data(sim, "stages")
  expect_identical(names(d), c("date", "series", "value"))
  expect_identical(levels(d$series), gsub("_", " ", ns$.sg_stage_names))
  expect_identical(nrow(d), 61L * 6L)
  expect_match(attr(d, "ylab"), "per km\\^2")
  g <- sg_plot_data(sim, "genotypes", per = "node")
  tot <- sg_totals(sim, by = "genotype", stages = "adults")
  expect_equal(g$value, tot$value)
  a <- sg_plot_data(cb, "alleles", replicates = 2:3)
  expect_identical(sort(unique(a$replicate)), 2:3)
  expect_identical(attr(a, "ylab"), "allele frequency")
  cs <- sg_plot_data(sim, "cases", period = "month", host_density = 10)
  expect_identical(unique(as.character(cs$series)), "cases per host")
  expect_equal(cs$value, sg_cases(sim, period = "month", host_density = 10)$cases_per_host)
  f <- sg_plot_data(simulate_cb(mod, init, n_rep = 1, seed = 2, keep = "totals"), "females")
  expect_false(any(c("egg laying", "fertile equivalent") %in% f$series)) # not available
})

test_that("every plot type builds for deterministic and stochastic runs", {
  skip_if_not_installed("ggplot2")
  for (s in list(sim, cb)) {
    for (w in whats) {
      g <- plot(s, w)
      expect_s3_class(g, "ggplot")
      built <- ggplot2::ggplot_build(g)
      expect_identical(nrow(built$data[[1]]), nrow(sg_plot_data(s, w)))
    }
  }
  g <- plot(cb, "genotypes", stages = "eggs", per = "node", replicates = 1)
  expect_identical(g$labels$subtitle, "1 replicate")
  expect_match(g$labels$title, "eggs")
})

# ColorBrewer palettes (as given by RColorBrewer, used by ggplot2 via scales)
dark2 <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E", "#E6AB02", "#A6761D", "#666666")
paired <- c(
  "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99", "#E31A1C",
  "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A", "#FFFF99", "#B15928"
)

test_that("colour-blind-friendly palettes, with ColorBrewer Paired for many series", {
  skip_if_not_installed("ggplot2")
  colours <- function(g) unique(ggplot2::ggplot_build(g)$data[[1]]$colour)
  oi <- colours(plot(sim, "stages")) # six stages
  expect_true(all(toupper(oi) %in% ns$.sg_okabe_ito))
  expect_length(oi, 6)
  br <- colours(plot(sim, "stages", palette = "brewer"))
  expect_true(all(toupper(br) %in% dark2))
  vi <- colours(plot(sim, "stages", palette = "viridis"))
  expect_length(vi, 6)
  expect_false(any(toupper(vi) %in% ns$.sg_okabe_ito))
  skip_if_not_installed("MGDrivE")
  # ten genotypes: more than eight colours, so ColorBrewer Paired is used
  mh <- sg_model(p, sg_cube(MGDrivE::cubeHomingDrive()), wx)
  equal <- stats::setNames(rep(1, 10), mh$cube$genotypes)
  sh <- simulate_ode(mh, init = sg_init(mh, genotypes = equal))
  for (pal in c("okabe_ito", "brewer")) {
    many <- colours(plot(sh, "genotypes", palette = pal))
    expect_length(many, 10)
    expect_true(all(toupper(many) %in% paired))
  }
  expect_error(plot(sim, palette = "rainbow"), "must be one of")
})

test_that("plots need ggplot2 and valid arguments", {
  local_mocked_bindings(.sg_has_ggplot2 = function() FALSE)
  expect_error(plot(sim), "install.packages\\(\"ggplot2\"\\)")
  expect_error(sg_plot_data(sim, "ages"), "must be one of")
  expect_error(sg_plot_data(list()), "simulation from")
})
