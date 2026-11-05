# Landscapes of several nodes (MODEL_SPEC section 8a, D-032).

land_forcing <- function(temp = 25, days = 10, p = sg_params()) {
  sg_forcing(sg_constant_weather(temp, days = days), p)
}

test_that("sg_landscape fills in the node table and shares one forcing", {
  wx <- land_forcing()
  land <- sg_landscape(data.frame(id = c("a", "b"), area = c(10, 20)), wx)
  expect_s3_class(land, "sg_landscape")
  expect_identical(land$n_nodes, 2L)
  nd <- land$nodes
  expect_named(nd, c(
    "id", "area", "x", "y", "role", "wound_density", "host_density",
    "forcing_series"
  ))
  expect_identical(nd$id, c("a", "b"))
  expect_identical(nd$area, c(10, 20))
  expect_identical(nd$role, c("study", "study"))
  expect_true(all(is.na(nd$x)) && all(is.na(nd$wound_density)))
  expect_length(land$forcing, 1L)
  expect_identical(nd$forcing_series, c(1L, 1L))
  expect_identical(as.data.frame(land), nd)
})

test_that("identical forcings are stored once, with an index per node", {
  warm <- land_forcing(27)
  cool <- land_forcing(20)
  warm2 <- land_forcing(27) # built separately, same content
  land <- sg_landscape(
    data.frame(id = c("s", "m", "n", "x"), area = 1),
    list(s = warm, m = warm2, n = cool, x = warm)
  )
  expect_length(land$forcing, 2L)
  expect_identical(land$nodes$forcing_series, c(1L, 1L, 2L, 1L))
  expect_identical(land$forcing[[2L]]$tmean, cool$tmean)
  # a list in node order without names
  land2 <- sg_landscape(data.frame(id = c("s", "n"), area = 1), list(cool, warm))
  expect_identical(land2$nodes$forcing_series, c(1L, 2L))
  # names are matched to ids, whatever their order
  land3 <- sg_landscape(data.frame(id = c("s", "n"), area = 1), list(n = warm, s = cool))
  expect_identical(land3$forcing[[land3$nodes$forcing_series[1L]]]$tmean, cool$tmean)
})

test_that("optional node columns are checked and kept", {
  nd <- data.frame(
    id = c("a", "b"), area = c(5, 7), x = c(0, 10), y = c(0, 0),
    role = c("source", "study"), wound_density = c(NA, 3), host_density = c(40, NA)
  )
  land <- sg_landscape(nd, land_forcing())
  expect_identical(land$nodes$role, c("source", "study"))
  expect_identical(land$nodes$wound_density, c(NA, 3))
  expect_identical(land$nodes$host_density, c(40, NA))
  expect_identical(land$nodes$x, c(0, 10))
})

test_that("bad node tables give clear errors", {
  wx <- land_forcing()
  expect_error(sg_landscape(list(id = "a", area = 1), wx), "must be a data frame")
  expect_error(
    sg_landscape(data.frame(id = character(0), area = numeric(0)), wx),
    "at least one row"
  )
  expect_error(sg_landscape(data.frame(id = "a"), wx), "lacks column")
  expect_error(sg_landscape(data.frame(id = "a", area = 1, aera = 2), wx), "Unknown column")
  expect_error(sg_landscape(data.frame(id = c("a", "a"), area = 1), wx), "unique; repeated: a")
  expect_error(sg_landscape(data.frame(id = "", area = 1), wx), "non-empty")
  expect_error(sg_landscape(data.frame(id = "a@1", area = 1), wx), "must not contain")
  expect_error(sg_landscape(data.frame(id = "a", area = 0), wx), "area` must be > 0")
  expect_error(sg_landscape(data.frame(id = "a", area = NA_real_), wx), "area` must be > 0")
  expect_error(
    sg_landscape(data.frame(id = "a", area = 1, wound_density = -1), wx),
    "wound_density` must be NA or >= 0"
  )
  expect_error(
    sg_landscape(data.frame(id = "a", area = 1, role = "sink"), wx),
    "\"study\" or \"source\""
  )
  expect_error(
    sg_landscape(data.frame(id = "a", area = 1, role = "source"), wx),
    "at least one study"
  )
  expect_error(sg_landscape(data.frame(id = "a", area = 1, x = "far"), wx), "x` must be numeric")
})

test_that("bad forcings give clear errors", {
  nd <- data.frame(id = c("a", "b"), area = 1)
  wx <- land_forcing()
  expect_error(sg_landscape(nd, sg_constant_weather(25, days = 10)), "sg_forcing object")
  expect_error(sg_landscape(nd, list(wx)), "1 elements but the landscape has 2")
  expect_error(sg_landscape(nd, list(a = wx, c = wx)), "unknown: c")
  expect_error(sg_landscape(nd, list(wx, land_forcing(days = 11))), "same dates")
  expect_error(
    sg_landscape(nd, list(wx, land_forcing(p = sg_params(dev_temp_min = 14)))),
    "same parameters"
  )
})

test_that("print shows nodes and forcing series", {
  land <- sg_landscape(
    data.frame(id = c("s", "n"), area = c(1, 2), role = c("source", "study")),
    list(land_forcing(27), land_forcing(20))
  )
  out <- capture.output(print(land))
  expect_match(out[1L], "2 nodes \\(1 study, 1 source\\), total area 3 km\\^2")
  expect_match(out[2L], "2 distinct series, 10 days")
})
