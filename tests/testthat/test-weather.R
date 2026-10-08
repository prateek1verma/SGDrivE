test_that("example weather is the complete McAllen 1960-1983 series", {
  wx <- sg_example_weather()
  expect_identical(nrow(wx), 8766L)
  expect_named(wx, c("date", "tmax", "tmin", "tmean", "rain", "source"))
  expect_s3_class(wx$date, "Date")
  expect_identical(range(wx$date), as.Date(c("1960-01-01", "1983-12-31")))
  n_src <- table(wx$source)[c("primary", "backup", "interpolated")]
  expect_identical(as.vector(n_src), c(8612L, 152L, 2L))
  expect_false(anyNA(wx[, c("tmax", "tmin", "tmean", "rain")]))
})

# ---- milestone check: cold index reproduces the published values ------------------

test_that("McAllen cold index reproduces Gutierrez et al. values (1971 = 6.67, 1979 = 9.80)", {
  ci <- sg_cold_index(sg_example_weather())
  val <- stats::setNames(ci$mu_cold, ci$year)
  expect_equal(unname(round(val["1971"], 2)), 6.67)
  expect_equal(unname(round(val["1979"], 2)), 9.80)
  # further values quoted in Gutierrez et al. (2019)
  expect_equal(unname(round(val[c("1961", "1976", "1977", "1978", "1980", "1981")], 1)),
    c(8.3, 10.4, 9.0, 10.3, 9.3, 7.4),
    tolerance = 0.11
  )
  expect_true(all(ci$complete))
})

test_that("cold index works on forcing objects and with the July-June period", {
  wx <- sg_example_weather()
  f <- sg_forcing(wx)
  expect_identical(sg_cold_index(f), sg_cold_index(wx))
  jj <- sg_cold_index(f, period = "jul_jun")
  expect_identical(jj$year, 1960:1984)
  expect_false(jj$complete[1]) # Jan-Jun 1960 only
  expect_false(jj$complete[nrow(jj)]) # Jul-Dec 1983 only
  expect_true(all(jj$complete[2:24]))
  d <- as.integer(format(f$dates, "%Y")) + (as.integer(format(f$dates, "%m")) >= 7)
  expect_equal(jj$mu_cold[jj$year == 1972], sum(f$p_A[d == 1972 & f$tmean < 27.2]))
})

# ---- reading and gap filling ---------------------------------------------------

pri <- sg_read_ghcn(test_path("fixtures", "ghcn_USC00415701_subset.csv"))
bak <- sg_read_ghcn(test_path("fixtures", "ghcn_USW00012959_subset.csv"))

test_that("GHCN files are read and converted from tenths", {
  expect_named(pri, c("date", "tmax", "tmin", "rain"))
  expect_identical(attr(pri, "station"), "USC00415701")
  expect_identical(pri$date[1], as.Date("1960-08-15"))
  expect_identical(pri$tmax[1], 33.9)
  expect_identical(pri$tmin[1], 23.3)
  expect_error(sg_read_ghcn("does-not-exist.csv"), "existing GHCN-Daily")
  tmp <- tempfile(fileext = ".csv")
  utils::write.csv(data.frame(a = 1), tmp, row.names = FALSE)
  expect_error(sg_read_ghcn(tmp), "does not look like a GHCN-Daily")
})

test_that("gap filling reproduces the bundled McAllen series (backup and interpolation)", {
  wx <- sg_example_weather()
  for (win in list(c("1960-08-15", "1960-10-31"), c("1975-03-15", "1975-05-31"))) {
    filled <- sg_fill_weather(pri, bak, start = win[1], end = win[2])
    ref <- wx[wx$date >= as.Date(win[1]) & wx$date <= as.Date(win[2]), ]
    expect_identical(filled$date, ref$date)
    expect_identical(filled$source, ref$source)
    for (v in c("tmax", "tmin", "rain", "tmean")) {
      expect_equal(filled[[v]], ref[[v]], tolerance = sg_tol_exact)
    }
  }
  # the windows exercise both backup filling and interpolation
  f75 <- sg_fill_weather(pri, bak, start = "1975-03-15", end = "1975-05-31")
  f60 <- sg_fill_weather(pri, bak, start = "1960-08-15", end = "1960-10-31")
  expect_true(any(f75$source == "backup"))
  expect_true(any(f60$source == "interpolated"))
})

test_that("gap filling steps behave as documented", {
  w <- data.frame(
    date = as.Date("2000-01-01") + c(0:2, 4:12),
    tmax = c(20, NA, 22, 23, NA, NA, 26, 27, 28, 29, 30, 31),
    tmin = 10, rain = 0
  )
  f <- sg_fill_weather(w)
  expect_identical(nrow(f), 13L) # missing day 2000-01-04 added
  expect_identical(f$tmax[2], 21) # isolated gap: neighbour mean
  expect_identical(f$source[2], "interpolated")
  expect_identical(f$tmax[4], 22.5) # added day, isolated: neighbour mean
  # days 6 and 7 form a two-day gap: mean of the surrounding +-5 days,
  # filled forwards (day 7 uses the new value of day 6)
  d6 <- mean(c(20, 21, 22, 22.5, 23, 26, 27, 28, 29))
  expect_equal(f$tmax[6], d6, tolerance = sg_tol_exact)
  expect_equal(f$tmax[7], mean(c(21, 22, 22.5, 23, d6, 26, 27, 28, 29, 30)),
    tolerance = sg_tol_exact
  )
  expect_false(anyNA(f$tmax))
  expect_error(sg_fill_weather(w[, 1:3]), "columns date, tmax, tmin and rain")
  expect_error(sg_fill_weather(w, start = "2000-02-01", end = "2000-01-01"), "start <= end")
  expect_error(sg_fill_weather(rbind(w, w)), "without NA or repeats")
})
