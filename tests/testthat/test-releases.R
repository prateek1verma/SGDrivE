test_that("adult male releases on day numbers", {
  rel <- sg_releases("S", stage = "adult", sex = "male", times = c(30, 37, 44), number = 1000)
  expect_s3_class(rel, "sg_releases")
  ev <- as.data.frame(rel)
  expect_identical(nrow(ev), 3L)
  expect_identical(ev$time, c(30, 37, 44))
  expect_identical(unique(ev$time_type), "day")
  expect_identical(unique(ev$sex), "male")
  expect_identical(unique(ev$stage), "adult")
  expect_identical(ev$number, rep(1000, 3))
  expect_identical(unique(ev$per), "node")
})

test_that("egg and pupa releases have no sex", {
  ev <- as.data.frame(sg_releases("S", stage = "pupa", times = 5, number = 10))
  expect_true(is.na(ev$sex))
  expect_identical(ev$stage, "pupa")
  eggs <- as.data.frame(sg_releases("S", stage = "egg", times = 1, number = 1))
  expect_identical(eggs$stage, "egg")
  expect_error(
    sg_releases("S", stage = "pupa", sex = "male", times = 5, number = 10),
    "applies only to adult"
  )
})

test_that("'both' splits each release into males and females (D-013)", {
  ev <- as.data.frame(sg_releases("S", "adult", sex = "both", times = c(1, 8), number = 100))
  expect_identical(nrow(ev), 4L)
  expect_identical(ev$sex, c("female", "male", "female", "male"))
  expect_identical(ev$number, c(50, 50, 50, 50))
  ev <- as.data.frame(sg_releases("S", "adult", sex = "both", times = 1, number = 100, split = 0.2))
  expect_identical(ev$number[ev$sex == "female"], 20)
  expect_identical(ev$number[ev$sex == "male"], 80)
})

test_that("number may vary by event; per can be km2", {
  rel <- sg_releases("S", "adult", "male", times = 1:3, number = c(1, 2, 3), per = "km2")
  ev <- as.data.frame(rel)
  expect_identical(ev$number, c(1, 2, 3))
  expect_identical(unique(ev$per), "km2")
  expect_error(sg_releases("S", "adult", "male", times = 1:3, number = 1:2), "length 1 or the same")
})

test_that("dates are accepted as Date or text", {
  a <- sg_releases("S", "adult", "male", times = as.Date("1971-06-01") + c(0, 7), number = 5)
  b <- sg_releases("S", "adult", "male", times = c("1971-06-01", "1971-06-08"), number = 5)
  expect_identical(a, b)
  ev <- as.data.frame(a)
  expect_s3_class(ev$time, "Date")
  expect_identical(unique(ev$time_type), "date")
})

test_that("schedules are combined with c() and sorted by time", {
  r1 <- sg_releases("S", "adult", "male", times = c(10, 30), number = 1)
  r2 <- sg_releases("S", "pupa", times = 20, number = 2)
  r <- c(r1, r2, NULL)
  ev <- as.data.frame(r)
  expect_identical(ev$time, c(10, 20, 30))
  expect_identical(ev$stage, c("adult", "pupa", "adult"))
  d <- sg_releases("S", "adult", "male", times = "1971-01-01", number = 1)
  expect_error(c(r1, d), "day numbers with releases given as dates")
  expect_error(c(r1, list()), "Only sg_releases")
})

test_that("invalid arguments give informative errors", {
  expect_error(sg_releases(c("S", "W"), "adult", "male", times = 1, number = 1), "single genotype")
  expect_error(sg_releases("S", "larva", times = 1, number = 1), "must be one of.*\"egg\"")
  expect_error(sg_releases("S", "adult", times = 1, number = 1), "need `sex`")
  expect_error(sg_releases("S", "adult", "males", times = 1, number = 1), "Did you mean 'male'")
  expect_error(sg_releases("S", "adult", "male", times = 0, number = 1), "whole numbers >= 1")
  expect_error(sg_releases("S", "adult", "male", times = 1.5, number = 1), "whole numbers >= 1")
  expect_error(sg_releases("S", "adult", "male", times = numeric(0), number = 1), "at least one")
  expect_error(sg_releases("S", "adult", "male", times = "1971-13-01", number = 1), "valid dates")
  expect_error(sg_releases("S", "adult", "male", times = 1, number = -1), ">= 0")
  expect_error(
    sg_releases("S", "adult", "male", times = 1, number = 1, per = "ha"),
    "must be one of"
  )
  expect_error(sg_releases("S", "adult", "both", times = 1, number = 1, split = 2), "split")
  expect_warning(sg_releases("S", "adult", "male", times = 1, number = 1, split = 0.3), "ignored")
})

test_that("print() summarises the schedule", {
  rel <- c(
    sg_releases("S", "adult", "male",
      times = seq(30, by = 7, length.out = 26),
      number = 1000, per = "km2"
    ),
    sg_releases("S", "pupa", times = 100, number = 50)
  )
  expect_snapshot(print(rel))
  both <- sg_releases("S", "adult", "both", times = c("1971-06-01", "1971-06-15"), number = 500)
  expect_snapshot(print(both))
})

# ---- resolving releases for a model (used by sg_model()) -------------------------

dates <- as.Date("1971-01-01") + 0:364

test_that("day-number releases resolve to days and counts per node", {
  rel <- sg_releases("S", "adult", "male", times = c(1, 365), number = 10, per = "km2")
  out <- ScrewDrivE:::.sg_resolve_releases(rel, dates, area = 2.5, genotypes = c("W", "S"))
  expect_identical(out$day, c(1L, 365L))
  expect_identical(out$count, c(25, 25))
  expect_named(out, c("day", "genotype", "stage", "sex", "count"))
})

test_that("date releases resolve relative to the first forcing date (D-015)", {
  rel <- sg_releases("S", "adult", "female", times = c("1971-01-01", "1971-02-01"), number = 3)
  out <- ScrewDrivE:::.sg_resolve_releases(rel, dates, area = 1, genotypes = c("W", "S"))
  expect_identical(out$day, c(1L, 32L))
  expect_identical(out$count, c(3, 3)) # per node: area not applied
})

test_that("resolution checks genotypes and the simulation period", {
  rel <- sg_releases("Z", "adult", "male", times = 1, number = 1)
  expect_error(ScrewDrivE:::.sg_resolve_releases(rel, dates, 1, c("W", "S")), "not in the cube: Z")
  late <- sg_releases("S", "adult", "male", times = 400, number = 1)
  expect_error(
    ScrewDrivE:::.sg_resolve_releases(late, dates, 1, c("W", "S")),
    "outside the simulation period.*day 400"
  )
  early <- sg_releases("S", "adult", "male", times = "1970-12-31", number = 1)
  expect_error(ScrewDrivE:::.sg_resolve_releases(early, dates, 1, c("W", "S")), "1970-12-31")
  expect_identical(nrow(ScrewDrivE:::.sg_resolve_releases(NULL, dates, 1, "W")), 0L)
  expect_error(ScrewDrivE:::.sg_resolve_releases(list(), dates, 1, "W"), "sg_releases object")
})
