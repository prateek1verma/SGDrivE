test_that("parameter table has the documented columns and unique names", {
  tab <- sg_param_table()
  expect_s3_class(tab, "data.frame")
  expect_named(tab, c(
    "name", "symbol", "group", "type", "default", "lower", "upper",
    "bounds", "choices", "allow_na", "units", "source", "spec", "description"
  ))
  expect_false(anyDuplicated(tab$name) > 0)
  expect_true(all(tab$group %in% SGDrivE:::.sg_param_groups))
  expect_true(all(tab$type %in% c("real", "integer", "choice")))
  expect_true(all(tab$bounds %in% c("[]", "(]", "[)", "()")))
  expect_true(all(tab$source %in% c("GP2014", "G2019", "MGDrivE2", "adapted")))
})

test_that("every parameter has units, a spec equation and a description", {
  tab <- sg_param_table()
  expect_true(all(nzchar(tab$units)))
  expect_true(all(grepl("^S[0-9]+(, S[0-9]+)*$", tab$spec)))
  expect_true(all(nchar(tab$description) > 10))
})

test_that("choice parameters list their choices and the default is one of them", {
  tab <- sg_param_table()
  ch <- tab[tab$type == "choice", ]
  expect_true(all(!is.na(ch$choices)))
  defaults <- SGDrivE:::.sg_param_defaults()
  for (i in seq_len(nrow(ch))) {
    expect_true(defaults[[ch$name[i]]] %in% strsplit(ch$choices[i], "|", fixed = TRUE)[[1]])
  }
})

test_that("table text is ASCII only (portability)", {
  tab <- sg_param_table()
  txt <- unlist(tab[vapply(tab, is.character, TRUE)])
  txt <- txt[!is.na(txt)]
  expect_true(all(!grepl("[^ -~]", txt)))
})

test_that("default sub-stage numbers do not exceed 10 (D-003)", {
  d <- SGDrivE:::.sg_param_defaults()
  expect_true(all(unlist(d[c("n_E", "n_L", "n_P", "n_A")]) <= 10L))
})

test_that("key defaults match the model specification", {
  d <- SGDrivE:::.sg_param_defaults()
  expect_identical(d$nu, log(2))
  expect_identical(d$h_M, 0)
  expect_identical(d$search_denominator, "searching_females")
  expect_identical(d$X_pre, 61)
  expect_identical(d$X_A, 202.8)
  expect_identical(d$R, 90)
})
