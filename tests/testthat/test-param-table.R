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
  expect_true(all(grepl("^[SH][0-9]+(, [SH][0-9]+)*$", tab$spec)))
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
  stage_n <- d[c(
    "egg_substages", "larva_substages", "pupa_substages", "adult_substages",
    "infestation_substages"
  )]
  expect_true(all(unlist(stage_n) <= 10L))
})

test_that("key defaults match the model specification", {
  d <- SGDrivE:::.sg_param_defaults()
  expect_identical(d$mating_rate, log(2))
  expect_identical(d$male_density_half_mating, 0)
  expect_identical(d$wound_competitors, "searching_females")
  expect_identical(d$dd_scale_preoviposition, 61)
  expect_identical(d$adult_lifespan_dd, 202.8)
  expect_identical(d$eggs_per_day_max, 90)
})

test_that("host-wound parameters are in the table (D-030)", {
  tab <- sg_param_table()
  host <- tab[tab$group == "hosts", ]
  expect_identical(host$name, c(
    "wound_model", "host_density", "wounded_fraction", "wound_healing_days",
    "infestation_days", "infestation_substages", "infested_wound_attractiveness",
    "treatment_rate"
  ))
  expect_identical(host$default[host$name == "wound_model"], "\"static\"")
  expect_true(all(host$allow_na[host$name %in% c("host_density", "wounded_fraction")]))
  d <- SGDrivE:::.sg_param_defaults()
  expect_true(is.na(d$host_density) && is.na(d$wounded_fraction))
  expect_identical(d$infestation_substages, 6L)
  expect_match(tab$choices[tab$name == "larval_crowding_mode"], "per_infested_host")
})
