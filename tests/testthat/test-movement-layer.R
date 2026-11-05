# Movement layers and their merge into groups (MODEL_SPEC section 8a, D-032).

rates3 <- list(male = 0.05, virgin = 0.1, mated = 0.2)
chain3 <- rbind(c(0, 1, 0), c(0.5, 0, 0.5), c(0, 1, 0))
rates_open <- list(male = 0.05, virgin = 0.1, mated = c(0.2, 0, 0.3))

test_that("a matrix dest becomes a sorted edge table", {
  l <- sg_movement_layer(rates3, chain3)
  expect_s3_class(l, "sg_movement_layer")
  expect_identical(l$clock, "physiological")
  expect_identical(l$n_nodes, 3L)
  expect_identical(l$edges$from, c(1L, 2L, 2L, 3L))
  expect_identical(l$edges$to, c(2L, 1L, 3L, 2L))
  expect_identical(l$edges$prob, c(1, 0.5, 0.5, 1))
  expect_identical(l$rate, list(male = 0.05, virgin = 0.1, mated = 0.2))
  # dimnames give node ids
  m <- chain3
  dimnames(m) <- list(c("a", "b", "c"), c("a", "b", "c"))
  l2 <- sg_movement_layer(c(male = 0.05, virgin = 0.1, mated = 0.2), m, clock = "calendar")
  expect_identical(l2$edges$from, c("a", "b", "b", "c"))
  expect_identical(l2$clock, "calendar")
})

test_that("an edge table dest is checked and sorted", {
  e <- data.frame(from = c("b", "a", "b"), to = c("c", "b", "a"), prob = c(0.3, 1, 0.5))
  l <- sg_movement_layer(rates3, e)
  expect_identical(l$edges$from, c("a", "b", "b"))
  expect_identical(l$edges$to, c("b", "a", "c"))
  expect_identical(l$edges$prob, c(1, 0.5, 0.3))
  expect_true(is.na(l$n_nodes))
  # zero-probability edges are dropped
  e0 <- rbind(e, data.frame(from = "c", to = "a", prob = 0))
  expect_identical(nrow(sg_movement_layer(rates3, e0)$edges), 3L)
})

test_that("bad rates give clear errors", {
  expect_error(sg_movement_layer(list(male = 1, virgin = 1), chain3), "missing: mated")
  expect_error(
    sg_movement_layer(list(male = 1, virgin = 1, mated = 1, maled = 1), chain3),
    "unknown: maled"
  )
  expect_error(sg_movement_layer(list(1, 1, 1), chain3), "named list")
  expect_error(
    sg_movement_layer(list(male = -1, virgin = 1, mated = 1), chain3),
    "male` must be finite numbers >= 0"
  )
  expect_error(sg_movement_layer(list(male = NA, virgin = 1, mated = 1), chain3), "finite numbers")
  expect_error(
    sg_movement_layer(list(male = c(1, 2), virgin = 1, mated = 1), chain3),
    "one per node \\(3\\); got 2"
  )
  expect_silent(sg_movement_layer(list(male = c(1, 2, 3), virgin = 1, mated = 1), chain3))
})

test_that("bad destinations give clear errors", {
  expect_error(sg_movement_layer(rates3, matrix(0, 2, 3)), "square matrix")
  expect_error(sg_movement_layer(rates3, diag(2)), "diagonal of `dest` must be 0")
  expect_error(sg_movement_layer(rates3, rbind(c(0, 1.2), c(1, 0))), "\\[0, 1\\]")
  expect_error(
    sg_movement_layer(rates3, rbind(c(0, 0.6, 0.6), c(1, 0, 0), c(1, 0, 0))),
    "sum to 1.2"
  )
  m <- chain3
  dimnames(m) <- list(c("a", "b", "c"), c("a", "c", "b"))
  expect_error(sg_movement_layer(rates3, m), "same node ids in the same order")
  expect_error(sg_movement_layer(rates3, data.frame(from = "a", to = "b")), "missing: prob")
  expect_error(sg_movement_layer(rates3, data.frame(from = "a", to = "a", prob = 1)), "to itself")
  expect_error(
    sg_movement_layer(rates3, data.frame(from = c("a", "a"), to = c("b", "b"), prob = 0.5)),
    "more than once: a -> b"
  )
  expect_error(
    sg_movement_layer(rates3, data.frame(from = c("a", "a"), to = c("b", "c"), prob = 0.6)),
    "sum to 1.2"
  )
  expect_error(sg_movement_layer(rates3, data.frame(from = 1.5, to = 2, prob = 1)), "whole numbers")
  expect_error(sg_movement_layer(rates3, "chain"), "numeric matrix or a data frame")
  expect_error(sg_movement_layer(rates3, chain3, clock = "wind"), "\"physiological\", \"calendar\"")
})

test_that("one layer gives one group with its own probabilities and CSR edges", {
  l <- sg_movement_layer(rates_open, chain3 * c(1, 1, 0.6))
  g <- .sg_movement_groups(list(l), c("s", "m", "n"))
  expect_length(g, 1L)
  g <- g[[1L]]
  expect_identical(g$clock, "physiological")
  expect_identical(unname(g$rate["mated", ]), c(0.2, 0, 0.3))
  expect_identical(g$row_ptr, c(1L, 2L, 4L, 5L))
  expect_identical(g$col, c(2L, 1L, 3L, 2L))
  # probabilities exactly as given (no rescaling) where the class moves
  expect_identical(unname(g$prob["male", ]), c(1, 0.5, 0.5, 0.6))
  # no movement of mated females from node 2: their probabilities are 0
  expect_identical(unname(g$prob["mated", ]), c(1, 0, 0, 0.6))
  # node 3 keeps 60%: 40% of its leavers leave the landscape
  expect_equal(unname(g$loss["male", ]), c(0, 0, 0.4), tolerance = 1e-15)
  expect_identical(unname(g$loss["mated", 2L]), 0)
})

test_that("closed rows that sum to 1 only up to round-off have no loss", {
  third <- 1 / 3
  m <- rbind(c(0, third, third, third), c(1, 0, 0, 0), c(1, 0, 0, 0), c(1, 0, 0, 0))
  g <- .sg_movement_groups(list(sg_movement_layer(rates3, m)), letters[1:4])[[1L]]
  expect_identical(unname(g$loss), matrix(0, 3, 4))
})

test_that("layers with the same clock are merged by rate-weighted destinations", {
  ids <- c("s", "m", "n")
  l1 <- sg_movement_layer(rates_open, chain3 * c(1, 1, 0.6))
  l2 <- sg_movement_layer(
    list(male = 0.01, virgin = 0.01, mated = 0.01),
    data.frame(from = c("s", "n"), to = c("n", "s"), prob = c(1, 0.5))
  )
  g <- .sg_movement_groups(list(l1, l2), ids)
  expect_length(g, 1L)
  g <- g[[1L]]
  expect_equal(unname(g$rate["mated", ]), c(0.21, 0.01, 0.31))
  # edges are the union: s->m, s->n, m->s, m->n, n->s, n->m
  expect_identical(g$col, c(2L, 3L, 1L, 3L, 1L, 2L))
  # node n, mated: 0.3 * 0.6 to m and 0.01 * 0.5 to s, out of 0.31
  expect_equal(unname(g$prob["mated", 5:6]), c(0.005, 0.18) / 0.31)
  expect_equal(unname(g$loss["mated", 3L]), 1 - 0.185 / 0.31)
  # node m, mated: only layer 2 moves them, and it has no edges from m: all lost
  expect_equal(unname(g$loss["mated", 2L]), 1)
  # flows by rate are kept: sum over layers of rate x prob equals rate^G x prob^G
  expect_equal(0.21 * unname(g$prob["mated", 2L]), 0.01 * 1)
})

test_that("layers with different clocks stay separate groups", {
  ids <- c("s", "m", "n")
  l1 <- sg_movement_layer(rates3, chain3)
  l2 <- sg_movement_layer(rates3, chain3, clock = "calendar")
  g <- .sg_movement_groups(list(l1, l2), ids)
  expect_identical(vapply(g, `[[`, "", "clock"), c("physiological", "calendar"))
})

test_that("layers are checked against the landscape's nodes", {
  ids <- c("a", "b")
  expect_error(
    .sg_movement_groups(list(sg_movement_layer(rates3, chain3)), ids),
    "has 3 nodes but the landscape has 2"
  )
  l <- sg_movement_layer(rates3, data.frame(from = "a", to = "z", prob = 1))
  expect_error(.sg_movement_groups(list(l), ids), "unknown node \"z\"")
  l <- sg_movement_layer(rates3, data.frame(from = 1, to = 3, prob = 1))
  expect_error(.sg_movement_groups(list(l), ids), "refers to node 3")
  l <- sg_movement_layer(
    list(male = c(1, 2, 3), virgin = 1, mated = 1),
    data.frame(from = "a", to = "b", prob = 1)
  )
  expect_error(.sg_movement_groups(list(l), ids), "one per node \\(2\\); got 3")
})

test_that("print summarises the layer", {
  l <- sg_movement_layer(rates_open, chain3 * c(1, 1, 0.6))
  out <- capture.output(print(l))
  expect_match(out[1L], "clock \"physiological\"; 4 edges from 3 nodes \\(1 with open edges\\)")
  expect_match(out[2L], "male 0.05, virgin 0.1, mated 0-0.3")
})
