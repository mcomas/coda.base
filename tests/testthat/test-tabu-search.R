test_that("partial tabu search records default neighbourhood", {
  set.seed(1)
  X <- matrix(rexp(100), ncol = 5)
  lI <- lapply(seq_len(ncol(X)), identity)

  res <- partial_pb_tabu_search(X, lI, iter = 10, tabu_size = 5)

  expect_true(is.matrix(res$balance))
  expect_named(
    res$neighbourhoods,
    c(
      "remove_active",
      "add_left",
      "add_right",
      "flip_side",
      "swap_zero",
      "swap_sides"
    )
  )
})

test_that("partial constrained search matches first constrained balance", {
  set.seed(9)
  X <- matrix(rexp(120), ncol = 6)

  partial <- partial_pb_constrained(X, constrained.criterion = "variance")
  constrained <- pb_basis(
    X,
    method = "constrained",
    constrained.criterion = "variance",
    ordering = FALSE
  )[, 1, drop = FALSE]

  expect_equal(
    abs(as.numeric(crossprod(partial$balance, constrained))),
    1,
    tolerance = 1e-8
  )
  expect_equal(
    partial$variance,
    as.numeric(var(coordinates(X, constrained))),
    tolerance = 1e-8
  )
  expect_equal(partial$constrained.criterion, "variance")
})

test_that("partial constrained search supports the angle criterion", {
  set.seed(14)
  X <- matrix(rexp(160), ncol = 8)
  lI <- list(c(1, 2), c(3, 4), 5, 6, c(7, 8))

  partial <- partial_pb_constrained(
    X,
    lI = lI,
    constrained.criterion = "angle"
  )

  expect_length(partial$balance_raw, length(lI))
  expect_true(any(partial$balance_raw < 0L))
  expect_true(any(partial$balance_raw > 0L))
  expect_true(is.finite(partial$variance))
  expect_equal(partial$constrained.criterion, "angle")
})

test_that("partial exact search agrees with exact principal balance without restriction", {
  set.seed(10)
  X <- matrix(rexp(120), ncol = 6)

  partial <- partial_pb_exact(X)
  exact <- pb_basis(X, method = "exact2", ordering = FALSE)[, 1, drop = FALSE]

  expect_equal(abs(as.numeric(crossprod(partial$balance, exact))), 1, tolerance = 1e-8)
  expect_equal(partial$variance, as.numeric(var(coordinates(X, exact))), tolerance = 1e-8)
})

test_that("partial exact search respects active-group restrictions", {
  set.seed(11)
  X <- matrix(rexp(100), ncol = 5)

  exact <- partial_pb_exact(X, min_parts = 3, max_parts = 3)

  expect_equal(sum(exact$balance_raw != 0L), 3)
  expect_equal(exact$min_parts, 3)
  expect_equal(exact$max_parts, 3)
  expect_true(any(exact$balance_raw < 0L))
  expect_true(any(exact$balance_raw > 0L))
})

test_that("partial exact search works on grouped parts", {
  set.seed(12)
  X <- matrix(rexp(120), ncol = 6)
  lI <- list(c(1, 2), c(3, 4), 5, 6)

  exact <- partial_pb_exact(X, lI = lI, min_parts = 2, max_parts = 2)

  expect_length(exact$balance_raw, length(lI))
  expect_equal(sum(exact$balance_raw != 0L), 2)
  expect_equal(nrow(exact$balance), length(lI))
})

test_that("partial exact search records the selected method", {
  set.seed(13)
  X <- matrix(rexp(140), ncol = 7)
  lI <- list(c(1, 2), 3, 4, c(5, 6), 7)

  exact <- partial_pb_exact(
    X,
    lI = lI,
    min_parts = 3,
    max_parts = 4,
    method = "restricted"
  )

  expect_equal(exact$method, "restricted")
  expect_gte(sum(exact$balance_raw != 0L), 3)
  expect_lte(sum(exact$balance_raw != 0L), 4)
})

test_that("partial tabu search neighbourhood accepts extra neighbourhoods", {
  set.seed(2)
  X <- matrix(rexp(120), ncol = 6)
  lI <- lapply(seq_len(ncol(X)), identity)

  res <- partial_pb_tabu_search(
    X,
    lI,
    iter = 10,
    tabu_size = 6,
    flip_side = TRUE,
    swap_zero = TRUE,
    swap_sides = TRUE
  )

  expect_length(res$balance_raw, ncol(X))
  expect_true(any(res$balance_raw < 0L))
  expect_true(any(res$balance_raw > 0L))
  expect_true(all(res$neighbourhoods))
})

test_that("partial tabu search records constrained initialisation criterion", {
  set.seed(15)
  X <- matrix(rexp(120), ncol = 6)
  lI <- lapply(seq_len(ncol(X)), identity)

  res <- partial_pb_tabu_search(
    X,
    lI,
    iter = 10,
    tabu_size = length(lI),
    constrained.criterion = "angle"
  )

  expect_equal(res$constrained.criterion, "angle")
})

test_that("partial tabu search accepts NULL grouping", {
  set.seed(16)
  X <- matrix(rexp(120), ncol = 6)

  res <- partial_pb_tabu_search(
    X,
    iter = 10
  )

  expect_length(res$lI, ncol(X))
  expect_equal(res$lI, lapply(seq_len(ncol(X)), identity))
  expect_length(res$balance_raw, ncol(X))
})

test_that("partial tabu search neighbourhood requires an active neighbourhood", {
  set.seed(3)
  X <- matrix(rexp(80), ncol = 4)
  lI <- lapply(seq_len(ncol(X)), identity)

  expect_error(
    partial_pb_tabu_search(
      X,
      lI,
      iter = 5,
      tabu_size = 4,
      remove_active = FALSE,
      add_left = FALSE,
      add_right = FALSE
    ),
    "At least one neighbourhood type must be active"
  )
})

test_that("partial tabu search neighbourhood limits active groups", {
  set.seed(4)
  X <- matrix(rexp(160), ncol = 8)
  lI <- list(c(1, 2, 3), c(4, 5), 6, 7, 8)

  res <- partial_pb_tabu_search(
    X,
    lI,
    iter = 20,
    tabu_size = length(lI),
    max_parts = 2,
    add_left = TRUE,
    add_right = TRUE,
    remove_active = TRUE,
    swap_zero = TRUE
  )

  expect_lte(sum(res$balance_raw != 0L), 2)
  expect_equal(res$max_parts, 2)
})

test_that("partial tabu search respects lower active-group restriction", {
  set.seed(7)
  X <- matrix(rexp(160), ncol = 8)
  lI <- lapply(seq_len(ncol(X)), identity)

  res <- partial_pb_tabu_search(
    X,
    lI,
    min_parts = 5,
    max_parts = 5,
    iter = 10,
    tabu_size = length(lI),
    remove_active = TRUE,
    add_left = FALSE,
    add_right = FALSE,
    flip_side = TRUE,
    swap_zero = TRUE,
    swap_sides = TRUE
  )

  expect_equal(sum(res$balance_raw != 0L), 5)
  expect_equal(res$min_parts, 5)
  expect_equal(res$max_parts, 5)
})

test_that("automatic initialisation starts on the max_parts boundary", {
  set.seed(6)
  X <- matrix(rexp(200), ncol = 10)
  lI <- lapply(seq_len(ncol(X)), identity)

  res <- partial_pb_tabu_search(
    X,
    lI,
    iter = 1,
    tabu_size = length(lI),
    max_parts = 5,
    remove_active = FALSE,
    add_left = FALSE,
    add_right = FALSE,
    flip_side = TRUE,
    swap_zero = TRUE,
    swap_sides = TRUE
  )

  expect_equal(sum(res$balance_raw != 0L), 5)
})

test_that("partial tabu search neighbourhood rejects oversized user initialisation", {
  set.seed(5)
  X <- matrix(rexp(120), ncol = 6)
  lI <- list(c(1, 2), 3, 4, 5, 6)

  expect_error(
    partial_pb_tabu_search(
      X,
      lI,
      iter = 10,
      tabu_size = length(lI),
      ini = c(-1, -1, -1, 1, 1),
      max_parts = 4
    ),
    "more active groups than max_parts"
  )
})

test_that("partial tabu search rejects undersized user initialisation", {
  set.seed(8)
  X <- matrix(rexp(120), ncol = 6)
  lI <- lapply(seq_len(ncol(X)), identity)

  expect_error(
    partial_pb_tabu_search(
      X,
      lI,
      min_parts = 4,
      max_parts = 6,
      iter = 10,
      tabu_size = length(lI),
      ini = c(-1, 1, 0, 0, 0, 0)
    ),
    "fewer active groups than min_parts"
  )
})
