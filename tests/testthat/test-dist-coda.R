library(testthat)
library(coda.base)

test_that("dist_coda computes Aitchison distance like dist", {
  X <- exp(matrix(seq_len(20), nrow = 5))

  d1 <- dist_coda(X, method = "aitchison")
  d2 <- dist_coda(X, method = "L2")

  expect_equal(as.vector(d1), as.vector(suppressWarnings(dist(X, method = "aitchison"))))
  expect_equal(attr(d1, "method"), "aitchison")
  expect_equal(as.vector(d2), as.vector(suppressWarnings(dist(X, method = "L2"))))
  expect_equal(attr(d2, "method"), "L2")
  expect_equal(as.vector(d1), as.vector(d2))
})

test_that("dist warns when using the legacy Aitchison extension", {
  X <- exp(matrix(seq_len(20), nrow = 5))

  expect_warning(
    d1 <- dist(X, method = "aitchison"),
    "Use dist_coda"
  )
  expect_equal(attr(d1, "method"), "aitchison")

  expect_warning(
    d2 <- dist(X, method = "L2"),
    "Use dist_coda"
  )
  expect_equal(attr(d2, "method"), "L2")
})

test_that("dist_coda computes L1 distances on expected log-ratio coordinates", {
  X <- exp(matrix(c(
    1, 2, 4,
    3, 5, 6,
    7, 8, 10,
    11, 13, 14
  ), ncol = 3, byrow = TRUE))

  expect_equal(
    as.vector(dist_coda(X, method = "L1-clr")),
    as.vector(stats::dist(coordinates(X, "clr"), method = "manhattan"))
  )

  expect_equal(
    as.vector(dist_coda(X, method = "L1-pw")),
    as.vector(stats::dist(coordinates(X, "pw"), method = "manhattan")) / (ncol(X) - 1)
  )
})

test_that("dist_coda L1 is scale invariant and validates positive data", {
  X <- exp(matrix(c(
    1, 2, 4,
    3, 5, 6,
    7, 8, 10
  ), ncol = 3, byrow = TRUE))
  scaled_X <- X * c(2, 10, 100)

  expect_equal(
    as.vector(dist_coda(X, method = "L1")),
    as.vector(dist_coda(scaled_X, method = "L1")),
    tolerance = 1e-12
  )

  expect_error(dist_coda(matrix(c(1, 0, 2, 3), nrow = 2), method = "L1"), "positive")
})
