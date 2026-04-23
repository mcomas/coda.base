library(testthat)
library(coda.base)

test_that("coda_replacement imputes zeros and missings and preserves matrix metadata", {
  X <- matrix(c(
    0.00, 0.30, 0.70,
    0.20,   NA, 0.80,
    0.40, 0.60, 0.00,
    0.25, 0.25, 0.50,
    0.10, 0.30, 0.60
  ), ncol = 3, byrow = TRUE)
  colnames(X) <- c("sand", "silt", "clay")
  rownames(X) <- paste0("s", seq_len(nrow(X)))

  out <- coda_replacement(X, DL = c(0.05, 0.05, 0.05), maxit = 20)

  expect_true(is.matrix(out))
  expect_equal(dim(out), dim(X))
  expect_equal(dimnames(out), dimnames(X))
  expect_true(all(is.finite(out)))
  expect_true(all(out > 0))
  expect_equal(unname(rowSums(out)), rep(1, nrow(X)), tolerance = 1e-8)
})

test_that("coda_replacement preserves data.frame subclasses when possible", {
  X <- data.frame(
    a = c(0.00, 0.20, 0.40, 0.25, 0.10),
    b = c(0.30, NA, 0.60, 0.25, 0.30),
    c = c(0.70, 0.80, 0.00, 0.50, 0.60),
    check.names = FALSE
  )
  class(X) <- c("tbl_df", "tbl", "data.frame")

  out <- coda_replacement(X, DL = c(0.05, 0.05, 0.05), maxit = 20)

  expect_equal(class(out), class(X))
  expect_equal(names(out), names(X))
  expect_equal(unname(rowSums(as.data.frame(out))), rep(1, nrow(X)), tolerance = 1e-8)
})

test_that("coda_replacement validates X and DL before calling C++", {
  expect_error(
    coda_replacement(data.frame(a = 1:3, b = letters[1:3])),
    "numeric"
  )

  expect_error(
    coda_replacement(matrix(c(0, 0, 1, 2, 3, 4), nrow = 2)),
    "positive observed value"
  )

  X <- matrix(c(0, 2, 3, 1, 0, 4, 5, 6, 0), nrow = 3)

  expect_error(coda_replacement(X, DL = c(1, 2)), "length ncol")
  expect_error(coda_replacement(X, DL = c(1, 0, 2)), "strictly positive")
  expect_error(coda_replacement(X, DL = c("a", "b", "c")), "numeric")
})

test_that("coda_replacement validates control parameters", {
  X <- matrix(c(0, 2, 3, 1, 0, 4, 5, 6, 0), nrow = 3)
  DL <- c(1, 1, 1)

  expect_error(coda_replacement(X, DL = DL, dl_prop = 0), "dl_prop")
  expect_error(coda_replacement(X, DL = DL, dl_prop = 1), "dl_prop")
  expect_error(coda_replacement(X, DL = DL, eps = 0), "eps")
  expect_error(coda_replacement(X, DL = DL, maxit = 0), "maxit")
  expect_error(coda_replacement(X, DL = DL, parameters = NA), "parameters")
  expect_error(coda_replacement(X, DL = DL, debug = NA), "debug")
})

test_that("random_composition_with_fixed_pb validates sd1", {
  expect_error(
    random_composition_with_fixed_pb(c(-1, 1, 0), n = 5, sd1 = NA),
    "sd1"
  )
  expect_error(
    random_composition_with_fixed_pb(c(-1, 1, 0), n = 5, sd1 = 0),
    "sd1"
  )
})
