test_that("perturbation works for vectors and matrices", {
  expect_equal(
    perturbation(c(1, 2, 3), c(1, 1, 2)),
    closure(c(1, 2, 6))
  )

  X <- matrix(
    c(1, 2, 3,
      4, 5, 6),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(c("r1", "r2"), c("a", "b", "c"))
  )

  out <- perturbation(X, c(1, 1, 2))

  expect_equal(out, closure(X * matrix(c(1, 1, 2), nrow = 2, ncol = 3, byrow = TRUE)))
  expect_equal(dimnames(out), dimnames(X))

  out_reverse <- perturbation(c(a = 1, b = 1, c = 2), X)

  expect_equal(out_reverse, closure(X * matrix(c(1, 1, 2), nrow = 2, ncol = 3, byrow = TRUE)))
  expect_equal(dimnames(out_reverse), dimnames(X))

  Y <- as.data.frame(X)
  out_df <- perturbation(c(a = 1, b = 1, c = 2), Y)

  expect_s3_class(out_df, "data.frame")
  expect_equal(out_df, as.data.frame(out_reverse))
  expect_equal(row.names(out_df), row.names(Y))
})

test_that("powering works for vectors and row-wise exponents", {
  expect_equal(
    powering(c(1, 2, 3), 2),
    closure(c(1, 4, 9))
  )

  expect_equal(
    powering(c(a = 1, b = 2, c = 3), c(1, 2)),
    closure(rbind(
      c(a = 1, b = 2, c = 3),
      c(a = 1, b = 4, c = 9)
    ))
  )

  X <- matrix(
    c(1, 2, 3,
      4, 5, 6),
    nrow = 2,
    byrow = TRUE
  )

  expect_equal(
    powering(X, c(1, 2)),
    closure(rbind(X[1, ], X[2, ] ^ 2))
  )
})
