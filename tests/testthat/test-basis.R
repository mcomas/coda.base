library(testthat)
library(coda.base)

test_that("basis constructors return matrices with expected dimensions", {
  D <- 5
  N <- D * 10
  X <- data.frame(matrix(exp(rnorm(N)), ncol = D))
  colnames(X) <- paste0("X", 1:D)

  ALR <- alr_basis(D)
  ILR <- ilr_basis(D)
  OLR <- olr_basis(D)
  CLR <- clr_basis(D)
  CDP <- cdp_basis(D)
  SBP <- sbp_basis(c(
    b1 = X1 ~ X2 + X3 + X4 + X5,
    b2 = X2 ~ X3 + X4 + X5,
    b3 = X3 ~ X4 + X5,
    b4 = X4 ~ X5
  ), data = X)
  PC <- pc_basis(X)
  PB <- pb_basis(X, method = "exact")
  PW <- pairwise_basis(D)

  expect_true(is.matrix(ALR))
  expect_true(is.matrix(ILR))
  expect_true(is.matrix(OLR))
  expect_true(is.matrix(CLR))
  expect_true(is.matrix(CDP))
  expect_true(is.matrix(SBP))
  expect_true(is.matrix(PC))
  expect_true(is.matrix(PB))

  expect_equal(dim(ALR), c(D, D - 1))
  expect_equal(dim(ILR), c(D, D - 1))
  expect_equal(dim(OLR), c(D, D - 1))
  expect_equal(dim(CLR), c(D, D))
  expect_equal(dim(CDP), c(D, D - 1))
  expect_equal(dim(SBP), c(D, D - 1))
  expect_equal(dim(PC), c(D, D - 1))
  expect_equal(dim(PB), c(D, D - 1))

  expect_equal(nrow(PW), D)
  expect_equal(ncol(PW), choose(D, 2))
})

test_that("matrix bases reconstruct the closed composition", {
  D <- 5
  N <- 40
  X <- matrix(exp(rnorm(N * D)), ncol = D)
  colnames(X) <- paste0("c", 1:D)
  Xc <- X / rowSums(X)

  bases <- list(
    alr_basis(D),
    ilr_basis(D),
    olr_basis(D),
    clr_basis(D),
    cdp_basis(D),
    sbp_basis(cbind(
      c( 1,  1,  1, -1,  -1),
      c(-1,  1,  1,  0,  0),
      c( 0,  0,  0,  1, -1),
      c( 0,  1, -1,  0,  0)
    )),
    pc_basis(X),
    pb_basis(X, method = "exact")
  )

  for (B in bases) {
    H <- coordinates(X, B)
    Xrec <- composition(H, B)
    expect_equal(Xrec, Xc, tolerance = 1e-8)
  }
})

test_that("exact2 principal balances agree with exact principal balances", {
  for (D in 3:7) {
    set.seed(100 + D)
    X <- matrix(exp(rnorm(30 * D)), ncol = D)

    B1 <- pb_basis(X, method = "exact", ordering = FALSE)
    B2 <- pb_basis(X, method = "exact2", ordering = FALSE)

    expect_equal(unname(abs(crossprod(B1, B2))), diag(D - 1), tolerance = 1e-8)

    v1 <- apply(coordinates(X, B1), 2, stats::var)
    v2 <- apply(coordinates(X, B2), 2, stats::var)
    expect_equal(unname(v2), unname(v1), tolerance = 1e-8)
  }
})


test_that("character bases agree with their matrix counterparts when applicable", {
  D <- 5
  N <- 30
  X <- matrix(exp(rnorm(N * D)), ncol = D)
  colnames(X) <- paste0("X", 1:D)

  expect_equal(coordinates(X, "alr"), coordinates(X, alr_basis(D)), tolerance = 1e-10)
  expect_equal(coordinates(X, "clr"), coordinates(X, clr_basis(D)), tolerance = 1e-10)
  expect_equal(coordinates(X, "ilr"), coordinates(X, ilr_basis(D)), tolerance = 1e-10)
  expect_equal(coordinates(X, "olr"), coordinates(X, olr_basis(D)), tolerance = 1e-10)
  expect_equal(coordinates(X, "cdp"), coordinates(X, cdp_basis(D)), tolerance = 1e-10)

  expect_equal(unname(composition(coordinates(X, "alr"), "alr")),
               unname(X / rowSums(X)), tolerance = 1e-8)
  expect_equal(unname(composition(coordinates(X, "clr"), "clr")),
               unname(X / rowSums(X)), tolerance = 1e-8)
  expect_equal(unname(composition(coordinates(X, "ilr"), "ilr")),
               unname(X / rowSums(X)), tolerance = 1e-8)
  expect_equal(unname(composition(coordinates(X, "olr"), "olr")),
               unname(X / rowSums(X)), tolerance = 1e-8)
})

test_that("vector, matrix and data.frame inputs are handled consistently", {
  x <- c(1, 2, 3, 4, 5)
  X <- rbind(x, 2 * x)
  Xdf <- as.data.frame(X)
  colnames(Xdf) <- paste0("X", 1:5)

  hx_vec <- coordinates(x, "ilr")
  hx_mat <- coordinates(matrix(x, nrow = 1), "ilr")[1, ]
  expect_equal(hx_vec, hx_mat, tolerance = 1e-10)

  h_df <- coordinates(Xdf, "clr")
  h_mat <- coordinates(as.matrix(Xdf), "clr")
  expect_equal(as.matrix(h_df), h_mat, tolerance = 1e-10)

  x_rec <- composition(hx_vec, "ilr")
  expect_equal(as.numeric(x_rec), x / sum(x), tolerance = 1e-8)
})

test_that("coord and comp aliases behave as expected", {
  X <- matrix(exp(rnorm(20)), ncol = 5)
  x <- X[1, ]

  expect_equal(coord(X, basis = "ilr"), coordinates(X, "ilr"), tolerance = 1e-10)
  expect_equal(comp(coordinates(X, "ilr"), "ilr"), composition(coordinates(X, "ilr"), "ilr"),
               tolerance = 1e-10)
})

test_that("column names of coordinates follow the basis", {
  D <- 5
  X <- matrix(exp(rnorm(20)), ncol = D)
  colnames(X) <- paste0("X", 1:D)

  B <- ilr_basis(D)
  colnames(B) <- paste0("custom", 1:ncol(B))

  H <- coordinates(X, B)
  expect_equal(colnames(H), colnames(B))

  H2 <- coordinates(X, "ilr")
  expect_equal(colnames(H2), paste0("ilr", 1:ncol(H2)))
})

test_that("rownames are preserved when possible", {
  X <- matrix(exp(rnorm(20)), ncol = 5)
  rownames(X) <- paste0("r", 1:nrow(X))

  H <- coordinates(X, "clr")
  Xrec <- composition(H, "clr")

  expect_equal(rownames(H), rownames(X))
  expect_equal(rownames(Xrec), rownames(X))
})

test_that("pairwise basis returns expected number of coordinates", {
  D <- 6
  X <- matrix(exp(rnorm(5 * D)), ncol = D)

  H <- coordinates(X, "pw")
  expect_equal(ncol(H), choose(D, 2))
})

test_that("invalid inputs raise errors", {
  X_bad <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(coordinates(X_bad, "ilr"), "numeric")

  H_bad <- data.frame(a = 1:3, b = letters[1:3])
  expect_error(composition(H_bad, "ilr"), "numeric")

  X <- matrix(exp(rnorm(20)), ncol = 5)
  expect_error(composition(coordinates(X, "ilr"), "pb"),
               "arg.*ilr.*olr.*alr.*clr")

})
