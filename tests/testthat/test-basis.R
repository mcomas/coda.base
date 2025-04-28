library(testthat)
library(coda.base)
test_that("basis_matrix", {

  D = 5
  N = D * 10
  X = data.frame(matrix(exp(rnorm(N)), ncol = D))

  ALR = alr_basis(D)
  ILR = ilr_basis(D)
  CLR = clr_basis(D)
  SBP = sbp_basis(c(b1=X1~X2+X3+X4+X5,
                    b2=X2~X3+X4+X5,
                    b3=X3~X4+X5,
                    b4=X4~X5), data = X)
  PC = pc_basis(X)
  PB = pb_basis(X, method = 'exact')

  expect_that(ALR, is_a('matrix'))
  expect_that(ILR, is_a('matrix'))
  expect_that(CLR, is_a('matrix'))
  expect_that(SBP, is_a('matrix'))
  expect_that(PC, is_a('matrix'))
  expect_that(PB, is_a('matrix'))

  expect_equal(sum(composition(coordinates(X, ALR), ALR) / X * rowSums(X)), N)
  expect_equal(sum(composition(coordinates(X, ILR), ILR) / X * rowSums(X)), N)
  expect_equal(sum(composition(coordinates(X, CLR), CLR) / X * rowSums(X)), N)
  expect_equal(sum(composition(coordinates(X, SBP), SBP) / X * rowSums(X)), N)
  expect_equal(sum(composition(coordinates(X, PC), PC) / X * rowSums(X)), N)
  expect_equal(sum(composition(coordinates(X, PB), PB) / X * rowSums(X)), N)

})

