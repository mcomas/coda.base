library(testthat)
library(coda.base)
test_that("basis_matrix", {
  X2 = data.frame(matrix(exp(rnorm(20)), ncol = 2))

  ALR = alr_basis(2)
  ILR = ilr_basis(2)
  SBP = sbp_basis(X1~X2, data = X2)
  PC = pc_basis(X2)
  PB = pb_basis(X2, method = 'exact')

  expect_that(ALR, is_a('matrix'))
  expect_that(ILR, is_a('matrix'))
  expect_that(SBP, is_a('matrix'))
  expect_that(PC, is_a('matrix'))
  expect_that(PB, is_a('matrix'))
})

test_that("coordinates_default", {
  X = matrix(exp(rnorm(200)), ncol = 10)

  ALR = coordinates(X, 'alr')
  expect_that(ALR, is_a('matrix'))

  ILR = coordinates(X, 'ilr')
  expect_that(ILR, is_a('matrix'))

  CLR = coordinates(X, 'clr')
  expect_that(CLR, is_a('matrix'))

  PC = coordinates(X, 'pc')
  expect_that(PC, is_a('matrix'))

  PB = coordinates(X, 'pb')
  expect_that(PB, is_a('matrix'))

  CDP = coordinates(X, 'cdp')
  expect_that(CDP, is_a('matrix'))

})
