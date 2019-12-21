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
