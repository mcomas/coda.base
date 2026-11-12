test_that("package defaults preserve existing user options", {
  old <- options(coda.base.basis = TRUE)
  on.exit(options(old), add = TRUE)

  coda.base:::.setOptions()
  expect_identical(getOption("coda.base.basis"), TRUE)

  options(coda.base.basis = NULL)
  coda.base:::.setOptions()
  expect_identical(getOption("coda.base.basis"), FALSE)
})
