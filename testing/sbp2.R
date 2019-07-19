D = 10
X = matrix(exp(rnorm(D*3)), ncol = D)

sbp_basis(sbp2_test_1(X), silent = TRUE)
pb_basis(X, method = 'exact')[,1, drop=FALSE]

var(coordinates(X, pb_basis(X, method = 'exact')[,1, drop=FALSE]))

library(microbenchmark)

microbenchmark(
  pb_basis(X, method = 'exact'),
  sbp2_test_1(X))

