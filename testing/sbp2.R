D = 5
n = 100
Y = matrix(exp(rnorm(D*n)), ncol = D)
x = rnorm(n)
sbp_basis(find_predictive_balance(Y,x,method=1), silent = TRUE)
sbp_basis(find_predictive_balance(Y,x,method=2), silent = TRUE)

microbenchmark(
  sbp_basis(find_principal_balance(Y), silent = TRUE),
  sbp_basis(find_predictive_balance(Y,x,method=1), silent = TRUE),
  sbp_basis(find_predictive_balance(Y,x,method=2), silent = TRUE))

pb_basis(X, method = 'exact')[,1, drop=FALSE]

x = matrix(rnorm(n), ncol=1)
sbp_basis(sbp2_test_2(X, x), silent = TRUE)

var(coordinates(X, pb_basis(X, method = 'exact')[,1, drop=FALSE]))

library(microbenchmark)

microbenchmark(
  pb_basis(X, method = 'exact'),
  sbp2_test_1(X))

