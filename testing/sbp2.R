D = 10
X = matrix(exp(rnorm(D*3)), ncol = D)

sbp2_test_1(X)
var(coordinates(X, pb_basis(X, method = 'exact')[,1, drop=FALSE]))

# library(microbenchmark)
#
# microbenchmark(
#   sbp2_test_1(X),
#   pb_basis(X, method = 'exact'))

