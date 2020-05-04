D = 50
n = 60
X = matrix(exp(rnorm(D*n)), ncol = D)
microbenchmark::microbenchmark(
  MASS::ginv(cov(coordinates(X, 'clr'))),
  ilr_basis(D) %*% chol2inv(chol(cov(coordinates(X)))) %*% t(ilr_basis(D)))
