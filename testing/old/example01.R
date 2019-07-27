set.seed(1)
library(Rcpp)
library(coda.base)
sourceCpp('src/principal_balances.cpp')
X=as.data.frame(matrix(exp(rnorm(10*100)), nrow=100, ncol=10))
M = cov(log(X))

## This example has the lowest variance precisly in the balance separating all parts.
##
B1 = find_PB(M, 10)
var(coordinates(X, B1)[,1])

B2 = find_PB2(M, 100, .Machine$integer.max)
var(coordinates(X, B2)[,1])

B3 = find_PB2(M, 100, 10)
var(coordinates(X, B3)[,1])
