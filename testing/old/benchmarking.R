library(microbenchmark)
library(coda.base)
library(compositions)
library(robCompositions)
library(easyCODA)

K = 10
X = as.data.frame(matrix(exp(rnorm(K*100)), nrow=100, ncol=K))
aX = acomp(X)
microbenchmark(
  coordinates(X, 'alr'),
  alr(aX),
  addLR(X),
  ALR(X), times = 1000
)

microbenchmark(
  coordinates(X, 'clr'),
  clr(aX),
  cenLR(X),
  CLR(X), times = 1000
)

microbenchmark(
  coordinates(X, 'ilr'),
  ilr(aX),
  isomLR(X, fast=T),
  pivotCoord(X, fast=T),
  PLR(X), times = 1000
)

B = sbp_basis(b1 = V1~V2,
              b2 = b1~V3,
              b3 = b2~V4,
              b4 = b3~V5,
              b5 = b4~V6,
              b6 = b5~V7,
              b7 = b6~V8,
              b8 = b7~V9,
              b9 = b8~V10, data = X)
microbenchmark(
  coordinates(X, B),
  ilr(aX, B),
  times = 1000
)

microbenchmark(
  variation_array(X, only_variation = TRUE),
  compositions::variation(aX),
  robCompositions::variation(X, robust=FALSE),times = 1000)
