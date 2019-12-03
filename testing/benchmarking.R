library(microbenchmark)
library(coda.base)
library(compositions)
library(robCompositions)
library(easyCODA)
library(RcppCoDA)

K = 10
N = 100
# ALR default transformation
mX = matrix(exp(rnorm(K*N)), nrow=N, ncol=K)
dX = as.data.frame(mX)
aX = acomp(mX)
B = alr_basis(K)
alr_results = microbenchmark(
  coda.base::coordinates(mX, 'alr'),
  coda.base::coordinates2(mX, 'alr'),
  compositions::alr(aX),
  robCompositions::addLR(mX),
  easyCODA::ALR(mX),
  RcppCoDA::alr(mX),
  times = 100
)
autoplot(alr_results)


clr_results = microbenchmark(
  coda.base::coordinates(mX, 'clr'),
  coda.base::coordinates2(mX, 'clr'),
  compositions::clr(aX),
  RcppCoDA::clr(mX),
  robCompositions::cenLR(mX),
  easyCODA::CLR(mX), times = 100
)
autoplot(clr_results)

ilr_results = microbenchmark(
  coda.base::coordinates(mX),
  coda.base::coordinates2(mX),
  compositions::ilr(aX),
  RcppCoDA::ilr(mX),
  # robCompositions::isomLR(mX, fast=T),
  robCompositions::pivotCoord(mX, fast=T),
  easyCODA::PLR(mX), times = 100
)
autoplot(ilr_results)

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
