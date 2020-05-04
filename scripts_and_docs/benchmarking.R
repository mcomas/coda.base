library(microbenchmark)
library(coda.base)
library(compositions)
library(robCompositions)
library(easyCODA)
library(RcppCoDA)

K = 1000
N = 100
# ALR default transformation
mX = matrix(exp(rnorm(K*N)), nrow=N, ncol=K)
dX = as.data.frame(mX)
aX = acomp(mX)
B = alr_basis(K)
alr_results = microbenchmark(
  coda.base::coordinates(mX, 'alr', basis_return = FALSE),
  compositions::alr(aX),
  robCompositions::addLR(mX),
  easyCODA::ALR(mX),
  RcppCoDA::alr(mX),
  times = 500
)
autoplot(alr_results)
(res <- summary(alr_results)[,c('expr', 'lq', 'median', 'uq')])[order(res$median),]

clr_results = microbenchmark(
  coda.base::coordinates(mX, 'clr', basis_return = FALSE),
  compositions::clr(aX),
  RcppCoDA::clr(mX),
  robCompositions::cenLR(mX),
  easyCODA::CLR(mX), times = 500
)
autoplot(clr_results)
(res <- summary(clr_results)[,c('expr', 'lq', 'median', 'uq')])[order(res$median),]

ilr_results = microbenchmark(
  coda.base::coordinates(mX),
  compositions::ilr(aX),
  RcppCoDA::ilr(mX),
  # robCompositions::isomLR(mX, fast=T),
  robCompositions::pivotCoord(mX, fast=T),
  easyCODA::PLR(mX), times = 500
)
autoplot(ilr_results)
(res <- summary(ilr_results)[,c('expr', 'lq', 'median', 'uq')])[order(res$median),]

B = sbp_basis(b1 = V1~V2,
              b2 = b1~V3,
              b3 = b2~V4,
              b4 = b3~V5,
              b5 = b4~V6,
              b6 = b5~V7,
              b7 = b6~V8,
              b8 = b7~V9,
              b9 = b8~V10, data = dX)
B.t = t(B)
mX.t = t(mX)
B_results = microbenchmark(
  coda.base::coordinates(mX, B),
  compositions::ilr(aX, B),
  RcppCoDA::ilr(mX.t, B.t),
  times = 500
)
autoplot(B_results)

microbenchmark(
  variation_array(X, only_variation = TRUE),
  compositions::variation(aX),
  robCompositions::variation(X, robust=FALSE),times = 1000)
