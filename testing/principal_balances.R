source('testing/principal_balances_functions.R')


set.seed(1)
X = as.data.frame(matrix(exp(rnorm(10*100)), nrow=10, ncol=100))
M = cov(log(X))


TESTS = 10
sims = replicate(TESTS, {
  I = 3
  J = 5
  ind = 1:nrow(M)
  L = sample(ind, I)
  R = sample(setdiff(ind, L), J)
  list(R = R, L = L)
}, simplify=FALSE)
sols_01 = lapply(sims, function(sim) improve_until_finished(M, sim$L, sim$R))

sims_02 = combn(ncol(M),2)
sols_02 = apply(sims_02[,sample(1:ncol(sims_02),TESTS)], 2, function(v){
  improve_until_finished(M, v[1], v[2])
})
max(sapply(sols_01, function(sol) sol$var))
max(sapply(sols_02, function(sol) sol$var))

pb_01 = gsi.PrinBal(acomp(X), "PBhclust")
pb_02 = gsi.PrinBal(acomp(X), "PBmaxvar")
pb_03 = gsi.PrinBal(acomp(X), "angprox")

data("Hydrochem")
x = acomp(Hydrochem[,c(6:10)])
(v1 = gsi.PrinBal(x, "PBhclust"))
(v2 = gsi.PrinBal(x, "PBmaxvar"))
(v3 = gsi.PrinBal(x, "PBangprox"))

M = cov(log(Hydrochem[,c(6:10)]))
TESTS = 10
sims = replicate(TESTS, {
  I = 2
  J = 1
  ind = 1:nrow(M)
  L = sample(ind, I)
  R = sample(setdiff(ind, L), J)
  list(R = R, L = L)
}, simplify=FALSE)
sols_01 = lapply(sims, function(sim) improve_until_finished(M, sim$L, sim$R))
#
#
# L0 = 1
# R0 = 2
# (sol <- build(L = L0, R = R0, O = setdiff(1:nrow(M), c(R0,L0)), M))
# (sol <- addL(6, sol))
# (sol <- addR(4, sol))
# (sol <- removeL(1, sol))
# (sol <- removeR(4, sol))
# (sol <- addR(5, sol))
#
# variances = function(sols) sapply(sols, function(sol) sol$var)
#
# L0 = 1
# R0 = 2
# (sol <- build(L = L0, R = R0, O = setdiff(1:nrow(M), c(R0,L0)), M))
#
#
# I = 3
# J = 5
# sims = replicate(100, {
#   ind = 1:nrow(M)
#   L = sample(ind, I)
#   R = sample(setdiff(ind, L), J)
#   list(R = R, L = L)
# }, simplify=FALSE)
# sols_01 = lapply(sims, function(sim) improve_until_finished(M, sim$L, sim$R))
# sols_02 = apply(combn(ncol(M),2), 2, function(v){
#   improve_until_finished(M, v[1], v[2])
# })
# max(sapply(sols_01, function(sol) sol$var))
# max(sapply(sols_02, function(sol) sol$var))
#
#
#
# #
# #
# #
# #
# #
# #
# #
# #
#
#
# K = ncol(X)
# R = c(1,2)
# L = c(3,4,5)
# O = 6:10
# nR = length(R)
# nL = length(L)
# k = nL + nR
# w = rep(0,K)
# w[R] = 1/nR
# w[L] = -1/nL
#
# library(microbenchmark)
# (nL*nR)/(nL+nR) * sum((w %*% t(w)) * M)
# (nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) - 2/(nR*nL) * sum(M[R,L]))
# (nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) ) - 2 * sum(M[R,L]) / (nL+nR)
# sR = (nL/nR) * sum(M[R,R])
# sL = (nR/nL) * sum(M[L,L])
# sM = - 2*sum(M[R,L])
# (sR + sL + sM) / (nL+nR)
#
# K = ncol(X)
# R = c(1,2,6)
# L = c(3,4,5)
# O = 7:10
# nR = length(R)
# nL = length(L)
# k = nL + nR
# w = rep(0,K)
# w[R] = 1/nR
# w[L] = -1/nL
#
# (nL*nR)/(nL+nR) * sum((w %*% t(w)) * M)
# (nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) - 2/(nR*nL) * sum(M[R,L]))
# ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR)
# sR = sR * nL/(nL-1)
# sL = sL * (nL-1)/nL + (nR/nL) * (2*sum(M[6,L]) - M[6,6])
# sM = sM - 2*sum(M[R,6])
# (sR + sL + sM) / (nL+nR)
#
#
