library(coda.base)
source('testing/principal_balances_functions_2b.R')


#set.seed(1)
K = 100
X = as.data.frame(matrix(exp(rnorm(100*K)), nrow=100, ncol=K))
M = cov(log(X))

#####
PB = list()
current_trees = list()
current_sols = list()
new_trees = list()
new_trees[[1]] = list(M=M,nodes=as.list(1:nrow(M)))
new_sols = lapply(new_trees, findPB, TESTS = 20)


for(i in 1:(nrow(M)-1)){
  current_trees = c(current_trees, new_trees)
  current_sols = c(current_sols, new_sols)
  BEST = which.max(sapply(current_sols, function(sols) sols$var))
  PB[[i]] = current_sols[[BEST]]
  new_trees = list()
  if(length(PB[[i]]$nodes) != length(c(PB[[i]]$L, PB[[i]]$R))){
    new_trees[[length(new_trees)+1]] = toptree(PB[[i]])
  }
  if(length(PB[[i]]$L) > 1){
    new_trees[[length(new_trees)+1]] = subtreeL(PB[[i]])
  }
  if(length(PB[[i]]$R) > 1){
    new_trees[[length(new_trees)+1]] = subtreeR(PB[[i]])
  }
  new_sols = lapply(new_trees, findPB)

  current_sols[[BEST]] = NULL
  current_trees[[BEST]] = NULL
}

x.var = sapply(PB, function(pb) pb$var)
plot(x.var)

solPB1
if(NCOL(SOL) == 2){

}

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

library(compositions)

data("Hydrochem")
X = Hydrochem[,c(6:10)]
x = acomp(X)
(v1 <- gsi.PrinBal(x, "PBhclust"))
(v2 <- gsi.PrinBal(x, "PBmaxvar"))
(v3 <- gsi.PrinBal(x, "PBangprox"))

M = cov(log(Hydrochem[,c(6:10)]))
SOL = findPB(M, TESTS = 10)
b1 = sprintf("%s ~ %s",
             paste(names(X)[SOL$L], collapse='+'),
             paste(names(X)[SOL$R], collapse='+'))

apply(coordinates(X, v1), 2, var)
apply(coordinates(X, v2), 2, var)
apply(coordinates(X, v3), 2, var)
apply(coordinates(X, sbp_basis(b1 = eval(parse(text=b1)), data=X)), 2, var)

list(M[SOL$L, SOL$L], M[SOL$L, SOL$L], "(M, amb restriccio [M i L] mateix costat")
pb = function(lM){
  lapply(lM, findPB)
}


#
#
L0 = 1
R0 = 2
(sol <- build(M, L = L0, R = R0, nodes = 1:nrow(M)))
(sol <- addL(6, sol))
(sol <- addL(4, sol))
(sol <- removeL(1, sol))
(sol <- removeL(6, sol))



(sol <- addR(4, sol))
(sol <- removeL(1, sol))
(sol <- removeR(4, sol))
(sol <- addR(5, sol))
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
