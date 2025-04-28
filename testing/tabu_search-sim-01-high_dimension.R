library(mvtnorm)
source('testing/tabu_search_functions_2.R')

set.seed(1)
d = 100
for(d in 5:50){
  H = rmvnorm(1000, rep(0, d))
  S = cov(H)
  EIG = eigen(S)

  PB1 = sample(c(-1,0,1), d+1, replace = TRUE)
  PB = sbp_basis(cbind(PB1), fill = TRUE)
  SBP_exact = sign(PB)

  X = composition(H %*% EIG$vectors, basis = PB)
  SBP_ts = pb_tabu_search(X, iter = floor(0.5 * ncol(X)))

  cat(sprintf("D:%d\tOptimum: %d\tStep: %d\n",
              d+1,
              all(SBP_ts == SBP_exact | SBP_ts == -1*SBP_exact),
              attr(SBP_ts, 'max_steps')))
}

#
# set.seed(1)
# d = 100
# H = rmvnorm(1000, rep(0, d))
# S = cov(H)
# EIG = eigen(S)
#
# set.seed(1)
# PB1 = sample(c(-1,0,1), d+1, replace = TRUE)
# B = sbp_basis(cbind(PB1), fill = TRUE)
#
# X = composition(H %*% EIG$vectors, basis = B)
#
# PB1_ts = pb_tabu_search(X, iter = ncol(X)/2)
#
# ## Millor trobat
# PB1_const = pb_basis(X, method = 'constrained')
# all(sign(PB1_const) == PB1) | all(sign(PB1_const) == -1 * PB1)
# all(PB1_ts == B) | all(PB1_ts == -1 * B)
#
#
#
#
# cbind(sign(PB1_ts$balance), PB1)
# var(X %*% pc1)
# sapply(1:30, function(i) var(X %*% a$vectors[,i]))
# a$vectors
#
#
# sbp_basis(sbp = cbind(c(1,-1,0,0)), fill = TRUE)
