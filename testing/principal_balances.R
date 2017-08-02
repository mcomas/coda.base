set.seed(1)
X = as.data.frame(matrix(exp(rnorm(10*100)), nrow=100, ncol=10))

lX = log(X)
M = cov(lX)

build = function(L, R, O, M){

  nR = length(R)
  nL = length(L)

  sR = (nL/nR) * sum(M[R,R])
  sL = (nR/nL) * sum(M[L,L])
  sM = - 2*sum(M[R,L])
  list('L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var' = (sL + sR + sM) / (nL+nR))
}
addL = function(sol, i){
  if(!i %in% sol$O){
    stop('Not in O set')
  }
  L = union(sol$L,i)
  R = sol$R
  O = setdiff(sol$O, i)

  nR = length(R)
  nL = length(L)

  sR = sol$sR * nL/(nL-1)
  sL = sol$sL * (nL-1)/nL + (nR/nL) * (2*sum(M[i,L]) - M[i,i])
  sM = sol$sM - 2*sum(M[R,i])
  list('L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var' = (sL + sR + sM) / (nL+nR))
}
addR = function(sol, i){
  if(!i %in% sol$O){
    stop('Not in O set')
  }
  L = sol$L
  R = union(sol$R, i)
  O = setdiff(sol$O, i)

  nR = length(R)
  nL = length(L)

  sR = sol$sR * (nR-1)/nR + (nL/nR) * (2*sum(M[i,R]) - M[i,i])
  sL = sol$sL * nR/(nR-1)
  sM = sol$sM - 2*sum(M[L,i])
  list('L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var' = (sL + sR + sM) / (nL+nR))
}
removeL = function(sol, i){
  if(!i %in% sol$L){
    stop('Not in K set')
  }
  L = setdiff(sol$L,i)
  R = sol$R
  O = union(sol$O, i)

  nR = length(R)
  nL = length(L)

  sR = sol$sR * (nL+1)/nL
  sL = sol$sL * nL/(nL+1) - (nR/(nL+1)) * (2*sum(M[i,sol$L]) - M[i,i])
  sM = sol$sM + 2*sum(M[R,i])
  list('L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var' = (sL + sR + sM) / (nL+nR))
}

sol1 = build(L = c(3,4,5), R = c(1,2), O = 6:10, M)
sol2 = addL(sol1, 6)
sol3 = removeL(sol2, 6)
sol1$var
sol3$var
addR(sol, 6)
remove(sol, )

K = ncol(X)
R = c(1,2)
L = c(3,4,5)
O = 6:10
nR = length(R)
nL = length(L)
k = nL + nR
w = rep(0,K)
w[R] = 1/nR
w[L] = -1/nL

library(microbenchmark)
(nL*nR)/(nL+nR) * sum((w %*% t(w)) * M)
(nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) - 2/(nR*nL) * sum(M[R,L]))
(nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) ) - 2 * sum(M[R,L]) / (nL+nR)
sR = (nL/nR) * sum(M[R,R])
sL = (nR/nL) * sum(M[L,L])
sM = - 2*sum(M[R,L])
(sR + sL + sM) / (nL+nR)

K = ncol(X)
R = c(1,2,6)
L = c(3,4,5)
O = 7:10
nR = length(R)
nL = length(L)
k = nL + nR
w = rep(0,K)
w[R] = 1/nR
w[L] = -1/nL

(nL*nR)/(nL+nR) * sum((w %*% t(w)) * M)
(nL*nR)/(nL+nR) * ( (1/nR)^2 * sum(M[R,R]) + (1/nL)^2 * sum(M[L,L]) - 2/(nR*nL) * sum(M[R,L]))
((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR)
sR = sR * nL/(nL-1)
sL = sL * (nL-1)/nL + (nR/nL) * (2*sum(M[6,L]) - M[6,6])
sM = sM - 2*sum(M[R,6])
(sR + sL + sM) / (nL+nR)


