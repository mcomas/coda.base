build = function(L, R, O, M, forbid = NULL){

  nR = length(R)
  nL = length(L)

  sR = (nL/nR) * sum(M[R,R])
  sL = (nR/nL) * sum(M[L,L])
  sM = - 2*sum(M[R,L])
  list('L' = L, 'R' = R, 'O' = O, forbid = forbid,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = 0)
}
addL = function(i, sol){
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
  list('L' = L, 'R' = R, 'O' = O, forbid = sol$forbid,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = sol$iter + 1)
}
removeL = function(i, sol){
  if(!i %in% sol$L){
    stop('Not in L set')
  }
  L = setdiff(sol$L,i)
  R = sol$R
  O = union(sol$O, i)

  nR = length(R)
  nL = length(L)

  sR = sol$sR * nL/(nL+1)
  sL = (sol$sL  - (nR/(nL+1)) * (2*sum(M[i,sol$L]) - M[i,i])) * (nL+1)/nL
  sM = sol$sM + 2*sum(M[R,i])
  list('L' = L, 'R' = R, 'O' = O, forbid = sol$forbid,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = sol$iter + 1)
}
removeR = function(i, sol){
  if(!i %in% sol$R){
    stop('Not in R set')
  }
  L = sol$L
  R = setdiff(sol$R,i)
  O = union(sol$O, i)

  nR = length(R)
  nL = length(L)

  sR = (sol$sR  - (nL/(nR+1)) * (2*sum(M[i,sol$R]) - M[i,i])) * (nR+1)/nR
  sL = sol$sL * nR/(nR+1)
  sM = sol$sM + 2*sum(M[L,i])
  list('L' = L, 'R' = R, 'O' = O, forbid = sol$forbid,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = sol$iter + 1)
}
addR = function(i, sol){
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
  list('L' = L, 'R' = R, 'O' = O, forbid = sol$forbid,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = ((nL/nR) * sum(M[R,R]) + (nR/nL) * sum(M[L,L]) - 2*sum(M[R,L])) / (nL+nR),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = sol$iter + 1)
}
variances = function(sols) sapply(sols, function(sol) sol$var)
improve = function(sol){
  bestSol = sol
  if(length(sol$O) > 0){
    sols.L = lapply(sol$O, addL, sol)
    sols.R = lapply(sol$O, addR, sol)
    v.addL = max(variances(sols.L))
    v.addR = max(variances(sols.R))
    if(v.addL > bestSol$var){
      bestSol = sols.L[[which.max(variances(sols.L))]]
    }
    if(v.addR > bestSol$var){
      bestSol = sols.R[[which.max(variances(sols.R))]]
    }
  }
  if(length(sol$L) > 1){
    sols = lapply(sol$L, removeL, sol)
    v.remove = max(variances(sols))
    if(v.remove > bestSol$var){
      bestSol = sols[[which.max(variances(sols))]]
    }
  }
  if(length(sol$R) > 1){
    sols = lapply(sol$R, removeR, sol)
    v.remove = max(variances(sols))
    if(v.remove > bestSol$var){
      bestSol = sols[[which.max(variances(sols))]]
    }
  }
  bestSol
}
improve_until_finished = function(M, L0, R0){
  sol.prev = list(var=-Inf)
  sol.next <- build(L = L0, R = R0, O = setdiff(1:nrow(M), c(R0,L0)), M)
  while(sol.prev$var != sol.next$var){
    sol.prev = sol.next
    sol.next = improve(sol.prev)
  }
  sol.next
}
findPB = function(M, TESTS = 10, I = 1, J = 2){
  if(nrow(M) == 2){
    return(list('L' = 1, 'R' = 2))
  }
  sims = replicate(TESTS, {
    ind = 1:nrow(M)
    L = sample(ind, I)
    R = sample(setdiff(ind, L), J)
    list(R = R, L = L)
  }, simplify=FALSE)
  sols_01 = lapply(sims, function(sim) improve_until_finished(M, sim$L, sim$R))
  sols_01[[which.max(variances(sols_01))]]
}

