eval_var = function(M, L, R, nodes){
  return(NULL)
  nL = sum(sapply(nodes[L], length))
  nR = sum(sapply(nodes[R], length))

  (nL*nR)/(nL+nR) * (1/nL^2 * sum(M[L,L]) + 1/nR^2 * sum(M[R,R]) - 2/(nL*nR) * sum(M[L,R]))
}

build = function(M, L, R, nodes = as.list(1:nrow(M))){

  nL = sum(sapply(nodes[L], length))
  nR = sum(sapply(nodes[R], length))

  sL = (nR/nL) * sum(M[L,L])
  sR = (nL/nR) * sum(M[R,R])
  sM = - 2*sum(M[R,L])

  O = setdiff(1:length(nodes), c(L,R))
  list('M' = M, 'nodes' = nodes,
       'L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = eval_var(M,L,R,nodes),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = 0)
}
addL = function(i, sol){
  M = sol$M

  nL = sum(sapply(sol$nodes[sol$L], length))
  nR = sum(sapply(sol$nodes[sol$R], length))
  ni = sum(sapply(sol$nodes[i], length))

  sL = sol$sL * nL/(nL+ni) + (nR/(nL+ni)) * (2*sum(M[i,sol$L]) + M[i,i])
  sR = sol$sR * (nL+ni)/nL
  sM = sol$sM - 2*sum(M[sol$R,i])

  L = union(sol$L,i)
  R = sol$R
  O = setdiff(1:length(sol$nodes), c(L,R))

  list('M' = M, 'nodes' = sol$nodes,
       'L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = eval_var(M,L,R,sol$nodes),
       'var' = (sL + sR + sM) / (nL+nR+ni), 'iter' = sol$iter + 1)
}
removeL = function(i, sol){
  M = sol$M

  L = setdiff(sol$L,i)
  R = sol$R
  O = setdiff(1:length(sol$nodes), c(L,R))

  nR = sum(sapply(sol$nodes[R], length))
  nL = sum(sapply(sol$nodes[L], length))
  ni = sum(sapply(sol$nodes[i], length))

  sL = (sol$sL - nR/(nL+ni) * (sum(M[i,i]) + 2 * sum(M[i,L]))) * (nL+ni) / nL
  sR = sol$sR * nL/(nL+ni)
  sM = sol$sM + 2*sum(M[R,i])

  list('M' = M, 'nodes' = sol$nodes,
       'L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = eval_var(M,L,R,sol$nodes),
       'var' = (sL + sR + sM) / (nL+nR), 'iter' = sol$iter + 1)
}
addR = function(i, sol){
  M = sol$M

  nL = sum(sapply(sol$nodes[sol$L], length))
  nR = sum(sapply(sol$nodes[sol$R], length))
  ni = sum(sapply(sol$nodes[i], length))

  sL = sol$sL * (nR+ni)/nR
  sR = sol$sR * nR/(nR+ni) + (nL/(nR+ni)) * (2*sum(M[i,sol$R]) + M[i,i])
  sM = sol$sM - 2*sum(M[sol$L,i])

  L = sol$L
  R = union(sol$R,i)
  O = setdiff(1:length(sol$nodes), c(L,R))

  list('M' = M, 'nodes' = sol$nodes,
       'L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = eval_var(M,L,R,sol$nodes),
       'var' = (sL + sR + sM) / (nL+nR+ni), 'iter' = sol$iter + 1)
}
removeR = function(i, sol){
  M = sol$M

  L = sol$L
  R = setdiff(sol$R,i)
  O = setdiff(1:length(sol$nodes), c(L,R))

  nL = sum(sapply(sol$nodes[L], length))
  nR = sum(sapply(sol$nodes[R], length))
  ni = sum(sapply(sol$nodes[i], length))

  sL = sol$sL * nR/(nR+ni)
  sR = (sol$sR - nL/(nR+ni) * (sum(M[i,i]) + 2 * sum(M[i,R]))) * (nR+ni) / nR
  sM = sol$sM + 2*sum(M[L,i])

  list('M' = M, 'nodes' = sol$nodes,
       'L' = L, 'R' = R, 'O' = O,
       'sL' = sL, 'sR' =  sR, 'sM' = sM,
       'var-check' = eval_var(M,L,R,sol$nodes),
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
improve_until_finished = function(sol){
  sol.prev = list(var=-Inf)
  sol.next <- sol

  while(sol.prev$var != sol.next$var){
    sol.prev = sol.next
    sol.next = improve(sol.prev)
  }
  sol.next
}
findPB = function(tree, TESTS = 10){
  M = tree$M
  nodes = tree$nodes
  if(length(nodes) == 2){
    nL = length(nodes[[1]])
    nR = length(nodes[[2]])
    return(list('M' = M, 'nodes' = nodes, 'L' = 1, 'R' = 2,
                'var' = (nR*nL)/(nR+nL) * ( sum(M[1,1])/nL^2 + sum(M[2,2])/nR^2 - 2*sum(M[1,2])/(nL*nR))))
  }
  sims = replicate(TESTS, {
    ind = 1:nrow(M)
    v = sample(ind, 2)
    L = v[1]
    R = v[2]
    O = setdiff(ind, v)
    sel = sample(c(0,1,2), length(O), replace = TRUE)
    L = c(L, O[sel == 1])
    R = c(R, O[sel == 2])
    list(R = R, L = L)
  }, simplify=FALSE)
  sols_01 = lapply(sims, function(sim) improve_until_finished(build(M, sim$L, sim$R, nodes)))
  sols_01[[which.max(variances(sols_01))]]
}

rmeans = function(x) if(is.vector(x)) mean(x) else rowMeans(x)
cmeans = function(x) if(is.vector(x)) mean(x) else colMeans(x)
rsums = function(x) if(is.vector(x)) sum(x) else rowSums(x)

toptree = function(sol){
  M = sol$M
  V = c(sol$L,sol$R)
  I = setdiff(1:nrow(M),V)
  list(
    'nodes' = c(lapply(I, function(i) sol$nodes[[i]]), list(unlist(sol$nodes[V]))),
    'M' = rbind(
      cbind(sol$M[I,I],rsums(sol$M[I,V])),
      c(rsums(sol$M[I,V]), sum(sol$M[V,V]))))
}
subtreeL = function(sol){
  list(
    'nodes' = lapply(sol$L, function(i) sol$nodes[[i]]),
    'M' = sol$M[sol$L,sol$L])
}
subtreeR = function(sol){
  list(
    'nodes' = lapply(sol$R, function(i) sol$nodes[[i]]),
    'M' = sol$M[sol$R,sol$R])
}
