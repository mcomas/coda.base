library(coda.base)
load('~/X.RData')
X0 = X

pb_constrained = function(X, variables = 1:ncol(X), constraints = NULL){
  if(!is.null(constraints)){
    for(constraint in constraints){
      X[,constraint] = exp(rowMeans(log(X[,constraint])))
    }
  }
  pb = matrix(0, nrow = ncol(X), ncol = 1)
  pb[variables,1] = balance_using_pc(X[,variables])[,1]
  pb
}

pb_subcomposition = function(X, variables = 1:ncol(X), constraints = NULL){
  # pb = pb_constrained(X, variables, constraints)
  Xc = X
  if(!is.null(constraints)){
    for(constraint in constraints){
      Xc[,constraint] = exp(rowMeans(log(Xc[,constraint])))
    }
  }
  pb = matrix(0, nrow = ncol(X), ncol = 1)
  pb[variables,1] = balance_using_pc(Xc[,variables])[,1]

  lpb = list(pb)
  if(sum(pb<0)>1 ){
    if(max(apply(Xc[,pb<0], 1, var))>0){
      lpb = c(lpb, Recall(X, variables = which(pb<0), constraints = constraints))
    }
  }
  if(sum(pb>0)>1 ){
    if(max(apply(Xc[,pb>0], 1, var))>0){
      lpb = c(lpb, Recall(X, variables = which(pb>0), constraints = constraints))
    }
  }
  sel = rep(F, ncol(X))
  sel[variables] = T
  if(sum((pb == 0) & sel) > 0){
    if(is.null(constraints)){
      constraints = list()
    }
    constraints[[1+length(constraints)]] = (pb != 0) & sel
    # print("Up with: ")
    # print(variables)
    # print(constraints)
    lpb = c(lpb, Recall(X, variables, constraints = constraints))
    # print("Search complete")
  }
  return(lpb)
}
X = matrix(rlnorm(10*15), ncol = 10)
l_pb = pb_subcomposition(X)
l_pb = l_pb[order(sapply(l_pb, function(b) var(coordinates(X, b))), decreasing = TRUE)]
CPB = Matrix::Matrix(sapply(l_pb, identity), sparse=TRUE)
## Primer extenem branca dreta i esquerra. Al final cap amunt.
EPB = pb_basis(X, method = 'exact')
apply(coordinates(X, CPB), 2, var)
apply(coordinates(X, EPB), 2, var)



lpb = list()
if(sum(pb1 == 0) > 0){
  lpb[[1+length(lpb)]] = pb_constrained(X, constraint = apply(sapply(PB, function(pb) pb != 0), 1, any))
}
if(sum(pb1<0)>1){
  lpb[[1+length(lpb)]] = pb_constrained(X, variables = which(pb1<0))
}
lpb_var = sapply(lpb, function(b) var(coordinates(X, b)))
PB[[length(PB) + 1]] = lpb[[which.max(lpb_var)]]



func = function(X){
  balance = balance_using_pc(X)
  Xc = X
  Xc[,balance!=0] = exp(rowMeans(log(Xc[,balance!=0])))
  return(list(balance, Xc))
}
res1 = func(X)
res2 = func(res1[[2]])
res3 = func(res2[[2]])

balance = balance_using_pc(Xc)


X1 = X
bal1 = balance_using_pc(X)

initial_states = function(X) list(
  list(X = X, variables = 1:ncol(X), super_balance = NULL)
)

eval_states = function(state, X){
  X = state$X
  if(!is.null(state$super_balance)){
    X[,state$super_balance] = exp(rowMeans(log(X[,state$super_balance])))
  }
  X = X[,state$variables]
  balance = balance_using_pc(X)
  list(
    balance = balance,
    variance = var(coordinates(X, balance)))
}

possible_state = function(state, balance){
  states = list()
  if(sum(balance == 0) > 0){
    Xc = state$X
    Xc[,balance!=0] = exp(rowMeans(log(Xc[,balance!=0])))
    super_state = list(X = Xc,
                       variables = state$variables,
                       super_balance = balance!=0)
    states[[length(states)+1]] = super_state
  }
  if(sum(balance<0) > 1){
    left_state = list(X = state$X[,balance<0],
                      variables = which(balance<0),
                      super_balance = NULL)
    states[[length(states)+1]] = left_state
  }
  if(sum(balance>0) > 1){
    right_state = list(X = state$X[,balance>0],
                       variables = which(balance>0),
                       super_balance = NULL)
    states[[length(states)+1]] = right_state
  }
  states
}

curr_states = initial_states(X)
var_states = lapply(curr_states, eval_states)
best_state = which.max(sapply(var_states, `[[`, 2))
curr_states[[best_state]]

curr_states2 = possible_state(curr_states[[best_state]], var_states[[best_state]]$balance)
lapply(curr_states2, eval_states)
eval_states(curr_states2[[1]])
eval_states(curr_states2[[2]])

func = function(X, variables = 1:ncol(X), super_balance = NULL){
  if(!is.null(super_balance)){
    X[,balance!=0] = exp(rowMeans(log(X[,balance!=0])))
  }
  balance = balance_using_pc(X[,variables])

  if(sum(balance != 0) > 0){
    func(X, variables, super_balance = balance != 0)
  }


  super_balance = NULL
  variables = which(balance < 0)

  super_balance = NULL
  variables = which(balance > 0)

}


pb_constrained = function(lX){
  lbalances = lapply(lX, balance_using_pc)
  variances = mapply(function(X,bal) var(coordinates(X,bal)), lX, lbalances)
  balance = lbalances[[which.max(variances)]]
  datasets = lX[-which.max(variances)]
  X = lX[[which.max(variances)]]
  if(sum(balance == 0) > 0){
    Xc = X
    Xc[,balance!=0] = exp(rowMeans(log(Xc[,balance!=0])))
    datasets = c(list(Xc))
  }
  if(sum(balance<0) > 1){
    Xl = X[,balance<0]
    if(!all(apply(Xl, 1, var) == 0)){
      datasets = c(datasets, list(Xl))
    }
  }
  if(sum(balance>0) > 1){
    Xr = X[,balance>0]
    if(!all(apply(Xr, 1, var) == 0)){
      datasets = c(datasets, list(Xr))
    }
  }
  list(balance, datasets)
}

res1 = pb_constrained(list(X))
res1[[1]]
res2 = pb_constrained(res1[[2]])
res2[[1]]
res3 = pb_constrained(res2[[2]])
res3[[1]]
res4 = pb_constrained(res3[[2]])
res4[[1]]

X2 = X1
X2[,bal1!=0] = exp(rowMeans(log(X2[,bal1!=0])))
b_center = balance_using_pc(X2)
if(sum(bal1<0) > 1){
  b_left = balance_using_pc(X2[,bal1<0])
}
if(sum(bal1>0) > 1){
  b_right = balance_using_pc(X2[,bal1>0])
}
