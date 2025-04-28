pb_subcomposition = function(X, variables = 1:ncol(X), constraints = NULL){
  Xc = X
  if(!is.null(constraints)){
    for(constraint in constraints){
      Xc[,constraint] = exp(rowMeans(log(Xc[,constraint])))
    }
  }
  pb = matrix(0, nrow = ncol(X), ncol = 1)
  pb[variables,1] = get_balance_using_pc(Xc[,variables])[,1]

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

constrained_pb = function(X){
  l_pb = pb_subcomposition(X)
  return(sapply(l_pb, identity))
}
