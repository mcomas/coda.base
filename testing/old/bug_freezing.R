library(coda.base)
#for(i in 1:1000){
i = 408
  print(i)
  set.seed(i)
  X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
  # Optimal variance obtained with Principal components
  ##apply(coordinates(X, 'pc'), 2, var)
  # Solution obtained using a hill climbing algorithm from pc approximation
  print('v1')
  v1 = apply(coordinates(X,pb_basis(X, method='lsearch')), 2, var)
  print('v2')
  # Solution obtained using a hill climbing algorithm using 10 restartings
  v2 = apply(coordinates(X,pb_basis(X, method='lsearch', rep=10)), 2, var)
  print('v3')
  v3 = apply(coordinates(X,pb_basis(X, method='exact')), 2, var)
  if(!isTRUE(all.equal(v1,v2)) & !isTRUE(all.equal(v2,v3))){
    print(i)
  }
#}
