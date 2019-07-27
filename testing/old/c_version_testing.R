library(microbenchmark)
library(Rcpp)
sourceCpp('src/pbalances.cpp')
source('testing/principal_balances_functions_2b.R')
set.seed(1)
X=as.data.frame(matrix(exp(rnorm(10*10)), nrow=10, ncol=10))
M = cov(log(X))
R_test_01 = function(M, L0, R0){
  sol = build(M, L = L0, R = R0)
  print(unlist(sol)[c('sL', 'sR', 'sM', 'var')])
  sol = addL(3,sol)
  print(unlist(sol)[c('sL', 'sR', 'sM', 'var')])
  top = toptree(sol)
  sol.top = build(top$M, L0, R0, top$nodes)
  print(unlist((sol.top)[c('sL', 'sR', 'sM', 'var')]))
  sol.top = addR(3, sol.top)
  print(unlist((sol.top)[c('sL', 'sR', 'sM', 'var')]))
  top = toptree(sol.top)
  sol.top = build(top$M, L0, R0, top$nodes)
  print(unlist((sol.top)[c('sL', 'sR', 'sM', 'var')]))
  sol.top = addL(5,sol.top)
  sol.top = addL(6,sol.top)
  sol.top = addR(3,sol.top)
  sol.top = addR(4,sol.top)
  print(unlist((sol.top)[c('sL', 'sR', 'sM', 'var')]))
  left = subtreeL(sol.top)
  sol.left = build(left$M, L0, R0, left$nodes)
  print(unlist((sol.left)[c('sL', 'sR', 'sM', 'var')]))
  right = subtreeR(sol.top)
  sol.right = build(right$M, L0, R0, right$nodes)
  print(unlist((sol.right)[c('sL', 'sR', 'sM', 'var')]))
}
R_test_02 = function(M, L0, R0){
  sol = build(M, L = L0, R = R0)
  #print(unlist(sol)[c('sL', 'sR', 'sM', 'var')])
  ls_sol = improve_until_finished(sol)
  #print(unlist(ls_sol)[c('sL', 'sR', 'sM', 'var')])
}
R_test_01(M, 1, 2)
c_test_01(M, 0, 1)

R_test_02(M, 1, 2)
c_test_02(M, 0, 1)

microbenchmark(
  R_test_02(M, 2, 3),
  c_test_02(M, 1, 2))

microbenchmark(
  R_test_01(M, 1, 2),
  c_test_01(M, 0, 1))



