set.seed(8)
## 8 interessant
library(coda.base)

.CONT = TRUE
while(.CONT){
  X = matrix(rlnorm(10*10), ncol = 10)
  B0 = pb_basis(X, method = 'constrained')
  B0_opt = pb_basis(X, method = 'exact')
  v0 = apply(coordinates(X, B0), 2, var)
  v0_opt = apply(coordinates(X, B0_opt), 2, var)
  if(abs(v0[1]-v0_opt[1]) > 0.001) .CONT = FALSE
}

var(coordinates(X, pb_basis(X, method = 'constrained'))[,1])
var(coordinates(X, pb_basis(X, method = 'exact'))[,1])

source('testing/tabu_search_functions.R')
a = pb_tabu_search(X, iter = 100, tabu_size = Inf)#ncol(X))
a$iter_best
a$tabu_size
plot(a$steps, type = 'l')
points(1+a$iter_best, var(coordinates(X, a$balance)), col=2)



#### TEST 2
set.seed(28)
X = matrix(rlnorm(200*500), ncol = 200)
B0 = pb_basis(X, method = 'constrained')
a = pb_tabu_search(X, iter = 100, tabu_size = 10)#ncol(X))
a$iter_best
a$tabu_size
plot(a$steps, type = 'l')
points(1+a$iter_best, var(coordinates(X, a$balance)), col=2)
