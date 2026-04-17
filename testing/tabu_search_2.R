set.seed(1)
library(coda.base)

.CONT = TRUE
D = 15
while(.CONT){
  X = matrix(rlnorm(10*D), ncol = D)
  B0 = pb_basis(X, method = 'constrained')
  B0_opt = pb_basis(X, method = 'exact')
  v0 = apply(coordinates(X, B0), 2, var)
  v0_opt = apply(coordinates(X, B0_opt), 2, var)
  if(abs(v0[1]-v0_opt[1]) > 0.001) .CONT = FALSE
}


var(coordinates(X, pb_basis(X, method = 'constrained'))[,1])
var(coordinates(X, pb_basis(X, method = 'exact'))[,1])
TS = partial_pb_tabu_search(X, lapply(1:ncol(X), identity),
                            iter = 10, tabu_size = 100, debug = TRUE)
TS$variance
TS$iter_best
TS$tabu_size

