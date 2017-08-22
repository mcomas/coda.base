set.seed(1)
library(Rcpp)
library(coda.base)
sourceCpp('src/principal_balances.cpp')
K = 1000
X=as.data.frame(matrix(exp(rnorm(K*K*2)), nrow=2*K, ncol=K))
M = cov(log(X))

## This example has the lowest variance precisly in the balance separating all parts.
##
H = coordinates(X, 'pc')
h.pc = attr(H, 'basis')[,1]
eps = 0.01
div = ifelse(h.pc > eps, 1,
             ifelse(h.pc < -eps, -1, 0))
h = eval(parse(text = sprintf("coordinates(X, sbp_basis(%s~%s, data=X))",
                              paste(names(div[div == -1]), collapse='+'),
                              paste(names(div[div == 1]), collapse='+'))))
var(h)

ORDERING = order(abs(h.pc), decreasing = T)
LR = sign(h.pc)
sapply(500:503, function(i){
  ORD = head(ORDERING, i)
  h = eval(parse(text = sprintf("coordinates(X, sbp_basis(%s~%s, data=X))",
                                paste(names(LR[LR == -1 & 1:1000 %in% ORD]), collapse='+'),
                                paste(names(div[div == 1 & 1:1000 %in% ORD]), collapse='+'))))
  var(h)
})



B0 = pb_basis(X, 'ward.D2')
var(coordinates(X, B0)[,1])

B1 = find_PB(M, 1)
var(coordinates(X, B1)[,1])

B4 = find_PB4(M)
var(coordinates(X, B4)[,1])

B3a = find_PB3(M, 1, K, .Machine$integer.max)
var(coordinates(X, B3a)[,1])

B3b = find_PB3(M, 10, K, .Machine$integer.max, 100)
var(coordinates(X, B3b)[,1])

B3c = find_PB3(M,
               steps = 1,
               random = 0,
               optim = .Machine$integer.max,
               k = 0)
var(coordinates(X, B3c)[,1])


B2 = find_PB2(M, 1000, .Machine$integer.max)
var(coordinates(X, B2)[,1])


