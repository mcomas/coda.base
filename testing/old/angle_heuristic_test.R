library(microbenchmark)
library(Rcpp)
library(coda.base)
sourceCpp('src/pbalances.cpp')
source('testing/principal_balances_functions_2b.R')
set.seed(1)

K = 7
X=as.data.frame(matrix(exp(rnorm(K*100)), nrow=100, ncol=K))
lX =  log(X)
PC = eigen(cov(lX- rowMeans(lX)))$vectors[,-K]
lambda = eigen(cov(lX- rowMeans(lX)))$values[-K]

M = cov(log(X))
# Extract the basis
SOL = findPB(list(M=M, nodes=as.list(1:nrow(M))))
b1 = sprintf("%s ~ %s",
             paste(names(X)[SOL$L], collapse='+'),
             paste(names(X)[SOL$R], collapse='+'))
PB1 = sbp_basis(b1 = eval(parse(text=b1)), data=X)[,1]
sapply(1:(K-1), function(i) sum(PB1*PC[,i]))

var(coordinates(X, sbp_basis(b1 = eval(parse(text=b1)), data=X))[,1])


var(coordinates(X, PC)[,1])
