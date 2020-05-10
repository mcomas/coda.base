load

load('testing/Example.RData')
Y = X[,1:10]
pc_basis(Y)
pb_basis(Y, method = 'exact')
source('testing/PBArticle.R')
fBalChip(Y)
data(Aar, package = 'compositions')
X = Aar[,c("Al2O3", "CaO", "Fe2O3t", "K2O", "MgO", "MnO", "Na2O", "P2O5", "SiO2", "TiO2")]

pc_basis(X)
pb_basis(X, method = 'exact')
t(fBalChip(X)$bal)
pb_basis(X, method = 'ward')

library(mvtnorm)
D = 300
n = 100
diag_ = rep(0.01, D-1)
diag_[1:10] = 0.9^(1:10)
Loadings = matrix(runif((D-1)*(D-1), min = -1, max = +1), nrow = D-1)
Scores = rmvnorm(n, sigma = diag(diag_))

X = composition(Scores %*% Loadings)
cumsum(apply(coordinates(X, pb_basis(X, method = 'ward.D')[,1:5]), 2, var))
cumsum(apply(coordinates(X, pb_basis(X, method = 'constrained')[,1:5]), 2, var))
cumsum(apply(coordinates(X, pb_basis(X, method = 'lsearch', rep = 1)[,1:5]), 2, var))
