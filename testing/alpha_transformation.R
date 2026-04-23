library(coda.base)
data(package='coda.base')
X = pollen[1:3] |> as.matrix()
dim(X)
a = 0.5
alpha_tranform = function(X, a, H = ilr_basis(ncol(X))){
  Xa = X^a
  Xa = Xa/rowSums(Xa)
  ( ncol(X) * Xa -1 ) %*% H / a
}
head(alpha_tranform(X, -1))
head(X)
head(coord(X))
head(coord(X^a)/a)
library(coda.count)
data(package='coda.count')
K = 100
X = coda.count::rdm(size = rep(K,10), alpha = c(0.2,0.3,0.4))
plot(alpha_tranform(X, 0.01))
