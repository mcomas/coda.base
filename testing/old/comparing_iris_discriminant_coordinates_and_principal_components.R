X = iris[,1:4]
K = ncol(X)
gr = iris$Species

library(fpc)
dcf <- discrcoord(X,gr)
plot(dcf$proj,col=gr)

nG = length(unique(gr))

X.split = split(X, gr)
grMeans = t(sapply(1:nG, function(i){
  x = X.split[[i]]
  colMeans(x)
}))
grCovs = lapply(1:nG, function(i){
  x = X.split[[i]]
  cov(x) * (nrow(x)-1)
})
W = Reduce(`+`, grCovs)
B = cov(X) * (nrow(X)-6)

library(mvtnorm)
POST = sapply(1:nG, function(i){
  dmvnorm(X, grMeans[i,], B)
})

library(coda.base)
PC.coord = coordinates(POST, 'pc')
plot(-PC.coord[,1], -PC.coord[,2], col=gr)
apply(PC.coord, 2, var)
attr(PC.coord, 'basis')

PB = pb_basis(POST, method = 'exact')
PB.coord = coordinates(POST, PB)
plot(PB.coord[,1], PB.coord[,2], col=gr)
PB

##
##
##

PB = pb_basis(POST, method='exact')
sign(PB)
PB.coord = coordinates(POST, basis = PB)
plot(PB.coord[,1:2], col=gr)
legend('bottomright', legend =levels(gr), pch=1, col=1:3)
apply(PB.coord, 2, var)
