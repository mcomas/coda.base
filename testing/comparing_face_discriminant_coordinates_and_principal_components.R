library(fpc)
set.seed(1)
face = rFace(600,dMoNo=2,dNoEy=0, p=8)
X = as.data.frame(face)
K = ncol(X)
gr = attr(face, 'grouping')


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

eig = eigen(chol2inv(chol(W)) %*% B)
#eig = eigen(solve(W) %*% B)
X2 = (face %*% eig$vectors)[, 1:2]
plot(-X2[,1], -X2[,2], col = gr)

library(mvtnorm)
POST = sapply(1:nG, function(i){
  dmvnorm(X, grMeans[i,], W)
})
logPOST = sapply(1:nG, function(i){
  dmvnorm(X, grMeans[i,], W, log = TRUE)
})
POST2 = exp(scale(logPOST, scale = FALSE))

library(coda.base)
PC.coord = coordinates(POST, 'pc')
plot(PC.coord[,1], PC.coord[,2], col=gr)
v = apply(PC.coord, 2, var)
cumsum(v) / sum(v)
attr(PC.coord, 'basis')

PC.coord2 = coordinates(POST2, 'pc')
plot(PC.coord2[,1], -PC.coord2[,2], col=gr)


PB1 = pb_basis(POST, method = 'exact')
PB1.coord = coordinates(POST, PB1)
plot(PB1.coord[,1], PB1.coord[,2], col=gr)
PB1

PB2 = find_principal_balance3_01(POST)
PB2.coord = coordinates(POST, PB2)
plot(PB2.coord[,1], PB2.coord[,2], col=gr)
PB2

##
##
##

PB = pb_basis(POST, method='exact')
sign(PB)
PB.coord = coordinates(POST, basis = PB)
plot(PB.coord[,1:2], col=gr)
legend('bottomright', legend =levels(gr), pch=1, col=1:3)

v = apply(PB.coord, 2, var)
vtotal = sum(v)
cumsum(v)/sum(v)


# pc1       pc2       pc3       pc4       pc5
# 0.7396288 0.9740265 0.9979777 0.9991101 1.0000000

# PRINCIPAL BALANCES
# x1        x2        x3        x4        x5
# 0.6711526 0.8300339 0.9834595 0.9988723 1.0000000

sign(PB)[,1:2]

table(gr)
M = sbp_basis(V2+V6~V1+V3+V4,
              V2~V5, data = X)
PB.coord = coordinates(POST, M)
plot(PB.coord, col = gr)
v = apply(PB.coord, 2, var)
cumsum(v)/vtotal


