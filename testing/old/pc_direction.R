library(magrittr)
D = 5
n = 100
set.seed(3)
Y = matrix(exp(rnorm(D*n)), ncol = D)
pb_basis(Y, method = 'exact')

M = cov(coordinates(Y, 'clr'))
v = svd(M, nu = 1, nv = 0)$u
b1 = rep(0, D)
b1[which.max(v)] = +1
b1[which.min(v)] = -1
suppressWarnings(b1 <- sbp_basis(cbind(b1)))
b1[,1]
var(coordinates(Y, b1))
b2 = M %*% cbind(b1)
b2 = b2 / sqrt(sum(b2^2))
# var(coordinates(Y, b2))
(I <- which.max(abs(b2-b1)))

b1[I,] = b2[I,]
suppressWarnings(b1 <- sbp_basis(sign(b1)))
b1[,1]
var(coordinates(Y, b1))
b2 = M %*% cbind(b1)
b2 = b2 / sqrt(sum(b2^2))
#var(coordinates(Y, b2))
(I <- which.max(abs(b2-b1)))


b1[I,] = 0
suppressWarnings(b1 <- sbp_basis(sign(b1)))
b1[,1]
var(coordinates(Y, b1))   ## stop here
b2 = M %*% cbind(b1)
b2 = b2 / sqrt(sum(b2^2))
#var(coordinates(Y, b2))
(I <- which.max(abs(b2-b1)))
