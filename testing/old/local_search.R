library(magrittr)
D = 5
n = 100
#set.seed(34)
Y = matrix(exp(rnorm(D*n)), ncol = D)
microbenchmark::microbenchmark(
  find_principal_balance(Y),
  pb_basis(Y, method = 'exact'), times = 1)

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
b2[,1]
var(coordinates(Y, b2))
(I <- which.max(abs(b2-b1)))

b1[I,] = b2[I,]
suppressWarnings(b1 <- sbp_basis(sign(b1)))
b1[,1]
b2 = M %*% cbind(b1)
b2 = b2 / sqrt(sum(b2^2))
b2[,1]
var(coordinates(Y, b2))
(I <- which.max(abs(b2-b1)))




pb_basis(Y, method = 'exact')[,1:2]
find_all_principal_balance(Y)



suppressWarnings(B <- cbind(find_testing(Y, v),
                            find_principal_balance(Y)) %>% sbp_basis())
apply(coordinates(Y, B), 2, var)

PC1 = find_principal_balance(Y)
candidate1 = find_principal_balance(Y[,which(PC1 == -1)])




b[,1]
which(b[,1]==1)-1
which(b[,1]==-1)-1
suppressWarnings(var(coordinates(Y, sbp_basis(b))))
b[5] = 0
suppressWarnings(var(coordinates(Y, sbp_basis(b))))

