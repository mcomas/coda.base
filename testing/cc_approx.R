D1 = 3
D2 = 4
n = 10
X = matrix(exp(rnorm(D1*n)),ncol=D1)
Y = matrix(exp(rnorm(D2*n)),ncol=D2)

S11 = cov(X)
S12 = cov(X, Y)
S22 = cov(Y)
M = chol2inv(chol(S11)) %*% S12 %*% chol2inv(chol(S22)) %*% t(S12)

eigen(M)$vectors / cancor(X,Y)$xcoef
