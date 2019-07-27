D1 = 3
D2 = 4
n = 10
X = matrix(exp(rnorm(D1*n)),ncol=D1)
Y = matrix(rnorm(D2*n),ncol=D2)

S11 = cov(X)
S12 = cov(X, Y)
S22 = cov(Y)
M = chol2inv(chol(S11)) %*% S12 %*% chol2inv(chol(S22)) %*% t(S12)

a = suppressWarnings(coda.base::sbp_basis(cbind(c(-1,+1,0))))
b = (chol2inv(chol(S22)) %*% t(S12) %*% a)
#b = b / sqrt(sum(b^2))
(t(a) %*% S12 %*% b)^2 / ((t(a) %*% S11 %*% a) * (t(b) %*% S22 %*% b))

eigen(M)$vectors / cancor(X,Y)$xcoef
