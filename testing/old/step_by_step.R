library(coda.base)
set.seed(2)
X = matrix(rlnorm(10*7), ncol = 7)
lX = log(X)
M = cov(lX)
sourceCpp('src/balance2_testing.cpp')
A5 = testing_05(X)
A6 = testing_06(X)

B1 = ilr_basis(7)
B = B1

D1 = B1 %*% eigen( t(B1) %*% M %*% B1, symmetric = TRUE)$vectors[,1,drop=FALSE]
D1_alt = B %*% eigen( t(B) %*% M %*% B, symmetric = TRUE)$vectors[,1,drop=FALSE]
D1 == D1_alt
abs(D1) - abs(A5[[2]][,1,drop=FALSE])  # Approx

Bal1 = A5[[1]][,1,drop=FALSE]
Bal1

H_restrict = t(B1) %*% cbind(Bal1)
B2 = qr.Q(qr(H_restrict, tol = 1e-25, LAPACK = FALSE), complete = TRUE)[,-1,drop=FALSE]
D2 = B1 %*% B2 %*% eigen( t(B1 %*% B2) %*% M %*% (B1 %*% B2))$vectors[,1,drop=FALSE]
B1 %*% B2 %*% eigen( cov(lX %*% (B1 %*% B2)) )$vectors[,1,drop=FALSE]

var(coordinates(X,D2))
var(coordinates(X, A6[[2]][,2,drop=FALSE]))

t(A6[[1]][,1,drop=FALSE]) %*% A6[[2]][,2,drop=FALSE]

Bal2 = A6[[1]][,2,drop=FALSE]
Bal2

H1_restrict = t(B1) %*% cbind(Bal1, Bal2)
H2_restrict = t(B1 %*% B2) %*% cbind(Bal2)


B3 = B1 %*% qr.Q(qr(H1_restrict), complete = TRUE)[,-(1:2),drop=FALSE]
D3 = B3 %*% eigen( t(B3) %*% M %*% (B3))$vectors[,1,drop=FALSE]

B3_alt = B1 %*% B2 %*% qr.Q(qr(H2_restrict), complete = TRUE)[,-1,drop=FALSE]
D3_alt = B3_alt %*% eigen( t(B3_alt) %*% M %*% (B3_alt))$vectors[,1,drop=FALSE]

var(coordinates(X, D3))
var(coordinates(X, D3_alt))
var(coordinates(X, A6[[2]][,3,drop=FALSE]))

cbind(A6[[1]][,3,drop=FALSE], A6[[2]][,3,drop=FALSE])
Bal3 = A6[[1]][,3,drop=FALSE]

H1_restrict = t(B1) %*% cbind(Bal1, Bal2, Bal3)

B4 = B1 %*% qr.Q(qr(H1_restrict), complete = TRUE)[,-(1:3),drop=FALSE]
D4 = B4 %*% eigen( t(B4) %*% M %*% (B4))$vectors[,1,drop=FALSE]

var(coordinates(X, D4))
var(coordinates(X, A6[[2]][,4,drop=FALSE]))
