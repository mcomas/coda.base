X = matrix(rlnorm(100*5), nrow=100)

PC1a = pc_basis(X)[,1,drop=FALSE]
PC2a = pc_basis(X)[,2,drop=FALSE]
PC3a = pc_basis(X)[,3,drop=FALSE]

B1 = ilr_basis(5)
H1 = coordinates(X, B1)
PC1_H1 = eigen(cov(H1))$vectors[,1,drop=FALSE]
PC1b = B1 %*% PC1_H1


B2 = qr.Q(qr(PC1_H1),complete=TRUE)[,-1,drop=FALSE]
H2 = coordinates(X, B1 %*% B2)
PC1_H2 = eigen(cov(H2))$vectors[,1,drop=FALSE]
PC2b = B1 %*% B2 %*% PC1_H2


B3 = qr.Q(qr(PC1_H2),complete=TRUE)[,-1,drop=FALSE]
H3 = coordinates(X, B1 %*% B2 %*% B3)
PC1_H3 = eigen(cov(H3))$vectors[,1,drop=FALSE]
PC3b = B1 %*% B2 %*% B3 %*% PC1_H3




