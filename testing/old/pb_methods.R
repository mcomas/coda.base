library(magrittr)
X = parliament2017[,-1]

pb_basis(X, method = 'exact')
source('testing/PBArticle.R')
fBalChip(X)
pb_basis(X, method = 'lsearch')

PC = pc_basis(X)
D = nrow(PC)
gamma1 = PC[,1]
alpha1 = rep(0, D)
imax = which.max(pc_basis(X)[,1])
imin = which.min(pc_basis(X)[,1])
alpha1[imax] = +1
alpha1[imin] = -1

gamma1[imax] = 0
gamma1[imin] = 0
ord = order(abs(gamma1), decreasing = TRUE)
alpha1_candidates = cbind(alpha1, sapply(1:(D-2), function(I){
  alpha1[ord[1:I]] = sign(gamma1[ord[1:I]])
  alpha1
})) %>% sbp_basis(silent = TRUE)

B1 = ilr_basis(8)

C1 = alpha1_candidates[,which.max(t(PC[,1,drop=FALSE]) %*% alpha1_candidates),drop=FALSE]
C1_H1 = t(B1) %*% C1


B2 = qr.Q(qr(C1_H1),complete=TRUE)[,-1,drop=FALSE]
H2 = coordinates(X, B1 %*% B2)
PC1_H2 = eigen(cov(H2))$vectors[,1,drop=FALSE]
gamma2 = B1 %*% B2 %*% PC1_H2

R2 = cbind(sign(C1), gamma2)
tapply(R2[,2], R2[,1] == 0, sum)
R2
R2[R2[,1]==+1,]
R2[R2[,1]==-1,]


alpha2 = rep(0, D)
imax = which.max(gamma2)
imin = which.min(gamma2)
alpha2[imax] = +1
alpha2[imin] = -1

gamma2[imax] = 0
gamma2[imin] = 0
ord = order(abs(gamma2), decreasing = TRUE)
alpha2_candidates = cbind(alpha2, sapply(1:(D-2), function(I){
  alpha2[ord[1:I]] = sign(gamma2[ord[1:I]])
  alpha2
})) %>% sbp_basis(silent = TRUE)
alpha2_candidates * as.vector(C1)
