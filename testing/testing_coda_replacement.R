library(coda.base)
gen_X_with_zeros_and_missings = function(n, d, missings = TRUE, zeros = TRUE,
                                         dl_par = 0.05, na_p = 0.15){
  gen_X = function(n, d){
    H = mvtnorm::rmvnorm(n, rep(0, d))
    S = cov(H)
    EIG = eigen(S)

    as.matrix(composition(H))
  }

  X0 = gen_X(n, d)
  X = X0
  below_dl = X < dl_par
  if(zeros){
    X[below_dl] = 0
  }
  if(missings){
    i2 = rbinom(length(X), 1, na_p) == 1
    X[i2] = NA
  }
  DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
  DL[!is.na(X) & below_dl] = dl_par
  list(X = X, DL = DL, X0 = X0)
}

set.seed(2)
N = 50
lR = gen_X_with_zeros_and_missings(N, 9, zeros = TRUE, na_p = 0.35, dl_par = 0.01)
X = lR$X
sum(!is.na(X) & X == 0)
DL = NULL
Bc = NULL
dl_prop = 0.65
xnz1  = coda_replacement(X, lR$DL)
library(zCompositions)
xnz2 = lrEMplus(cbind(X,1), dl = rep(0.01, 1+ncol(X)), ini.cov = "multRepl")
xnz2 = unname(as.matrix(xnz2[,1:ncol(X)]))

xnz3 = lrSVDplus(cbind(X,1), dl = rep(0.01, 1+ncol(X)), ini.cov = "multRepl")
xnz3 = unname(as.matrix(xnz3[,1:ncol(X)]))

library(coda.plot)
clr_biplot(as.data.frame(rbind(xnz1,
                               xnz2,
                               xnz3,
                               lR$X0)),
           group = rep(c('coda.base', 'lrEMplus', 'lrSVDplus', 'Original'), each = N),
           biplot_type = 'form')

library(vegan)
dist_mat0 = coordinates(lR$X0)
dist_mat1 = coordinates(xnz1)
dist_mat2 = coordinates(xnz2)
dist_mat3 = coordinates(xnz3)
v0 = as.vector(dist_mat0)
v1 = as.vector(dist_mat1)
v2 = as.vector(dist_mat2)
v3 = as.vector(dist_mat3)
sqrt(sum((v1 - v0)^2) / sum(v0^2))
sqrt(sum((v2 - v0)^2) / sum(v0^2))
sqrt(sum((v3 - v0)^2) / sum(v0^2))

