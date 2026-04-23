library(tidyverse)
library(coda.base)

load("~/Data/GEMAS/gemas.RData")
source('testing/coda_lasso_functions.R')
y = gemas$pH
X = as.matrix(gemas[,2:12])


## coda-lasso
m1_coda = lasso_l1_coda_genlasso(X, y)
p1_coda = plot_coda_lasso(m1_coda)

m2_coda = lasso_l1_coda_genlasso_2(X, y)
p2_coda = plot_coda_lasso(m2_coda)

m3_coda = lasso_l1_coda_glmnet(X, y, lambda.min.ratio = 1e-10)
p3_coda = plot_coda_lasso(m3_coda)

p1_coda
p2_coda
p3_coda

## Tibshirani
Bpw = pairwise_basis(X)
m1_pw = glmnet(x = coordinates(X, Bpw), y = y, standardize = FALSE)

NW = 1000 * nrow(X)
X_aug = rbind(cbind(1,log(X)),
              t(c(0, rep(1/ncol(X), ncol(X)))))
y_aug = c(y, 0)
m2_pw = glmnet(x = X_aug, y = y_aug, intercept = FALSE,
            weights = c(rep(1, nrow(X)), NW),
            penalty.factor = c(0, rep(1, ncol(X))),
            standardize = FALSE)

dplot1_pw = do.call(
  rbind,
  lapply(seq_along(m1_pw$lambda),
         function(i){
           data.frame(
             lambda = m1_pw$lambda[i],
             var = colnames(X),
             coef = coordinates(composition(as.vector(m1_pw$beta[,i]), Bpw), 'clr')
           )
         }))
dplot2_pw = do.call(
  rbind,
  lapply(seq_along(m2_pw$lambda),
         function(i){
           data.frame(
             lambda = m2_pw$lambda[i],
             var = colnames(X),
             coef = m2_pw$beta[-1,i]
           )
         }))

plot_coda_lasso(dplot1_pw, m = FALSE)
plot_coda_lasso(dplot2_pw, m = FALSE)




B = ilr_basis(X)
m3 = genlasso(X = cbind(1, coord(X, 'clr') %*% B), y = y, D = cbind(0, t(pairwise_basis(X)) %*% B))
dplot3 = do.call(
  rbind,
  lapply(seq_along(m3$lambda),
         function(i){
           data.frame(
             lambda = m3$lambda[i],
             var = colnames(X),
             coef =  as.vector(B %*% m3$beta[-c(1),i])
           )
         }))
plot_coda_lasso(dplot3, m = FALSE)


m1 = genlasso(X = cbind(1, coord(X, 'clr')), y = y, D = cbind(0, t(pairwise_basis(X))))

dplot1 = do.call(
  rbind,
  lapply(seq_along(m1$lambda),
         function(i){
           data.frame(
             lambda = m1$lambda[i],
             var = colnames(X),
             coef = as.vector(m1$beta[-c(1),i])
           )
         }))

p1 = plot_coda_lasso(dplot1, m = FALSE)


Bpw = pairwise_basis(X)
m3 = glmnet(x = coordinates(X, Bpw), y = y)



dplot3 = do.call(
  rbind,
  lapply(seq_along(m3$lambda),
         function(i){
           data.frame(
             lambda = m3$lambda[i],
             var = colnames(X),
             coef = coordinates(composition(m3$beta[,i], Bpw), 'clr')
           )
         })
)
m3$beta
coordinates(composition(coef(m3, s = 0.2)[-1,1], Bpw), 'clr')



plot_coda_lasso(dplot3, m = FALSE)


m3 = genlasso(X = cbind(1, coord(X, 'clr')), y = y, D = cbind(0, diag(ncol(X))))
m4 = genlasso(X = cbind(1, coord(X, 'clr') %*% B), y = y, D = cbind(0, B))
dplot3 = do.call(
  rbind,
  lapply(seq_along(m3$lambda),
         function(i){
           data.frame(
             lambda = m3$lambda[i],
             var = colnames(X),
             coef = as.vector(m3$beta[-c(1),i])
           )
         }))
dplot4 = do.call(
  rbind,
  lapply(seq_along(m4$lambda),
         function(i){
           data.frame(
             lambda = m4$lambda[i],
             var = colnames(X),
             coef =  as.vector(B %*% m4$beta[-c(1),i])
           )
         }))
plot_coda_lasso(dplot3, m = FALSE)
plot_coda_lasso(dplot4, m = FALSE)
