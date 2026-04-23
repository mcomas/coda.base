library(coda.base)
library(ggplot2)
library(patchwork)
library(ggnuplot)

X = as.matrix(waste[,5:9])
y = waste$floating_population

source('testing/coda_lasso_functions.R')

d_l1_coda_coefs_genlasso = lasso_l1_coda_genlasso(X, y)
P1 = plot_coda_lasso(d_l1_coda_coefs_genlasso)

d_l1_coda_coefs_glmnet = lasso_l1_coda_glmnet(X, y)
P2 = plot_coda_lasso(d_l1_coda_coefs_glmnet)

P1 | P2

y_bin = waste$floating_population_cat == "+"
d_binomial = lasso_l1_coda_glmnet(X, y_bin, family = 'binomial')
plot_coda_lasso(d_binomial)
