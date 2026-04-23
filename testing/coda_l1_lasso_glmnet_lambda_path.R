library(glmnet)
library(genlasso)
library(ggplot2)
library(coda.base)
library(patchwork)

X = as.matrix(waste[,5:9])
y = waste$floating_population

LX = log(X)
K = ncol(X)
n = nrow(X)

#
# fit = glm(y~0+cbind(1,LX-rowMeans(LX)))
# b1 = as.vector(tidyr::replace_na(coef(fit)[-1], 0))
# b1 = b1 - mean(b1)
# logLik(fit)  # 'log Lik.' 17.54899 (df=6)
#
fit = glm(y~0+cbind(1,LX %*% ilr_basis(LX)))
b0 = as.vector(ilr_basis(LX) %*% coef(fit)[-1])
logLik(fit)  # 'log Lik.' 17.54899 (df=6)

B = ilr_basis(LX)
H = LX %*% B

X1 = cbind(1, H, 0)
D = cbind(0, B, 1)
fit1 = genlasso(y, X1, D, svd = TRUE, maxsteps = 1e7, rtol = 1e-15, btol = 1e-15)
as.vector(B %*% coef(fit1, lambda = 0)$beta[2:5,])
b0

b1 = coef(fit1, lambda = 0.08)$beta
m = b1[6]
clr_beta_genlasso = as.vector(B %*% b1[2:5])



CLR = LX - rowMeans(LX)

fit2 = glmnet(
  CLR, y,
  intercept = TRUE,
  standardize = FALSE,
  thresh = 1e-20,
  maxit = 1e7
)



#######
lp_genlasso = function(fit){
  dplot = do.call(rbind, lapply(seq_along(fit$lambda),
                                function(i){
                                  data.frame(
                                    lambda = fit$lambda[i],
                                    var = colnames(LX),
                                    coef = as.vector(fit$beta[-c(1, K+1),i]  %*% t(B)),
                                    m = as.vector(fit$beta[K+1, i])
                                  )
                                }))
  p1 = ggplot(data=dplot) +
    geom_path(aes(x = lambda, y = coef, col = var))  +
    theme_minimal()
  p2 = ggplot(data=dplot) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_path(aes(x = lambda, y = m))  +
    theme_minimal()
  p1  / p2 +
    plot_layout(design = "AAAA
                          AAAA
                          AAAA
                          AAAA
                          BBBB")
}
lp_glmnet = function(fit){
  dplot = do.call(rbind, lapply(seq_along(fit$lambda),
                                function(i){
                                  data.frame(
                                    lambda = fit$lambda[i] * fit$nobs,
                                    var = colnames(LX),
                                    coef = as.vector(fit$beta[,i] - mean(fit$beta[,i])),
                                    m = mean(fit$beta[,i])
                                  )
                                }))
  p1 = ggplot(data=dplot) +
    geom_path(aes(x = lambda, y = coef, col = var))  +
    theme_minimal()
  p2 = ggplot(data=dplot) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_path(aes(x = lambda, y = m))  +
    theme_minimal()
  p1  / p2 +
    plot_layout(design = "AAAA
                          AAAA
                          AAAA
                          AAAA
                          BBBB")
}
lp_genlasso(fit1) | lp_glmnet(fit2)

y = 1*(waste$floating_population_cat == "+")
lp_glmnet(glmnet(
  CLR, y, family = 'binomial',
  intercept = TRUE,
  standardize = FALSE,
  thresh = 1e-20,
  maxit = 1e7
))
