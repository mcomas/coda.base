lasso_l1_coda_genlasso = function(X, y){
  require(genlasso)
  LX = log(X)
  K = ncol(X)
  n = nrow(X)

  B = ilr_basis(LX)
  H = LX %*% B

  X1 = cbind(1, H, 0)
  D = cbind(0, B, 1)
  fit = genlasso(
    y, X1, D,
    svd = TRUE,
    maxsteps = 1e7,
    eps = 1e-10,
    rtol = 1e-15,
    btol = 1e-15)

  do.call(rbind,
          lapply(seq_along(fit$lambda),
                 function(i){
                   data.frame(
                     lambda = fit$lambda[i],
                     var = colnames(LX),
                     coef = as.vector(fit$beta[-c(1, K+1),i]  %*% t(B)),
                     m = as.vector(fit$beta[K+1, i])
                   )
                 }))


}
lasso_l1_coda_genlasso_2 = function(X, y){
  require(genlasso)
  LX = log(X)
  K = ncol(X)
  n = nrow(X)

  CLR = LX - rowMeans(LX)

  X1 = cbind(1, CLR, 0)
  D = cbind(0, diag(K), 1)
  fit = genlasso(
    y, X1, D,
    svd = TRUE,
    maxsteps = 1e7,
    rtol = 1e-15,
    btol = 1e-15)

  do.call(rbind,
          lapply(seq_along(fit$lambda),
                 function(i){
                   data.frame(
                     lambda = fit$lambda[i],
                     var = colnames(LX),
                     coef = as.vector(fit$beta[-c(1, K+2),i]),
                     m = as.vector(fit$beta[K+2, i])
                   )
                 }))


}
lasso_l1_coda_glmnet = function(X, y, ...){
  require(glmnet)
  LX = log(X)
  K = ncol(X)
  n = nrow(X)

  CLR = LX - rowMeans(LX)

  fit = glmnet(
    CLR, y,
    intercept = TRUE,
    standardize = FALSE,
    thresh = 1e-20,
    maxit = 1e7,
    ...
  )

  do.call(rbind,
          lapply(seq_along(fit$lambda),
                 function(i){
                   data.frame(
                     lambda = fit$lambda[i] * fit$nobs,
                     var = colnames(LX),
                     coef = as.vector(fit$beta[,i] - mean(fit$beta[,i])),
                     m = mean(fit$beta[,i])
                   )
                 }))
}

plot_coda_lasso = function(dplot, m = TRUE){
  require(ggplot2)
  require(patchwork)
  require(ggnuplot)
  p1 = ggplot(data=dplot) +
    geom_path(aes(x = lambda, y = coef, col = var))  +
    scale_x_continuous(trans = 'log') +
    labs(y='clr(b)', col = '') +
    theme_gnuplot()
  if(m == FALSE) return(p1)
  p2 = ggplot(data=dplot) +
    geom_hline(yintercept = 0, linetype = 'dotted') +
    geom_path(aes(x = lambda, y = m))  +
    scale_x_continuous(trans = 'log') +
    labs(y='Asymmetry', col = '') +
    theme_gnuplot()

  p1  / p2 +
    plot_layout(design = "AAAA
                          AAAA
                          AAAA
                          AAAA
                          BBBB")

}
