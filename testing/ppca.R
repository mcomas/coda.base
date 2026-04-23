X0 = as.matrix(iris[,1:4])
## Latent Factor Model
#
# p(x|mu, W, psi), W orthonormal, psi diagonal
# where x ~ N(mu, psi + W W')
#
## Probabilistic PCA
#
# p(x|mu, W, psi), W orthonormal, psi = sigma^2 Identity
# where x ~ N(mu, psi + W W')


fit1 <- psych::fa(X0, nfactors = 2, rotate = "varimax", fm="ml")
W1 = fit1$loadings
W1 = matrix(W1, nrow = nrow(W1))
psi1 = diag(fit1$uniquenesses)
mu = colMeans(X0)

library(mvtnorm)
loglik = dmvnorm(X0, mu, psi1 + W1 %*% t(W1), log = TRUE)
sum(loglik)
loglik_eps = dmvnorm(X0, mu, 0.01 * diag(4) + psi1 + W1 %*% t(W1), log = TRUE)
sum(loglik_eps)



fit2 = pcaMethods::ppca(X0, nPcs = 2, seed = 2, maxIterations = 1e5)
W2 = fit2@loadings
W2
fit2@varLimit
t(W2) %*% W2
em_ppca = function(X, L){
  i = 1
  I_L = diag(L)

  sigma_i = I_L + solve(t(W) %*% solve(psi) %*% W)
  m = sapply(1:nrow(X), function(i)
    sigma_i %*% (t(W) %*% solve(psi) %*% (X[i,] - mu)))
  b = rbind(m,1)
  W = t(b %*% X) %*% solve( rbind(
    cbind(nrow(X) * sigma_i, rowSums(m) ),
    cbind(t(rowSums(m)), 1)) )
  psi = diag((t(X) - W %*% b) %*% X) / nrow(X)


}

L = 2
threshold = 1e-05
maxIterations = 10000


# pcaMethods -> adapted to Kevin P. Murphy

pca_imp = function(X0){

  #set.seed(1)
  #########
  N = nrow(X0)
  D = ncol(X0)
  Obs <- !is.na(X0)
  hidden <- which(is.na(X0))
  missing <- length(hidden)
  if (missing) {
    X0[hidden] <- 0
  }

  ## inicialitza amb L valors agafats a l'atzar
  #W <- t(X0[sample(N)[1:L], , drop = FALSE])
  ## agafa valors a l'atzar seguint normal
  #W <- matrix(rnorm(W), nrow(W), ncol(W), dimnames = labels(W))
  W = matrix(rnorm(L*D), D, L, dimnames = list(colnames(X0),NULL))

  WtW <- t(W) %*% W
  X <- X0 %*% W %*% solve(WtW)   # E-step
  recon <- X %*% t(W)
  recon[hidden] <- 0
  ss <- sum(sum((recon - X0)^2))/(N * D - missing) ## error quadràtic
  count <- 1
  old <- Inf
  while (count > 0) {
    Sx <- solve(diag(L) + WtW/ss)
    ss_old <- ss
    if (missing) {
      proj <- X %*% t(W)
      X0[hidden] <- proj[hidden]
    }
    X <- X0 %*% W %*% Sx/ss   # equiv. Sx/ss = solve(WtW + diag(L)*ss)
    # First, since ss diag(L) is diagonal, the inversion in the e-step can be performed efficiently using the
    # matrix inversion lemma:
    # solve{WtW + ss diag(L)} = (diag(L)/ss - W solve(diag(L) + WtW/ss) t(W)/ss^2).
    # Second, since we
    # are only taking the trace of the matrix in the m-step, we do not need to compute the full sample
    # covariance XtX but instead can compute only the variance along each coordinate

    SumXtX <- t(X) %*% X
    W <- (t(X0) %*% X) %*% solve((SumXtX + N * Sx))   # M-step
    WtW <- t(W) %*% W
    ss <- (sum(sum((W %*% t(X) - t(X0))^2)) + N * sum(sum(WtW %*%
                                                            Sx)) + missing * ss_old)/(N * D)
    objective <- N * (D * log(ss) + sum(diag(Sx)) - log(det(Sx))) +
      sum(diag(SumXtX)) - missing * log(ss_old)
    rel_ch <- abs(1 - objective/old)
    old <- objective
    count <- count + 1
    if (rel_ch < threshold & count > 5) {
      print(count)
      count <- 0
    } else if (count > maxIterations) {
      count <- 0
      warning("stopped after max iterations, but rel_ch was > threshold")
    }
  }
  W <- orth(W)
  evs <- eigen(cov(Matrix %*% W))
  vals <- evs[[1]]
  vecs <- evs[[2]]
  W <- W %*% vecs
  X <- Matrix %*% W
  R2cum <- rep(NA, L)
  TSS <- sum(Matrix^2, na.rm = TRUE)
  for (i in 1:ncol(W)) {
    difference <- Matrix - (X[, 1:i, drop = FALSE] %*% t(W[,
                                                           1:i, drop = FALSE]))
    R2cum[i] <- 1 - (sum(difference^2, na.rm = TRUE)/TSS)
  }
  res <- new("pcaRes")
  res@scores <- X
  res@loadings <- W
  res@R2cum <- R2cum
  res@method <- "ppca"
  return(res)
}
