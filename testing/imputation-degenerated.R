library(coda.base)
load("~/Research/COUNT-3WATER/dades.RData")
source('R/imputation-functions.R')
species = L_taxonomy$species
X = as.matrix(species[,-1])
X[X==0] = NA
X = matrix(1:24, ncol = 8)
X[sample(1:24, 6)] = NA
X



impute_coda_1 = function(X, eps = 1e-4){

  Bc = c_conditional_obasis(t(X), is.na(t(X)))

  Xclr = matrix(nrow = nrow(X), ncol = ncol(X))
  for(i in 1:nrow(X)){
    xi = X[i,]
    I1 = !is.na(xi)
    I0 = is.na(xi)
    d1 = sum(I1)-1
    D0 = sum(I0)
    H1 = Bc[1:d1,I1,i] %*% log(xi[I1])
    Xclr[i,] = c(H1,rep(0, D0)) %*% Bc[,,i]
  }
  SVD = svd(Xclr)
  dim(SVD$u %*% diag(SVD$d) %*% t(SVD$v))
  PC = prcomp(Xclr, center = FALSE)

  B0 = PC$rotation[,PC$sdev>1e-8]
  L_B0c = list()
  for(i in 1:nrow(X)){
    xi = X[i,]
    I1 = !is.na(xi)
    I0 = is.na(xi)
    d1 = sum(I1)-1
    D0 = sum(I0)

    L_B0c[[i]] = Bc[1:d1,,i] %*% B0
  }

  H0 = Xclr %*% B0

  MU0 = colMeans(H0)
  SIGMA0 = cov(H0)

  # logLik_current = sum(mvtnorm::dmvnorm(H, MU_h, SIGMA_h, log = TRUE))
  logLik_current = mLogLik(X, MU0, SIGMA0, B0c)
  CONT = TRUE
  it = 0
  EPS = 1e-6 * diag(d+1)
  while(CONT){
    it = it + 1
    Lm = lapply(1:nrow(X), function(i){
      # i=10  no missing
      x = X[i,]
      is_na = is.na(x)
      n0 = sum(is_na)
      if(n0 >= d){
        m1 = MU
        m2 = MU %*% t(MU) + SIGMA
      }else{
        if(n0 == 0){
          lx = log(x)
          m1 = lx - mean(lx)
          m2 = m1 %*% t(m1)
        }else{
          Bi = Bc[,,i]
          B1 = Bi[+(1:n0),,drop=FALSE]
          B2 = Bi[-(1:n0),,drop=FALSE]
          h2 = B2[,!is_na,drop=FALSE] %*% log(x[!is_na])

          mu1 = B1 %*% MU
          mu2 = B2 %*% MU
          sigma11 = B1 %*% SIGMA %*% t(B1)
          sigma22 = B2 %*% SIGMA %*% t(B2)
          sigma12 = B1 %*% SIGMA %*% t(B2)
          sigma21 = B2 %*% SIGMA %*% t(B1)

          h1_m1 = mu1 + sigma12 %*% MASS::ginv(sigma22) %*% (h2 - mu2)
          h1_m2 = h1_m1 %*% t(h1_m1) + sigma11 - sigma12 %*% MASS::ginv(sigma22) %*% sigma21


          m1 = t(Bi) %*% c(h1_m1, h2)
          m2 = t(Bi) %*% rbind(
            cbind(h1_m2, h1_m1 %*% t(h2)),
            cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi
        }
      }
      list('m1' = m1,
           'm2' = m2)
    })
    MU    = Reduce(`+`, lapply(Lm, nth, 1)) / length(Lm)
    SIGMA = Reduce(`+`, lapply(Lm, nth, 2)) / length(Lm) - MU %*% t(MU)

    MU_h = t(B) %*% MU
    SIGMA_h = t(B) %*% SIGMA %*% B
    H = t(do.call(cbind, lapply(Lm, nth, 1))) %*% B
    # logLik_next = sum(mvtnorm::dmvnorm(H, MU_h, SIGMA_h, log = TRUE))
    logLik_next = mLogLik(X, MU_h, SIGMA_h)

    cat(sprintf("Iteration: %d, LogLik: %0.3f\n", it, logLik_next))
    if( abs(logLik_next - logLik_current)/abs(logLik_current) < eps ) CONT = FALSE

    logLik_current = logLik_next
  }

  Ximp = composition(t(do.call(cbind, lapply(Lm, nth, 1))), 'clr')
  K = suppressWarnings(apply(X/Ximp, 1, min, na.rm = TRUE))
  K.gmean = exp(mean(log(K[is.finite(K)])))
  K[!is.finite(K)] = K.gmean
  Ximp * K

}
