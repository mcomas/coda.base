library(coda.base)

X = as.matrix(iris[,1:4])
set.seed(6)
i = rbinom(length(X), 1, 0.05) == 1
i = i & (1:length(X) > nrow(X))
X[i] = NA
X = X[apply(is.na(X), 1, sum) < 3,]
X
mLogLik = function(X, mu, sigma, B = ilr_basis(ncol(X))){
  MU = B %*% mu
  SIGMA = B %*% sigma %*% t(B)
  Bc = c_conditional_obasis(t(X), !is.na(t(X)))
  d = ncol(X) - 1
  lprobs = sapply(1:nrow(X), function(i){
    x = X[i,]
    is_na = is.na(x)
    n0 = sum(is_na)
    lprob = 0
    if(n0 < d){
      if(n0 == 0){
        lprob = mvtnorm::dmvnorm(as.vector( t(B) %*% log(x)), mu, sigma, log = TRUE)
      }else{
        Bi = Bc[,,i]
        B2 = Bi[-(1:n0),,drop=FALSE]
        h2 = as.vector(B2[,!is_na,drop=FALSE] %*% log(x[!is_na]))

        mu2 = B2 %*% MU
        sigma2 = B2 %*% SIGMA %*% t(B2)
        lprob = mvtnorm::dmvnorm(h2, mu2, sigma2, log = TRUE)
      }
    }
    lprob
  })
  sum(lprobs)
}
mLogLik(X, rep(0, 3), diag(3))


impute_coda_1 = function(X, B = ilr_basis(ncol(X)), eps = 1e-4){
  nth = function(x, i){
    x[[i]]
  }
  Bc = c_conditional_obasis(t(X), !is.na(t(X)))
  d = ncol(X) - 1
  MU = rep(0, d+1)
  SIGMA = diag(d+1)

  Xnz = X
  Xnz[is.na(Xnz)] = mean(X, na.rm=TRUE)

  MU_h = t(B) %*% MU
  SIGMA_h = t(B) %*% SIGMA %*% B

  H = log(Xnz) %*% B
  # logLik_current = sum(mvtnorm::dmvnorm(H, MU_h, SIGMA_h, log = TRUE))
  logLik_current = mLogLik(X, MU_h, SIGMA_h)
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


# mean(fitted(lm(ilr3~., data = H)))
# mean(lm.fit(x = cbind(1,H[,-d]), y = H[,d])$fitted.values)

impute_coda_2 = function(X, eps = 1e-4){
  iNA = is.na(X)

  m = colMeans(X, na.rm = TRUE)
  Xt = t(X)

  Xi = t(apply(ifelse(is.na(Xt), m, Xt), 2, prop.table))

  D = ncol(X)
  d = D-1

  ####
  H = coordinates(Xi)
  mom = cov.wt(H)
  logLik_current = sum(mvtnorm::dmvnorm(H, mom$center, mom$cov, log = TRUE))
  CONT = TRUE
  it = 0
  while(CONT){
    it = it + 1
    for(I in 1:D){

      Xi[,c(I,D)] = Xi[,c(D,I)]
      H = coordinates(Xi)

      fit_ = lm.fit(x = cbind(1,H[,-d]), y = H[,d]) #lm(H[,d]~.,data=data.frame(H[,-d])) #
      H[iNA[,I],d] = fit_$fitted.values[iNA[,I]]     #fitted(fit_)[iNA[,I]] #

      Xi = composition(H)
      Xi[,c(I,D)] = Xi[,c(D,I)]
    }

    H = coordinates(Xi)
    mom = cov.wt(H)
    logLik_next = sum(mvtnorm::dmvnorm(H, mom$center, mom$cov, log = TRUE))
    cat(sprintf("Iteration: %d, complete-LogLik: %0.3f\n", it, logLik_next))
    if( abs(logLik_next - logLik_current)/abs(logLik_current) < eps ) CONT = FALSE

    logLik_current = logLik_next
  }
  K = suppressWarnings(apply(X/Xi, 1, min, na.rm = TRUE))
  K.gmean = exp(mean(log(K[is.finite(K)])))
  K[!is.finite(K)] = K.gmean
  cat(sprintf("Iterations: %d\n", it))
  Xi * K
}
impute_coda_3 = function(X, eps = 1e-4){
  iNA = is.na(X)

  m = colMeans(X, na.rm = TRUE)
  Xt = t(X)

  lXi = log(t(apply(ifelse(is.na(Xt), m, Xt), 2, prop.table)))
  lXi = lXi - rowMeans(lXi)

  D = ncol(X)
  d = D-1

  B = ilr_basis(D)
  ####
  H = lXi %*% B
  mom = cov.wt(H)
  logLik_current = sum(mvtnorm::dmvnorm(H, mom$center, mom$cov, log = TRUE))
  CONT = TRUE
  it = 0
  while(CONT){
    it = it + 1
    for(I in 1:D){

      lXi[,c(I,D)] = lXi[,c(D,I)]
      H = lXi %*% B

      fit_ = lm.fit(x = cbind(1,H[,-d]), y = H[,d])  #lm(H[,d]~.,data=data.frame(H[,-d])) #
      H[iNA[,I],d] = fit_$fitted.values[iNA[,I]]     #fitted(fit_)[iNA[,I]] #

      lXi = H %*% t(B)
      lXi[,c(I,D)] = lXi[,c(D,I)]
    }

    H = lXi %*% B
    mom = cov.wt(H)
    logLik_next = sum(mvtnorm::dmvnorm(H, mom$center, mom$cov, log = TRUE))
    cat(sprintf("Iteration: %d, complete-LogLik: %0.3f\n", it, logLik_next))
    if( abs(logLik_next - logLik_current)/abs(logLik_current) < eps ) CONT = FALSE

    logLik_current = logLik_next
  }

  Xi = composition(lXi, 'clr')
  K = suppressWarnings(apply(X/Xi, 1, min, na.rm = TRUE))
  K.gmean = exp(mean(log(K[is.finite(K)])))
  K[!is.finite(K)] = K.gmean
  cat(sprintf("Iterations: %d\n", it))
  Xi * K
}
impute_coda_4 = function(X, B = ilr_basis(ncol(X)), eps = 1e-4){
  nth = function(x, i){
    x[[i]]
  }
  Bc = c_conditional_obasis(t(X), !is.na(t(X)))
  d = ncol(X) - 1
  MU = rep(0, d+1)
  SIGMA = diag(d+1)

  Xnz = X
  Xnz[is.na(Xnz)] = mean(X, na.rm=TRUE)

  MU_h = t(B) %*% MU
  SIGMA_h = t(B) %*% SIGMA %*% B

  H = log(Xnz) %*% B
  # logLik_current = sum(mvtnorm::dmvnorm(H, MU_h, SIGMA_h, log = TRUE))
  logLik_current = mLogLik(X, MU_h, SIGMA_h)

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
        m1 = t(B) %*% MU
        m2 = m1 %*% t(m1) + t(B) %*% SIGMA %*% B
      }else{
        if(n0 == 0){
          lx = log(x)
          m1 = t(B) %*% lx
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


          m1 = t(B) %*% t(Bi) %*% c(h1_m1, h2)
          m2 = t(B) %*% t(Bi) %*% rbind(
            cbind(h1_m2, h1_m1 %*% t(h2)),
            cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi %*% B
        }
      }
      list('m1' = m1,
           'm2' = m2)
    })
    MU_h    = Reduce(`+`, lapply(Lm, nth, 1)) / length(Lm)
    SIGMA_h = Reduce(`+`, lapply(Lm, nth, 2)) / length(Lm) - MU_h %*% t(MU_h)

    MU = B %*% MU_h
    SIGMA = B %*% SIGMA_h %*% t(B)
    H = t(do.call(cbind, lapply(Lm, nth, 1)))
    # logLik_next = sum(mvtnorm::dmvnorm(H, MU_h, SIGMA_h, log = TRUE))
    logLik_next = mLogLik(X, MU_h, SIGMA_h)

    cat(sprintf("Iteration: %d, LogLik: %0.3f\n", it, logLik_next))
    if( abs(logLik_next - logLik_current)/abs(logLik_current) < eps ) CONT = FALSE

    logLik_current = logLik_next
  }

  Ximp = composition(t(do.call(cbind, lapply(Lm, nth, 1))), B)
  K = suppressWarnings(apply(X/Ximp, 1, min, na.rm = TRUE))
  K.gmean = exp(mean(log(K[is.finite(K)])))
  K[!is.finite(K)] = K.gmean
  cat(sprintf("Iterations: %d\n", it))
  Ximp * K

}

iris.coda_1 = impute_coda_1(X, eps = 1e-4)
iris.coda_2 = impute_coda_2(X, eps = 1e-4)
iris.coda_3 = impute_coda_3(X, eps = 1e-4)
iris.coda_4 = impute_coda_4(X, eps = 1e-4)
library(zCompositions)
iris.zcomp = zCompositions::lrEM(X, label = NA, imp.missing = TRUE, ini.cov = 'complete.obs')
iris.robcomp = robCompositions::impCoda(X, method = 'lm', eps = 1e-4)

logLikNormImputation = function(X, Ximp){
  H = coordinates(Ximp)
  M = cov.wt(H, method = 'ML')
  mLogLik(X, M$center, M$cov)
}
logLikNormImputation(X, iris.coda_1)
logLikNormImputation(X, iris.coda_2)
logLikNormImputation(X, iris.coda_3)
logLikNormImputation(X, iris.coda_4)
logLikNormImputation(X, iris.zcomp)
logLikNormImputation(X, iris.robcomp$xImp)
sapply(list('iris' = iris[,1:4],
            'robComp' = iris.robcomp$xImp,
            'EM' = iris.coda_1,
            'zComp' = iris.zcomp),
       function(x) sum(eigen(cov(coordinates(x)))$values))

load("~/Research/COUNT-3WATER/dades.RData")
species = L_taxonomy$species
X = as.matrix(species[,-1])
X[X==0] = NA
X
Xna = impute_coda_1(X, eps = 0.01)
