library(coda.base)
#
# nth = function(x, i) x[[i]]
# X = as.matrix(iris[,1:4])
# set.seed(6)
# i = rbinom(length(X), 1, 0.25) == 1
# X[i] = NA
# X


coda_impute = function(X, debug = FALSE, parameters = FALSE, eps = 1e-4){
  source('testing/imputation_version2_functions.R')
  lX = t(unname(as.matrix(log(X))))
  lX.is_na = is.na(lX)
  n.na = colSums(lX.is_na)

  N = ncol(lX)
  D = nrow(lX)
  d = D - 1
  G = diag(D) - 1/D * matrix(1, D, D)

  MU = G %*% rowMeans(lX, na.rm = TRUE)
  strict_positive = function(X){
    eig = eigen(X)
    eig$vectors %*% diag(pmax(1e-6, eig$values - min(eig$values))) %*% t(eig$vectors)
  }
  SIGMA = strict_positive(G %*% (cov(t(lX), use = "pairwise.complete.obs") * (N-1)/N) %*% G)
  # logLik_current = mLogLik(X, MU[,1], SIGMA)
  # cat(sprintf("Iteration: %d, LogLik: %0.3f\n", 0, logLik_current))

  Bc = conditional_obasis(!is.na(t(X)))
  CONT = TRUE
  it = 0
  while(CONT){
    it = it + 1
    #Lm = lapply(1:nrow(X), function(i){
    M1 = matrix(0, ncol = N, nrow = D)
    M2 = matrix(0, ncol = D, nrow = D)
    for(i in 1:N){
      lx = lX[,i,drop=FALSE]

      if(n.na[i] >= d){ # all relations missing
        m1 = MU
        m2 = m1 %*% t(m1) + SIGMA
      }else{
        if(n.na[i] == 0){
          m1 = lx - mean(lx)
          m2 = m1 %*% t(m1)
        }else{
          is_na = is.na(lx)
          Bi = Bc[,,i]
          h2 = Bi[-(1:n.na[i]),!is_na,drop=FALSE] %*% lx[!is_na,,drop=FALSE]
          mu_i = Bi %*% MU
          sigma_i = Bi %*% SIGMA %*% t(Bi)
          aug_ilr = rbind(cbind(sigma_i, mu_i),
                          cbind(t(mu_i), 1))
          swept = SWEEP(aug_ilr, (1+n.na[i]):d)

          h1_m1 = crossprod(swept[(1+n.na[i]):D, 1:n.na[i], drop=FALSE], rbind(h2, 1))
          h1_m2 = h1_m1 %*% t(h1_m1) + swept[1:n.na[i], 1:n.na[i], drop=FALSE]

          m1 = t(Bi) %*% rbind(h1_m1, h2)
          m2 = t(Bi) %*% rbind(
            cbind(h1_m2, h1_m1 %*% t(h2)),
            cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi

        }
      }
      M1[,i] = m1
      M2 = M2 + m2
    }
    MU_new    = matrix(rowMeans(M1), ncol = 1)
    SIGMA = strict_positive(M2 / N - MU_new %*% t(MU_new))
    if(debug){
      logLik_current = mLogLik(X, MU_new[,1], SIGMA)
      cat(sprintf("Iteration: %d, LogLik: %0.3f\n", it, logLik_current))
      #print(as.vector(MU))
    }
    if( max(abs(MU_new - MU)/abs(MU)) < eps ) CONT = FALSE
    MU = MU_new
  }
  MU = MU_new
  if(parameters) return(list(clr_mu = MU[,1], clr_sigma = SIGMA, clr_h = t(M1)))
  t(apply(exp(M1), 2, prop.table))
}
coda_impute_with_zeros = function(X, DL = NULL, dl_prop = 0.65,
                                  debug = FALSE, parameters = FALSE, eps = 1e-4){
  tX = t(unname(as.matrix(X)))
  if(is.null(DL)){
    DL = apply(tX, 1, function(x) min(x[!is.na(x) & x > 0]))
  }
  if(is.vector(DL)){
    tDL = matrix(DL, nrow = nrow(tX), ncol = ncol(tX))
  }
  if(is.data.frame(DL) | is.matrix(DL)){
    tDL = t(as.matrix(DL))
  }
  if(!exists('tDL')){
    stop("Detection limit parameter (DL) must be a vector, a matrix or unset")
  }
  lX = log(tX)
  lX.is_na = is.na(lX)
  lX.is_zero = !is.na(lX) & is.infinite(lX)

  lX[lX.is_zero] = log(tDL[lX.is_zero]) + log(dl_prop)

  n.na = colSums(lX.is_na)
  n.zero = colSums(lX.is_zero)

  n.na_zero = n.na + pmax(n.zero - 1, 0)

  N = ncol(lX)
  D = nrow(lX)
  d = D - 1
  G = diag(D) - 1/D * matrix(1, D, D)

  MU = G %*% rowMeans(lX, na.rm = TRUE)
  strict_positive = function(X){
    eig = eigen(X)
    eig$vectors %*% diag(pmax(1e-6, eig$values - min(eig$values))) %*% t(eig$vectors)
  }
  SIGMA = strict_positive(G %*% (cov(t(lX), use = "pairwise.complete.obs") * (N-1)/N) %*% G)
  # logLik_current = mLogLik(X, MU[,1], SIGMA)
  # cat(sprintf("Iteration: %d, LogLik: %0.3f\n", 0, logLik_current))

  Bc = zero_na_conditional_obasis(tX)
  Bc_inv = array(0, dim = dim(Bc))
  for(i in 1:N){
    Bc_inv[,,i] = t(pinv(Bc[,,i])) #t(MASS::ginv(Bc[,,i]))
  }
  CONT = TRUE
  it = 0
  while(CONT){
    it = it + 1
    #Lm = lapply(1:nrow(X), function(i){
    M1 = matrix(0, ncol = N, nrow = D)
    M2 = matrix(0, ncol = D, nrow = D)
    for(i in 1:N){
      lx = lX[,i,drop=FALSE]

      if(n.na_zero[i] >= d){ # all relations missing
        m1 = MU
        m2 = m1 %*% t(m1) + SIGMA
      }else{
        if(n.na_zero[i] == 0){
          m1 = lx - mean(lx)
          m2 = m1 %*% t(m1)
        }else{
          is_na = is.na(lx)
          Bi = Bc[,,i]
          Bi_inv = Bc_inv[,,i]
          # Bi[,!is_na,drop=FALSE] %*% lx[!is_na,,drop=FALSE]
          # Bi[-(1:n.na_zero[i]),!is_na,drop=FALSE] %*% lx[!is_na,,drop=FALSE]
          h2 = Bi[-(1:n.na_zero[i]),!is_na,drop=FALSE] %*% lx[!is_na,,drop=FALSE]
          mu_i = Bi %*% MU
          sigma_i = Bi %*% SIGMA %*% t(Bi)
          aug_ilr = rbind(cbind(sigma_i, mu_i),
                          cbind(t(mu_i), 1))
          swept = SWEEP(aug_ilr, (1+n.na_zero[i]):d)

          h1_m1 = crossprod(swept[(1+n.na_zero[i]):D, 1:n.na_zero[i], drop=FALSE], rbind(h2, 1))
          h1_m2 = h1_m1 %*% t(h1_m1) + swept[1:n.na_zero[i], 1:n.na_zero[i], drop=FALSE]

          m1 = t(Bi_inv) %*% rbind(h1_m1, h2)
          m2 = t(Bi_inv) %*% rbind(
            cbind(h1_m2, h1_m1 %*% t(h2)),
            cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi_inv

        }
      }
      M1[,i] = m1
      M2 = M2 + m2
    }
    MU_new    = matrix(rowMeans(M1), ncol = 1)
    SIGMA = strict_positive(M2 / N - MU_new %*% t(MU_new))
    if(debug){
      # logLik_current = mLogLik(X, MU_new[,1], SIGMA)
      logLik_zr_current = mLogLik_zr(X, DL, dl_prop, MU_new[,1], SIGMA, Bc)
      cat(sprintf("Iteration: %d, LogLik: %0.3f \n", it, logLik_zr_current))
      #print(as.vector(MU))
    }
    if( max(abs(MU_new - MU)/abs(MU)) < eps ) CONT = FALSE
    MU = MU_new
  }
  MU = MU_new
  if(parameters) return(list(clr_mu = MU[,1], clr_sigma = SIGMA, clr_h = t(M1)))
  t(apply(exp(M1), 2, prop.table))
}
# P = coda_impute(X, debug = TRUE)
# limp = coda_impute(X, debug = TRUE, parameters = TRUE)
# mLogLik(X, limp$clr_mu, limp$clr_sigma)
# X / composition(limp$clr_h, 'clr')
