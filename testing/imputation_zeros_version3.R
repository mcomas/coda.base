library(coda.base)

# nth = function(x, i) x[[i]]
# X = as.matrix(iris[,1:4])
# i1 = rbinom(length(X), 1, 0.15) == 1
# i2 = rbinom(length(X), 1, 0.15) == 1
# X[i1] = 0
# X[i2] = NA
# X
# DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
# DL[!is.na(X) & X == 0] = rlnorm(length(X), meanlog = log(0.01))[!is.na(X) & X == 0]
# DL
gen_X_with_zeros_and_missings = function(d, dl_par = 0.05, na_p = 0.15){
  gen_X = function(d = 5){
    H = mvtnorm::rmvnorm(300, rep(0, d))
    S = cov(H)
    EIG = eigen(S)

    as.matrix(composition(H %*% EIG$vectors, 'cdp'))
  }

  X0 = gen_X(d)
  X = X0
  below_dl = X < dl_par
  sum(below_dl)
  X[below_dl] = 0

  i2 = rbinom(length(X), 1, na_p) == 1
  X[i2] = NA

  DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
  DL[!is.na(X) & below_dl] = dl_par
  list(X = X, DL = DL, X0 = X0)
}
set.seed(0)
l_gen = gen_X_with_zeros_and_missings(5)
X = l_gen$X
DL = l_gen$DL

source('testing/imputation_version2_functions.R')


coda_zero_replace = function(X, DL, dl_prop = 0.65, eps = 1e-4, debug = FALSE, parameters = FALSE){
  tX = t(X)
  tDL = t(DL)

  N = ncol(tX)
  D = nrow(tX)
  d = D - 1

  tX_na = is.na(tX)
  tX_zero = !tX_na & tX == 0
  n_na = colSums(tX_na)
  n_zero = colSums(tX_zero)
  n_positive = D - n_na - n_zero


  tXna_imp = tX
  tXna_imp[tX_zero] = tDL[tX_zero] * dl_prop

  lX = unname(as.matrix(log(tXna_imp)))
  lDL = log(tDL)

  Bc = zero_na_conditional_obasis(tX)
  Bc_inv = array(0, dim = dim(Bc))
  for(i in 1:N){
    Bc_inv[,,i] = t(MASS::ginv(Bc[,,i]))
  }

  G = diag(D) - 1/D * matrix(1, D, D)
  MU = G %*% rowMeans(lX, na.rm = TRUE)
  strict_positive = function(X){
    eig = eigen(X)
    eig$vectors %*% diag(pmax(1e-6, eig$values - min(eig$values))) %*% t(eig$vectors)
  }
  SIGMA = strict_positive(G %*% (cov(t(lX), use = "pairwise.complete.obs") * (N-1)/N) %*% G)

  CONT = TRUE
  it = 0
  while(CONT){
    it = it + 1
    M1 = matrix(0, ncol = N, nrow = D)
    M2 = matrix(0, ncol = D, nrow = D)
    for(i in 1:N){
      lx = lX[,i,drop=FALSE]

      if(D - n_positive[i] > 0){ # Missing or zero
        if(n_zero[i] == 0){ # Only missing
          if(n_na[i] >= d){ # All relations missing
            m1 = MU
            m2 = m1 %*% t(m1) + SIGMA
          }else{
            is_na = tX_na[,i]
            Bi = Bc[,,i]
            Bi_inv = Bc_inv[,,i]
            h2 = Bi[-(1:n_na[i]),!is_na,drop=FALSE] %*% lx[!is_na,,drop=FALSE]
            mu_i = Bi %*% MU
            sigma_i = Bi %*% SIGMA %*% t(Bi)

            aug_ilr = rbind(cbind(sigma_i, mu_i),
                            cbind(t(mu_i), 1))
            # Conditioning on last coordinates (associated to positive parts)
            swept = SWEEP(aug_ilr, (1+n_na[i]):d)

            h1_m1 = crossprod(swept[(1+n_na[i]):D, 1:n_na[i], drop=FALSE], rbind(h2, 1))
            h1_m2 = h1_m1 %*% t(h1_m1) + swept[1:n_na[i], 1:n_na[i], drop=FALSE]

            m1 = t(Bi_inv) %*% rbind(h1_m1, h2)
            m2 = t(Bi_inv) %*% rbind(
              cbind(h1_m2, h1_m1 %*% t(h2)),
              cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi_inv
          }
        }else{ # With zeros
          if(n_zero[i] == 1 & n_na[i] == 0){ # fully determined
            m1 = lx - mean(lx)
            m2 = m1 %*% t(m1)
          }else{
            is_na = tX_na[,i]
            is_zero = tX_zero[,i]
            Bi = Bc[,,i]
            Bi_inv = Bc_inv[,,i]
            b = rep(0, D)
            b[is_zero] = lx[is_zero]
            b[!is_na & !is_zero] = lx[!is_na & !is_zero]
            h2_dl = crossprod(Bi[n_na[i] + n_zero[i],], b)
            if(n_positive[i] > 1){
              h2_pos = Bi[(1+n_na[i] + n_zero[i]):d,!is_na & !is_zero,drop=FALSE] %*% lx[!is_na & !is_zero,,drop=FALSE]
              h2 = rbind(h2_dl, h2_pos)
            }else{
              h2 = h2_dl
            }

            mu_i = Bi %*% MU
            sigma_i = Bi %*% SIGMA %*% t(Bi)

            aug_ilr = rbind(cbind(sigma_i, mu_i),
                            cbind(t(mu_i), 1))

            # Conditioning on last coordinates (associated to positive parts
            # and relationa positive and zero parts)
            swept = SWEEP(aug_ilr, (n_na[i]+n_zero[i]):d)

            h1_m1 = crossprod(swept[(n_na[i]+n_zero[i]):D, 1:(n_na[i]+n_zero[i]-1), drop=FALSE], rbind(h2, 1))
            h1_m2 = h1_m1 %*% t(h1_m1) + swept[1:(n_na[i]+n_zero[i]-1), 1:(n_na[i]+n_zero[i]-1), drop=FALSE]

            m1 = t(Bi_inv) %*% rbind(h1_m1, h2)
            m2 = t(Bi_inv) %*% rbind(
              cbind(h1_m2, h1_m1 %*% t(h2)),
              cbind(h2 %*% t(h1_m1), h2 %*% t(h2)))  %*% Bi_inv
          }
        }
      }else{ # Complete case
        m1 = lx - mean(lx)
        m2 = m1 %*% t(m1)
      }
      M1[,i] = m1
      M2 = M2 + m2
    }
    MU_next    = matrix(rowMeans(M1), ncol = 1)
    SIGMA = M2 / N - MU_next %*% t(MU_next)
    #print(as.vector(MU_next))
    if( max(abs(MU-MU_next)) < eps) CONT = FALSE
    MU = MU_next
    # logLik_next = mLogLik(X, MU[,1], SIGMA)
    # if(debug) cat(sprintf("Iteration: %d, LogLik: %0.3f\n", it, logLik_next))
    # if( abs(logLik_next - logLik_current)/abs(logLik_current) < EPS ) CONT = FALSE
    #logLik_current = logLik_next
  }
  if(parameters) return(list(clr_mu = MU[,1], clr_sigma = SIGMA, clr_h = t(M1)))
  t(apply(exp(M1), 2, prop.table))
}
Xr = X
Xr[!is.na(X) & X == 0] = -DL[!is.na(X) & X == 0]
P = coda_zero_replace(X,DL)
head(Xr)
head(P  * apply(X/P, 1, max, na.rm=TRUE))
