library(coda.base)

nth = function(x, i) x[[i]]
# X = as.matrix(iris[,1:4])
# i1 = rbinom(length(X), 1, 0.15) == 1
# i2 = rbinom(length(X), 1, 0.15) == 1
# X[i1] = 0
# X[i2] = NA
# X
# DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
# DL[!is.na(X) & X == 0] = rlnorm(length(X), meanlog = log(0.01))[!is.na(X) & X == 0]
# DL
gen_X = function(d = 5){
  H = mvtnorm::rmvnorm(300, rep(0, d))
  S = cov(H)
  EIG = eigen(S)

  PB1 = sample(c(-1,0,1), d+1, replace = TRUE)
  PB = sbp_basis(cbind(PB1), fill = TRUE)
  SBP_exact = sign(PB)

  as.matrix(composition(H %*% EIG$vectors, basis = PB))
}
set.seed(11)
X0 = gen_X()
X = X0
dl_par = 0.05
below_dl = X < dl_par
sum(below_dl)
X[below_dl] = 0

i2 = rbinom(length(X), 1, 0.15) == 1
X[i2] = NA

DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
DL[!is.na(X) & below_dl] = dl_par


source("testing/imputation_version2.R")
source('testing/imputation_version2_functions.R')
coda_zero_replace = function(X, DL, debug = FALSE, parameters = FALSE){
  Xna = X
  Xna[X == 0] = NA
  l_imp = coda_impute(Xna, debug = debug, parameters = TRUE)
  Xna_imp = exp(l_imp$clr_h) * apply(X / exp(l_imp$clr_h), 1, max, na.rm=TRUE)




  lX = t(unname(as.matrix(log(Xna_imp))))
  lDL = t(log(DL))
  lX.is_zero = is.finite(lDL)

  Bc = c_conditional_obasis(!lX.is_zero)

  n.zero = colSums(lX.is_zero)


  N = ncol(lX)
  D = nrow(lX)
  d = D - 1
  G = diag(D) - 1/D * matrix(1, D, D)

  MU = cbind(l_imp$clr_mu)
  SIGMA = l_imp$clr_sigma

  CONT = TRUE
  it = 0
  EPS = 1e-6
  while(CONT){
    it = it + 1
    M1 = matrix(0, ncol = N, nrow = D)
    M2 = matrix(0, ncol = D, nrow = D)
    for(i in 1:N){
      lx = lX[,i,drop=FALSE]

      if(n.zero[i] == D){ # all relations missing
        m1 = MU
        m2 = m1 %*% t(m1) + SIGMA
      }else{
        if(n.zero[i] == 0){
          m1 = lx - mean(lx)
          m2 = m1 %*% t(m1)
        }else{
          is_zero = lX.is_zero[,i]
          Bi = Bc[,,i]

          h2 = Bi %*% lx

          mu_i = Bi %*% MU
          sigma_i = Bi %*% SIGMA %*% t(Bi)

          aug_ilr = rbind(cbind(sigma_i, mu_i),
                          cbind(t(mu_i), 1))
          swept = SWEEP(aug_ilr, (1:d)[-n.zero[i]])

          m = crossprod(swept[(1:D)[-n.zero[i]], n.zero[i], drop=FALSE],
                        rbind(h2[-n.zero[i], , drop=FALSE], 1))
          v = swept[n.zero[i], n.zero[i], drop=FALSE]

          ldl = lDL[,i]
          ldl[!is_zero] = lx[!is_zero]

          b = crossprod(ldl, Bi[n.zero[i],])

          h1_m1 = e_truncnorm(b, m , sqrt(v))
          h1_m2 = h1_m1 %*% t(h1_m1) + v_truncnorm(b, m , sqrt(v))

          h2[n.zero[i],] = h1_m1
          h2_m2 = h2 %*% t(h2)
          h2_m2[n.zero[i],n.zero[i]] = h1_m2

          m1 = t(Bi) %*% h2
          m2 = t(Bi) %*% h2_m2  %*% Bi

        }
      }
      M1[,i] = m1
      M2 = M2 + m2
    }
    MU_next    = matrix(rowMeans(M1), ncol = 1)
    SIGMA = M2 / N - MU_next %*% t(MU_next)
    #print(as.vector(MU_next))
    if( max(abs(MU-MU_next)) < EPS) CONT = FALSE
    MU = MU_next
    # logLik_next = mLogLik(X, MU[,1], SIGMA)
    # if(debug) cat(sprintf("Iteration: %d, LogLik: %0.3f\n", it, logLik_next))
    # if( abs(logLik_next - logLik_current)/abs(logLik_current) < EPS ) CONT = FALSE
    #logLik_current = logLik_next
  }
  if(parameters) return(list(clr_mu = MU[,1], clr_sigma = SIGMA, clr_h = t(M1)))
  t(apply(exp(M1), 2, prop.table))
}
P = coda_zero_replace(X, DL)
X_replaced = P * apply(X / P, 1, max, na.rm=TRUE)
i_zeros = which(!is.na(X) & X == 0)
delta.prop = X_replaced[i_zeros]/DL[i_zeros]

summary(delta.prop)
summary(X0[below_dl]/dl_par)

X <- matrix(c(26.91,8.08,12.59,31.58,6.45,14.39,
              39.73,26.20,0.00,15.22,6.80,12.05,
              10.76,31.36,7.10,12.74,31.34,6.70,
              10.85,46.40,31.89,10.86,0.00,0.00,
              7.57,11.35,30.24,6.39,13.65,30.80,
              38.09,7.62,23.68,9.70,20.91,0.00,
              27.67,7.15,13.05,32.04,6.54,13.55,
              44.41,15.04,7.95,0.00,10.82,21.78,
              11.50,30.33,6.85,13.92,30.82,6.58,
              19.04,42.59,0.00,38.37,0.00,0.00),byrow=TRUE,ncol=6)
X
X_lrSVD = zCompositions::lrSVD(X,label=0,dl=rep(1,6))
X_lrDA <- zCompositions::lrDA(X,label=0,dl=rep(1,6),ini.cov="multRepl",n.iters=150)
head(X_lrSVD)
head(X_lrDA)
X <- matrix(c(26.91,8.08,12.59,31.58,6.45,14.39,
              39.73,41.42,0.00,NA,6.80,12.05,
              NA,35.13,7.96,14.28,35.12,7.51,
              10.85,46.40,31.89,10.86,0.00,0.00,
              10.85,16.27,NA,9.16,19.57,44.15,
              38.09,7.62,23.68,9.70,20.91,0.00,
              NA,9.89,18.04,44.30,9.04,18.73,
              44.41,15.04,7.95,0.00,10.82,21.78,
              11.50,30.33,6.85,13.92,30.82,6.58,
              19.04,42.59,0.00,38.37,0.00,0.00),byrow=TRUE,ncol=6)
X_lrEMplus <- zCompositions::lrEMplus(X,dl=rep(1,6),ini.cov="multRepl")
DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
DL[!is.na(X) & X == 0] = 1
X_replaced = coda_zero_replace(X, DL, debug = TRUE)
head(X_lrEMplus)
head(X_replaced*100)
