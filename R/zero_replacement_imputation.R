#' Replacement of Missing Values and Below-Detection Zeros in Compositional Data
#'
#' @description
#' Performs imputation (replacement) of missing values and/or values below the detection limit (BDL) in compositional datasets,
#' using a parametric approach based on the multivariate normal distribution.
#' This function is designed to prepare compositional data for subsequent log-ratio transformations.
#'
#' @param X A compositional dataset: numeric matrix or data frame where rows represent observations and columns represent parts.
#' @param DL An optional matrix or vector of detection limits. If \code{NULL}, the minimum non-zero value in each column of \code{X} is used.
#' @param dl_prop A numeric value between 0 and 1, used for initialization in the EM algorithm (default is 0.65).
#' @param eps A small positive value controlling the convergence criterion for the EM algorithm (default is \code{1e-4}).
#' @param debug Logical. If \code{TRUE}, displays intermediate steps and diagnostic messages (default is \code{FALSE}).
#' @param parameters Logical. If \code{TRUE}, returns additional output including estimated multivariate normal parameters (default is \code{FALSE}).
#'
#' @return
#' If \code{parameters = FALSE}, returns a numeric matrix with imputed values.
#' If \code{parameters = TRUE}, returns a list with two components:
#' \describe{
#'   \item{X_imp}{The imputed compositional data matrix.}
#'   \item{info}{A list containing information about the EM algorithm parameters and convergence diagnostics.}
#' }
#'
#' @details
#' - Missing values are imputed based on a multivariate normal model.
#' - Zeros are treated as censored values and replaced accordingly.
#' - The EM algorithm iteratively estimates the missing parts and model parameters.
#'
#' @examples
#' # Simulate compositional data with zeros
#' set.seed(123)
#' X <- abs(matrix(rnorm(100), ncol = 5))
#' X[sample(length(X), 10)] <- 0  # Introduce some zeros
#'
#' # Apply replacement
#' X_imp <- coda_replacement(X)
#' summary(X_imp)
#'
#' @seealso
#' \code{\link{cenLR}}, \code{\link{multRepl}} for alternative imputation methods.
#'
#' @export
coda_replacement = function(X, DL = NULL, dl_prop = 0.65, eps = 1e-4,
                            debug = FALSE, parameters = FALSE){

  tX = t(as.matrix(X))
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
      if(n_positive[i] == 0){
        m1 = MU
        m2 = m1 %*% t(m1) + SIGMA
      }else{
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

zero_na_conditional_obasis = function(tX){
  D = nrow(tX)

  tX_na = is.na(tX)
  tX_zero = !tX_na & tX == 0
  n_na = colSums(tX_na)
  n_zero = colSums(tX_zero)
  n_positive = D - n_na - n_zero
  l_B = list(NULL)
  for(nparts in 2:D){
    l_B[[nparts]] = t(ilr_basis(nparts))
  }
  B = array(0, dim = c(D-1, D, ncol(tX)))
  for(i in 1:ncol(tX)){
    if(n_positive[i] == 0){
      B[,,i] = l_B[[D]]
      next
    }
    I = c(which(tX_na[,i]), which(tX_zero[,i]), which(!tX_na[,i] & !tX_zero[,i]))
    if(n_na[i] > 0){
      if(n_na[i] > 1){
        B[1:(n_na[i]-1), I[1:n_na[i]],i] = l_B[[n_na[i]]]
      }
      w = sqrt((n_na[i] * n_positive[i]) / (n_na[i] + n_positive[i]))
      B[n_na[i], I[1:n_na[i]], i] = 1/n_na[i] * w
      B[n_na[i], I[(n_na[i]+n_zero[i]+1):D], i] = -1/n_positive[i] * w
    }
    if(n_zero[i] > 0){
      if(n_zero[i] > 1){
        B[n_na[i] + 1:(n_zero[i]-1), I[(n_na[i] + 1):(n_na[i]+n_zero[i])], i] = l_B[[n_zero[i]]]
      }
      w = sqrt((n_zero[i] * n_positive[i]) / (n_zero[i] + n_positive[i]))
      B[n_na[i]+n_zero[i], I[n_na[i] + 1:n_zero[i]], i] = 1/n_zero[i] * w
      B[n_na[i]+n_zero[i], I[(n_na[i]+n_zero[i]+1):D], i] = -1/n_positive[i] * w
    }
    if(n_positive[i] > 1){
      B[n_na[i]+n_zero[i]+1:(n_positive[i]-1), I[(n_na[i]+n_zero[i]+1):D], i] = l_B[[n_positive[i]]]
    }
  }
  B
}
SWEEP = function(A, K){
  for(k in K){
    D = A[k,k]
    A[k,] = A[k,] / D
    for(i in which(1:nrow(A)!=k)){
      B = A[i,k]
      A[i,] = A[i,] - B * A[k,]
      A[i,k] = -B/D
    }
    A[k,k] = 1/D
  }
  A
}
#' @export
gen_X_with_zeros_and_missings = function(n, d, missings = TRUE, zeros = TRUE,
                                         dl_par = 0.05, na_p = 0.15){
  gen_X = function(n, d){
    H = mvtnorm::rmvnorm(n, rep(0, d))
    S = cov(H)
    EIG = eigen(S)

    as.matrix(composition(H %*% EIG$vectors, 'cdp'))
  }

  X0 = gen_X(n, d)
  X = X0
  below_dl = X < dl_par
  if(zeros){
    X[below_dl] = 0
  }
  if(missings){
    i2 = rbinom(length(X), 1, na_p) == 1
    X[i2] = NA
  }
  DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
  DL[!is.na(X) & below_dl] = dl_par
  list(X = X, DL = DL, X0 = X0)
}
