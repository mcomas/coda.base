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
#' @export
coda_replacement = function(X, DL = NULL, dl_prop = 0.65,
                            eps = 1e-4, parameters = FALSE,
                            debug = FALSE){

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

mLogLik_zr = function(X, DL, dl_prop = 0.65, mu, sigma, Bc = NULL){
  if(is.null(DL)){
    DL = apply(X, 2, function(x) min(x[!is.na(x) & x > 0]))
  }
  if(is.vector(DL)){
    mDL = matrix(DL, ncol = ncol(X), nrow = nrow(X), byrow = TRUE)
  }else{
    mDL = DL
  }
  if(is.null(Bc)){
    Bc = zero_na_conditional_obasis(t(X)) #c_conditional_obasis(!is.na(t(X)))
  }
  Bc_inv = array(0, dim = dim(Bc))
  for(i in 1:nrow(X)){
    Bc_inv[,,i] = t(pinv(Bc[,,i])) #t(MASS::ginv(Bc[,,i]))
  }
  D = ncol(X)
  d = D - 1
  lprobs = sapply(1:nrow(X), function(i){
    x = X[i,]
    Bc_i = Bc[,,i]

    is_na = is.na(x)
    is_zero = !is_na & x == 0
    x[is_zero] = mDL[i,is_zero] * dl_prop

    lx = log(x)

    n.na_zero = sum(is_na) + pmax(sum(is_zero) - 1, 0)


    lprob = 0
    if(n.na_zero < d){
      if(n.na_zero == 0){
        h_obs = as.vector(Bc_i[,!is_na,drop=FALSE] %*% lx[!is_na])
        mu_obs = as.vector(Bc_i[,,drop=FALSE] %*% mu)
        sigma_obs = Bc_i[,,drop=FALSE] %*% sigma %*% t(Bc_i[,,drop=FALSE])
      }else{
        h_obs = as.vector(Bc_i[-(1:n.na_zero),!is_na,drop=FALSE] %*% lx[!is_na])
        mu_obs = as.vector(Bc_i[-(1:n.na_zero),,drop=FALSE] %*% mu)
        sigma_obs = Bc_i[-(1:n.na_zero),,drop=FALSE] %*% sigma %*% t(Bc_i[-(1:n.na_zero),,drop=FALSE])
      }
      lprob = dmvnorm(h_obs, mu_obs, sigma_obs, log = TRUE)
    }
    lprob
  })
  sum(lprobs)
}


# gen_X_with_zeros_and_missings = function(n, d, missings = TRUE, zeros = TRUE,
#                                          dl_par = 0.05, na_p = 0.15){
#   gen_X = function(n, d){
#     H = mvtnorm::rmvnorm(n, rep(0, d))
#     S = cov(H)
#     EIG = eigen(S)
#
#     as.matrix(composition(H %*% EIG$vectors))
#   }
#
#   X0 = gen_X(n, d)
#   X = X0
#   below_dl = X < dl_par
#   if(zeros){
#     X[below_dl] = 0
#   }
#   if(missings){
#     i2 = rbinom(length(X), 1, na_p) == 1
#     X[i2] = NA
#   }
#   DL = matrix(0, nrow = nrow(X), ncol = ncol(X))
#   DL[!is.na(X) & below_dl] = dl_par
#   list(X = X, DL = DL, X0 = X0)
# }
