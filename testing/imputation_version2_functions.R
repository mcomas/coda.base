mLogLik = function(X, mu_clr, sigma_clr){
  B = ilr_basis(ncol(X))
  mu =  mu_clr %*% B
  sigma = t(B) %*% sigma_clr %*% B

  Bc = c_conditional_obasis(!is.na(t(X)))
  d = ncol(X) - 1
  lprobs = sapply(1:nrow(X), function(i){
    x = X[i,]
    is_na = is.na(x)
    n0 = sum(is_na)
    lprob = 0
    if(n0 < d){
      if(n0 == 0){
        lprob = mvtnorm::dmvnorm(log(x) %*% B, mu, sigma, log = TRUE)
      }else{
        Bi = Bc[,,i]
        B2 = Bi[-(1:n0),,drop=FALSE]
        h2 = as.vector(B2[,!is_na,drop=FALSE] %*% log(x[!is_na]))

        mu2 =  mu_clr %*% t(B2)
        sigma2 = B2 %*% sigma_clr %*% t(B2)
        lprob = mvtnorm::dmvnorm(h2, mu2, sigma2, log = TRUE)
      }
    }
    lprob
  })
  sum(lprobs)
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

e_truncnorm = function(b, m, s){
  beta = (b - m) / s
  phi_beta = dnorm(beta, log = TRUE)
  Phi_beta = pnorm(beta, log.p = TRUE)

  m - s * exp(phi_beta - Phi_beta)
}

v_truncnorm = function(b, m, s){
  beta = (b - m) / s
  phi_beta = dnorm(beta, log = TRUE)
  Phi_beta = pnorm(beta, log.p = TRUE)

  r = exp(phi_beta - Phi_beta)
  s * s * (1.0 - r * ( beta - r))
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
