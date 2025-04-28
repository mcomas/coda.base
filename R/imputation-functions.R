

mLogLik = function(X, mu, sigma, B = ilr_basis(ncol(X))){
  MU = mu %*% t(B)
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
