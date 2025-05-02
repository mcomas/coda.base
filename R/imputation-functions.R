
dmvnorm = function (x, mean = rep(0, p), sigma = diag(p), log = FALSE,
          checkSymmetry = TRUE) {
  if (is.vector(x))
    x <- matrix(x, ncol = length(x))
  p <- ncol(x)
  if (!missing(mean)) {
    if (!is.null(dim(mean)))
      dim(mean) <- NULL
    if (length(mean) != p)
      stop("x and mean have non-conforming size")
  }
  if (!missing(sigma)) {
    if (p != ncol(sigma))
      stop("x and sigma have non-conforming size")
    if (checkSymmetry && !isSymmetric(sigma, tol = sqrt(.Machine$double.eps),
                                      check.attributes = FALSE))
      stop("sigma must be a symmetric matrix")
  }
  dec <- tryCatch(base::chol(sigma), error = function(e) e)
  if (inherits(dec, "error")) {
    x.is.mu <- colSums(t(x) != mean) == 0
    logretval <- rep.int(-Inf, nrow(x))
    logretval[x.is.mu] <- Inf
  }
  else {
    tmp <- backsolve(dec, t(x) - mean, transpose = TRUE)
    rss <- colSums(tmp^2)
    logretval <- -sum(log(diag(dec))) - 0.5 * p * log(2 *
                                                        pi) - 0.5 * rss
  }
  names(logretval) <- rownames(x)
  if (log)
    logretval
  else exp(logretval)
}
mLogLik = function(X, mu, sigma, B = ilr_basis(ncol(X))){
  MU = mu %*% t(B)
  SIGMA = B %*% sigma %*% t(B)
  Bc = c_conditional_obasis(!is.na(t(X)))
  d = ncol(X) - 1
  lprobs = sapply(1:nrow(X), function(i){
    x = X[i,]
    is_na = is.na(x)
    n0 = sum(is_na)
    lprob = 0
    if(n0 < d){
      if(n0 == 0){
        lprob = dmvnorm(as.vector( t(B) %*% log(x)), mu, sigma, log = TRUE)
      }else{
        Bi = Bc[,,i]
        B2 = Bi[-(1:n0),,drop=FALSE]
        h2 = as.vector(B2[,!is_na,drop=FALSE] %*% log(x[!is_na]))

        mu2 = B2 %*% MU
        sigma2 = B2 %*% SIGMA %*% t(B2)
        lprob = dmvnorm(h2, mu2, sigma2, log = TRUE)
      }
    }
    lprob
  })
  sum(lprobs)
}
