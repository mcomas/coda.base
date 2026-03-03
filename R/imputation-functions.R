
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
mLogLik <- function(X, mu, sigma, B = ilr_basis(ncol(X))) {
  X <- as.matrix(X)
  n <- nrow(X)
  D <- ncol(X)
  d <- D - 1

  MU <- as.vector(mu %*% t(B))
  SIGMA <- B %*% sigma %*% t(B)

  obs <- !is.na(X)
  pattern_id <- apply(obs, 1, paste0, collapse = "")
  pattern_levels <- unique(pattern_id)

  Bc <- c_conditional_obasis(t(obs))

  lprobs <- numeric(n)

  for (pat in pattern_levels) {
    idx <- which(pattern_id == pat)
    i <- idx[1]

    is_obs <- obs[i, ]
    n0 <- sum(!is_obs)

    if (n0 >= d) {
      lprobs[idx] <- 0
      next
    }

    if (n0 == 0) {
      H <- log(X[idx, , drop = FALSE]) %*% B
      lprobs[idx] <- dmvnorm(H, mean = mu, sigma = sigma, log = TRUE)
    } else {
      Bi <- Bc[, , i]
      B2 <- Bi[-seq_len(n0), , drop = FALSE]

      mu2 <- as.vector(B2 %*% MU)
      sigma2 <- B2 %*% SIGMA %*% t(B2)

      H2 <- log(X[idx, is_obs, drop = FALSE]) %*% t(B2[, is_obs, drop = FALSE])

      lprobs[idx] <- dmvnorm(H2, mean = mu2, sigma = sigma2, log = TRUE)
    }
  }

  sum(lprobs)
}
