getDim <- function(X) {
  if (is.vector(X)) length(X) else NCOL(X)
}

#' Closure operation for compositional data
#'
#' Applies the closure operation to a numeric vector, matrix or data frame so
#' that each composition sums to a prescribed constant \code{k}.
#'
#' If \code{X} is:
#' \itemize{
#'   \item a vector, the returned vector sums to \code{k};
#'   \item a matrix or data frame, closure is applied row-wise, and each row
#'   sums to the corresponding value of \code{k}.
#' }
#'
#' The argument \code{k} may be:
#' \itemize{
#'   \item a single positive number, recycled to all rows;
#'   \item a numeric vector of length \code{nrow(X)}, specifying a different
#'   closure constant for each row.
#' }
#'
#' @param X A numeric vector, matrix, data frame, or an object coercible to one
#'   of these. For matrices and data frames, rows are interpreted as
#'   compositions.
#' @param k A numeric vector of length 1 or length \code{nrow(X)}. Must contain
#'   strictly positive values.
#'
#' @return
#' If \code{X} is a vector, a numeric vector of the same length.
#'
#' If \code{X} is a matrix, a numeric matrix with the same dimensions,
#' dimnames, and row-wise sums equal to \code{k}.
#'
#' If \code{X} is a data frame, a data frame with the same row and column names,
#' and row-wise sums equal to \code{k}.
#'
#' @details
#' For a composition \eqn{x = (x_1, \dots, x_D)} with positive sum,
#' the closure to constant \eqn{k} is
#' \deqn{C(x) = k \frac{x}{\sum_{j=1}^D x_j}.}
#'
#' This function requires all entries of \code{X} to be finite and
#' non-negative, and every row sum (or the vector sum) must be strictly
#' positive.
#'
#' @examples
#' closure(c(2, 3, 5))
#' closure(c(2, 3, 5), k = 100)
#'
#' X <- matrix(c(1, 1, 2,
#'               2, 3, 5), nrow = 2, byrow = TRUE)
#' closure(X)
#' closure(X, k = c(1, 100))
#'
#' df <- data.frame(a = c(1, 2), b = c(1, 3), c = c(2, 5))
#' closure(df, k = 10)
#'
#' @export
closure <- function(X, k = 1) {

  if (is.atomic(X) && is.null(dim(X))) {
    if (!is.numeric(X)) {
      stop("'X' must be numeric.")
    }
    if (length(k) != 1L) {
      stop("If 'X' is a vector, 'k' must have length 1.")
    }
    if (!is.numeric(k) || !is.finite(k) || k <= 0) {
      stop("'k' must be a strictly positive finite number.")
    }
    if (any(!is.finite(X))) {
      stop("'X' must contain only finite values.")
    }
    if (any(X < 0)) {
      stop("'X' must contain non-negative values.")
    }

    s <- sum(X)
    if (s <= 0) {
      stop("The sum of 'X' must be strictly positive.")
    }

    return(k * X / s)
  }

  is_df <- is.data.frame(X)

  if (is_df) {
    X_mat <- as.matrix(X)
  } else if (is.matrix(X)) {
    X_mat <- X
  } else {
    X_mat <- tryCatch(as.matrix(X), error = function(e) NULL)
    if (is.null(X_mat)) {
      stop("'X' must be a numeric vector, matrix, data.frame, or coercible to matrix.")
    }
  }

  if (!is.numeric(X_mat)) {
    stop("'X' must be numeric.")
  }
  if (any(!is.finite(X_mat))) {
    stop("'X' must contain only finite values.")
  }
  if (any(X_mat < 0)) {
    stop("'X' must contain non-negative values.")
  }

  n <- nrow(X_mat)
  if (is.null(n)) {
    stop("Could not determine 'nrow(X)'.")
  }

  if (!(length(k) %in% c(1L, n))) {
    stop("'k' must have length 1 or length nrow(X).")
  }
  if (!is.numeric(k) || any(!is.finite(k)) || any(k <= 0)) {
    stop("'k' must contain strictly positive finite values.")
  }

  if (length(k) == 1L) {
    k <- rep.int(k, n)
  }

  rs <- rowSums(X_mat)
  if (any(rs <= 0)) {
    stop("All row sums of 'X' must be strictly positive.")
  }

  out <- X_mat / rs
  out <- out * k

  if (is_df) {
    out <- as.data.frame(out, stringsAsFactors = FALSE)
    names(out) <- names(X)
    rownames(out) <- rownames(X)
  }

  out
}

#' Variation array is returned.
#'
#' @param X Compositional dataset
#' @param include_means if TRUE logratio means are included in the lower-left triangle
#' @param ml_covariance if TRUE Maximum-likelihood estimation of the covariance for the multivariate normal distribution is used (dividing the scatter matrix by n instead of n-1)
#' @return variation array matrix
#' @examples
#' set.seed(1)
#' X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
#' variation_array(X)
#' variation_array(X, include_means = TRUE)
#' @export
variation_array <- function(X, include_means = FALSE, ml_covariance = FALSE) {
  var_arr <- c_variation_array(
    as.matrix(X),
    as.logical(include_means),
    as.logical(ml_covariance)
  )

  if (!is.null(colnames(X))) {
    colnames(var_arr) <- colnames(X)
    rownames(var_arr) <- colnames(X)
  }

  var_arr
}

#' Distance Matrix Computation (including Aitchison distance)
#'
#' Compute a distance matrix for compositional data, including the Aitchison
#' distance as an extension of \code{\link[stats]{dist}}.
#'
#' @param x A data matrix whose rows are compositions.
#' @param method The distance measure to be used. This must be one of
#'   \code{"aitchison"}, \code{"euclidean"}, \code{"maximum"},
#'   \code{"manhattan"}, \code{"canberra"}, \code{"binary"}, or
#'   \code{"minkowski"}. Any unambiguous abbreviation can be given.
#' @param ... Additional arguments passed to \code{\link[stats]{dist}}.
#'
#' @return An object of class \code{"dist"}.
#'
#' @seealso \code{\link[stats]{dist}}
#'
#' @examples
#' X <- exp(matrix(rnorm(10 * 50), ncol = 50, nrow = 10))
#'
#' (d <- dist(X, method = "aitchison"))
#' plot(hclust(d))
#'
#' # In contrast to Euclidean distance
#' dist(rbind(c(1, 1, 1), c(100, 100, 100)), method = "euc")
#'
#' # Using Aitchison distance, only relative information is of importance
#' dist(rbind(c(1, 1, 1), c(100, 100, 100)), method = "ait")
#'
#' @export
dist = function(x, method = "euclidean", ...) {
  methods_available <- c(
    "aitchison", "euclidean", "maximum",
    "manhattan", "canberra", "binary", "minkowski"
  )

  imethod <- pmatch(method, methods_available)

  if (is.na(imethod)) {
    stop(
      "'method' must be one of: ",
      paste(methods_available, collapse = ", ")
    )
  }

  is_aitchison <- (imethod == 1)

  if (is_aitchison) {
    x <- coordinates(x)
    method <- "euclidean"
  } else {
    method <- methods_available[imethod]
  }

  d <- stats::dist(x, method = method, ...)

  if (is_aitchison) {
    attr(d, "method") <- "aitchison"
  }

  d
}

#' Geometric Mean
#'
#' Generic function for the (trimmed) geometric mean.
#'
#' @param x A nonnegative vector.
#' @param zero.rm a logical value indicating whether zero values should be stripped
#' before the computation proceeds.
#' @param trim the fraction (0 to 0.5) of observations to be trimmed from each
#' end of x before the mean is computed. Values of trim outside that range are
#' taken as the nearest endpoint.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @seealso \code{\link{center}}
#' @export
gmean <- function(x, zero.rm = FALSE, trim = 0, na.rm = FALSE) {
  if (any(x < 0, na.rm = TRUE)) {
    stop("Negative values")
  }

  if (zero.rm) {
    x <- x[x != 0]
  }

  lmean <- mean(log(x), trim = trim, na.rm = na.rm)

  if (is.infinite(lmean) && lmean < 0) {
    return(0)
  }

  exp(lmean)
}

#' Dataset center
#'
#' Generic function to calculate the center of a compositional dataset
#'
#' @param X compositional dataset
#' @param zero.rm a logical value indicating whether zero values should be stripped
#' before the computation proceeds.
#' @param na.rm	a logical value indicating whether NA values should be stripped
#' before the computation proceeds.
#' @examples
#' X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
#' g = rep(c('a','b','c','d'), 25)
#' center(X)
#' (by_g <- by(X, g, center))
#' center(t(simplify2array(by_g)))
#' @export
center <- function(X, zero.rm = FALSE, na.rm = FALSE) {
  C <- apply(X, 2, gmean, zero.rm = zero.rm, na.rm = na.rm)
  setNames(C / sum(C), colnames(X))
}


.fill_partition_list <- function(left, right, D) {
  if (right <= left) {
    return(list())
  }

  new_row <- rep(0, D)

  if (right - left == 1) {
    new_row[left] <- 1
    new_row[right] <- -1
    return(list(new_row))
  }

  middle <- left + (0.5 + right - left) / 2

  new_row[left:floor(middle)] <- 1
  new_row[ceiling(middle):right] <- -1

  c(
    list(new_row),
    Recall(left, floor(middle), D),
    Recall(ceiling(middle), right, D)
  )
}
# Function used to build balanced partitions. default basis in CoDaPack
fill_partition <- function(D) {
  do.call(rbind, .fill_partition_list(1, D, D))
}

fill_sbp <- function(iRES) {
  D <- nrow(iRES)
  K <- ncol(iRES)
  target_K <- D - 1

  if (K == 0) {
    return(sign(ilr_basis(D)))
  }

  # Add one balance for completely free components
  free_comp <- rowSums(iRES != 0) == 0
  if (any(free_comp)) {
    iRES <- cbind(iRES, 1 - 2 * free_comp)
    K <- ncol(iRES)
  }

  # Choose root balance and, if needed, force a split of its free components
  iROOT <- which.max(colSums(iRES != 0))
  free_comp <- iRES[, iROOT] == 0
  if (any(free_comp)) {
    iRES <- cbind(iRES, 1 - 2 * free_comp)
    iROOT <- ncol(iRES)
    K <- ncol(iRES)
  }

  iBAL <- which(iRES[, iROOT] > 0)

  BAL <- matrix(0, D, target_K)
  BAL[, seq_len(K)] <- iRES

  # Fill positive branch
  pBAL <- iRES[iBAL, -iROOT, drop = FALSE]
  pBAL <- pBAL[, colSums(pBAL != 0) > 0, drop = FALSE]
  pDIM <- length(iBAL) - 1
  pK <- ncol(pBAL)

  if (pK < pDIM) {
    pBAL_full <- Recall(pBAL)
    BAL[iBAL, (K + 1):(K + pDIM - pK)] <- pBAL_full[, (pK + 1):pDIM, drop = FALSE]
  }

  # Fill negative branch
  nBAL <- iRES[-iBAL, -iROOT, drop = FALSE]
  nBAL <- nBAL[, colSums(nBAL != 0) > 0, drop = FALSE]
  nDIM <- nrow(nBAL) - 1
  nK <- ncol(nBAL)

  if (nK < nDIM) {
    nBAL_full <- Recall(nBAL)
    BAL[-iBAL, (K + 1 + pDIM - pK):(K + pDIM + nDIM - pK)] <-
      nBAL_full[, (nK + 1):nDIM, drop = FALSE]
  }

  BAL
}

#' CoDaPack's default binary partition
#'
#' Compute the default binary partition used in CoDaPack's software
#'
#' @param ncomp number of parts
#' @return matrix
#' @examples
#' cdp_partition(4)
#' @export
cdp_partition <- function(ncomp) {
  if (length(ncomp) != 1 || !is.numeric(ncomp) || ncomp < 2) {
    stop("'ncomp' must be a single integer greater than or equal to 2.")
  }

  ncomp <- as.integer(ncomp)
  unname(t(fill_partition(ncomp)))
}

#' Conditional orthonormal basis
#'
#' Compute orthonormal ilr bases adapted to row-wise conditioning patterns.
#'
#' Each row of `X` defines one conditioning pattern on the parts of a
#' composition. According to `scheme`, the parts are split into ordered blocks:
#'
#' \itemize{
#'   \item `"zero"`: parts equal to `0` and parts with strictly positive values,
#'   \item `"zero_na"`: missing values (`NA`), zeros, and strictly positive values.
#' }
#'
#' For each row, the function constructs an orthonormal basis of the clr-plane
#' preserving the block structure induced by the selected scheme.
#'
#' Under `scheme = "zero"`, if a row contains `nz` zeros, then:
#' \itemize{
#'   \item the first `nz - 1` coordinates describe the internal log-ratio
#'   structure of the zero block,
#'   \item the coordinate `nz` describes the balance between the zero block and
#'   the positive block,
#'   \item the remaining coordinates describe the internal log-ratio structure
#'   of the positive block.
#' }
#'
#' Under `scheme = "zero_na"`, the blocks are ordered as:
#' \itemize{
#'   \item missing values (`NA`),
#'   \item zeros,
#'   \item strictly positive values.
#' }
#'
#' In this case:
#' \itemize{
#'   \item the first coordinates describe the internal structure of the `NA` block,
#'   \item the next coordinate contrasts the `NA` block with the positive block,
#'   \item the following coordinates describe the internal structure of the zero block,
#'   \item the next coordinate contrasts the zero block with the positive block,
#'   \item the remaining coordinates describe the internal structure of the
#'   positive block.
#' }
#'
#' @param X A numeric matrix or data frame with one observation or conditioning
#'   pattern per row and one part per column.
#' @param scheme Character string indicating the conditioning scheme. Possible
#'   values are `"zero"` and `"zero_na"`. Default is `"zero"`.
#'
#' @return A three-dimensional array of dimension `(D - 1, D, nrow(X))`, where
#'   `D` is the number of parts. Each slice contains one orthonormal ilr basis.
#'
#' @examples
#' C <- rbind(
#'   c(0, 0, 1, 1, 0),
#'   c(0, 1, 0, 1, 0)
#' )
#'
#' conditional_obasis(C)
#'
#' X <- rbind(
#'   c(1, NA, 0, 2),
#'   c(NA, 3, 0, 4),
#'   c(1, 2, 3, 4)
#' )
#'
#' conditional_obasis(X, scheme = "zero_na")
#'
#' @export
conditional_obasis <- function(X, scheme = c("zero", "zero_na")) {
  scheme <- match.arg(scheme)

  if (is.data.frame(X)) {
    X <- as.matrix(X)
  }

  if (!is.matrix(X)) {
    stop("'X' must be a matrix or data.frame.", call. = FALSE)
  }

  if (!is.numeric(X)) {
    stop("'X' must be numeric.", call. = FALSE)
  }

  if (ncol(X) < 2) {
    stop("'X' must contain at least two parts.", call. = FALSE)
  }

  out <- switch(
    scheme,
    zero = c_conditional_obasis(t(X)),
    zero_na = c_zero_na_conditional_obasis(t(X))
  )

  dn <- dimnames(X)
  dimnames(out) <- list(NULL, dn[[2]], dn[[1]])
  out
}

#' Generate compositional data with zeros and missing values
#'
#' @description
#' Simulate compositional data and optionally introduce structural zeros
#' (interpreted as values below a detection limit) and missing values.
#'
#' The function first generates a compositional data set `X0`, then creates a
#' modified version `X` by:
#' \itemize{
#'   \item replacing values below `dl_par` by zero, if `zeros = TRUE`,
#'   \item introducing missing values at random, if `missings = TRUE`.
#' }
#'
#' A matrix of detection limits `DL` is also returned. It contains `dl_par` in
#' the positions that were censored to zero, and `0` elsewhere.
#'
#' @param n Number of observations.
#' @param d Dimension of the latent coordinate space used to generate the
#'   compositions.
#' @param missings Logical; if `TRUE`, introduce missing values at random.
#' @param zeros Logical; if `TRUE`, replace values below `dl_par` by zero.
#' @param dl_par Detection-limit threshold used to generate zeros.
#' @param na_p Probability that any entry is replaced by `NA` when
#'   `missings = TRUE`.
#'
#' @return A list with three components:
#' \describe{
#'   \item{X}{The generated compositional data set with simulated zeros and/or missing values.}
#'   \item{DL}{A matrix of detection limits, with `dl_par` in censored positions and `0` elsewhere.}
#'   \item{X0}{The original simulated compositional data set before introducing zeros or missing values.}
#' }
#'
#' @details
#' Compositions are generated from multivariate normal coordinates and mapped to
#' the simplex through `composition()`. The eigenvector rotation is included to
#' induce a non-trivial covariance structure in the generated coordinates.
#'
#' Missing values are introduced completely at random, independently for each
#' cell, with probability `na_p`.
#'
#' @examples
#' set.seed(123)
#' sim <- gen_coda_with_zeros_and_missings(100, 4)
#'
#' str(sim)
#' summary(sim$X0)
#' summary(sim$X)
#' table(sim$X == 0, useNA = "ifany")
#'
#' @export
gen_coda_with_zeros_and_missings <- function(n, d,
                                          missings = TRUE,
                                          zeros = TRUE,
                                          dl_par = 0.05,
                                          na_p = 0.15) {
  if (!is.numeric(n) || length(n) != 1 || is.na(n) || n < 1) {
    stop("'n' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(d) || length(d) != 1 || is.na(d) || d < 1) {
    stop("'d' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(dl_par) || length(dl_par) != 1 || is.na(dl_par) || dl_par <= 0) {
    stop("'dl_par' must be a positive number.", call. = FALSE)
  }
  if (!is.numeric(na_p) || length(na_p) != 1 || is.na(na_p) || na_p < 0 || na_p > 1) {
    stop("'na_p' must be a number between 0 and 1.", call. = FALSE)
  }

  n <- as.integer(n)
  d <- as.integer(d)

  gen_X <- function(n, d) {
    H <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
    S <- stats::cov(H)
    EIG <- eigen(S)

    as.matrix(composition(H %*% EIG$vectors))
  }

  X0 <- gen_X(n, d)
  X <- X0

  below_dl <- X < dl_par
  zero_mask <- matrix(FALSE, nrow = nrow(X), ncol = ncol(X))

  if (zeros) {
    X[below_dl] <- 0
    zero_mask <- below_dl
  }

  if (missings) {
    miss_mask <- matrix(stats::rbinom(length(X), 1, na_p) == 1,
                        nrow = nrow(X), ncol = ncol(X))
    X[miss_mask] <- NA
  }

  DL <- matrix(0, nrow = nrow(X), ncol = ncol(X))
  DL[zero_mask] <- dl_par

  list(X = X, DL = DL, X0 = X0)
}
