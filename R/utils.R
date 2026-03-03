getDim <- function(X) {
  if (is.vector(X)) length(X) else NCOL(X)
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
#' Compute orthonormal ilr bases associated with conditioning patterns on the
#' parts of a composition.
#'
#' Each column of `C` defines one conditioning pattern. For a given column, the
#' ilr basis is constructed by separating the parts marked with `0` from the
#' parts marked with a positive value.
#'
#' If a column contains `nz` zeros, then:
#' \itemize{
#'   \item the first `nz - 1` coordinates describe the internal log-ratio
#'   structure of the parts marked with `0`,
#'   \item the coordinate `nz` describes the balance between the block of parts
#'   marked with `0` and the block of parts marked with positive values,
#'   \item the remaining coordinates describe the internal log-ratio structure
#'   of the parts marked with positive values.
#' }
#'
#' Thus, each basis preserves the split defined by the conditioning pattern and
#' completes it to an orthonormal basis of the clr-plane.
#'
#' @param C A numeric matrix with one conditioning pattern per column. Rows
#'   correspond to parts. For each column, entries equal to `0` define one block
#'   and positive entries define the complementary block.
#'
#' @return A three-dimensional array of dimension `(D - 1, D, ncol(C))`, where
#'   `D` is the number of parts. Each slice contains one orthonormal ilr basis.
#'
#' @examples
#' C <- cbind(
#'   c(0, 0, 1, 1, 0),
#'   c(0, 1, 0, 1, 0)
#' )
#'
#' conditional_obasis(C)
#'
#' @export
conditional_obasis <- c_conditional_obasis
