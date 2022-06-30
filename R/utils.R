getDim = function(X) ifelse(is.vector(X), length(X), NCOL(X))

#' Variation array is returned.
#'
#' @param X Compositional dataset
#' @param only_variation if TRUE only the variation matrix is calculated
#' @return variation array matrix
#' @examples
#' set.seed(1)
#' X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
#' variation_array(X)
#' variation_array(X, only_variation = TRUE)
#' @export
variation_array = function(X, only_variation = FALSE){
  var_arr = c_variation_array(as.matrix(X), as.logical(only_variation))
  if(!is.null(colnames(X))) colnames(var_arr) = rownames(var_arr) = colnames(X)
  var_arr
}

#' Distance Matrix Computation (including Aitchison distance)
#'
#' This function overwrites \code{\link[stats]{dist}} function to contain Aitchison distance between
#' compositions.
#'
#' @param x compositions
#' method
#' @param method the distance measure to be used. This must be one of "aitchison", "euclidean", "maximum",
#' "manhattan", "canberra", "binary" or "minkowski". Any unambiguous substring can be given.
#' @param ... arguments passed to \code{\link[stats]{dist}} function
#' @return \code{dist} returns an object of class "dist".
#' @seealso See functions \code{\link[stats]{dist}}.
#' @examples
#' X = exp(matrix(rnorm(10*50), ncol=50, nrow=10))
#'
#' (d <- dist(X, method = 'aitchison'))
#' plot(hclust(d))
#'
#' # In contrast to Euclidean distance
#' dist(rbind(c(1,1,1), c(100, 100, 100)), method = 'euc') # method = 'euclidean'
#' # using Aitchison distance, only relative information is of importance
#' dist(rbind(c(1,1,1), c(100, 100, 100)), method = 'ait') # method = 'aitchison'
#'
#' @export
dist = function(x, method = 'euclidean', ...){
  METHODS <- c("aitchison", "euclidean", "maximum",
               "manhattan", "canberra",  "binary", "minkowski")
  imethod <- pmatch(method, METHODS)
  .coda = FALSE
  if (!is.na(imethod) & imethod == 1) {
    .coda = TRUE
    x = coordinates(x)
    method = 'euclidean'
  }
  adist = stats::dist(x, method = method, ...)
  if(.coda){
    attr(adist, 'method') = 'aitchison'
  }
  adist
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
gmean = function(x, zero.rm = FALSE, trim = 0, na.rm = FALSE){
  if(any(x < 0)) stop('Negative values')

  if(zero.rm) x = x[x != 0]

  lmean = mean(log(x), trim = trim, na.rm = na.rm)

  if(lmean == -Inf) return(0)
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
center = function(X, zero.rm = FALSE, na.rm = FALSE){
  C = apply(X, 2, gmean, zero.rm = zero.rm, na.rm = na.rm)
  setNames(C/sum(C), colnames(X))
}

fillPartition = function(partition, row, left, right){
  new_row = rep(0, ncol(partition))
  if(right - left <= 0){
    return(partition)
  }
  if(right - left == 1){
    new_row[left] = 1
    new_row[right] = -1
    if(row == 0){
      partition = rbind(new_row)
    }else{
      partition = rbind(partition, new_row)
    }
    return(partition)
  }
  middle = left + (0.5 + right - left)/2
  new_row[left:floor(middle)] = 1
  new_row[ceiling(middle):right] = -1
  if(row == 0){
    partition = rbind(new_row)
  }else{
    partition = rbind(partition, new_row)
  }
  partition = fillPartition(partition, nrow(partition), left, floor(middle))
  partition = fillPartition(partition, nrow(partition), ceiling(middle), right)
  return(partition)
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
cdp_partition = function(ncomp) unname(t(fillPartition(matrix(0, nrow = 1, ncol = ncomp), 0, 1, ncomp)))
