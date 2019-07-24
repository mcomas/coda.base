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
  X = as.matrix(X)
  c_variation_array(X, as.logical(only_variation))
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

#
# pc_basis = function(X){
#   lX =  log(X)
#   pr = stats::princomp(lX - rowMeans(lX))
#   pr$loadings[,-NCOL(X)]
# }

#' @title Get coordinates from compositions w.r.t. an specific basis
#'
#' @description
#' Calculate the coordinates of a composition with respect a given basis
#'
#' @details
#' \code{coordinates} function calculates the coordinates of a compositiona w.r.t. a given basis. `basis` parameter is
#' used to set the basis, it can be either a matrix defining the log-contrasts in columns or a string defining some well-known
#' log-contrast: 'alr' 'clr', 'ilr', 'pc', 'pb' and 'cdp', for the additive log-ratio, centered log-ratio, isometric log-ratio,
#' clr principal components, clr principal balances or default's CoDaPack balances respectively.
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @param basis basis used to calculate the coordinates. \code{basis} can be either a string or a matrix.
#' Accepted values for strings are: 'ilr' (default), 'clr', 'alr', 'pc', 'pb' and 'cdp'. If \code{basis} is a matrix, it is expected
#' to have log-ratio basis given in columns.
#' @param label name given to the coordinates
#' @param sparse_basis Is the given matrix basis sparse? If TRUE calculation are carried
#' taking into an account sparsity (default `FALSE`)
#'
#' @return
#' Coordinates of composition \code{X} with respect the given \code{basis}.
#'
#' @seealso See functions \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}
#' to define different compositional basis.
#' See function \code{\link{composition}} to obtain details on how to calculate
#' a compositions from given coordinates.
#' @examples
#' coordinates(c(1,2,3,4,5))
#' # basis is shown if 'coda.base.basis' option is set to TRUE
#' options('coda.base.basis' = TRUE)
#' coordinates(c(1,2,3,4,5))
#' # Setting sparse_basi to TRUE can improve performance if log-ratio basis is sparse.
#' N = 100
#' K = 1000
#' X = matrix(exp(rnorm(N*K)), nrow=N, ncol=K)
#' system.time(coordinates(X, alr_basis(K), sparse_basis = FALSE))
#' system.time(coordinates(X, alr_basis(K), sparse_basis = TRUE))
#' system.time(coordinates(X, 'alr'))
#' @export
coordinates = function(X, basis = 'ilr', label = NULL, sparse_basis = FALSE){
  class_type = class(X)
  is_vector = is.atomic(X) & !is.list(X) & !is.matrix(X)
  is_data_frame = inherits(X, 'data.frame')
  RAW = X
  if(is.list(basis) & !is.data.frame(basis)){
    basis = do.call('sbp_basis', args = c(basis, list(data = X)), envir = as.environment('package:coda.base'))
  }
  if(is_vector){
    class_type = 'double'
    RAW = matrix(X, nrow=1)
  }
  if(is_data_frame){
    RAW = as.matrix(X)
  }
  non_compositional = rowSums(is.na(RAW) | RAW <= 0)
  if(sum(non_compositional) > 0){
    warning("Some observations are not compositional (either missing or non-strictly positive). They are returned as missing values.",
            call. = FALSE)
  }
  sel_compositional = non_compositional == 0
  if(is.character(basis)){
    if(is.null(label)){
      label = basis
    }
    dim = ncol(RAW)
    RAW.coda = matrix(RAW[sel_compositional, ], ncol = dim)
    coord.dim = dim - 1
    if(basis == 'ilr'){
      basis = ilr_basis(dim)
      COORD.coda = coordinates_basis(RAW.coda, ilr_basis(dim), sparse = FALSE)
    }else{
      if(basis == 'alr'){
        basis = alr_basis(dim)
        COORD.coda = coordinates_alr(RAW.coda, 0)
      }else{
        if(basis == 'clr'){
          coord.dim = ncol(RAW.coda)
          basis = clr_basis(dim)
          COORD.coda = coordinates_basis(RAW.coda, clr_basis(dim), sparse = FALSE)
        }else{
          if(basis == 'pc'){
            lRAW =  log(RAW.coda)
            SVD = svd(scale(lRAW - rowMeans(lRAW), scale = FALSE))
            basis = SVD$v[,-dim]
            COORD.coda = coordinates(RAW.coda, basis = basis, label = 'pc')
          }else{
            if(basis == 'pb'){
              if(ncol(RAW.coda) > 15){
                message("Number of columns is high for 'pb' method. Depending on your system this computation can take too long. Consider using an approximate method. Consult 'pb_basis()' function for more details.")
              }
              basis = pb_basis(RAW.coda, method = 'exact')
              COORD.coda = coordinates_basis(RAW.coda, basis, sparse = FALSE)
            }else{
              if(basis == 'cdp'){
                basis = sbp_basis(cdp_partition(dim))
                COORD.coda = coordinates(RAW.coda, basis = basis, label = 'cdp')
              }else{
                stop(sprintf('Basis %d not recognized'))
              }
            }
          }
        }
      }
    }
    COORD = matrix(NA_real_, ncol = coord.dim, nrow = nrow(RAW))
    COORD[sel_compositional,] = COORD.coda
  }else{
    if(is.matrix(basis)){
      if(is.null(label)){
        label = 'x'
      }
      if(max(colSums(basis)) > .Machine$double.eps^0.5){
        warning("Supplied basis matrix is not a log-contrast.")
      }
      dim = nrow(basis)
      coord.dim = ncol(basis)
      RAW.coda = matrix(RAW[sel_compositional, ], ncol = dim)
      COORD.coda = coordinates_basis(RAW.coda, basis, sparse_basis)
      COORD = matrix(NA_real_, ncol = coord.dim, nrow = nrow(RAW))
      COORD[sel_compositional,] = COORD.coda
    }else{
      stop(sprintf('Basis need to be either an string or a matrix'))
    }
  }

  colnames(COORD) = sprintf(sprintf('%s%%0%dd',label, 1+floor(log(ncol(COORD), 10))),1:ncol(COORD))
  if(is_vector){
    COORD = COORD[1,]
    names(COORD) = sprintf(sprintf('%s%%0%dd', label, 1+floor(log(length(COORD), 10))),1:length(COORD))
  }
  if(is_data_frame){
    COORD = as.data.frame(COORD)
  }
  class(COORD) = class_type
  attr(COORD, 'basis') = basis
  suppressWarnings(row.names(COORD) <- row.names(X))
  set.coda(COORD)
}

#' Get composition from coordinates w.r.t.  an specific basis
#'
#' Calculate a composition from coordinates with respect a given basis
#'
#' @param H coordinates of a composition. Either a matrix, a data.frame or a vector
#' @param basis basis used to calculate the coordinates
#' @param label name given to the coordinates
#' @param sparse_basis Is the given matrix basis sparse? If TRUE calculation are carried
#' taking into an account sparsity (default `FALSE`)
#' @return coordinates with respect the given basis
#' @seealso See functions \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}
#' to define different compositional basis.
#' See function \code{\link{coordinates}} to obtain details on how to calculate
#' coordinates of a given composition.
#' @export
composition = function(H, basis = NULL, label = 'x', sparse_basis = FALSE){
  class_type = class(H)
  if(is.null(basis) & "basis" %in% names(attributes(H))){
    basis = attr(H, 'basis')
  }
  if(is.null(basis)){
    basis = 'ilr'
  }
  is_vector = is.atomic(H) & !is.list(H) & !is.matrix(H)
  is_data_frame = inherits(H, 'data.frame')

  COORD = H
  if(is_vector){
    class_type = 'double'
    COORD = matrix(H, nrow=1)
  }
  if(is_data_frame){
    COORD = as.matrix(H)
  }
  if(is.character(basis)){
    dim = ncol(COORD)+1
    if(basis == 'ilr'){
      basis = ilr_basis(dim)
      RAW = exp(COORD %*% t(basis))
    }else{
      if(basis == 'alr'){
        basis = alr_basis(dim)
        RAW = cbind(exp(COORD), 1)
      }else{
        if(basis == 'clr'){
          dim = ncol(COORD)
          basis = clr_basis(dim)
          RAW = exp(COORD)
        }else{
          if(basis == 'cdp'){
            basis = sbp_basis(cdp_partition(dim))
            RAW = composition(COORD, basis = basis)
          }else{
            stop(sprintf('Basis %d not recognized'))
          }
        }
      }
    }
  }else{
    if(is.matrix(basis)){
      RAW = exp(COORD %*% pinv(basis))
    }else{
      stop(sprintf('Basis need to be either an string or a matrix'))
    }
  }
  RAW = RAW / rowSums(RAW)
  colnames(RAW) = sprintf(sprintf('%s%%0%dd',label, 1+floor(log(ncol(RAW), 10))),1:ncol(RAW))
  if(is_vector){
    RAW = RAW[1,]
    names(RAW) = sprintf(sprintf('%s%%0%dd', label, 1+floor(log(length(RAW), 10))),1:length(RAW))
  }
  if(is_data_frame){
    RAW = as.data.frame(RAW)
  }
  class(RAW) = setdiff(class_type, 'coda')
  suppressWarnings(row.names(RAW) <- row.names(H))
  #attr(RAW, 'basis') = basis
  RAW
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
