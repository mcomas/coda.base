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

alr = function(X, basis_return){
  COORD = alr_coordinates(X, ncol(X))
  colnames(COORD) = paste0('alr', 1:ncol(COORD))
  if(basis_return){
    B = alr_basis(ncol(X))
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

ilr = function(X, basis_return){
  COORD = ilr_coordinates(X)
  colnames(COORD) = paste0('ilr', 1:ncol(COORD))
  if(basis_return){
    B = ilr_basis(ncol(X))
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

clr = function(X, basis_return){
  COORD = clr_coordinates(X)
  colnames(COORD) = paste0('clr', 1:ncol(COORD))
  if(basis_return){
    B = clr_basis(ncol(X))
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

pc = function(X, basis_return){
  B = as.matrix(ilr_basis(ncol(X)))
  B = B %*% svd(scale(log(X) %*% B, scale=FALSE))$v
  COORD = matrix_coordinates(X, B)
  colnames(COORD) = paste0('pc', 1:ncol(B))
  if(basis_return){
    colnames(B) = paste0('pc', 1:ncol(B))
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

cdp = function(X, basis_return){
  B = cdp_basis(ncol(X))
  COORD = sparse_coordinates(X, B)
  colnames(COORD) = colnames(B)
  if(basis_return){
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

pb = function(X, basis_return){
  B = pb_basis(X, method = 'exact')
  COORD = sparse_coordinates(X, B)
  colnames(COORD) = colnames(B)
  if(basis_return){
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

pw = function(X, basis_return){
  B = pairwise_basis(ncol(X))
  COORD = sparse_coordinates(X, B)
  colnames(COORD) = colnames(B)
  if(basis_return){
    if(!is.null(colnames(X))){
      rownames(B) = colnames(X)
    }
    attr(COORD, 'basis') = B
  }
  COORD
}

#' @title Get coordinates from compositions w.r.t. an specific basis
#'
#' @description
#' Calculate the coordinates of a composition with respect a given basis
#'
#' @details
#' \code{coordinates} function calculates the coordinates of a compositiona w.r.t. a given basis. `basis` parameter is
#' used to set the basis, it can be either a matrix defining the log-contrasts in columns or a string defining some well-known
#' log-contrast: 'alr' 'clr', 'ilr', 'pw', 'pc', 'pb' and 'cdp', for the additive log-ratio, centered log-ratio, isometric log-ratio,
#' pairwise log-ratio, clr principal components, clr principal balances or default's CoDaPack balances respectively.
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @param basis basis used to calculate the coordinates. \code{basis} can be either a string or a matrix.
#' Accepted values for strings are: 'ilr' (default), 'clr', 'alr', 'pw', 'pc', 'pb' and 'cdp'. If \code{basis} is a matrix, it is expected
#' to have log-ratio basis given in columns.
#' @param basis_return Should the basis be returned as attribute? (default: \code{TRUE})
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
#' h = coordinates(c(1,2,3,4,5))
#' basis(h)
#' # basis is shown if 'coda.base.basis' option is set to TRUE
#' options('coda.base.basis' = TRUE)
#' coordinates(c(1,2,3,4,5))
#' # Default transformation can improve performance.
#' N = 100
#' K = 1000
#' X = matrix(exp(rnorm(N*K)), nrow=N, ncol=K)
#' system.time(coordinates(X, alr_basis(K)))
#' system.time(coordinates(X, 'alr'))
#' @export
coordinates = function(X, basis = 'ilr', basis_return = TRUE){
  if(is.matrix(X)){
    if(is.character(basis)){   # default's basis with characters
      COORD = get(basis)(X, basis_return)
      #sprintf(sprintf('%s%%0%dd', basis, 1+floor(log(ncol(COORD), 10))),1:ncol(COORD))
    }else{                      # matrix basis
      if('dgCMatrix' %in% class(basis)){
        COORD = sparse_coordinates(X, basis)
      }else{
        COORD = matrix_coordinates(X, basis)
      }
      if(basis_return) attr(COORD, 'basis') = basis
      if(!is.null(colnames(basis))){
        colnames(COORD) = colnames(basis)
      }else{
        colnames(COORD) = sprintf("h%d", 1:ncol(COORD))
      }
    }
  }else{
    if(is.atomic(X) & !is.list(X)){ # vector
      COORD = Recall(matrix(X, nrow = 1), basis, basis_return)
      B = attr(COORD, 'basis')
      COORD = COORD[1,]
      attr(COORD, 'basis') = B
    }else{
      class_type = class(X)
      if(inherits(X, 'data.frame')){
        if(!all(sapply(X, is.numeric))) stop("All parts must be numeric", call. = FALSE)
        mCOORD = Recall(as.matrix(X), basis, basis_return)
        COORD = as.data.frame(mCOORD)
        attr(COORD, 'basis') = attr(mCOORD, 'basis')
      }

      class(COORD) = class_type
    }

  }
  suppressWarnings(row.names(COORD) <- row.names(X))
  set.coda(COORD)
}

#' @rdname coordinates
#' @export
coord = function(..., basis = 'ilr'){
  largs = list(...)
  # cat("---\n")
  # str(largs)
  # cat("---\n")
  lpars = as.list(substitute(largs))
  inum = sapply(1:length(lpars), function(i) is.numeric(lpars[[i]]) & !is.matrix(lpars[[i]]) & length(lpars[[i]] > 1))
  if(sum(inum) == 0){
    if(length(lpars) == 1 & (is.matrix(lpars[[1]]) | is.data.frame(lpars[[1]]) | is.vector(lpars[[1]]))){
      coordinates(lpars[[1]], basis = basis)
    }else{
      stop("Please specify second argument", call. = FALSE)
    }

  }else{
    if(sum(inum) == 1){
      stop("Compositions should have at leat two parts", call. = FALSE)
    }
    if(1 < sum(inum)  & sum(inum) < length(lpars)){
      stop("All components should be numeric", call. = FALSE)
    }
    coordinates(cbind(...), basis = basis)
  }
  # print(inum)
  # print(lpars)
  # cat("---\n")
  # coordinates(a)
}

#' #' Define a composition passing parts one by one.
#' #'
#' #' @param ... compositional dataset. Composition column names
#' #' @param basis basis used to calculate the coordinates. \code{basis} can be either a string or a matrix.
#' #' Accepted values for strings are: 'ilr' (default), 'clr', 'alr', 'pw', 'pc', 'pb' and 'cdp'. If \code{basis} is a matrix, it is expected
#' #' to have log-ratio basis given in columns.
#' #' @export
#' coda = function(..., basis = 'ilr'){
#'   X = coordinates(cbind(...), basis = basis)
#'   X
#' }

coordinates_sparse = function(X, B, sparse){
  coordinates_basis(X,B,sparse)
}

#' Get composition from coordinates w.r.t.  an specific basis
#'
#' Calculate a composition from coordinates with respect a given basis
#'
#' @param H coordinates of a composition. Either a matrix, a data.frame or a vector
#' @param basis basis used to calculate the coordinates
#' @return coordinates with respect the given basis
#' @seealso See functions \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}
#' to define different compositional basis.
#' See function \code{\link{coordinates}} to obtain details on how to calculate
#' coordinates of a given composition.
#' @export
composition = function(H, basis = NULL){
  rnames = rownames(H)

  class_type = class(H)
  if(is.null(basis) & "basis" %in% names(attributes(H))){
    basis = attr(H, 'basis')
  }
  if(is.null(basis)){
    stop("Basis is not defined", call. = FALSE)
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
      basis = as.matrix(ilr_basis(dim))
      RAW = exp(COORD %*% t(basis))
    }else{
      if(basis == 'alr'){
        RAW = cbind(exp(COORD), 1)
      }else{
        if(basis == 'clr'){
          dim = ncol(COORD)
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
    if(is.matrix(basis) | 'dgCMatrix' %in% class(basis)){
      RAW = exp(COORD %*% pinv(as.matrix(basis)))
    }else{
      stop(sprintf('Basis need to be either an string or a matrix'))
    }
  }
  RAW = RAW / rowSums(RAW)
  cnames = sprintf('c%d', 1:ncol(RAW))
  if(!is.character(basis) & !is.null(colnames(basis))){
    cnames = rownames(basis)
  }
  colnames(RAW) = cnames
  if(is_vector){
    RAW = RAW[1,]
  }
  if(is_data_frame){
    RAW = as.data.frame(RAW)
  }
  class(RAW) = setdiff(class_type, 'coda')
  suppressWarnings(row.names(RAW) <- row.names(H))
  #attr(RAW, 'basis') = basis
  RAW
}

#' @rdname composition
#' @export
comp = composition

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
