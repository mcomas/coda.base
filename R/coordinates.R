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
#' @param ... components of the compositional data
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
#' # Default ilr given by ilr_basis(5) is given
#' coordinates(1:5)
#' B = ilr_basis(5)
#' coordinates(1:5, B)
#' @export
coordinates = function(X, basis = 'ilr'){
  if(is.matrix(X)){
    if(!is.numeric(X)){
      stop("Composition must be numeric", call. = FALSE)
    }
    if(is.character(basis)){   # default's basis with characters
      COORD = get(basis)(X)
    }else{                      # matrix basis
      COORD = coordinates.matrix(X, basis)
      colnames(COORD) = colnames(basis)
      if(is.null(colnames(COORD))){
        colnames(COORD) = sprintf("h%d", 1:ncol(COORD))
      }
    }
  }else{
    if(is.atomic(X) & !is.list(X)){ # vector
      COORD = Recall(matrix(X, nrow=1), basis)
      COORD = COORD[1,]
    }else{
      class_type = class(X)
      if(inherits(X, 'data.frame')){
        if(!all(sapply(X, is.numeric))) stop("All parts must be numeric", call. = FALSE)
        mCOORD = Recall(as.matrix(X), basis)
        COORD = as.data.frame(mCOORD)
      }
      class(COORD) = class_type
    }
  }
  suppressWarnings(row.names(COORD) <- row.names(X))
  COORD
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
composition = function(H, basis = 'ilr'){
  if(is.matrix(H)){
    D = ncol(H) + 1
    if(!is.numeric(H)){
      stop("Coordinates must be numeric", call. = FALSE)
    }
    P = NULL
    if(is.character(basis)){   # default's basis with characters
      if(basis == 'ilr' | basis == 'olr') P = exp(H %*% t(ilr_basis(D)))
      if(basis == 'alr') P = cbind(exp(H), 1)
      if(basis == 'clr') P = exp(H)
      if(is.null(P)){
        stop(sprintf("Basis '%s' not recognized", basis))
      }
    }else{
      if(is.matrix(basis)){
        P = exp(as.matrix(H) %*% pinv(basis))
      }else{
        stop(sprintf('Basis need to be either a string or a matrix'))
      }
    }
    P = P / rowSums(P)
    colnames(P) = rownames(basis)
  }else{
    if(is.atomic(H) & !is.list(H)){ # vector
      P = Recall(matrix(H, nrow=1), basis)[1,]
    }else{
      class_type = class(H)
      if(inherits(H, 'data.frame')){
        if(!all(sapply(H, is.numeric))) stop("All parts must be numeric", call. = FALSE)
        P = as.data.frame(Recall(as.matrix(H), basis))
      }
      class(P) = class_type
    }
  }
  suppressWarnings(row.names(P) <- row.names(H))
  P
}



alr = function(X){
  COORD = alr_coordinates(X, ncol(X))
  colnames(COORD) = paste0('alr', 1:ncol(COORD))
  COORD
}

ilr = function(X){
  COORD = coordinates(X, ilr_basis(ncol(X)))
  colnames(COORD) = paste0('ilr', 1:ncol(COORD))
  COORD
}

olr = function(X){
  COORD = coordinates(X, ilr_basis(ncol(X)))
  colnames(COORD) = paste0('olr', 1:ncol(COORD))
  COORD
}

clr = function(X){
  COORD = clr_coordinates(X)
  colnames(COORD) = paste0('clr', 1:ncol(COORD))
  COORD
}

pc = function(X){
  B = ilr_basis(ncol(X))
  B = B %*% svd(scale(log(X) %*% B, scale=FALSE))$v
  COORD = matrix_coordinates(X, B)
  colnames(COORD) = paste0('pc', 1:ncol(B))
  COORD
}

cdp = function(X){
  B = cdp_basis(ncol(X))
  COORD = matrix_coordinates(X, B)
  colnames(COORD) = colnames(B)
  COORD
}
pb = function(X){
  B = pb_basis(X, method = 'exact')
  COORD = matrix_coordinates(X, B)
  colnames(COORD) = colnames(B)
  COORD
}

pw = function(X){
  B = pairwise_basis(ncol(X))
  COORD = matrix_coordinates(X, B)
  colnames(COORD) = colnames(B)
  COORD
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
      if(is.character(lpars[[2]])){
        coordinates(lpars[[1]], basis = lpars[[2]])
      }else{
        stop("Please specify second argument", call. = FALSE)
      }
    }

  }else{
    if(sum(inum) == 1){
      stop("Compositions should have at least two parts. If you want to calculate the coordinates of one composition, please use function coordinates().", call. = FALSE)
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

#' @rdname coordinates
#' @export
alr_c = function(X){
  coordinates(X, 'alr')
}

#' @rdname coordinates
#' @export
clr_c = function(X){
  coordinates(X, 'clr')
}

#' @rdname coordinates
#' @export
ilr_c = function(X){
  coordinates(X, 'ilr')
}

#' @rdname coordinates
#' @export
olr_c = function(X){
  coordinates(X, 'olr')
}

#' @rdname composition
#' @export
comp = composition

coordinates.matrix = function(X, basis){
  if(inherits(basis, "sparseMatrix")){
    sparse_coordinates(X, basis)
  }else{
    matrix_coordinates(X, basis)
  }
}
