getDim = function(X) ifelse(is.vector(X), length(X), NCOL(X))

coda_basis = function(){

}

#' Default isometric log-ratio basis
#'
#' @param dim number of components
#' @return matrix
#' @examples
#' ilr_basis(5)
#' @export
ilr_basis = function(dim){
  .Call('coda_base_ilr_basis_default', PACKAGE = 'coda.base', dim)
}

#' Default additive log-ratio basis
#'
#' @param dim number of components
#' @return matrix
#' @examples
#' alr_basis(5)
#' @export
alr_basis = function(dim, denominator = dim, numerator = which(denominator != 1:dim)){
  res = .Call('coda_base_alr_basis_default', PACKAGE = 'coda.base', dim)
  res = cbind(res, 0)
  if(dim != denominator){
    res[c(denominator, dim),] = res[c(dim, denominator),]
    res[,c(denominator, dim)] = res[,c(dim, denominator)]
  }
  res[,numerator]
}

#' Default centered log-ratio basis
#'
#' @param dim number of components
#' @return matrix
#' @examples
#' clr_basis(5)
#' @export
clr_basis = function(dim){
  .Call('coda_base_clr_basis_default', PACKAGE = 'coda.base', dim)
}

#' Principal components orthonormal basis
#'
#' @param X pc_basis
#' @return matrix
#' @examples
#' X = data.frame(a=exp(rnorm(10)), b=exp(rnorm(10)),
#'                c=exp(rnorm(10)), d=exp(rnorm(10)),
#'                e=exp(rnorm(10)), f=exp(rnorm(10)))
#' pc_basis(X)
#' @export
pc_basis = function(X){
  lX =  log(X)
  pr = princomp(lX - rowMeans(lX))
  pr$loadings[,-NCOL(X)]
}

#' Define basis using binary partitions.
#'
#' @param X composition from where to extract parts names
#' @param ... balances to consider
#' @param silent inform about orthgonality
#' @return matrix
#' @examples
#' X = data.frame(a=1:2, b=2:3, c=4:5, d=5:6, e=10:11, f=100:101, g=1:2)
#' sbp_basis(b1 = a~b+c+d+e+f+g,
#'           b2 = b~c+d+e+f+g,
#'           b3 = c~d+e+f+g,
#'           b4 = d~e+f+g,
#'           b5 = e~f+g,
#'           b6 = f~g, data = X)
#' sbp_basis(b1 = a~b,
#'          b2 = b1~c,
#'          b3 = b2~d,
#'          b4 = b3~e,
#'          b5 = b4~f,
#'          b6 = b5~g, data = X)
#' # A non-orthogonal basis can also be calculated.
#' sbp_basis(b1 = a+b+c~e+f+g,
#'           b2 = d~a+b+c,
#'           b3 = d~e+g,
#'           b4 = a~e+b,
#'           b5 = b~f,
#'           b6 = c~g, data = X)
#' @export
sbp_basis = function(..., data, silent=F){
  if (!is.data.frame(data) && !is.environment(data) && !is.null(attr(data, "class")))
    data <- as.data.frame(data)
  else if (is.array(data))
    stop("'data' must be a data.frame, not a matrix or an array")

  sbp = list(...)
  if(!all(unlist(lapply(sbp, all.vars)) %in% c(names(data), names(sbp)))){
    stop("Balances should be columns of 'data'")
  }
  nms = setdiff(names(sbp), "")
  if(length(nms) > 0){
    substitutions = lapply(sbp, all.vars)
    substitutions = substitutions[nms]
    while(!all(is.na(substitutions)) &
          !all(unlist(substitutions) %in% names(data))){
      for(nm in nms){
        substitutions = lapply(substitutions, function(subs){
          I = match(nm, subs)
          if(!is.na(I)){
            c(subs[-I], substitutions[[nm]])
          }else{
            subs
          }
        })
      }
    }
  }

  sbp_split = function(part){
    RIGHT = attr(terms(part), 'term.labels')
    LEFT = setdiff(all.vars(part), RIGHT)
    if(length(nms) > 0){
      for(nm in nms){
        I = match(nm, RIGHT)
        if(!is.na(I)){
          RIGHT = c(RIGHT[-I], substitutions[[nm]])
        }
        I = match(nm, LEFT)
        if(!is.na(I)){
          LEFT = c(LEFT[-I], substitutions[[nm]])
        }
      }
    }
    list(LEFT, RIGHT)
  }

  sbp_clean = lapply(sbp, sbp_split)
  RES = sapply(sbp_clean, function(balance){
    I1 = length(balance[[1]])
    I2 = length(balance[[2]])
    l = +1/I1 * sqrt(I1*I2/(I1+I2))
    r = -1/I2 * sqrt(I1*I2/(I1+I2))
    bal = setNames(rep(0, length(names(X))), names(X))
    bal[balance[[1]]] = bal[balance[[1]]] + l
    bal[balance[[2]]] = bal[balance[[2]]] + r
    bal
  })
  if(!silent){
    if(qr(RES)$rank != NCOL(X)-1){
      warning('Given partition is not a basis')
    }else{
      Z = t(RES) %*% RES
      if( !all( Z - diag(diag(Z), nrow=NROW(Z), ncol=NCOL(Z)) < 10e-10 ) ){
        warning('Given basis is not orthogonal')
      }else{
        if(!all( Z - diag(1, nrow=NROW(Z), ncol=NCOL(Z)) < 10e-10 )){
          warning('Given basis is not orthonormal')
        }
      }
    }
  }
  RES
}

#'
#' coordinates with respect an specific basis
#'
#' Calculate the coordinates with respect a given basis
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @param basis ilr base used to obtain coordinates
#' @param label name given to the coordinates
#' @param sparse_basis Is the given matrix basis sparse? If TRUE calculation are carried
#' taking into an account sparsity (default `FALSE`)
#' @return coordinates with respect the given basis
#' @seealso See functions \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}, \code{\link{pc_basis}}
#' to define different compositional basis.
#' @examples
#' coordinates(c(1,2,3,4,5))
#' # Setting sparse_basi to TRUE can improve performance if log-ratio basis is sparse.
#' N = 100
#' K = 1000
#' X = matrix(exp(rnorm(N*K)), nrow=N, ncol=K)
#' system.time(coordinates(X, alr_basis(K), sparse_basis = FALSE))
#' system.time(coordinates(X, alr_basis(K), sparse_basis = TRUE))
#' system.time(coordinates(X, 'alr', sparse_basis = TRUE))
#'
#' @export
coordinates = function(X, basis = 'ilr', label = 'x', sparse_basis = FALSE){
  class_type = class(X)
  is_vector = is.vector(X)
  is_data_frame = is.data.frame(X)
  RAW = X
  if(is_vector){
    class_type = 'double'
    RAW = matrix(X, nrow=1)
  }
  if(is_data_frame){
    RAW = as.matrix(X)
  }
  if(is.character(basis)){
    dim = ncol(RAW)
    if(basis == 'ilr'){
      basis = ilr_basis(dim)
      COORD = .Call('coda_base_coordinates_basis', PACKAGE = 'coda.base', RAW, ilr_basis(dim), sparse = FALSE)
    }else{
      if(basis == 'alr'){
        basis = alr_basis(dim)
        COORD = .Call('coda_base_coordinates_alr', PACKAGE = 'coda.base', RAW, 0)
      }else{
        if(basis == 'clr'){
          basis = clr_basis(dim)
          COORD = .Call('coda_base_coordinates_basis', PACKAGE = 'coda.base', RAW, clr_basis(dim), sparse = FALSE)
        }else{
          if(basis == 'pc'){
            lRAW =  log(RAW)
            pr = princomp(lRAW - rowMeans(lRAW))
            basis = pr$loadings[,-dim]
            COORD = pr$scores[,-dim]
          }else{
            stop(sprintf('Basis %d not recognized'))
          }
        }
      }
    }
  }else{
    if(is.matrix(basis)){
      COORD = .Call('coda_base_coordinates_basis', PACKAGE = 'coda.base', RAW, basis, sparse_basis)
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
  COORD
}
