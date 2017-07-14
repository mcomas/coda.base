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
#' sbp_basis(X,
#'           b1 = a~b+c+d+e+f+g,
#'           b2 = b~c+d+e+f+g,
#'           b3 = c~d+e+f+g,
#'           b4 = d~e+f+g,
#'           b5 = e~f+g,
#'           b6 = f~g)
#' @export
sbp_basis = function(X, ..., silent=F){
  sbp = list(...)

  nms = setdiff(names(sbp), "")
  if(length(nms) > 0){
    substitutions = lapply(sbp, all.vars)
    substitutions = substitutions[nms]
    while(!all(is.na(substitutions)) &
          !all(unlist(substitutions) %in% names(X))){
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
#' Calculate the coordinates with respect a fiven basis
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @param basis ilr base used to obtain coordinates
#' @param label name given to the coordinates
#' @return coordinates with respect the given basis
#' @examples
#' coordinates(c(1,2,3,4,5))
#' @export
coordinates = function(X, basis = ilr_basis(getDim(X)), label = 'x'){
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
  COORD = .Call('coda_base_coordinates', PACKAGE = 'coda.base', RAW, basis)

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
