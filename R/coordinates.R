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

#' Build an isometric log-ratio basis
#'
#' @param dim number of components
#' @return matrix
#' @references
#' Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and Barceló-Vidal C. (2003).
#' \emph{Isometric logratio transformations for compositional data analysis}.
#' Mathematical Geology, \strong{35}(3) 279-300
#' @examples
#' ilr_basis(5)
#' @export
ilr_basis = function(dim){
  ilr_basis_default(dim)
}

#' Additive log-ratio basis
#'
#' Compute the transformation matrix to express a composition using the oblique additive log-ratio
#' coordinates.
#'
#' @param dim number of parts
#' @param denominator part used as denominator (default behaviour is to use last part)
#' @param numerator parts to be used as numerator. By default all except the denominator parts are chosen following original order.
#' @return matrix
#' @examples
#' alr_basis(5)
#' # Third part is used as denominator
#' alr_basis(5, 3)
#' # Third part is used as denominator, and
#' # other parts are rearranged
#' alr_basis(5, 3, c(1,5,2,4))
#' @references
#' Aitchison, J. (1986)
#' \emph{The Statistical Analysis of Compositional Data}.
#' Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). 416p.
#' @export
alr_basis = function(dim, denominator = dim, numerator = which(denominator != 1:dim)){
  res = alr_basis_default(dim)
  res = cbind(res, 0)
  if(dim != denominator){
    res[c(denominator, dim),] = res[c(dim, denominator),]
    res[,c(denominator, dim)] = res[,c(dim, denominator)]
  }
  res[,numerator]
}

#' Centered log-ratio basis
#'
#' Compute the transformation matrix to express a composition using
#' the linearly dependant centered log-ratio coordinates.
#'
#' @param dim number of parts
#' @return matrix
#' @references
#' Aitchison, J. (1986)
#' \emph{The Statistical Analysis of Compositional Data}.
#' Monographs on Statistics and Applied Probability. Chapman & Hall Ltd., London (UK). 416p.
#' @examples
#' (B <- clr_basis(5))
#' # CLR coordinates are linearly dependant coordinates.
#' (clr_coordinates <- coordinates(c(1,2,3,4,5), B))
#' # The sum of all coordinates equal to zero
#' sum(clr_coordinates) < 1e-15
#' @export
clr_basis = function(dim){
  clr_basis_default(dim)
}

#
# pc_basis = function(X){
#   lX =  log(X)
#   pr = stats::princomp(lX - rowMeans(lX))
#   pr$loadings[,-NCOL(X)]
# }

#' Build an \code{\link{ilr_basis}} using a sequential binary partition or
#' or a generic coordinate system based on balances.
#'
#' @param ... balances to consider
#' @param data composition from where name parts are extracted
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
    RIGHT = attr(stats::terms(part), 'term.labels')
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
    bal = stats::setNames(rep(0, length(names(data))), names(data))
    bal[balance[[1]]] = bal[balance[[1]]] + l
    bal[balance[[2]]] = bal[balance[[2]]] + r
    bal
  })
  if(!silent){
    if(qr(RES)$rank != NCOL(data)-1){
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

#' Isometric log-ratio basis based on Principal Balances.
#'
#' Different approximations to approximate the principal balances of a compositional dataset.
#'
#' @param X compositional dataset
#' @param method method to be used with Principal Balances. Methods available are: 'lsearch' or
#' method to be passed to hclust function (for example `ward.D` or `ward.D2` to use Ward method).
#' @param rep Number of restartings to be used with the local search algorithm. If zero is supplied
#' (default), one local search is performed using an starting point close to the principal component
#' solution.
#' @param ordering should the principal balances found be returned ordered? (first column, first
#' principal balance and so on)
#' @param ... parameters passed to hclust function
#' @return matrix
#' @references
#' Pawlowsky-Glahn, V., Egozcue, J.J., Tolosana-Delgado R. (2011).
#' \emph{Principal balances}.
#' in proceeding of the 4th International Workshop on Compositional Data Analysis (CODAWORK'11) (available online at \url{http://www-ma3.upc.edu/users/ortego/codawork11-Proceedings/Admin/Files/FilePaper/p55.pdf})
#' @examples
#' X = matrix(exp(rnorm(20*100)), nrow=100, ncol=20)
#' # Optimal variance obtained with Principal components
#' head(apply(coordinates(X, 'pc'), 2, var))
#'# Solution obtained using a hill climbing algorithm from pc approximation
#'head(apply(coordinates(X,pb_basis(X, method='lsearch')), 2, var))
#'# Solution obtained using a hill climbing algorithm using 10 restartings
#'head(apply(coordinates(X,pb_basis(X, method='lsearch', rep=10)), 2, var))
#' # Solution obtained using Ward method
#' head(apply(coordinates(X,pb_basis(X, method='ward.D2')), 2, var))
#' # Solution obtained using Old Ward function (in R versions <= 3.0.3)
#' head(apply(coordinates(X,pb_basis(X, method='ward.D')), 2, var))
#' @export
pb_basis = function(X, method, rep = 0, ordering = TRUE, ...){
  X = as.matrix(X)
  if(method == 'lsearch'){
    if(rep == 0){
      B = find_PB_pc_local_search(X)
    }else{
      B = find_PB_rnd_local_search(stats::cov(log(X)), rep=rep)
    }

  }else{
    hh = stats::hclust(stats::as.dist(variation_array(X, only_variation = TRUE)), method=method, ...)
    bin = hh$merge
    df = as.data.frame(X)
    names(df) = paste0('P.', 1:NCOL(df))
    nms = paste0('P',gsub('-','.', bin))
    dim(nms) = dim(bin)
    sbp = apply(nms, 1, paste, collapse='~')
    id = seq_along(sbp)
    sbp.exp = paste(sprintf("%s = %s ~ %s", paste0('P', id), nms[,1], nms[,2]),
                    collapse=', ')
    B = eval(parse(text = sprintf("sbp_basis(%s,data=df)", sbp.exp)))[,rev(id)]
  }
  if(ordering){
    B = B[,order(apply(coordinates(X, B), 2, stats::var), decreasing = TRUE)]
  }
  B
}

#' coordinates with respect an specific basis
#'
#' Calculate the coordinates of a composition with respect a given basis
#'
#' @param X compositional dataset. Either a matrix, a data.frame or a vector
#' @param basis basis used to calculate the coordinates. \code{basis} can be either a string or a matrix.
#' Accepted values for strings are: 'ilr' (default), 'clr', 'alr' and 'pc'. If \code{basis} is a matrix, it is expected
#' to have log-ratio basis given in columns.
#' @param label name given to the coordinates
#' @param sparse_basis Is the given matrix basis sparse? If TRUE calculation are carried
#' taking into an account sparsity (default `FALSE`)
#' @return coordinates with respect the given basis
#' @seealso See functions \code{\link{ilr_basis}}, \code{\link{alr_basis}},
#' \code{\link{clr_basis}}, \code{\link{sbp_basis}}
#' to define different compositional basis.
#' See function \code{\link{composition}} to obtain details on how to calculate
#' a compositions from given coordinates.
#' @examples
#' coordinates(c(1,2,3,4,5))
#' # Setting sparse_basi to TRUE can improve performance if log-ratio basis is sparse.
#' N = 100
#' K = 1000
#' X = matrix(exp(rnorm(N*K)), nrow=N, ncol=K)
#' system.time(coordinates(X, alr_basis(K), sparse_basis = FALSE))
#' system.time(coordinates(X, alr_basis(K), sparse_basis = TRUE))
#' system.time(coordinates(X, 'alr', sparse_basis = TRUE))
#' @export
coordinates = function(X, basis = 'ilr', label = 'x', sparse_basis = FALSE){
  class_type = class(X)
  is_vector = is.atomic(X) & !is.list(X) & !is.matrix(X)
  is_data_frame = inherits(X, 'data.frame')
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
      COORD = coordinates_basis(RAW, ilr_basis(dim), sparse = FALSE)
    }else{
      if(basis == 'alr'){
        basis = alr_basis(dim)
        COORD = coordinates_alr(RAW, 0)
      }else{
        if(basis == 'clr'){
          basis = clr_basis(dim)
          COORD = coordinates_basis(RAW, clr_basis(dim), sparse = FALSE)
        }else{
          if(basis == 'pc'){
            lRAW =  log(RAW)
            pr = stats::princomp(lRAW - rowMeans(lRAW))
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
      COORD = coordinates_basis(RAW, basis, sparse_basis)
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

#' coordinates with respect an specific basis
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
  }else{
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
          stop(sprintf('Basis %d not recognized'))
        }
      }
    }
  }else{
    if(is.matrix(basis)){
      RAW = exp(COORD %*% MASS::ginv(basis))
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
  class(RAW) = class_type
  #attr(RAW, 'basis') = basis
  RAW
}

#' Distance Matrix Computation (including Aitchison distance)
#'
#' This function overwrite \code{\link[stats]{dist}} function to contain Aitchison distance between
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
