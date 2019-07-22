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

#' Default Isometric log-ratio basis
#'
#' Build an isometric log-ratio basis for a composition with k+1 parts
#' \deqn{h_i = \sqrt{\frac{i}{i+1}} \log\frac{\sqrt[i]{\prod_{j=1}^i x_j}}{x_{i+1}}}{%
#' h[i] = \sqrt(i/(i+1)) ( log(x[1] \ldots x[i])/i - log(x[i+1]) )}
#' for \eqn{i in 1\ldots k}.
#'
#'Modifying parameter type (pivot or cdp) other ilr basis can be generated
#'
#' @param dim number of components
#' @param type if different than `pivot` (pivot balances) or `cdp` (codapack balances) default balances are returned, which computes a triangular Helmert matrix as defined by Egozcue et al., 2013.
#' @return matrix
#' @references
#' Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G. and Barceló-Vidal C. (2003).
#' \emph{Isometric logratio transformations for compositional data analysis}.
#' Mathematical Geology, \strong{35}(3) 279-300
#' @examples
#' ilr_basis(5)
#' @export
ilr_basis = function(dim, type = 'default'){
  B = ilr_basis_default(dim)
  if(type == 'pivot'){
    return((-B)[,ncol(B):1, drop = FALSE][nrow(B):1,])
  }
  if(type == 'cdp'){
    return(sbp_basis(cdp_partition(dim)))
  }
  B
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
    res[c(denominator, dim),] = res[c(dim, denominator),, drop = FALSE]
    res[,c(denominator, dim)] = res[,c(dim, denominator), drop = FALSE]
  }
  res[,numerator, drop = FALSE]
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

#' Isometric log-ratio basis based on Balances
#' Build an \code{\link{ilr_basis}} using a sequential binary partition or
#' a generic coordinate system based on balances.
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
sbp_basis = function(..., data = NULL, silent=F){
  sbp = list(...)
  if(is.null(data) & is.matrix(sbp[[1]])){
    P = t(sbp[[1]])
    df = as.data.frame(matrix(1, ncol(P), nrow = 1))
    str_to_frm = function(vec){
      frm = paste(stats::aggregate(nm ~ vec, subset(data.frame(nm = paste0('`',names(df), '`'), vec = -1 * vec,
                                                               stringsAsFactors = FALSE), vec != 0),
                                   FUN = paste, collapse= ' + ')[['nm']], collapse=' ~ ')
      stats::as.formula(frm)
    }
    return(do.call('sbp_basis', c(apply(P, 1, str_to_frm), list(silent=silent), list(data=df)))) #, envir = as.environment('package:coda.base')
  }

  if (!is.data.frame(data) && !is.environment(data) && ( (is.matrix(data) && !is.null(colnames(data))) | !is.null(attr(data, "class"))))
    data <- as.data.frame(data)
  else if (is.array(data))
    stop("'data' must be a data.frame or a matrix with column names")


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

#' Isometric log-ratio basis based on Principal Components.
#'
#' Different approximations to approximate the principal balances of a compositional dataset.
#'
#' @param X compositional dataset
#' @return matrix
#'
#' @export
pc_basis = function(X){
  X = as.matrix(X)
  lX =  log(X)
  SVD = svd(scale(lX - rowMeans(lX), scale = FALSE))
  B = SVD$v[,-ncol(X), drop = FALSE]
  B
}

#' Isometric log-ratio basis based on Principal Balances.
#'
#' Exact method to calculate the principal balances of a compositional dataset. Different methods to approximate the principal balances of a compositional dataset are also included.
#'
#' @param X compositional dataset
#' @param method method to be used with Principal Balances. Methods available are: 'exact', 'lsearch' or
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
#' set.seed(1)
#' X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
#' # Optimal variance obtained with Principal components
#' (v1 <- apply(coordinates(X, 'pc'), 2, var))
#' # Optimal variance obtained with Principal balances
#' (v2 <- apply(coordinates(X,pb_basis(X, method='exact')), 2, var))
#' # Solution obtained using a hill climbing algorithm from pc approximation
#' apply(coordinates(X,pb_basis(X, method='lsearch')), 2, var)
#' # Solution obtained using a hill climbing algorithm using 10 restartings
#' apply(coordinates(X,pb_basis(X, method='lsearch', rep=10)), 2, var)
#' # Solution obtained using Ward method
#' (v3 <- apply(coordinates(X,pb_basis(X, method='ward.D2')), 2, var))
#' # Solution obtained using Old Ward function (in R versions <= 3.0.3)
#' apply(coordinates(X,pb_basis(X, method='ward.D')), 2, var)
#' # Plotting the variances
#' barplot(rbind(v1,v2,v3), beside = TRUE,
#'         legend = c('Principal Components','PB (Exact method)','PB (Ward approximation)'),
#'         names = paste0('Comp.', 1:4), args.legend = list(cex = 0.8), ylab = 'Variance')
#'
#' @export
pb_basis = function(X, method, rep = 0, ordering = TRUE, ...){
  X = as.matrix(X)
  if(!(all(X > 0))){
    stop("All components must be strictly positive.", call. = FALSE)
  }
  if(method %in% c('lsearch', 'exact')){
    if(method == 'exact'){
      B = find_PB(X)
    }
    if(method == 'lsearch'){
      if(rep == 0){
        B = find_PB_pc_local_search(X)
      }else{
        B = find_PB_rnd_local_search(stats::cov(log(X)), rep=rep)
      }
    }
  }else{
    # Passing arguments to hclust function
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
    B = eval(parse(text = sprintf("sbp_basis(%s,data=df)", sbp.exp)))[,rev(id), drop = FALSE]
  }
  if(ordering){
    B = B[,order(apply(coordinates(X, B), 2, stats::var), decreasing = TRUE), drop = FALSE]
  }
  B
}

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
