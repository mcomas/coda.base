#' @title Coordinates basis
#'
#' @description
#' Obtain coordinates basis
#' @param H coordinates for which basis should be shown
#' @return basis used to create coordinates H
#' @export
basis = function(H){
  if(is.null(attr(H, 'basis'))) return(message('No basis specified'))
  attr(H, 'basis')
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
  if(type == 'cdp'){
    return(cdp_basis_(dim))
  }
  B = ilr_basis_default(dim)
  if(type == 'pivot'){
    return((-B)[,ncol(B):1, drop = FALSE][nrow(B):1,])
  }
  colnames(B) = sprintf("ilr%d", 1:ncol(B))
  rownames(B) = sprintf("c%d", 1:nrow(B))
  B
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
  B = clr_basis_default(dim)
  colnames(B) = sprintf("clr%d", 1:ncol(B))
  rownames(B) = sprintf("c%d", 1:nrow(B))
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
  B = res[,numerator, drop = FALSE]
  colnames(B) = sprintf("alr%d", 1:ncol(B))
  rownames(B) = sprintf("c%d", 1:nrow(B))
  B
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
  parts = colnames(X)
  if(is.null(parts)){
    parts = paste0('c', 1:nrow(B))
  }
  rownames(B) = parts
  colnames(B) = paste0('pc', 1:ncol(B))
  B
}

#' Isometric log-ratio basis based on canonical correlations
#'
#'
#' @param Y compositional dataset
#' @param X explanatory dataset
#' @return matrix
#'
#' @export
cc_basis = function(Y, X){
  Y = as.matrix(Y)
  X = cbind(X)
  B = ilr_basis(ncol(Y))
  cc = stats::cancor(coordinates(Y), X)
  B = B %*% cc$xcoef
  parts = colnames(Y)
  if(is.null(parts)){
    parts = paste0('c', 1:nrow(B))
  }
  rownames(B) = parts
  colnames(B) = paste0('cc', 1:ncol(B))
  B
}

#' Balance generated from the first canonical correlation component
#'
#'
#' @param Y compositional dataset
#' @param X explanatory dataset
#' @return matrix
#'
#' @export
cbalance_approx = function(Y,X){
  Y = as.matrix(Y)
  X = cbind(X)
  B = ilr_basis(ncol(Y))
  cc1 = B %*% stats::cancor(coordinates(Y), X)$xcoef[,1,drop=F]
  ord = order(abs(cc1))
  cb1_ = sign(cc1)
  cb1 = cb1_
  cor1 = abs(suppressWarnings(stats::cancor(coordinates(Y,sbp_basis(cb1_)), X)$cor))
  for(i in 1:(ncol(Y)-2)){
    cb1_[ord[i]] = 0
    cor1_ = abs(suppressWarnings(stats::cancor(coordinates(Y,sbp_basis(cb1_)), X)$cor))
    if(cor1_ > cor1){
      cb1 = cb1_
      cor1 = cor1_
    }
  }
  suppressWarnings(sbp_basis(cb1))
}

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
    return(do.call('sbp_basis', c(apply(P, 1, str_to_frm), list(data=df,
                                                                silent = silent)))) #, envir = as.environment('package:coda.base')
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

#' Isometric log-ratio basis based on Principal Balances.
#'
#' Exact method to calculate the principal balances of a compositional dataset. Different methods to approximate the principal balances of a compositional dataset are also included.
#'
#' @param X compositional dataset
#' @param method method to be used with Principal Balances. Methods available are: 'exact', 'constrained' or 'cluster'.
#' @param constrained.complete_up When searching up, should the algorithm try to find possible siblings for the current balance (TRUE) or build a parent directly forcing current balance to be part of the next balance (default: FALSE). While the first is more exhaustive and given better results the second is faster and can be used with highe dimensional datasets.
#' @param cluster.method Method to be used with the hclust function (default: `ward.D2`) or any other method available in  hclust function
#' @param ordering should the principal balances found be returned ordered? (first column, first
#' principal balance and so on)
#' @param ... parameters passed to hclust function
#' @return matrix
#' @references
#' Martín-Fernández, J.A., Pawlowsky-Glahn, V., Egozcue, J.J., Tolosana-Delgado R. (2018).
#' Advances in Principal Balances for Compositional Data.
#' \emph{Mathematical Geosciencies}, 50, 273-298.
#' @examples
#' set.seed(1)
#' X = matrix(exp(rnorm(5*100)), nrow=100, ncol=5)
#'
#' # Optimal variance obtained with Principal components
#' (v1 <- apply(coordinates(X, 'pc'), 2, var))
#' # Optimal variance obtained with Principal balances
#' (v2 <- apply(coordinates(X,pb_basis(X, method='exact')), 2, var))
#' # Solution obtained using constrained method
#' (v3 <- apply(coordinates(X,pb_basis(X, method='constrained')), 2, var))
#' # Solution obtained using Ward method
#' (v4 <- apply(coordinates(X,pb_basis(X, method='cluster')), 2, var))
#'
#' # Plotting the variances
#' barplot(rbind(v1,v2,v3,v4), beside = TRUE, ylim = c(0,2),
#'         legend = c('Principal Components','PB (Exact method)',
#'                    'PB (Constrained)','PB (Ward approximation)'),
#'         names = paste0('Comp.', 1:4), args.legend = list(cex = 0.8), ylab = 'Variance')
#'
#' @export
pb_basis = function(X, method, constrained.complete_up = FALSE, cluster.method = 'ward.D2',
                    ordering = TRUE, ...){
  X = as.matrix(X)
  if(!(all(X > 0))){
    stop("All components must be strictly positive.", call. = FALSE)
  }
  if(method %in% c('constrained', 'exact')){
    if(method == 'exact'){
      M = 'PB'
      B = find_PB(X)
    }
    if(method == 'constrained'){
      M = 'CS'
      # B = t(fBalChip(X)$bal)
      B = find_PB_using_pc_recursively_forcing_parents(X)
    }
    if(method == 'constrained2'){
      M = 'CS'
      B = find_PB_using_pc(X)
    }
    # if(method == 'lsearch'){
    #   if(rep == 0){
    #     B = find_PB_pc_local_search(X)
    #   }else{
    #     B = find_PB_rnd_local_search(stats::cov(log(X)), rep=rep)
    #   }
    # }
  }else if(method == 'cluster'){
    M = 'CL'
    # Passing arguments to hclust function
    hh = stats::hclust(stats::as.dist(variation_array(X, only_variation = TRUE)), method=cluster.method, ...)
    B = matrix(0, ncol = nrow(hh$merge), nrow = ncol(X))
    for(i in 1:nrow(hh$merge)){
      if(hh$merge[i,1] < 0 & hh$merge[i,2] < 0){
        B[-hh$merge[i,],i] = c(-1,+1)
      }else{
        if(hh$merge[i,1] > 0){
          B[B[,hh$merge[i,1]] != 0,i] = -1
        }else{
          B[-hh$merge[i,1],i] = -1
        }
        if(hh$merge[i,2] > 0){
          B[B[,hh$merge[i,2]] != 0,i] = +1
        }else{
          B[-hh$merge[i,2],i] = +1
        }
      }
    }
    B = sbp_basis(B[,nrow(hh$merge):1, drop = FALSE])
    # bin = hh$merge
    # df = as.data.frame(X)
    # names(df) = paste0('P.', 1:NCOL(df))
    # nms = paste0('P',gsub('-','.', bin))
    # dim(nms) = dim(bin)
    # sbp = apply(nms, 1, paste, collapse='~')
    # id = seq_along(sbp)
    # sbp.exp = paste(sprintf("%s = %s ~ %s", paste0('P', id), nms[,1], nms[,2]),
    #                 collapse=', ')
    # B = eval(parse(text = sprintf("sbp_basis(%s,data=df)", sbp.exp)))[,rev(id), drop = FALSE]
  } else{
    stop(sprintf("Method %s does not exist", method))
  }
  if(ordering){
    B = B[,order(apply(coordinates(X, B, basis_return = FALSE), 2, stats::var), decreasing = TRUE), drop = FALSE]
  }
  parts = colnames(X)
  if(is.null(parts)){
    parts = paste0('c', 1:nrow(B))
  }
  rownames(B) = parts
  colnames(B) = paste0('pb', 1:ncol(B))
  B
}

#' Isometric log-ratio basis based on Balances.
#'
#' The function return default balances used in CoDaPack software.
#'
#' @param dim dimension to build the ILR basis based on balanced balances
#' @return matrix
#' @export
cdp_basis = function(dim){
  B = cdp_basis_(dim)
  B
}

cdp_basis_ = function(dim, wR = 1:ceiling(dim/2), wL = ceiling(dim/2) + 1:floor(dim/2)){
  R = length(wR)
  L = length(wL)
  D = R + L
  v = rep(0, dim)
  v[wR] = +sqrt(L/R/D)
  v[wL] = -sqrt(R/L/D)
  if(R == 1 & L == 1){
    return(v)
  }
  if(R == 1){
    return(cbind(v,
                 Recall(dim, wR = wL[1:ceiling(L/2)], wL = wL[ceiling(L/2) + 1:floor(L/2)])))
  }
  if(L == 1){
    return(cbind(v,
                 Recall(dim, wR = wR[1:ceiling(R/2)], wL = wR[ceiling(R/2) + 1:floor(R/2)])))
  }
  cbind(v,
        Recall(dim, wR = wR[1:ceiling(R/2)], wL = wR[ceiling(R/2) + 1:floor(R/2)]),
        Recall(dim, wR = wL[1:ceiling(L/2)], wL = wL[ceiling(L/2) + 1:floor(L/2)]))
}
