#' Isometric/Orthonormal Log-Ratio Basis for Log-Transformed Compositions
#'
#' Builds an isometric log-ratio (ilr) basis for a composition with \code{k+1} parts, also called orthonormal log-ratio (olr) basis.
#'
#' The basis vectors are constructed as:
#' \deqn{h_i = \sqrt{\frac{i}{i+1}} \log\frac{\sqrt[i]{\prod_{j=1}^i x_j}}{x_{i+1}}}{%
#' h[i] = sqrt(i/(i+1)) * ( log(x[1] * ... * x[i]) / i - log(x[i+1]) )}
#' for \eqn{i = 1, \ldots, k}.
#'
#' Setting the \code{type} parameter to \code{"pivot"} (pivot balances) or \code{"cdp"} (codapack balances) allows generating alternative ilr/olr bases.
#'
#' @param dim An integer indicating the number of components.
#'            If a dataframe or matrix is provided, the number of components is inferred from the number of columns. If a character vector specifying the names of the parts is provided the number of component is its length.
#' @param type Character string specifying the type of basis to generate.
#'             Options are \code{"pivot"}, \code{"cdp"}. Any other option will return the Helmert basis defined by Egozcue et al., 2013..
#' @return A matrix representing the orthonormal basis.
#' @references
#' Egozcue, J.J., Pawlowsky-Glahn, V., Mateu-Figueras, G., & Barceló-Vidal, C. (2003).
#' \emph{Isometric logratio transformations for compositional data analysis}.
#' Mathematical Geology, \strong{35}(3), 279–300.
#' @examples
#' ilr_basis(5)
#' ilr_basis(alimentation[,1:9])
#' @export
ilr_basis = function(dim, type = 'default'){
  parts = colnames(dim)
  if(is.character(dim)){
    parts = dim
    dim = length(dim)
  }
  D = check_dim(dim)
  if(type == 'cdp'){
    B = cdp_basis_(D)
  }else{
    B = ilr_basis_default(D)
    if(type == 'pivot'){
      B = (-B)[,ncol(B):1, drop = FALSE][nrow(B):1,]
    }
  }
  rownames(B) = sprintf("c%d", 1:nrow(B))
  if(!is.null(parts)) rownames(B) = parts
  colnames(B) = sprintf("ilr%d", 1:ncol(B))
  B
}

#' @rdname ilr_basis
#' @export
olr_basis = function(dim, type = 'default'){
  B = ilr_basis(dim, type)
  colnames(B) = sprintf("olr%d", 1:ncol(B))
  B
}

cdp_basis = function(dim){
  parts = colnames(dim)
  if(is.character(dim)){
    parts = dim
    dim = length(dim)
  }
  D = check_dim(dim)
  B = cdp_basis_(D)
  rownames(B) = paste0("c", 1:D)
  if(!is.null(parts)) rownames(B) = parts
  colnames(B) = paste0("ilr", 1:ncol(B))
  B
}

#' Centered log-ratio basis
#'
#' Compute the transformation matrix to express a composition using
#' the linearly dependant centered log-ratio coordinates.
#'
#' @param dim An integer indicating the number of components.
#'            If a dataframe or matrix is provided, the number of components is inferred from the number of columns. If a character vector specifying the names of the parts is provided the number of component is its length.
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
  parts = colnames(dim)
  if(is.character(dim)){
    parts = dim
    dim = length(dim)
  }
  D = check_dim(dim)
  B = clr_basis_default(D)
  colnames(B) = sprintf("clr%d", 1:ncol(B))
  rownames(B) = sprintf("c%d", 1:nrow(B))
  if(!is.null(parts)) rownames(B) = parts
  B
}


#' Additive log-ratio basis
#'
#' Compute the transformation matrix to express a composition using the oblique additive log-ratio
#' coordinates.
#'
#' @param dim An integer indicating the number of components.
#'            If a dataframe or matrix is provided, the number of components is inferred from the number of columns. If a character vector specifying the names of the parts is provided the number of component is its length.
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
alr_basis = function(dim, denominator = NULL, numerator = NULL){
  parts = colnames(dim)
  if(is.character(dim)){
    parts = dim
    dim = length(dim)
  }
  D = check_dim(dim)
  if(is.null(denominator)) denominator = D
  if(is.null(numerator)) numerator = which(1:D != denominator)
  res = alr_basis_default(D)
  res = cbind(res, 0)
  if(D != denominator){
    res[c(denominator, D),] = res[c(D, denominator),, drop = FALSE]
    res[,c(denominator, D)] = res[,c(D, denominator), drop = FALSE]
  }
  B = res[,numerator, drop = FALSE]
  colnames(B) = sprintf("alr%d", 1:ncol(B))
  rownames(B) = sprintf("c%d", 1:nrow(B))
  if(!is.null(parts)) rownames(B) = parts
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
  B = ilr_basis(ncol(X))
  B = B %*% svd(scale(log(as.matrix(X)) %*% B, scale=FALSE))$v

  parts = colnames(X)
  if(is.null(parts)){
    parts = paste0('c', 1:nrow(B))
  }
  rownames(B) = parts
  colnames(B) = paste0('pc', 1:ncol(B))
  as.matrix(B)
}

#' Isometric Log-Ratio Basis Based on Canonical Correlations
#'
#' Constructs an isometric log-ratio (ilr) basis for a compositional dataset,
#' optimized with respect to canonical correlations with an explanatory dataset.
#'
#' @param Y A compositional dataset (matrix or data frame).
#' @param X An explanatory dataset (matrix or data frame).
#' @return A matrix representing the isometric log-ratio basis.
#' @export
cc_basis = function(Y, X){
  Y = as.matrix(Y)
  X = cbind(X)
  B = ilr_basis_default(ncol(Y))
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



#' Isometric log-ratio basis based on Balances
#'
#' Build an \code{\link{ilr_basis}} using a sequential binary partition or
#' a generic coordinate system based on balances.
#'
#' @param sbp parts to consider in the numerator and the denominator. Can be
#' defined either using a list of formulas setting parts (see examples) or using
#' a matrix where each column define a balance. Positive values are parts in
#' the numerator, negative values are parts in the denominator, zeros are parts
#' not used to build the balance.
#' @param data composition from where name parts are extracted
#' @param fill should the balances be completed to become an orthonormal basis?
#'  if the given balances are not orthonormal, the function will complete the
#'  balance to become a basis.
#' @param silent inform about orthogonality
#' @return matrix
#' @examples
#' X = data.frame(a=1:2, b=2:3, c=4:5, d=5:6, e=10:11, f=100:101, g=1:2)
#' sbp_basis(list(b1 = a~b+c+d+e+f+g,
#'                b2 = b~c+d+e+f+g,
#'                b3 = c~d+e+f+g,
#'                b4 = d~e+f+g,
#'                b5 = e~f+g,
#'                b6 = f~g), data = X)
#' sbp_basis(list(b1 = a~b,
#'                b2 = b1~c,
#'                b3 = b2~d,
#'                b4 = b3~e,
#'                b5 = b4~f,
#'                b6 = b5~g), data = X)
#' # A non-orthogonal basis can also be calculated.
#' sbp_basis(list(b1 = a+b+c~e+f+g,
#'               b2 = d~a+b+c,
#'               b3 = d~e+g,
#'               b4 = a~e+b,
#'               b5 = b~f,
#'               b6 = c~g), data = X)
#' @export
sbp_basis = function(sbp, data = NULL, fill = FALSE, silent=FALSE){
  if(is.null(data) & is.matrix(sbp)){
    # P = t(sbp)
    df = as.data.frame(matrix(1, nrow(sbp), nrow = 1))
    if(!is.null(rownames(sbp))){
      colnames(df) = rownames(sbp)
    }
    str_to_frm = function(vec){
      frm = paste(stats::aggregate(nm ~ vec, subset(data.frame(nm = paste0('`',names(df), '`'), vec = -1 * vec,
                                                               stringsAsFactors = FALSE), vec != 0),
                                   FUN = paste, collapse= ' + ')[['nm']], collapse=' ~ ')
      stats::as.formula(frm)
    }
    return(sbp_basis(apply(sbp, 2, str_to_frm),
                     data = df,
                     fill = fill,
                     silent = silent))
    # return(do.call('sbp_basis', c(apply(P, 1, str_to_frm), list(data=df,
                                                                # fill = fill,
                                                                # silent = silent)))) #, envir = as.environment('package:coda.base')
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
  if(fill){
    return(Recall(fill_sbp(sign(RES))))
  }
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
pb_basis = function(X, method, cluster.method = 'ward.D2',
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
      B = constrained_pb(as.matrix(X))
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
    hh = stats::hclust(stats::as.dist(variation_array(X)), method=cluster.method, ...)
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
  } else{
    stop(sprintf("Method %s does not exist", method))
  }
  if(ordering){
    B = B[,order(apply(coordinates(X, B), 2, stats::var),
                 decreasing = TRUE),
          drop = FALSE]
  }
  parts = colnames(X)
  if(is.null(parts)){
    parts = paste0('c', 1:nrow(B))
  }
  rownames(B) = parts
  colnames(B) = paste0('pb', 1:ncol(B))
  B
}



#' Pairwise log-ratio generator system
#'
#' The function returns all combinations of pairs of log-ratios.
#'
#' @param dim An integer indicating the number of components.
#'            If a dataframe or matrix is provided, the number of components is inferred from the number of columns. If a character vector specifying the names of the parts is provided the number of component is its length.
#' @return matrix
#' @export
pairwise_basis = function(dim){
  parts = colnames(dim)
  if(is.character(dim)){
    parts = dim
    dim = length(dim)
  }
  D = check_dim(dim)
  I = utils::combn(D,2)
  B = apply(I, 2, function(i){
    b = rep(0, D)
    b[i] = c(1,-1)
    b
  })
  colnames(B) = paste0('pw', apply(I, 2, paste, collapse = '_'))
  rownames(B) = paste0("c", 1:D)
  if(!is.null(parts)) rownames(B) = parts
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

check_dim = function(dim){
  if(!is.null(ncol(dim))){
    dim = ncol(dim)
  }
  if(!is.numeric(dim)){
    stop("Dimension should be a number", call. = FALSE)
  }
  if(!dim>1){
    stop("Dimension should be at least 2", call. = FALSE)
  }
  dim
}
