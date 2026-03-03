#' Plot a balance with node labels under horizontal branches
#'
#' @param B Balance basis matrix
#' @param data Optional compositional data used to compute balance summaries
#' @param main Plot title
#' @param summary_fun Optional function applied to each balance coordinate vector.
#'   It must take a numeric vector and return a character string.
#' @param cex_node Character expansion for node labels
#' @param offset_node Vertical offset below the horizontal branch, relative to max height
#' @param ... Further arguments passed to plot
#' @return Invisibly returns a data.frame with node coordinates and labels
#'
#' @examples
#' X = waste[,5:9]
#' B = pb_basis(X, method = 'exact')
#'
#' plot_balance(B)
#'
#' plot_balance(B, data = X,
#'              summary_fun = function(x){
#'                q = quantile(x, probs = c(0.25, 0.5, 0.75))
#'                sprintf("%0.2f [%0.2f-%0.2f]", q[2], q[1], q[3])
#'              })
#'
#' @export
plot_balance = function(B, data = NULL,
                        main = "Balance dendrogram",
                        summary_fun = NULL,
                        cex_node = 0.9,
                        offset_node = 0.05,
                        ...){

  if(!is.matrix(B)){
    stop("'B' must be a matrix.")
  }
  if(is.null(rownames(B))){
    stop("'B' must have row names corresponding to parts.")
  }
  if(!is.null(summary_fun) && !is.function(summary_fun)){
    stop("'summary_fun' must be NULL or a function.")
  }

  hc = hclust_dendrogram(B)
  coords = hclust_segment_coordinates(hc)

  node_labels = build_balance_labels(
    B = B,
    data = data,
    summary_fun = summary_fun
  )

  plot(
    hc,
    main = main,
    sub = "",
    xlab = "",
    ylab = "",
    axes = FALSE,
    ...
  )

  yoff = offset_node * max(c(1, hc$height))

  text(
    x = coords$x_mid,
    y = pmax(coords$y - yoff, 0),
    labels = node_labels,
    cex = cex_node,
    pos = 1
  )

  invisible(
    data.frame(
      node = seq_along(node_labels),
      x_left = coords$x_left,
      x_right = coords$x_right,
      x_mid = coords$x_mid,
      y = coords$y,
      label = node_labels,
      stringsAsFactors = FALSE
    )
  )
}


# Build labels for internal nodes
build_balance_labels = function(B, data = NULL, summary_fun = NULL){

  bal_names = colnames(B)
  if(is.null(bal_names)){
    bal_names = paste0("b", seq_len(ncol(B)))
  }

  if(is.null(data) || is.null(summary_fun)){
    return(bal_names)
  }

  X = as.data.frame(data)

  if(!all(rownames(B) %in% colnames(X))){
    stop("All row names of 'B' must be present as column names in 'data'.")
  }

  X = X[, rownames(B), drop = FALSE]
  X = as.matrix(X)

  if(!is.numeric(X)){
    stop("'data' must contain numeric values.")
  }
  if(any(!is.finite(X))){
    stop("'data' contains non-finite values.")
  }
  if(any(X <= 0)){
    stop("Balances require strictly positive data.")
  }

  Z = log(X) %*% B

  stat_txt = vapply(seq_len(ncol(B)), function(j){
    txt = summary_fun(Z[, j])
    if(length(txt) != 1 || !is.character(txt)){
      stop("'summary_fun' must return a character string of length 1.")
    }
    txt
  }, character(1))

  paste0(bal_names, "\n", stat_txt)
}


# Coordinates tied to each internal horizontal segment
hclust_segment_coordinates = function(hc){

  n = length(hc$order)
  merge = hc$merge
  height = hc$height

  leaf_x = numeric(n)
  leaf_x[hc$order] = seq_len(n)

  node_center = numeric(nrow(merge))

  x_left  = numeric(nrow(merge))
  x_right = numeric(nrow(merge))
  x_mid   = numeric(nrow(merge))
  y       = height

  get_center = function(idx){
    if(idx < 0){
      leaf_x[-idx]
    } else {
      node_center[idx]
    }
  }

  for(i in seq_len(nrow(merge))){
    left  = merge[i, 1]
    right = merge[i, 2]

    xl = get_center(left)
    xr = get_center(right)

    x_left[i] = xl
    x_right[i] = xr
    x_mid[i] = (xl + xr) / 2
    node_center[i] = x_mid[i]
  }

  list(
    x_left = x_left,
    x_right = x_right,
    x_mid = x_mid,
    y = y
  )
}


hclust_dendrogram = function(B){

  if(!is.matrix(B)){
    stop("'B' must be a matrix.")
  }
  if(is.null(rownames(B))){
    stop("'B' must have row names.")
  }

  MERGE = matrix(0, nrow = ncol(B), ncol = 2)
  ORD = order(colSums(B != 0))

  for(i in seq_len(ncol(B))){

    neg_idx = which(B[, ORD[i]] < 0)
    pos_idx = which(B[, ORD[i]] > 0)

    if(length(neg_idx) == 1){
      MERGE[i, 1] = -neg_idx
    } else {
      cand = which(sapply(seq_len(i - 1), function(j){
        all(B[neg_idx, ORD[j]] * B[neg_idx, ORD[i]] != 0)
      }))
      if(length(cand) == 0){
        stop("Could not reconstruct left merge structure from 'B'.")
      }
      MERGE[i, 1] = utils::tail(cand, 1)
    }

    if(length(pos_idx) == 1){
      MERGE[i, 2] = -pos_idx
    } else {
      cand = which(sapply(seq_len(i - 1), function(j){
        all(B[pos_idx, ORD[j]] * B[pos_idx, ORD[i]] != 0)
      }))
      if(length(cand) == 0){
        stop("Could not reconstruct right merge structure from 'B'.")
      }
      MERGE[i, 2] = utils::tail(cand, 1)
    }
  }

  left = function(pair){
    B_ = sign(B)[pair, , drop = FALSE]
    active = apply(B_, 2, function(x) all(x != 0))
    which.min(rowSums(B_[, active, drop = FALSE]))
  }

  ORDER = seq_len(nrow(B))
  for(i in seq_len(nrow(B))){
    if(i < nrow(B)){
      for(j in (i + 1):nrow(B)){
        x = c(ORDER[i], ORDER[j])
        ileft_ = left(x)
        ORDER[i] = x[ileft_]
        ORDER[j] = x[3 - ileft_]
      }
    }
  }

  HEIGHT = rep(0, ncol(B))
  for(i in seq_len(nrow(MERGE))){
    l_ = 1
    r_ = 1
    if(MERGE[i, 1] > 0){
      l_ = HEIGHT[MERGE[i, 1]] + 1
    }
    if(MERGE[i, 2] > 0){
      r_ = HEIGHT[MERGE[i, 2]] + 1
    }
    HEIGHT[i] = max(l_, r_)
  }

  structure(
    list(
      merge  = MERGE,
      height = HEIGHT,
      order  = ORDER,
      labels = rownames(B)
    ),
    class = "hclust"
  )
}
