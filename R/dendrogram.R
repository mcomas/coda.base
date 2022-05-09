hclust_dendrogram = function(B){
  MERGE = matrix(0, nrow = ncol(B), ncol = 2)
  ORD = order(colSums(B != 0))
  for(i in 1:ncol(B)){
    if(sum(B[,ORD[i]] < 0) == 1){
      MERGE[i,1] = -which(B[,ORD[i]] < 0)
    }else{
      MERGE[i,1] = which(sapply(1:(i-1), function(j) all(B[which(B[,ORD[i]] < 0), ORD[j]] * B[which(B[,ORD[i]] < 0), ORD[i]] != 0)))
    }
    if(sum(B[,ORD[i]] > 0) == 1){
      MERGE[i,2] = -which(B[,ORD[i]] > 0)
    }else{
      MERGE[i,2] = which(sapply(1:(i-1), function(j) all(B[which(B[,ORD[i]] > 0), ORD[j]] * B[which(B[,ORD[i]] > 0), ORD[i]] != 0)))
    }
  }
  left = function(pair){
    B_ = sign(B)[pair,]
    which.min(rowSums(B_[,apply(B_, 2, function(x) all(x != 0)), drop=FALSE]))
  }
  ORDER = 1:nrow(B)
  for(i in 1:nrow(B)){
    if(i != nrow(B)){
      for(j in (i+1):nrow(B)){
        x = c(ORDER[i],ORDER[j])
        ileft_ = left(x)
        ORDER[i] = x[ileft_]
        ORDER[j] = x[3-ileft_]
      }
    }
  }
  HEIGHT = rep(0, ncol(B))
  for(i in 1:nrow(MERGE)){
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
    list(merge = MERGE, height = HEIGHT, order = ORDER, labels = rownames(B)), class = 'hclust'
  )
}

#' Plot a balance
#'
#' @param B Balance to plot
#' @param data (Optional) Data used to calculate the statistics associated to a balance
#' @param main Plot title
#' @param ... further arguments passed to plot
#' @return Balance plot
#'
#' @export
plot_balance = function(B, data = NULL, main = 'Balance dendrogram', ...){
  if(is.null(data)){
    B = as.matrix(B)
    plot(as.dendrogram(hclust_dendrogram(B)), main = main, ylab = "", axes = FALSE, ...)
  }else{

  }

}


#######
