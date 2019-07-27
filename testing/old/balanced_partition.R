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

cdp_partition = function(ncomp) fillPartition(matrix(0, nrow = 1, ncol = ncomp), 0, 1, ncomp)

cdp_partition(2)
cdp_partition(3)
cdp_partition(4)
