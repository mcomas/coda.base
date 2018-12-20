public static int[][] defaultPartition(int n_components){
  int partition[][] = new int[n_components-1][n_components];
  fillPartition(partition, 0, 0, n_components);
  return partition;
}

defaultPartition = function(n_comp){

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

fillPartition(matrix(0, nrow = 1, ncol = 2), 0, 1, 2)
fillPartition(matrix(0, nrow = 1, ncol = 3), 0, 1, 3)
fillPartition(matrix(0, nrow = 1, ncol = 4), 0, 1, 4)
fillPartition(matrix(0, nrow = 1, ncol = 5), 0, 1, 5)

partition = fillPartition(partition, 1, 1, n_comp)
