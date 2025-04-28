SWEEP = function(A, K){
  for(k in K){
    D = A[k,k]
    A[k,] = A[k,] / D
    for(i in which(1:nrow(A)!=k)){
      B = A[i,k]
      A[i,] = A[i,] - B * A[k,]
      A[i,k] = -B/D
    }
    A[k,k] = 1/D
  }
  A
}

Y = as.matrix(iris[,4,drop=FALSE])
X = as.matrix(iris[,1:3])

A = rbind(
  cbind(t(X) %*% X, t(X) %*% Y),
  cbind(t(Y) %*% X, t(Y) %*% Y))

SWEEP(SWEEP(A, 1), 2)

S = cov(iris[1:4])
SWEEP(S, 1)

library(mvtnorm)  # Per a mostreig multivariant (opcional)

# Definim una matriu de covariància per a 3 variables
Sigma <- matrix(c(
  4,  2,  1,
  2,  3,  1,
  1,  1,  2
), nrow = 3, byrow = TRUE)
colnames(Sigma) <- rownames(Sigma) <- c("Y1", "Y2", "Y3")
Sigma

# Augmentem amb mitjanes fictícies (totes zero)
mu <- c(0, 0, 0)
aug <- cbind(Sigma, mu)
aug <- rbind(aug, c(mu, 1))
colnames(aug) <- rownames(aug) <- c("Y1", "Y2", "Y3", "Mean")

sweep_matrix <- Sigma
sweep_var <- c(2, 3)  # Y2 i Y3
# Sweep manual (amb la funció sweep definida a mà)
sweep_op <- function(S, k) {
  for (j in k) {
    d <- S[j, j]
    for (i in 1:nrow(S)) {
      for (l in 1:ncol(S)) {
        if (i != j && l != j) {
          S[i, l] <- S[i, l] - S[i, j] * S[j, l] / d
        }
      }
    }
    for (i in 1:nrow(S)) {
      if (i != j) S[i, j] <- S[i, j] / d
    }
    for (l in 1:ncol(S)) {
      if (l != j) S[j, l] <- S[j, l] / d
    }
    S[j, j] <- -1 / d
  }
  return(S)
}

sweep_op(aug, c(2,3))
SWEEP(aug, c(2, 3))


