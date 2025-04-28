# Generate orthonormal coordinates (from different ilr systems) representing information equivalent to alr coordinates
# (called backward pivot coordinates (bpc) in paper Hron et al., 2021)

# It is an adaptation of the idea of pivot coordinates

# Example data
dat <- matrix(c(404.0555,738.4495,236.9774,57.90058,2.616169,
         522.9290,573.0399,253.1222,87.79391,3.115526), ncol = 5, byrow = TRUE,
         dimnames=list(NULL,c("SLEEP","SB","LPA","MPA","VPA")))

bpc <- function(x){

  gm <- function (x){
    if (!is.numeric(x))
      stop("x has to be a vector of class numeric")
    if (any(na.omit(x == 0)))
      0
    else exp(mean(log(unclass(x)[is.finite(x) & x > 0])))
  }

  x.bpc <- matrix(NA, nrow = nrow(x), ncol = ncol(x) - 1)
  D <- ncol(x)
  for (i in 1:ncol(x.bpc)) {
    x.bpc[, i] <- sqrt(i/(i+1))*log(apply(as.matrix(x[,1:i]), 1, gm)/(x[, i+1]))
  }
  return(x.bpc)
}

namx <- colnames(dat)
D <- ncol(dat)

# bpc with the ratio x[1]/x[D] represented by first coordinate of the ilr system (first column)
dat.bpc <- bpc(cbind(dat[,1],dat[,D],dat[,namx[-c(1,D)]]))
dat.bpc[,1]

# bpc with the ratio x[2]/x[D] represented by first coordinate of the ilr system (first column)
dat.bpc <- bpc(cbind(dat[,2],dat[,D],dat[,namx[-c(2,D)]]))
dat.bpc[,1]

# And so on for x[3]/x[D] ... x[D-1]/x[D]

coordinates(cbind(dat[,1],dat[,D],dat[,namx[-c(1,D)]]), sbp_basis(alr_basis(D)))
