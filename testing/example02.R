set.seed(1)
library(Rcpp)
library(coda.base)
sourceCpp('src/principal_balances.cpp')
K=10
X=as.data.frame(matrix(exp(rnorm(2*K^2)), nrow=2*K, ncol=K))
lX = log(X)
M = cov(lX)

## This example has the lowest variance precisly in the balance separating all parts.
##

##
library(compositions)

hide = function(f) suppressMessages(invisible(capture.output(f)))

hide(b0 <- find_PB_rnd_local_search(M, 10))
max(apply(coordinates(X, b0), 2, var))
hide(b1 <- gsi.PrinBal(acomp(X), method = 'PBhclust'))
max(apply(coordinates(X, b1), 2, var))
hide(b2 <- gsi.PrinBal(acomp(X), method = 'PBmaxvar'))
max(apply(coordinates(X, b2), 2, var))
hide(b3 <- gsi.PrinBal(as.data.frame(X), method = 'PBangprox'))
max(apply(coordinates(X, b3), 2, var))


library(microbenchmark)
microbenchmark(
  find_PB_rnd_local_search(M, rep=10),
  gsi.PrinBal(acomp(X), method = 'PBhclust'),
  gsi.PrinBal(acomp(X), method = 'PBmaxvar'),
  gsi.PrinBal(as.data.frame(X), method = 'PBangprox'))
# Unit: microseconds
#                                                expr        min         lq        mean
#                                find_PB(M, rep = 10)    483.733    559.595    740.1333
#                               gsi.PrinBal(acomp(X))   1530.630   1924.890   2395.1544
#          gsi.PrinBal(acomp(X), method = "PBmaxvar") 784356.609 833910.183 876705.3402
# gsi.PrinBal(as.data.frame(X), method = "PBangprox")  51218.971  64074.697  96679.9755
#
#     median          uq         max neval
#    619.838    937.5665    1546.248   100
#   2115.884   2771.8675    6658.016   100
# 866365.108 897012.5055 1309996.216   100
#  69534.091  83199.3105  468865.299   100


