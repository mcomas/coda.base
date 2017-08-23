set.seed(1)
library(Rcpp)
library(coda.base)
sourceCpp('src/principal_balances.cpp')
K = 5
X=as.data.frame(matrix(exp(rnorm(K*K*2)), nrow=2*K, ncol=K))
M = cov(log(X))
Xclr = coordinates(X, 'clr')
pc1 = eigen(cov(Xclr))$vectors[,1]

pL = exp(-pc1)/sum(exp(-pc1))
pR = exp(pc1)/sum(exp(pc1))


l_p = pL / ( pL + pR )
r_p = pR / ( pL + pR )

# Firstly, random location
# Secondly, removing elements optimizing variance
# adding elements either to left or right
l_p
r_p


pc1

pL[-(1:2)]
pR[-(1:2)]
