library(coda.base)
library(magrittr)
D = 10
n = 100
#set.seed(1)
Y = matrix(exp(rnorm(D*n)), ncol = D)
ev = function(B) apply(coordinates(Y, B), 2, var)
find_principal_balance_01(Y) %>% ev()
find_principal_balance2_01(Y) %>% ev()
X = matrix(rnorm((D-1)*n), ncol = D-1)

waste=read.table("testing/old/waste.csv", header=TRUE, sep=";", dec=",")

Y = as.matrix(waste[,6:10])
X = cbind(waste$floating_population)
sign(find_canonical_balance_01(Y, X))
