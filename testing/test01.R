suppressMessages(library(coda.base))
library(magrittr)

cat("\n\n\n# Two different approaches to find principal balances.\n\n")
D = 9
n = 100
#set.seed(1)
Y = matrix(exp(rnorm(D*n)), ncol = D)
ev = function(B) apply(coordinates(Y, B), 2, var)
find_principal_balance_01(Y) %>% ev()
find_principal_balance2_01(Y) %>% ev()
X = matrix(rnorm((D-1)*n), ncol = D-1)

cat("\n\n\n# Canonical balance example using waste dataset.\n\n")
waste=read.table("testing/old/waste.csv", header=TRUE, sep=";", dec=",")
Y = as.matrix(waste[,6:10])
X = cbind(waste$floating_population)
sign(find_canonical_balance_01(Y, X))
