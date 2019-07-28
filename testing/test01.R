suppressMessages(library(coda.base))
library(magrittr)

cat("\n\n\n# Two different approaches to find principal balances.\n\n")
D = 5
n = 100
#set.seed(1)
Y = matrix(exp(rnorm(D*n)), ncol = D)
ev = function(B) apply(coordinates(Y, B), 2, var)
find_principal_balance_01(Y) %>% ev()
find_principal_balance2_01(Y) %>% ev()
find_principal_balance3_01(Y) %>% ev()

B2 = find_principal_balance2_01(Y)[,1:3,drop=F]
t(B2) %*% cov(coordinates(Y, 'clr')) %*% B2

B3 = find_principal_balance3_01(Y)[,1:3,drop=F]
t(B3) %*% cov(coordinates(Y, 'clr')) %*% B3

cor(coordinates(Y, B2))
cor(coordinates(Y, B3))

det(cov(coordinates(Y, find_principal_balance2_01(Y))[,1:2]))
det(cov(coordinates(Y, find_principal_balance3_01(Y))[,1:2]))

cor(coordinates(Y, find_principal_balance2_01(Y)))[1:2,1:2]
cor(coordinates(Y, find_principal_balance3_01(Y)))[1:2,1:2]

cat("\n\n\n# Canonical balance example using waste dataset.\n\n")
waste=read.table("testing/old/waste.csv", header=TRUE, sep=";", dec=",")
Y = as.matrix(waste[,6:10])
X = cbind(waste$floating_population)
sign(find_canonical_balance_01(Y, X))

sign(find_principal_balance2_01(Y))
sign(find_principal_balance3_01(Y))
