suppressMessages(library(coda.base))
library(magrittr)

cat("\n\n\n# Approaches to find principal balances.\n\n")
D = 5
n = 100

set.seed(80)
Y = matrix(exp(rnorm(D*n)), ncol = D)
ev = function(B) apply(coordinates(Y, B), 2, var)
cat("## Two equivalent approaches. First faster, second clearer.\n")
(PB1 <- find_principal_balance_01(Y)) %>% ev()
(PB2 <- find_principal_balance2_01(Y)) %>% ev()
cat("## Balances are build optimising the total variance.\n")
(PB3 <- find_principal_balance3_01(Y)) %>% ev()


det(t(PB2[,1:3, drop=F]) %*% cov(coordinates(Y, 'clr')) %*% PB2[,1:3, drop=F])
det(t(PB3[,1:3, drop=F]) %*% cov(coordinates(Y, 'clr')) %*% PB3[,1:3, drop=F])


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
