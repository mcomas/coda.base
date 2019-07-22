library(coda.base)

Y = as.matrix(iris[,1:4])

m0 = lm(coordinates(Y)~1)
m0.res = composition(m0$residuals)

# Partint de la idea que els principal balances aplicats
# a un model lineal constant és equivalent a fer-ho sobre
# la composició.
B0 = pb_basis(m0.res, method = 'exact')
B0.alt = pb_basis(Y, method='exact')
colnames(B0) = colnames(B0.alt) = c('PB1', 'PB2', 'PB3')
rownames(B0) = rownames(B0.alt) = colnames(Y)
B0
B0.alt

# Podem extendra-ho utilizant una covariable (o més).
x = iris[,5]
m1 = lm(coordinates(Y)~x)
m1.res = composition(m1$residuals)

B1 = pb_basis(m1.res, method = 'exact')
colnames(B1) = c('PB1', 'PB2', 'PB3')
rownames(B1) = colnames(Y)
B1
