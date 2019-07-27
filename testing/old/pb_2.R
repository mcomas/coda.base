
#setwd("c:/f/compositions")
waste=read.table("testing/waste.csv", header=TRUE, sep=";", dec=",")
names(waste)

Y = as.matrix(waste[,6:10])

m0 = lm(coordinates(Y)~waste$floating_population)

m0.res = composition(m0$residuals)
m0.fitted = composition(m0$fitted.values)

B0fitted= pb_basis(m0.fitted, method = 'exact')

# Seguint amb l'exemple del waste.
# Crec que el que es veu a
B0fitted
# és una reinterpretació dels coeficients de la regressió
B0fitted[,1]
# dient que, en mitjana, les dues primeres parts es
# perturben en direcció contraria a les altres parts.

# De fet, si enlloc de principal balances agafem els principal
# components. Tenim que és exactament el mateix:
v1 = coordinates(composition(coef(m0)[2,]), 'clr')
v1/sqrt(sum(v1^2))   # Faig el vector unitari per comparar
# De fet,
v2 = pc_basis(m0.fitted)[,1]
v2 / sqrt(sum(v2^2)) # Faig el vector unitari per comparar


# Sobre lo que necessitem més observacions que parts, no ho veig clar.
# Els principal balances no tenen restricció en quan a relació
# d'observacions i parts
X = matrix(exp(rnorm(5*10)), ncol = 10, nrow=5)
pb_basis(X, method='exact')[,1]
pb_basis(X, method='ward.D2')[,1]
