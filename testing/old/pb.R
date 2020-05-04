library(coda.base)

#agafo aquest arxiu perqu? la composici? ?s explicativa i hi ha una variable dependent
#poblacio_flotant del municipi com % del cens. Les dades son les del lab de regressi? amb composici? explicativa.
#x1 fins x5 s?n les fraccions de recollida de brossa:
#x1 non recyclable (grey waste container in Catalonia)
#x2 glass (bottles and jars of any colour: green waste container)
#x3 light containers (plastic packaging, cans and tetra packs: yellow container)
#x4 paper and cardboard (blue container)
#x5 biodegradable waste (brown container)


#setwd("c:/f/compositions")
waste=read.table("testing/waste.csv", header=TRUE, sep=";", dec=",")
names(waste)

Y = as.matrix(waste[,6:10])

m0 = lm(coordinates(Y)~waste$floating_population)

m0.res = composition(m0$residuals)
m0.fitted = composition(m0$fitted.values)



# Partint de la idea que els principal balances aplicats

# a un model lineal constant ?s equivalent a fer-ho sobre

# la composici?.

B0 = pb_basis(m0.res, method = 'exact')
B0fitted= pb_basis(m0.fitted, method = 'exact')
B0.alt = pb_basis(Y, method='exact')


colnames(B0) = colnames(B0.alt) = colnames(B0fitted) = c('PB1', 'PB2', 'PB3', 'PB4')

rownames(B0) = rownames(B0.alt) = rownames(B0fitted) = colnames(Y)

B0
B0fitted
B0.alt

#he afegit al que tenies el mateix sobre els valors previstos en lloc de sobre els residus



# Podem extendra-ho utilizant una covariable (o m?s).

#x = iris[,5]

#m1 = lm(coordinates(Y)~x)

#m1.res = composition(m1$residuals)




#B1 = pb_basis(m1.res, method = 'exact')

#colnames(B1) = c('PB1', 'PB2', 'PB3')

#rownames(B1) = colnames(Y)

#B1

#ara es tracta d'explicar la variable amb els PB creats, per exemple els 1 o 2 primers, en ordre
#de vari?ncia i no necess?riament en ordre d'extracci?

summary(lm(waste$floating_population~coordinates(Y, basis = B0)[,1]))
summary(lm(waste$floating_population~coordinates(Y, basis = B0fitted)[,1]))

b = sbp_basis(find_predictive_balance(Y, waste$floating_population), silent = TRUE)
summary(lm(waste$floating_population~coordinates(Y, basis = b)[,1]))


optimal_balance=sqrt(2/3)*log(waste$x2/(sqrt(waste$x3*waste$x4)))
summary(lm(waste$floating_population~optimal_balance))

m0 = lm(coordinates(Y)~waste$floating_population)
H = coordinates(Y, 'ilr')
cc = cancor(H, cbind(waste$floating_population))
cc_H = (H %*% cc$xcoef)[,1]
summary(lm(waste$floating_population~cc_H))

CC_CLR = ilr_basis(5) %*% cc$xcoef
summary(lm(waste$floating_population~coordinates(Y, CC_CLR[,1, drop=FALSE])))



clr1 = matrix(coordinates(composition(coef(m0)[2,]), 'clr'), ncol=1)
summary(lm(waste$floating_population~coordinates(Y, basis = clr1)[,1]))


summary(lm(waste$floating_population~coordinates(Y, basis = B0.alt)[,1]))

summary(lm(waste$floating_population~coordinates(Y, basis = B0)[,1:2]))
summary(lm(waste$floating_population~coordinates(Y, basis = B0fitted)[,1:2]))
summary(lm(waste$floating_population~coordinates(Y, basis = B0.alt)[,1:2]))

#em surt la vari?ncia explicada m?s gran quan es fa amb els fitted que no pas quan es fa amb els residus.
#en tot cas la gr?cia del selbal ?s que tamb? funciona si hi ha m?s parts que individus. Aix? que fem requerei m?s parts que individus


