library(coda.base)

X = matrix(c(1,2,3,
             4,5,6,
             7,8,9), byrow = TRUE, ncol = 3)
colnames(X) = c('c1', 'c2', 'c3')

## Definir coordenades a partir d'una matriu de logcontrastos qualsevol
B = matrix(c(-1, 0,-1,
              1, 1,-1,
              0,-1, 2), byrow = TRUE, ncol = 3)
colSums(B)
coordinates(X, basis = B)

## Definir alr coordinates
coordinates(X, basis = 'alr')
B = alr_basis(3)
coordinates(X, basis = B, label = 'alr')   # Equivalent a anterior
B = alr_basis(3, denominator = 1)
coordinates(X, basis = B, label = 'alr')
B = alr_basis(3, denominator = 1, numerator = 3:2)
coordinates(X, basis = B, label = 'alr')

## Definir clr coordinates
coordinates(X, basis = 'clr')
B = clr_basis(3)
coordinates(X, basis = B, label = 'clr')   # Equivalent a anterior

## Definir ilr coordinates
coordinates(X, basis = 'ilr')
B = ilr_basis(3)
coordinates(X, basis = B, label = 'ilr')   # Equivalent a anterior

coordinates(X, basis = 'pb')
B = pb_basis(X, method = 'exact')
coordinates(X, basis = B, label = 'pb')

coordinates(X, basis = 'pc')
B = pc_basis(X)
coordinates(X, basis = B, label = 'pc')

## A partir de balanços
coordinates(X, list(b1 = c1+c2~c3,
                    b2 = c1~c2))
B = sbp_basis(b1 = c1+c2~c3,
              b2 = c1~c2, data = X)
coordinates(X, basis = B)
B = sbp_basis(matrix(c(-1,-1, 1,
                       -1, 1, 0), byrow = TRUE, ncol = 3))
coordinates(X, basis = B)

B = sbp_basis(matrix(c(1, 1,-1,
                       1,-1, 0), byrow = TRUE, ncol = 3))
coordinates(c(10.10, 1.40, 0.50), basis = B)


