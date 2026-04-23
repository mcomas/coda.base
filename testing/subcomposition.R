library(coda.base)

X = as.matrix(alimentation[,1:9])

I = combn(9, 3, simplify = FALSE)
i.max = I[[which.max(sapply(I, \(i) sum(diag(cov(coordinates(X[,i], 'clr'))))))]]
colnames(X)[i.max]

sum(diag(cov(coordinates(X[,1:3], 'clr'))))
H = coordinates(X, 'clr')
colnames(H) = colnames(X)
SVD = svd(base::scale(H, scale = FALSE))
order(diag(cov()), decreasing = TRUE)[1:3]
sign(pb_basis(X, method = 'exact')[,1])
library(ggbiplot)
ggbiplot(prcomp(H))

RAYS = SVD$v %*% diag(SVD$d)
convhulln(RAYS[,1:4])

