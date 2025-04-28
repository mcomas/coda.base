library(coda.base)

X = as.matrix(iris[,1:4])
set.seed(6)
i = rbinom(length(X), 1, 0.3) == 1
# i = i & (1:length(X) > nrow(X))
X[i] = NA
X = X[apply(is.na(X), 1, sum) < 3,]
X1 = X
X2 = X/rowSums(X, na.rm = TRUE)

library(mice)
set.seed(1)
X1.imp = mice(X1, method='norm', printFlag = FALSE)
set.seed(1)
X2.imp = mice(X2, method='norm', printFlag = FALSE)
summary(complete(X1.imp, 1))
summary(complete(X2.imp, 1))


set.seed(1)
B = pairwise_basis(ncol(X))
H = coordinates(X, B)
n.na = colSums(is.na(H))
n.na
sort(n.na)
ind = c(6, 3, 5)
mi = mice(H[,ind])
X.imp = composition(complete(mi, 1), B[,ind])
summary(X.imp)
