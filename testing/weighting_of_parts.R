library(coda.base)
library(coda.count)
set.seed(1)
n = 100
X = rlrnormal(n, c(0,0,0), diag(3))
b = rnorm(3)
y = coordinates(X) %*% b + rnorm(n, 0, 0.1)

X[,2] = X[,2]^2
dat = data.frame(as.matrix(coordinates(X)))
dat$y = y |> as.vector()

summary(m <- lm(y~., dat))

x = c(1,2,3)
p = c(3,2,1)
coordinates(x, 'clr')
log(x) - weighted.mean(log(x),p)

x
