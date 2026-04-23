set.seed(123)
set_zeros = function(X, p = 0.05){
  tX = t(X)
  DL = apply(tX, 1, function(x) quantile(x[!is.na(x) & x > 0], p))
  tX[tX < DL] = 0
  t(tX)
}
library(coda.count)
alpha = fit_dm(X)

X = alimentation[,1:9] |>
  as.matrix() |>
  set_zeros(p = 0.1)

sum(X == 0)

Xr = coda_replacement(X)

tX = t(X)
DL = apply(tX, 1, function(x) min(x[!is.na(x) & x > 0]))


X = t(tX)
sum(X == 0)

X[sample(length(X), 10)] <- 0  # Introduce some zeros
# X[sample(length(X), 10)] <- NA  # Introduce some NAs
tX = t(X)

