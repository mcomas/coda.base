library(tidyverse)
library(coda.base)
get_exact_diff = function(D){
  .CONT = TRUE
  while(.CONT){
    X = matrix(rlnorm(10*D), ncol = D)
    B0 = pb_basis(X, method = 'constrained')
    B0_opt = pb_basis(X, method = 'exact')
    v0 = apply(coordinates(X, B0), 2, var)
    v0_opt = apply(coordinates(X, B0_opt), 2, var)
    if(abs(v0[1]-v0_opt[1]) > 0.001) .CONT = FALSE
  }
  X
}
get_ts_diff = function(D){
  .CONT = TRUE
  while(.CONT){
    X = matrix(rlnorm(10*D), ncol = D)
    B0 = pb_basis(X, method = 'constrained')
    TS = partial_pb_tabu_search(X, lapply(1:ncol(X), identity), iter = 100, tabu_size = 100)
    v0 = apply(coordinates(X, B0), 2, var)
    if(abs(v0[1]-TS$variance) > 0.001) .CONT = FALSE
  }
  X
}


dplot = map_df(1:5000, function(i){
  set.seed(i)
  X = get_ts_diff(10)
  TS = partial_pb_tabu_search(X, lapply(1:ncol(X), identity), iter = 500, tabu_size = 500)
  tibble(
    seed = i,
    var_constrained = var(coordinates(X, pb_basis(X, method = 'constrained')[,1,drop=F])),
    var_tabu_search = TS$variance,
    iter_best = TS$iter_best
  )
})
ggplot(data = dplot) +
  geom_bar(aes(x = iter_best))
slice_max(dplot, iter_best)

set.seed(56)
X = get_ts_diff(10)
partial_pb_tabu_search(X, lapply(1:ncol(X), identity),
                       iter = 100, tabu_size = 100)$variance
partial_pb_tabu_search(X, lapply(1:ncol(X), identity),
                       iter = 100, tabu_size = 100)$variance


TXT = capture.output(
  TS <- partial_pb_tabu_search(X, lapply(1:ncol(X), identity),
                               iter = 100, tabu_size = 100, debug = TRUE))
TXT
vars <- TXT |>
  grep(pattern = "current variance:", value = TRUE) |>
  sub(".*current variance: ", "", x = _) |>
  as.numeric()

plot(vars, type = 'l')

set.seed(4182)
X = get_ts_diff(10)


TXT = capture.output(
  TS <- partial_pb_tabu_search(X, lapply(1:ncol(X), identity),
                               iter = 500, tabu_size = 500, debug = TRUE))

vars <- TXT |>
  grep(pattern = "current variance:", value = TRUE) |>
  sub(".*current variance: ", "", x = _) |>
  as.numeric()

plot(vars, type = 'l')
