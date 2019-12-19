library(coda.base)
load('testing/testingX.RData')

(A <- .Call('_coda_base_find_PB', PACKAGE = 'coda.base', X))
apply(coordinates(X,A), 2, var)
(B <- pb_basis(X, 'lsearch', rep=10))
apply(coordinates(X, B), 2, var)

# library(microbenchmark)
# microbenchmark(
#   .Call('_coda_base_find_PB', PACKAGE = 'coda.base', X),
#   pb_basis(X, 'lsearch'))
