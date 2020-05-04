# install.packages('devtools')
library(devtools)

# install_github('mcomas/coda.count')
library(coda.count)

download.file('https://github.com/cran/zCompositions/raw/master/data/Pigs.rda', 'Pigs.rda')
load('Pigs.rda')

XZ = as.matrix(Pigs)

# Dirichlet-multinomial fitting

pars.dm = fit_dm(XZ)
ALPHA = pars.dm[,1]
print( sprintf("alpha = (%s)", paste(sprintf("%0.4f", ALPHA), collapse= ', ')))

clo = function(X) X/rowSums(X)
# Posteriori replacement
X.dm <- clo( t(t(XZ) + ALPHA) )
round(X.dm, 4)

H.dm = ilr_coordinates(X.dm)
MU = colMeans(H.dm)
SIGMA = cov(H.dm)


# Logratio-normal-multinomial fitting (using dirichlet-multinomial starting point)

## Normal random vector
set.seed(2)
N = 10000
d = ncol(Pigs) - 1
Z = matrix(rnorm(N * d), ncol = d)

fit.lnm = c_lrnm_fit_mc(X = XZ, mu0 = MU, sigma0 = SIGMA, Z = Z, tol = 0.01, em_max_steps = 100)
X.lnm = fit.lnm[[3]]
H.lnm = ilr_coordinates(X.lnm)

# Using pseudo-random sequences
library(randtoolbox)
Z.pseudo = sobol(N, dim = ncol(XZ) - 1, normal = TRUE, init = TRUE, scrambling = 1, seed = 2)

fit.lnm.pseudo = c_lrnm_fit_mc(X = XZ, mu0 = MU, sigma0 = SIGMA, Z = Z.pseudo, tol = 0.01, em_max_steps = 100)
X.lnm.pseudo = fit.lnm.pseudo[[3]]
H.lnm.pseudo = ilr_coordinates(X.lnm.pseudo)
