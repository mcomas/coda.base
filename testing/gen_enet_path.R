Rcpp::sourceCpp("src/admm_gen_enet_path.cpp")

set.seed(1)
n <- 100; p <- 20
X <- matrix(rnorm(n*p), n, p)
beta_true <- c(rep(2, 3), rep(0, p-3))
y <- X %*% beta_true + rnorm(n)

# Exemple D: fused lasso 1D sobre coeficients (penalitza diferències consecutives)
D <- diag(p) - rbind(0, diag(p-1))
D <- D[-1, , drop = FALSE]  # (p-1) x p

# Path de lambdes (de gran a petit)
lmax <- 50
lambdas <- exp(seq(log(lmax), log(1e-2), length.out = 50))

fit <- admm_gen_enet_path(y = as.vector(y), X = X, D = D,
                          lambdas = lambdas, gamma = 1.0,
                          rho = 1.0, maxit = 5000)

dim(fit$beta)  # p x length(lambdas)


#####################
Rcpp::sourceCpp("src/fista_dual_gen_enet_path.cpp")

set.seed(1)
n <- 120; p <- 30
X <- matrix(rnorm(n*p), n, p)
beta0 <- c(rep(1.5, 4), rep(0, p-4))
y <- as.vector(X %*% beta0 + rnorm(n))

# D: fused lasso 1D sobre coeficients (p-1 x p)
D <- diag(p)[-1, ] - diag(p)[-p, ]
D <- D[-1, , drop = FALSE]

lambdas <- exp(seq(log(50), log(1e-2), length.out = 60))

fit_pg <- fista_dual_gen_enet_path(y, X, D, lambdas, gamma = 1.0,
                                   maxit = 10000, tol = 1e-15,
                                   powerit = 80, verbose = FALSE)

fit <- admm_gen_enet_path(y = as.vector(y), X = X, D = D, reltol = 1e-15,
                          lambdas = lambdas, gamma = 1.0,
                          rho = 1.0, maxit = 5000)

dim(fit_pg$beta)   # p x length(lambdas)
