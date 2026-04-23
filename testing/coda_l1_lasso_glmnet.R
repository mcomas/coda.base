library(glmnet)
library(genlasso)
library(coda.base)

X = as.matrix(waste[,5:9])[1:15,]
y = waste$floating_population[1:15]

LX = log(X)
K = ncol(X)
n = nrow(X)

coda_l1_genlasso = function(LX, y, lambda){
  B = ilr_basis(LX)
  H = LX %*% B
  X1 = cbind(1, H, 0)
  D = cbind(0, B, 1)
  m1 = genlasso(y, X1, D, svd = TRUE)
  b1 = coef(m1, lambda = 0.08)$beta
  clr_beta = as.vector(B %*% b1[2:5])
  m = b1[6]

  glmnet:::glmnet.fit(
    x = CLR, y = y,
    weights = weights,
    offset = off,
    lambda = LAMBDA,
    intercept = TRUE
  )
}
beta_coda_l1_genlasso = function(LX, y, ...){
  B = ilr_basis(LX)
  H = LX %*% B
  X1 = cbind(1, H, 0)
  D = cbind(0, B, 1)
  m1 = genlasso(y, X1, D, svd = TRUE)
  do.call(rbind, lapply(seq_along(m1$lambda),
                        function(i){
                          data.frame(
                            lambda = m1$lambda[i],
                            var = colnames(LX),
                            coef = as.vector(m1$beta[1+1:ncol(H),i]  %*% t(B))
                          )
                        }))
}
beta_coda_l1_glmnet_m = function(LX, y, lambda = NULL, weights = NULL,
                                 m_interval = NULL, ...) {
  # LX: log(X) (n x K) o directament X (n x K) si tu ja li fas log abans
  # y : vector (n)
  # lambda: vector de lambdes (opcional). Si NULL, glmnet genera el camí.
  # weights: pesos (opcional)
  # m_interval: interval per optimise (opcional)

  stopifnot(is.matrix(LX) || is.data.frame(LX))
  LX <- as.matrix(LX)
  y  <- as.numeric(y)

  n <- nrow(LX)
  K <- ncol(LX)

  if (is.null(weights)) weights <- rep(1, n)

  CLR <- LX - rowMeans(LX)

  # Camí de lambdes (sobre H) per tenir una graella raonable
  fit_path <- glmnet(CLR, y, intercept = TRUE, standardize = FALSE,
                     weights = weights, lambda = lambda, ...)

  lambdas <- fit_path$lambda

  # Interval per m (si no el dones, el fem “prou ample”)
  if (is.null(m_interval)) {
    # escala orientativa via OLS sobre H (si petés, fallback)
    scale0 <- tryCatch({
      b_ols <- coef(lm.fit(x = cbind(1, CLR), y = y))[-1]
      max(abs(b_ols), na.rm = TRUE)
    }, error = function(e) 1)

    if (!is.finite(scale0) || scale0 == 0) scale0 <- 1
    m_interval <- c(-2 * scale0, 2 * scale0)
  }

  # Funció per avaluar l'objectiu a un m concret i un lambda concret
  f_obj <- function(m, LAMBDA) {
    # offset = CLR %*% (m * 1_q)
    off <- as.vector(CLR %*% rep(m, ncol(CLR)))

    # Ajusta gamma amb lasso: penalitza |gamma|
    # on y ~ b0 + CLR %*% gamma + offset
    gfit <- glmnet:::glmnet.fit(
      x = CLR, y = y,
      weights = weights,
      offset = off,
      lambda = LAMBDA,
      intercept = TRUE
    )

    b0    <- as.numeric(gfit$a0)
    gamma <- as.vector(gfit$beta)

    b <- gamma + m
    yhat <- b0 + as.vector(CLR %*% b)

    0.5 * mean((y - yhat)^2) + LAMBDA * sum(abs(b - m)) # = LAMBDA * sum(|gamma|)
  }

  # Loop per cada lambda: optimitza m i guarda coeficients
  out <- lapply(lambdas, function(LAMBDA) {
    m_opt <- optimise(f_obj, interval = m_interval, LAMBDA = LAMBDA)$minimum

    off <- as.vector(CLR %*% rep(m_opt, ncol(CLR)))
    gfit <- glmnet:::glmnet.fit(
      x = CLR, y = y,
      weights = weights,
      offset = off,
      lambda = LAMBDA,
      intercept = TRUE
    )

    gamma <- as.vector(gfit$beta)
    b_clr <- as.vector(gamma + mean(gamma))

    data.frame(
      lambda = LAMBDA,
      m = m_opt,
      var = colnames(LX),
      coef = b_clr
    )
  })

  do.call(rbind, out)
}

# beta_coda_l1_glmnet = function(LX, y, ...){
#
#   CLR = LX - rowMeans(LX)
#
#   m_clr = glmnet(CLR, y, intercept = TRUE, standardize = FALSE)
#
#   do.call(rbind, lapply(m_clr$lambda[1:10], function(LAMBDA){
#     f_opt = function(m, LAMBDA){
#       om = -as.vector(CLR %*% rep(m, ncol(LX)))
#
#       m2 = glmnet:::glmnet.fit(CLR,
#                                y,
#                                weights = rep(1, nrow(CLR)),
#                                offset = om,
#                                lambda = LAMBDA, intercept = TRUE)
#
#
#       b0 = as.numeric(m2$a0)
#       gamma = as.vector(m2$beta)
#
#       beta = gamma - m
#       if(m2$jerr<0)return(Inf)
#       yhat = b0 + as.vector(CLR %*% beta)
#
#       0.5 * mean( (y - yhat)^2 ) + LAMBDA * sum(abs(beta + m))
#     }
#
#     rng <- range(as.vector(m_clr$beta), finite = TRUE)
#     m_optim = optimise(f_opt, interval = rng, LAMBDA = LAMBDA)$minimum
#
#     om = -as.vector(CLR %*% rep(m_optim, ncol(LX)))
#     m2 = glmnet:::glmnet.fit(CLR,
#                              y,
#                              weights = rep(1, nrow(CLR)),
#                              offset = om,
#                              lambda = LAMBDA,
#                              intercept = TRUE)
#
#
#     gamma = as.vector(m2$beta)
#     beta = gamma - mean(gamma)
#
#     data.frame(
#       lambda = LAMBDA,
#       m = m_optim,
#       var = colnames(LX),
#       coef = beta
#     )
#   }))
# }

dplot1 = beta_coda_l1_genlasso(LX, y)
dplot2 = beta_coda_l1_glmnet_m(LX, y)
# dplot2 = beta_coda_l1_glmnet(LX, y)
X_LIMITS = range(dplot1$lambda)
Y_LIMITS = range(c(dplot1$coef, dplot2$coef))
FACTOR = max(dplot1$lambda)/max(dplot2$lambda)
library(ggplot2)
library(patchwork)
p1 = ggplot(data = dplot1) +
  geom_path(aes(x = lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS) +
  scale_y_continuous(limits = Y_LIMITS) +
  labs(title = 'genlasso')
p1
p2 = ggplot(data = dplot2) +
  geom_path(aes(x = FACTOR * lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS)+
  scale_y_continuous(limits = Y_LIMITS) +
  labs(title = 'glmnet')
p1 / p2

  #   m = 0.01
  #
  #   f_opt = function(m){
  #     om = -as.vector(LX %*% rep(m, ncol(LX)))
  #     LAMBDA = 0.0001
  #
  #     CLR = LX - rowMeans(LX)
  #     m2 = glmnet(CLR, y, offset = om, lambda = LAMBDA, intercept = TRUE, standardize = FALSE)
  #     b0_gamma = as.vector(coef(m2))
  #
  #     b0 = b0_gamma[1]
  #     gamma = b0_gamma[-1]
  #
  #     beta = gamma - m
  #
  #     yhat = b0 + as.vector(CLR %*% beta)
  #
  #     0.5 * mean( (y - yhat)^2 ) + LAMBDA * sum(abs(beta + m))
  #   }
  #   m = optimise(f_opt, c(-3,3))$minimum
  #   om = -as.vector(LX %*% rep(m, ncol(LX)))
  #   LAMBDA = 0.0001
  #
  #   m2 = glmnet(CLR, y, offset = om, lambda = LAMBDA, intercept = TRUE, standardize = FALSE)
  #   b0_gamma = as.vector(coef(m2))
  #
  #   b0 = b0_gamma[1]
  #   gamma = b0_gamma[-1]
  #
  #   beta = gamma - mean(gamma)
  #   beta
  #
  #
  #   coef(m2)
  #   sum(abs(coef(m2)[-1]+0.00274646375))
  #   do.call(rbind, lapply(seq_along(m2$lambda),
  #                         function(i){
  #                           data.frame(
  #                             lambda = m2$lambda[i],
  #                             var = colnames(LX),
  #                             coef = as.vector(m2$beta[-1,i])
  #                           )
  #                         }))
  #   # lambda_new = 1/sum(w_aug) * (n / sum(w_aug))  * lambda
  #   # coef(m2, s = lambda_new)[-(1:2)]
  # }

  #
  # beta_coda_l1_glmnet(LX, y, lambda = 0.05)
  #
  # sapply(seq(0, 0.1, 0.001), function(l) beta_coda_l1_glmnet(LX, y, lambda = l)) |> t()
  #
  #
  #
  # coef(m2, s = lam_best)[-(1:2)]
  #
  # b1_clr
  # b2_clr
  # plot(m2, xvar = 'lambda')
  #
  #
  #
  # n <- nrow(LX); p <- ncol(LX)
  #
  # # Disseny amb intercept explícit (1a columna)
  # Xtilde <- cbind(1, LX)
  #
  # # Pseudo-observació per imposar sum(u)=0 sense afectar intercept
  # x_star <- c(0, rep(1, p))  # 0 a l'intercept, 1 als coeficients de log(X)
  # y_star <- 0
  #
  # # Pes enorme
  # w_star <- 1e6
  #
  # X_aug <- rbind(Xtilde, x_star)
  # y_aug <- c(y, y_star)
  # w_aug <- c(rep(1, n), w_star)
  #
  # # Penalització: NO penalitzar intercept (1a col), sí la resta
  # pf <- c(0, rep(1, p))
  #
  # fit <- glmnet(
  #   x = X_aug, y = y_aug,
  #   weights = w_aug,
  #   intercept = FALSE,      # IMPORTANT: intercept ja és a Xtilde
  #   standardize = FALSE,
  #   penalty.factor = pf
  # )
  #
  # # Coeficients per un lambda concret
  # b <- as.numeric(coef(fit, s = 0.01/n)[-c(1:2)])  # treu l'intercept explícit
  #
  # sum(b)  # hauria de sortir molt a prop de 0
  # b
  #
  #
  #
  #
  # B = ilr_basis(X)    # 5 x 4
  # H = coordinates(X, B)
  # H = log(X) %*% B
  #
  #
  # m2 = glmnet(y = y, x = H %*% t(B), standardize = FALSE)
  #
  # coef(m1, lambda = 1)$beta[-1] %*% t(B)
  # coef(m2, s = 1/nrow(X))[-1,1] |> (function(x) x-mean(x))()
  #
