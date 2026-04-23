library(glmnet)
library(genlasso)
library(coda.base)

X = as.matrix(waste[,5:9])
y = waste$floating_population

LX = log(X)
K = ncol(X)
n = nrow(X)

beta_pairwise_l1_glmnet = function(LX, y, ...){
  X_aug = rbind(c(0, rep(1, ncol(LX))), cbind(1, LX))
  y_aug = c(0, y)
  w = c(10*nrow(LX), rep(1, nrow(LX)))
  m1 = glmnet(X_aug, y_aug, weights = w, intercept = FALSE)

  B = ilr_basis(LX)
  H = LX %*% B
  b0_norm = sum( as.vector(B %*% coef(lm(y~H))[-1])^2 )
  b0_norm_weigh = sum( (coef(m1, s = 0)[-c(1:2)])^2 )
  FACTOR = b0_norm / b0_norm_weigh
  do.call(rbind, lapply(seq_along(m1$lambda),
         function(i){
           data.frame(
             lambda = m1$lambda[i],
             var = colnames(LX),
             coef = FACTOR * as.vector(m1$beta[-1,i] - mean(m1$beta[-1,i]))
           )
         }))
}
beta_pairwise_l1_glmnet_pairs = function(LX, y, ...){
  B = pairwise_basis(LX)
  H = LX %*% B
  m1 = glmnet(H, y, intercept = TRUE)
  CLR_BETA = coordinates(composition(as.matrix(t(m1$beta)), basis = pairwise_basis(LX)), 'clr')
  do.call(rbind, lapply(seq_along(m1$lambda),
                        function(i){
                          data.frame(
                            lambda = m1$lambda[i],
                            var = colnames(LX),
                            coef = CLR_BETA[i,]
                          )
                        }))
}
dplot1 = beta_pairwise_l1_glmnet(LX, y)
dplot2 = beta_pairwise_l1_glmnet_pairs(LX, y)

X_LIMITS = range(dplot1$lambda)
Y_LIMITS = range(c(dplot1$coef, dplot2$coef))
FACTOR = max(dplot1$lambda) / max(dplot2$lambda)
library(ggplot2)
library(patchwork)
p1 = ggplot(data = dplot1) +
  geom_path(aes(x = lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS) +
  scale_y_continuous(limits = Y_LIMITS)

p2 = ggplot(data = dplot2) +
  geom_path(aes(x = FACTOR * lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS) +
  scale_y_continuous(limits = Y_LIMITS)

p1 / p2

#
# beta_clr_l1_glmnet(LX, y, lambda = 0.05)
#
# sapply(seq(0, 0.1, 0.001), function(l) beta_clr_l1_glmnet(LX, y, lambda = l)) |> t()
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
