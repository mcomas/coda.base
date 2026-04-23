library(glmnet)
library(genlasso)
library(coda.base)

X = as.matrix(waste[,5:9])[1:15,]
y = waste$floating_population[1:15]

LX = log(X)
K = ncol(X)
n = nrow(X)

beta_clr_l1_genlasso = function(LX, y, ...){
  B = ilr_basis(X)
  H = LX %*% B
  X1 = cbind(1, H)
  D = cbind(0, B)
  m1 = genlasso(y, X1, D, svd = TRUE, ...)
  do.call(rbind, lapply(seq_along(m1$lambda),
         function(i){
           data.frame(
             lambda = m1$lambda[i],
             var = colnames(LX),
             coef = as.vector(m1$beta[-1,i]  %*% t(B))
           )
         }))
}
beta_clr_l1_glmnet = function(LX, y, ...){
  CLR = LX - rowMeans(LX)
  X_aug = rbind(
    c(0, rep(1, K)),
    cbind(1, CLR))
  y_aug = c(0, y)
  w_aug = c(10*n, rep(1, n))
  pf = c(0, rep(1, K))
  m2 = glmnet(X_aug, y_aug, weights = w_aug, penalty.factor = pf,
              intercept = FALSE, standardize = FALSE, ...)
  do.call(rbind, lapply(seq_along(m2$lambda),
                        function(i){
                          data.frame(
                            lambda = m2$lambda[i],
                            var = colnames(LX),
                            coef = as.vector(m2$beta[-1,i])
                          )
                        }))
}
beta_clr_l1_glmnet_no_sum_0 = function(LX, y, ...){
  CLR = LX - rowMeans(LX)

  m2 = glmnet(CLR, y, intercept = TRUE, standardize = FALSE, ...)
  do.call(rbind, lapply(seq_along(m2$lambda),
                        function(i){
                          data.frame(
                            lambda = m2$lambda[i],
                            var = colnames(LX),
                            coef = as.vector(m2$beta[,i]-mean(m2$beta[,i]))
                          )
                        }))
}
dplot1 = beta_clr_l1_genlasso(LX, y)
dplot2 = beta_clr_l1_glmnet(LX, y)
dplot3 = beta_clr_l1_glmnet_no_sum_0(LX, y)
X_LIMITS = range(dplot1$lambda)
Y_LIMITS = range(c(dplot1$coef, dplot2$coef, dplot3$coef))
FACTOR1 = max(dplot1$lambda)/max(dplot2$lambda)
FACTOR2 = max(dplot1$lambda)/max(dplot3$lambda)
print(FACTOR1)
library(ggplot2)
library(patchwork)
p1 = ggplot(data = dplot1) +
  geom_path(aes(x = lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS) +
  scale_y_continuous(limits = Y_LIMITS) +
  labs(title = 'genlasso')
p2 = ggplot(data = dplot2) +
  geom_path(aes(x = FACTOR1 * lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS)+
  scale_y_continuous(limits = Y_LIMITS) +
  labs(title = 'glmnet')
p3 = ggplot(data = dplot3) +
  geom_path(aes(x = FACTOR2 * lambda, y = coef, col = var)) +
  scale_x_continuous(limits = X_LIMITS)+
  scale_y_continuous(limits = Y_LIMITS) +
  labs(title = 'glmnet (no 0)')
p1 / p2 / p3

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
