library(Rcpp)
sourceCpp("testing/genlasso_path.cpp")

library(coda.base)

summary(y <- waste$floating_population)
summary(X <- waste[,5:9])

B = ilr_basis(X)
H = coordinates(X, B)

LAMBDAS = seq(0, 5, length= 100)
sol_clr_l1 = generalized_lasso_path_admm_cpp(as.matrix(H), y, B, LAMBDAS, intercept = TRUE,
                                             tol = 1e-10, max_iter = 10000)

library(tidyverse)
dplot_clr_l1 = as_tibble(t(sol_clr_l1$coef), .name_repair = ~ c('b0', colnames(B))) %>%
  select(-b0) %>%
  composition(basis = B) %>%
  pivot_longer(everything()) %>%
  mutate(lambda = rep(LAMBDAS, each = nrow(B)))

p1 = ggplot(data=dplot_clr_l1) +
  geom_hline(yintercept = 1/nrow(B), linetype = 'dotted') +
  geom_line(aes(x = lambda, y = value, col = name))

sol_coda_l1 = generalized_lasso_path_admm_cpp(cbind(as.matrix(H),0), y, cbind(B,-1), LAMBDAS,
                                              intercept = TRUE, tol = 1e-10, max_iter = 10000)

dplot_coda_l1 = as_tibble(t(sol_coda_l1$coef), .name_repair = ~ c('b0', colnames(B), 'm')) %>%
  select(-m, -b0) %>%
  composition(basis = B) %>%
  pivot_longer(everything()) %>%
  mutate(lambda = rep(LAMBDAS, each = nrow(B)))

p2 = ggplot(data=dplot_coda_l1) +
  geom_hline(yintercept = 1/nrow(B), linetype = 'dotted') +
  geom_line(aes(x = lambda, y = value, col = name))

library(patchwork)
p1 / p2
