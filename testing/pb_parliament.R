library(dplyr)
library(coda.base)

data('catalan_elections')

X2015 = select(parl2015_com, jxsi:cup)
X2017 = select(parl2017_com, jxcat:cup)

pb2015 = pb_basis(X2015, method = 'exact')
rownames(pb2015) = names(X2015)
pb2017 = pb_basis(X2017, method = 'exact')
rownames(pb2017) = names(X2017)

pb2015
pb2017


X.parl = bind_rows(
  X2015 %>%
    mutate(year = '2015'),
  X2017 %>%
    mutate(jxsi = jxcat + erc) %>%
    select_(.dots = names(X2015)) %>%
    mutate(year = '2017'))

library(broom)
pb.period = X.parl %>%
  group_by(year) %>%
  do(pb = pb_basis(.[,1:6], method = 'exact')) %>%
  tidy(pb) %>%
  mutate(
    vars = names(X.parl[,1:6])
  ) %>%
  select(year, vars, starts_with('X'))
pb.period
