library(dplyr)
library(coda.base)

data('catalan_elections')

by_prov = function(YEAR){
  parliament_elections %>%
    subset(year == YEAR) %>%
    left_join(municipality_information %>% select(mun, com), by = 'mun') %>%
    group_by(com, party) %>%
    summarise(
      votes = sum(votes, na.rm = TRUE)
    ) %>%
    spread(party, votes) %>%
    ungroup()
}

prinbal = function(.parl){
  X = .parl %>% select(-prov, -com, -other)
  pb = pb_basis(X, method = 'exact')
  rownames(pb) = names(X)
  colnames(pb) = paste0('pb', 1:ncol(pb))
  pb
}

by_prov(1980)
by_prov(1984)
by_prov(1988)
by_prov(1992)
by_prov(1995)
by_prov(1999)
by_prov(2003)
by_prov(2006)
by_prov(2010)
by_prov(2012)
by_prov(2015)
by_prov(2017)

prinbal(by_prov(1999))
prinbal(parl2010_com)
prinbal(parl2012_com)
prinbal(parl2015_com)
prinbal(parl2017_com)


