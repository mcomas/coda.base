library(dplyr)
library(tidyr)
load('testing/catalan_elections.RData')
parliament2017 = parliament_elections %>%
  left_join(municipality_information) %>%
  subset(year == 2017) %>%
  select(com, party, votes) %>%
  group_by(com, party) %>%
  summarise(votes = sum(votes)) %>%
  spread(party, votes) %>%
  ungroup() %>%
  select(com, cs, jxcat, erc, psc, catsp, cup, pp, other)

save(parliament2017, file='data/catalan_elections_2017.RData',
     compress = 'xz')
