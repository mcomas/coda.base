library(dplyr)
library(tidyr)
load('testing/emex.RData')
load('testing/parliament2015.RData')

mun_renaming = function(.data) .data  %>%
  mutate(
    mun = gsub('^Brunyola$', 'Brunyola i Sant Martí Sapresa', mun),
    mun = gsub('^Calonge$', 'Calonge i Sant Antoni', mun)
  )
emex.data[emex.data$com == 'Alt Camp', 'prov'] = 'Tarragona'
emex.data[emex.data$prov == 'LLeida', 'prov'] = 'Lleida'

data_prov = emex.data %>% select(prov, com, mun) %>%
  full_join(parlament2015 %>%
              select(mun = name, lon, lat), by = 'mun') %>%
  mun_renaming()
data_prov %>% subset(!complete.cases(data_prov))
data_prov_pos = data_prov %>% subset(complete.cases(data_prov))


parl2015_mun = data_prov_pos %>%
  full_join(parlament2015 %>%
              select(mun = name,
                     jxsi, psc, pp, catsp, cs, cup, other,
                     pop) %>% mun_renaming(), by = 'mun') %>%
  ungroup()
parl2015_mun %>% subset(!complete.cases(parl2015_mun))
parl2015 = parl2015_mun %>%
  subset(complete.cases(parl2015_mun))

parl2015_com = parl2015 %>%
  gather(party, num, jxsi:other) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, jxsi, psc, pp, catsp, cs, cup, other) %>%
  ungroup()


library(readr)
parl2017 = read_delim('testing/t6011mun.csv', skip = 5, delim = ";",
               col_types = 'cciiiiiiiii_') %>%
  select(mun = Literal, cs = `C's`, jxcat = JUNTSxCAT, erc = `ERC-Cat Sí`, psc = PSC,
         catsp = `Cat Comú-Podem`, cup = CUP, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  group_by(mun) %>%
  gather(party, votes, cs:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl2017_mun = full_join(data_prov_pos, parl2017, by = 'mun')
subset(parl2017_mun, !complete.cases(parl2017_mun))

parl2017 = subset(parl2017_mun, complete.cases(parl2017_mun))

parl2017_com = parl2017 %>%
  gather(party, num, catsp:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, jxcat, erc, psc, pp, catsp, cs, cup, other) %>%
  ungroup()

save(parl2015, parl2017,
     parl2015_com, parl2017_com, file = 'data/catalan_elections.RData')
