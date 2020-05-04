library(dplyr)
library(tidyr)
load('testing/emex.RData')
load('testing/parliament2015.RData')

mun_renaming = function(.data) .data  %>%
  mutate(
    mun = gsub('^Brunyola$', 'Brunyola i Sant Martí Sapresa', mun),
    mun = gsub('^Calonge$', 'Calonge i Sant Antoni', mun),
    mun = gsub('^Palmerola$', 'Llosses, les', mun)
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
                     jxsi, psc, pp, catsp, cs, cup, other) %>% mun_renaming(), by = 'mun') %>%
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

parl2012 = read_delim('testing/t6011mun201200.csv', skip = 5, delim = ";",
                      col_types = 'cciiiiiiiii_') %>%
  select(mun = Literal, cs = `C's`, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, cup = CUP, pp = PP, other = `Altres candidatures`) %>%
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

parl2012_mun = full_join(data_prov_pos, parl2012, by = 'mun')
subset(parl2012_mun, !complete.cases(parl2012_mun))

parl2012 = subset(parl2012_mun, complete.cases(parl2012_mun))

parl2012_com = parl2012 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, cs, cup, other) %>%
  ungroup()

parl2010 = read_delim('testing/t6011mun201000.csv', skip = 6, delim = ";",
                      col_types = 'cciiiiii_ii_') %>%
  select(mun = Literal, cs = `C's`, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
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

parl2010_mun = full_join(data_prov_pos, parl2010, by = 'mun')
subset(parl2010_mun, !complete.cases(parl2010_mun))

parl2010 = subset(parl2010_mun, complete.cases(parl2010_mun))

parl2010_com = parl2010 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, cs, other) %>%
  ungroup()

parl2006 = read_delim('testing/t6011mun200600.csv', skip = 6, delim = ";",
                      col_types = 'cciiiiii_ii_') %>%
  select(mun = Literal, cs = `C's`, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
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

parl2006_mun = full_join(data_prov_pos, parl2006, by = 'mun')
subset(parl2006_mun, !complete.cases(parl2006_mun))

parl2006 = subset(parl2006_mun, complete.cases(parl2006_mun))

parl2006_com = parl2006 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, cs, other) %>%
  ungroup()

parl2003 = read_delim('testing/t6011mun200300.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl2003_mun = full_join(data_prov_pos, parl2003, by = 'mun')
subset(parl2003_mun, !complete.cases(parl2003_mun))

parl2003 = subset(parl2003_mun, complete.cases(parl2003_mun))

parl2003_com = parl2003 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1999 = read_delim('testing/t6011mun199900.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1999_mun = full_join(data_prov_pos, parl1999, by = 'mun')
subset(parl1999_mun, !complete.cases(parl1999_mun))

parl1999 = subset(parl1999_mun, complete.cases(parl1999_mun))

parl1999_com = parl1999 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1995 = read_delim('testing/t6011mun199500.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1995_mun = full_join(data_prov_pos, parl1995, by = 'mun')
subset(parl1995_mun, !complete.cases(parl1995_mun))

parl1995 = subset(parl1995_mun, complete.cases(parl1995_mun))

parl1995_com = parl1995 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1992 = read_delim('testing/t6011mun199200.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  mun_renaming() %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1992_mun = full_join(data_prov_pos, parl1992, by = 'mun')
subset(parl1992_mun, !complete.cases(parl1992_mun))

parl1992 = subset(parl1992_mun, complete.cases(parl1992_mun))

parl1992_com = parl1992 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1988 = read_delim('testing/t6011mun198800.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  mun_renaming() %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1988_mun = full_join(data_prov_pos, parl1988, by = 'mun')
subset(parl1988_mun, !complete.cases(parl1988_mun))

parl1988 = subset(parl1988_mun, complete.cases(parl1988_mun))

parl1988_com = parl1988 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1984 = read_delim('testing/t6011mun198400.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  mun_renaming() %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1984_mun = full_join(data_prov_pos, parl1984, by = 'mun')
subset(parl1984_mun, !complete.cases(parl1984_mun))

parl1984 = subset(parl1984_mun, complete.cases(parl1984_mun))

parl1984_com = parl1984 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

parl1980 = read_delim('testing/t6011mun198000.csv', skip = 6, delim = ";",
                      col_types = 'cciiiii__ii_') %>%
  select(mun = Literal, ciu = CiU, erc = ERC, psc = PSC,
         ic = IC, pp = PP, other = `Altres candidatures`) %>%
  subset(mun != 'Catalunya') %>%
  mutate(
    mun = gsub("Residents absents ", "", mun),
    mun = gsub("Residents absents ", "", mun)) %>%
  mun_renaming() %>%
  group_by(mun) %>%
  gather(party, votes, ciu:other) %>%
  group_by(mun, party) %>%
  summarise(votes = sum(votes, na.rm=TRUE)) %>%
  spread(party, votes) %>%
  ungroup()

parl1980_mun = full_join(data_prov_pos, parl1980, by = 'mun')
subset(parl1980_mun, !complete.cases(parl1980_mun))

parl1980 = subset(parl1980_mun, complete.cases(parl1984_mun))

parl1980_com = parl1980 %>%
  gather(party, num, ciu:psc) %>%
  group_by(prov, com, party) %>%
  summarise(
    num = sum(num, na.rm = TRUE)
  ) %>%
  spread(party, num) %>%
  select(prov, com, ciu, erc, psc, pp, ic, other) %>%
  ungroup()

get_votes = function(.data, year){
  .data %>%
    select(-prov,-com, -lon, -lat) %>%
    gather(party, votes, -mun) %>%
    mutate(year = year)
}

parliament_elections = bind_rows(
  parl1980 %>% get_votes(1980),
  parl1984 %>% get_votes(1984),
  parl1988 %>% get_votes(1988),
  parl1992 %>% get_votes(1992),
  parl1995 %>% get_votes(1995),
  parl1999 %>% get_votes(1999),
  parl2003 %>% get_votes(2003),
  parl2006 %>% get_votes(2006),
  parl2010 %>% get_votes(2010),
  parl2012 %>% get_votes(2012),
  parl2015 %>% get_votes(2015),
  parl2017 %>% get_votes(2017)) %>%
  select(year, mun, party, votes) %>%
  arrange(desc(year), mun, party, votes)
municipality_information = data_prov_pos

save(municipality_information, parliament_elections,
     file = 'data/catalan_elections.RData')
