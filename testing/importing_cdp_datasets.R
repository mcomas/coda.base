library(coda.base)
alimentation = read_cdp('~/Software/CoDaPack/data/alimentation.cdp')
alimentation

arctic_lake = read_cdp('~/Software/CoDaPack/data/arctic_lake.cdp')
save(arctic_lake, file = 'data/arctic_lake.RData')

blood_mn = read_cdp('~/Software/CoDaPack/data/bloodMN.cdp')
blood_mn

bmi_activity = read_cdp('~/Software/CoDaPack/data/BMIPhisActi.cdp')
bmi_activity

eurostat_employment = read_cdp('~/Software/CoDaPack/data/eurostat_employment.cdp')
eurostat_employment = eurostat_employment[,1:17]
eurostat_employment

foraminiferals = read_cdp('~/Software/CoDaPack/data/foraminiferalsNA.cdp')
foraminiferals_na = foraminiferals[[2]]

house_expend = read_cdp('~/Software/CoDaPack/data/houseexpend.cdp')$hexpcomplete
house_expend

household_budget = read_cdp('~/Software/CoDaPack/data/householdbudget.cdp')
save(household_budget, file = 'data/household_budget.RData')

mammals_milk =read_cdp('~/Software/CoDaPack/data/mammalsmilk.cdp')
mammals_milk

petrafm = read_cdp('~/Software/CoDaPack/data/petrafm.cdp')
petrafm = petrafm$Petrafm
petrafm

pottery = read_cdp('~/Software/CoDaPack/data/pottery.cdp')
save(pottery, file='data/pottery.RData')

statistitian_time = read_cdp('~/Software/CoDaPack/data/statisticiantimebudget.cdp')
statistitian_time

waste = read_cdp('~/Software/CoDaPack/data/waste.cdp')
waste

trondelag = read_cdp('~/Software/CoDaPack/data/trondelagO.cdp')
trondelag

fu = read_cdp('~/Software/CoDaPack/data/FU.cdp')$FU
fu

library(readxl)
read_excel_df = function(file) as.data.frame(read_excel(file))
bacteria = read_excel_df('~/Software/CoDaPack/data/bacteria.xlsx')
bacteria

foraminiferals = read_excel_df('~/Software/CoDaPack/data/foraminiferal.xls')

milk_cows = read_excel_df('~/Software/CoDaPack/data/milkcows.xls')
milk_cows

montana = read_excel_df('~/Software/CoDaPack/data/montana.xls')
montana

serprot = read_excel_df('~/Software/CoDaPack/data/serprot.xls')
serprot

weibo_hotels = read_excel_df('~/Software/CoDaPack/data/weibo_hotels.xls')
weibo_hotels

library(readr)
read_delim_df = function(file) as.data.frame(read_delim(file))

pollen = read_delim_df('~/Software/CoDaPack/data/pollen.txt')
pollen
