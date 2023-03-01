library(coda.base)
library(readxl)
l_df = list()
# A.1 Chemical compositions of Roman-British pottery
pottery = read_cdp('~/Software/CoDaPack/data/pottery.cdp') |>
  transform(ident = as.integer(ident))
l_df[['pottery']] = pottery

# A.2 Arctic lake sediments at different depths
arctic_lake = read_cdp('~/Software/CoDaPack/data/arctic_lake.cdp')
l_df[['arctic_lake']] = arctic_lake

# A.3 Household budget patterns
household_budget = read_cdp('~/Software/CoDaPack/data/householdbudget.cdp') |>
  transform(Id = as.integer(Id))
l_df[['household_budget']] = household_budget

# A.4 Milk composition study
milk_cows = read_excel('~/Software/CoDaPack/data/milkcows.xls')
str(milk_cows)
l_df[['milk_cows']] = milk_cows

# A.5 A statistician’s time budget
statistitian_time = read_cdp('~/Software/CoDaPack/data/statisticiantimebudget.cdp') |>
  transform(Day = as.integer(Day))
str(statistitian_time)
l_df[['statistitian_time']] = statistitian_time

# A.6 The MN blood system
blood_mn = read_cdp('~/Software/CoDaPack/data/bloodMN.cdp')
blood_mn$popul_code = NULL
str(blood_mn)
l_df[['blood_mn']] = blood_mn

# A.7 Mammal’s milk
mammals_milk =read_cdp('~/Software/CoDaPack/data/mammalsmilk.cdp')
mammals_milk$code = NULL
str(mammals_milk)
l_df[['mammals_milk']] = mammals_milk

# A.8 Calc-alkaline and tholeiitic volcanic rocks
petrafm = read_cdp('~/Software/CoDaPack/data/petrafm.cdp')
petrafm = petrafm$Petrafm
petrafm$row.names = NULL
str(petrafm)
l_df[['petrafm']] = petrafm

# A.9 Concentration of minor elements in carbon ashes
read_excel_df = function(file) as.data.frame(read_excel(file))
montana = read_excel_df('~/Software/CoDaPack/data/montana.xls')
montana$type = NULL
montana$ident = as.integer(montana$ident)
str(montana)
l_df[['montana']] = montana

# A.10 Paleocological compositions
foraminiferals = read_cdp('~/Software/CoDaPack/data/foraminiferalsNA.cdp')
foraminiferals = foraminiferals[[1]]
foraminiferals$glob_obesa.dl = NULL
foraminiferals$glob_triloba.dl = NULL
foraminiferals$zeros = NULL
foraminiferals$code = as.integer(foraminiferals$code)
str(foraminiferals)
l_df[['foraminiferals']] = foraminiferals

# A.11 Pollen composition in fossils
library(readr)
read_delim_df = function(file) as.data.frame(read_delim(file))
pollen = read_delim_df('~/Software/CoDaPack/data/pollen.txt')
str(pollen)
l_df[['pollen']] = pollen

# A.12 Food consumption in European countries
alimentation = read_cdp('~/Software/CoDaPack/data/alimentation.cdp')
alimentation$ident = NULL
str(alimentation)
l_df[['alimentation']] = alimentation

# A.13 Household expenditures
house_expend = read_cdp('~/Software/CoDaPack/data/houseexpend.cdp')$hexpcomplete
str(house_expend)
l_df[['house_expend']] = house_expend

# A.14 Serum proteins
serprot = read_excel_df('~/Software/CoDaPack/data/serprot.xls')
serprot$case = as.integer(serprot$case)
serprot$disease = as.integer(serprot$disease)
serprot$group =  as.integer(serprot$group)
str(serprot)
l_df[['serprot']] = serprot

# A.15 Physical activity and body mass index
bmi_activity = read_cdp('~/Software/CoDaPack/data/BMIPhisActi.cdp')
bmi_activity$zBMI.dl = NULL
str(bmi_activity)
l_df[['bmi_activity']] = bmi_activity

# A.16 Hotel posts in social media
weibo_hotels = read_excel_df('~/Software/CoDaPack/data/weibo_hotels.xls')
str(weibo_hotels)
l_df[['weibo_hotels']] = weibo_hotels

# A.17 The waste composition in Catalonia
waste = read_cdp('~/Software/CoDaPack/data/waste.cdp')
str(waste)
l_df[['waste']] = waste

# A.18 Employment distribution in EUROSTAT countries
eurostat_employment = read_cdp('~/Software/CoDaPack/data/eurostat_employment.cdp')
eurostat_employment = eurostat_employment[,1:17]
str(eurostat_employment)
l_df[['eurostat_employment']] = eurostat_employment

for(i in seq_along(l_df)){
  write_csv(l_df[[i]], file = sprintf("testing/csv_files/%s.csv", names(l_df)[i]))
  eval(parse(text = sprintf("save(%s, file = 'data/%s.RData')", names(l_df)[i], names(l_df)[i])))
}
#
# arctic_lake = read_cdp('~/Software/CoDaPack/data/arctic_lake.cdp')
# save(arctic_lake, file = 'data/arctic_lake.RData')
#
# blood_mn = read_cdp('~/Software/CoDaPack/data/bloodMN.cdp')
# blood_mn
#
# bmi_activity = read_cdp('~/Software/CoDaPack/data/BMIPhisActi.cdp')
# bmi_activity
#
# eurostat_employment = read_cdp('~/Software/CoDaPack/data/eurostat_employment.cdp')
# eurostat_employment = eurostat_employment[,1:17]
# eurostat_employment
#
# foraminiferals = read_cdp('~/Software/CoDaPack/data/foraminiferalsNA.cdp')
# foraminiferals_na = foraminiferals[[2]]
#
# house_expend = read_cdp('~/Software/CoDaPack/data/houseexpend.cdp')$hexpcomplete
# house_expend
#
# household_budget = read_cdp('~/Software/CoDaPack/data/householdbudget.cdp')
# save(household_budget, file = 'data/household_budget.RData')
#
# mammals_milk =read_cdp('~/Software/CoDaPack/data/mammalsmilk.cdp')
# mammals_milk
#
# petrafm = read_cdp('~/Software/CoDaPack/data/petrafm.cdp')
# petrafm = petrafm$Petrafm
# petrafm
#
# pottery = read_cdp('~/Software/CoDaPack/data/pottery.cdp')
# save(pottery, file='data/pottery.RData')
#
# statistitian_time = read_cdp('~/Software/CoDaPack/data/statisticiantimebudget.cdp')
# statistitian_time
#
# waste = read_cdp('~/Software/CoDaPack/data/waste.cdp')
# waste
#
# trondelag = read_cdp('~/Software/CoDaPack/data/trondelagO.cdp')
# trondelag
#
# fu = read_cdp('~/Software/CoDaPack/data/FU.cdp')$FU
# fu
#
# library(readxl)
# read_excel_df = function(file) as.data.frame(read_excel(file))
# bacteria = read_excel_df('~/Software/CoDaPack/data/bacteria.xlsx')
# bacteria
#
# foraminiferals = read_excel_df('~/Software/CoDaPack/data/foraminiferal.xls')
#
# milk_cows = read_excel_df('~/Software/CoDaPack/data/milkcows.xls')
# milk_cows
#
# montana = read_excel_df('~/Software/CoDaPack/data/montana.xls')
# montana
#
# serprot = read_excel_df('~/Software/CoDaPack/data/serprot.xls')
# serprot
#
# weibo_hotels = read_excel_df('~/Software/CoDaPack/data/weibo_hotels.xls')
# weibo_hotels
#
# library(readr)
# read_delim_df = function(file) as.data.frame(read_delim(file))
#
# pollen = read_delim_df('~/Software/CoDaPack/data/pollen.txt')
# pollen
