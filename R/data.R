#' Results of catalan parliament elections in 2017 by regions.
#' @name parliament2017
#' @format A data frame with 42 rows and 9 variables:
#' \describe{
#' \item{com}{Region}
#' \item{cs}{Votes to Ciutadans party}
#' \item{jxcat}{Votes to Junts per Catalunya party}
#' \item{erc}{Votes to Esquerra republicana de Catalunya party}
#' \item{psc}{Votes to Partit socialista de Catalunya party}
#' \item{catsp}{Votes to Catalunya si que es pot party}
#' \item{cup}{Votes to Candidatura d'unitat popular party}
#' \item{pp}{Votes to Partit popular party}
#' \item{other}{Votes to other parties}
#' }
#' @source \url{https://www.idescat.cat/tema/elecc}
"parliament2017"

#' Chemical compositions of Romano-British pottery
#' @description
#' The pottery data set consists of data pertaining to the chemical composition of 45
#' specimens of Romano-British pottery. The method used to generate these data is
#' atomic absorption spectophotometry, and readings for nine oxides
#' (Al2O3, Fe2O3, MgO, CaO, Na2O, K2O, TiO2 , MnO, BaO) are provided. These samples come
#' from five different kiln sites.
#' @name pottery
"pottery"

#' Arctic lake sediments at different depths
#' @description
#' The arctic lake data set records the [sand, silt, clay] compositions of 39 sediment
# samples at different water depths in an Arctic lake.
#' @name arctic_lake
"arctic_lake"

#' Household budget patterns
#' @description
#' In a sample survey of single persons living alone in rented accommodation, twenty
#' men and twenty women were randomly selected and asked to record over a period of
#' one month their expenditures on the following four mutually exclusive and
#' exhaustive commodity groups:
#' * Hous: Housing, including fuel and light.
#' * Food: Foodstuffs, including alcohol and tobacco.
#' * Serv: Services, including transport and vehicles.
#' * Other: Other goods, including clothing, footwear and durable goods.
#' @name household_budget
"household_budget"

#' Food consumption in European countries
#'
#' The alimentation data set contains the percentages of consumption of several
#' types of food in 25 European countries during the 80s. The categories are:
#' * RM: red meat (pork, veal, beef),
#' * WM: white meat (chicken),
#' * E: eggs,
#' * M: milk,
#' * F: fish,
#' * C: cereals,
#' * S: starch (potatoes),
#' * N: nuts, and
#' * FV: fruits and vegetables.
#'
#' Moreover, the dataset contains a categorical variable that
#' shows if the country is from the North or a Southern Mediterranean country. In
#' addition, the countries are classified as Eastern European or as Western European.
#' @name alimentation
"alimentation"

#' The MN blood system
#' @description
#'
#' In humans the main blood group systems are the ABO system, the Rh system and
#' the MN system. The MN blood system is a system of blood antigens also related
#' to proteins of the red blood cell plasma membrane. The inheritance pattern of the
#' MN blood system is autosomal with codominance, a type of lack of dominance in
#' which the heterozygous manifests a phenotype totally distinct from the homozygous.
#' The possible phenotypical forms are three blood types: type M blood, type
#' N blood and type MN blood. The frequencies of M, N and MN blood types vary
#' widely depending on the ethnic population. However, the Hardy-Weinberg
#' principle states that allele and genotype frequencies in a population will
#' remain constant from generation to generation in the absence of other
#' evolutionary influences. This implies that, in the long run, it holds that
#'
#' \eqn{\frac{x_{MM}x_{NN}}{x_{MN}} = \frac{1}{4}}{ (x_MM x_NN) / x_MN = 1/4}
#'
#'where xM M and xN N are the genotype relative frequencies of MM and NN
#'homozygotes, respectively, and xM N is the genotype relative frequency
#'of MN heterozygotes. This principle was named after G.H. Hardy and W.
#' Weinberg demonstrated it mathematically.
#' @name blood_mn
"blood_mn"

#' Physical activity and body mass index
#'
#' @description
#' The `bmi_activity` data set records the proportion of daily time spent to sleep
#' (sleep), sedentary behaviour (sedent), light physical activity (Lpa), moderate
#' physical activity (Mpa) and vigorous physical activity (Vpa) measured on a small
#' population of 393 children. Moreover the standardized body mass index (zBMI) of
#' each child was also registered.
#'
#' This data set was used in the example of the article (Dumuid et al. 2019) to examine the
#' expected differences in zBMI for reallocations of daily time between sleep, physical
#' activity and sedentary behaviour. Because the original data is confidential, the
#' data set BMIPhisActi includes simulated data that mimics the main features of the
#' original data.
#'
#' @references
#' D. Dumuid, Z. Pedisic, T.E. Stanford, J.A. Martín-Fernández, K. Hron, C.
#' Maher, L.K. Lewis and T.S. Olds, \emph{The Compositional Isotemporal Sub-
#'  stitution Model: a Method for Estimating Changes in a Health Outcome
#' for Reallocation of Time between Sleep, Sedentary Behaviour, and Physical
#' Activity}. Statistical Methods in Medical Research \strong{28}(3) (2019), 846–857
#' @name bmi_activity
"bmi_activity"

#' Employment distribution in EUROSTAT countries
#'
#' @description
#'
#' According to the three–sector theory, as a country’s economy develops, employment
#' shifts from the primary sector (raw material extraction: farming, hunting, fishing,
#' mining) to the secondary sector (industry, energy and construction) and finally to
#' the tertiary sector (services). Thus, a country’s employment distribution can be
#' used as a predictor of economic wealth.
#'
#' The `eurostat_employment` data set contains EUROSTAT data on employment
#' aggregated for both sexes, and all ages distributed by economic activity
#' (classification 1983-2008, NACE Rev. 1.1) in 2008 for the 29 EUROSTAT member
#' countries, thus reflecting reality just before the 2008 financial crisis.
#' Country codes in alphabetical order according to the country name in its
#' own language are: Belgium (BE), Cyprus (CY), Czechia (CZ), Denmark (DK),
#' Deutchland–Germany (DE), Eesti–Estonia (EE), Eire–Ireland (IE),
#' España–Spain (ES), France (FR), Hellas-Greece (GR), Hrvatska–Croatia (HR),
#' Iceland (IS), Italy (IT), Latvia (LV), Lithuania (LT), Luxembourg (LU),
#' Macedonia (MK), Magyarország-Hungary (HU), Malta (MT), Netherlands (NL),
#' Norway (NO), Österreich–Austria (AT), Portugal (PT), Romania (RO),
#' Slovakia (SK), Suomi–Finland (FI), Switzerland (CH), Turkey (TR),
#' United Kingdom (GB).
#'
#' A key related variable is the logarithm of gross domestic product per person in
#' EUR at current prices (“logGDP”). For the purposes of exploratory data analyses
#' it has also been categorised as a binary variable indicating values higher or lower
#' than the median (“Binary GDP”). The employment composition (D = 11) is:
#'
#' * Primary sector (agriculture, hunting, forestry, fishing, mining, quarrying)
#' * Manufacturing
#' * Energy (electricity, gas and water supply)
#' * Construction
#' * Trade repair transport (wholesale and retail trade, repair, transport,
#' storage, communications)
#' * Hotels restaurants
#' * Financial intermediation
#' * Real estate (real estate, renting and business activities)
#' * Educ admin defense soc sec (education, public administration, defence,
#' social security)
#' * Health social work
#' * Other services (other community, social and personal service activities)
#' @name eurostat_employment
"eurostat_employment"

#' Paleocological compositions
#'
#' @description
#' The foraminiferal data set (Aitchison, 1986) is a typical example of
#' paleocological data. It contains compositions of 4 different fossils
#' (Neogloboquadrina atlantica, Neogloboquadrina pachyderma, Globorotalia
#' obesa, and Globigerinoides triloba) at 30 different depths. Due to the
#' rounded zeros present in the data set we will apply some zero replacement
#' techniques to impute these values in advance. After data preprocessing,
#' the analysis that should be undertaken is the association between
#' the composition and the depth.
#' @name foraminiferals
"foraminiferals"

#' Household expenditures
#'
#' @description
#' From Eurostat (the European Union’s statistical information service) the
#' houseexpend data set records the composition on proportions of mean
#' consumption expenditure of households expenditures on 12 domestic year
#' costs in 27 states of the European Union. Some values in the data set are
#' rounded zeros. In addition the data set contains the gross domestic
#' product (GDP05) and (GDP14) in years 2005 and 2014, respectively. An
#' interesting analysis is the potential association between expenditures
#' compositions and GDP. Once a linear regression model is established,
#' predictions can be provided.
#' @name house_expend
"house_expend"

#' Mammal’s milk
#'
#' @description
#' The mammalsmilk data set contains the percentages of five constituents (W: water,
#' P: protein, F: fat, L: lactose, and A: ash) of the milk of 24 mammals. The data are
#' taken from [Har75].
#' @name mammals_milk
"mammals_milk"

#' Milk composition study
#'
#' @description
#' In an attempt to improve the quality of cow milk, milk from each of thirty cows
#' was assessed by dietary composition before and after a strictly controlled dietary
#' and hormonal regime over a period of eight weeks. Although seasonal variations in
#' milk quality might have been regarded as negligible over this period, it was decided
#' to have a control group of thirty cows kept under the same conditions but on a
#' regular established regime. The sixty cows were of course allocated to control and
#' treatment groups at random. The `milk_cows` data set provides the complete set of
#' before and after milk compositions for the sixty cows, showing the protein (pr),
#' milk fat (mf), carbohydrate (ch), calcium (Ca), sodium (Na) and potassium (K)
#' proportions by weight of total dietary content.
#' @name milk_cows
"milk_cows"

#' Concentration of minor elements in carbon ashes
#'
#' @description
#' The montana data set consists of 229 samples of the concentration (in ppm) of
#' minor elements [Cr, Cu, Hg, U, V] in carbon ashes from the Fort Union
#' formation (Montana, USA), side of the Powder River Basin. The formation is
#' mostly Palaeocene in age, and the coal is the result of deposition in
#' conditions ranging from fluvial to lacustrine. All samples were taken from
#' the same seam at different sites over an area of 430 km by 300 km, which
#' implies that on average, the sampling spacing is 24 km. Using the spatial
#' coordinates of the data, a semivariogram analysis was conducted for each
#' chemical element in order to check for a potential spatial dependence
#' structure in the data (not shown here). No spatial dependence patterns
#' were observed for any component, which allowed us to assume an independence
#' of the chemical samples at different locations.
#'
#' The aforementioned chemical components actually represent a fully observed
#' subcomposition of a much larger chemical composition. The five elements are
#' not closed to a constant sum. Note that, as the samples are expressed in
#' parts per million and all concentrations were originally measured, a residual
#' element could be defined to fill up the gap to 10^6.
#' @name montana
"montana"

#' Calc-alkaline and tholeiitic volcanic rocks
#'
#' @description
#'
#' This petrafm data set is formed by 100 classified volcanic rock samples from
#' Ontario (Canada). The three parts are:
#'
#' \eqn{[A: {Na}_2 O + K_2 O; F: Fe O + 0.8998 Fe_2 O_3 ; M: Mg O]}{[A: Na2 O + K2 O; F: FeO + 0.8998 · Fe2 O3 ; M: MgO]}
#'
#' Rocks from the calc-alkaline magma series (25) can be well distinguished from
#' samples from the tholeiitic magma series (75) on an AFM diagram.
#' @name petrafm
"petrafm"

#' Pollen composition in fossils
#'
#' @description
#' The pollen data set is formed by 30 fossil pollen samples from three different
#' locations (recorded in variable group) . The samples were analysed and the 3-part
#' composition [pinus, abies, quercus] was measured.
#' @name pollen
"pollen"

#' Serum proteins
#'
#' @description
#' The `serprot` data set records the percentages of the four serum proteins
#' from the blood samples of 30 patients. Fourteen patients have one
#' disease (1) and sixteen are known to have another different disease (2).
#' The 4-compositions are formed by the proteins [albumin, pre-albumin,
#' globulin A, globulin B].
#' @name serprot
"serprot"

#' A statistician’s time budget
#'
#' @description
#' Time budgets –how a day or a period of work is divided up into different
#' activities have become a popular source of data in psychology and
#' sociology. To illustrate such problems we consider six daily activities
#' undertaken by an academic statistician: teaching (T); consultation (C);
#' administration (A); research (R); other wakeful activities (O); and sleep (S).
#'
#' The `statistician_time` data set records the daily time (in hours) devoted
#' to each activity, recorded on each of 20 days, selected randomly from working days
#' in alternate weeks so as to avoid possible carry-over effects such as a short-sleep
#' day being compensated by make-up sleep on the succeeding day. The six activities
#' may be divided into two categories: 'work' comprising activities T, C, A, and R, and
#' 'leisure', comprising activities O and S. Our analysis may then be directed towards
#' the work pattern consisting of the relative times spent in the four work activities,
#' the leisure pattern, and the division of the day into work time and leisure time.
#' Two obvious questions are as follows. To what extent, if any, do the patterns of
#' work and of leisure depend on the times allocated to these major divisions of the
#' day? Is the ratio of sleep to other wakeful activities dependent on the times spent
#' in the various work activities?
#' @name statistitian_time
"statistitian_time"

#' The waste composition in Catalonia
#'
#' The actual population residing in a municipality of Catalonia is composed by the
#' census count and the so-called floating population (tourists, seasonal visitors,
#' hostel students, short-time employees, and the like). Since actual population
#' combines long and short term residents it is convenient to express it as
#' equivalent full-time residents. Floating population may be positive if the +
#' municipality is receiving more short term residents than it is sending elsewhere,
#' or negative if the opposite holds (expressed as a percentage above –if positive–
#' or below –if negative–  the census count). The waste data set includes this
#' information in the variable floating population. Floating population has a
#' large impact on solid waste generation and thus waste can be used to predict
#' floating population which is a hard to estimate demographic variable. This
#' case study was presented in
#'
#' Tourists and census population do not generate the same volume of waste and have
#' different consumption and recycling patterns (waste composition). The Catalan
#' Statistical Institute (IDESCAT) publishes official floating population data for all
#' municipalities in Catalonia (Spain) above 5000 census habitants. The composition
#' of urban solid waste is classified into D = 5 parts:
#' * x1 : non recyclable (grey waste container in Catalonia),
#' * x2 : glass (bottles and jars of any colour: green waste container),
#' * x3 : light containers (plastic packaging, cans and tetra packs: yellow container),
#' * x4 : paper and cardboard (blue container), and
#' * x5 : biodegradable waste (brown container).
#'
#' @references
#' G. Coenders, J.A.Martín-Fernández and B. Ferrer-Rosell, \emph{When relative and
#' absolute information matter: compositional predictor with a total in generalized
#' linear models}. Statistical Modelling \strong{17}(6) (2017), 494–512.
#' @name waste
"waste"

#' Hotel posts in social media
#'
#' @description
#' The `weibo_hotels` data set aims at comparing the use of Weibo (Facebook
#' equivalent in China) in hospitality e-marketing between small and medium
#' accommodation establishments (private hostels, small hotels) and big and
#' well-established business (such as international hotel chains or large hotels)
#' in China. The 50 latest posts of the Weibo pages of each hotel (n = 10) are
#' content-analyzed and coded regarding the count of posts featuring information
#'  on a 4-part composition [facilities, food, events, promotions]. Hotels were
#'  coded as large “L” or small “S” in the hotel size categorical variable.
#' @name weibo_hotels
"weibo_hotels"

#' @title Chemical Composition of Volcanic Rocks from Kilauea Iki
#' @description This dataset contains the chemical composition of volcanic rocks sampled from the lava lake at Kilauea Iki (Hawaii). The data represents major oxide concentrations in fractional form.
#'
#' @format A data frame with 17 observations and 11 variables:
#' \describe{
#'   \item{SiO2}{Silicon dioxide (fraction)}
#'   \item{TiO2}{Titanium dioxide (fraction)}
#'   \item{Al2O3}{Aluminium oxide (fraction)}
#'   \item{Fe2O3}{Ferric oxide (fraction)}
#'   \item{FeO}{Ferrous oxide (fraction)}
#'   \item{MnO}{Manganese oxide (fraction)}
#'   \item{MgO}{Magnesium oxide (fraction)}
#'   \item{CaO}{Calcium oxide (fraction)}
#'   \item{Na2O}{Sodium oxide (fraction)}
#'   \item{K2O}{Potassium oxide (fraction)}
#'   \item{P2O5}{Phosphorus pentoxide (fraction)}
#' }
#'
#' @details The variability in the oxide concentrations is attributed to magnesic olivine fractionation, starting from a single magmatic mass as suggested by Richter & Moore (1966).
#'
#' @source Richter, D.H., & Moore, J.G. (1966). Petrology of Kilauea Iki lava lake, Hawaii. *Geological Survey Professional Paper* 537-B.
#' @name kilauea_iki
"kilauea_iki"
