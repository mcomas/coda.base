#' Catalan Parliament election results in 2017 by region
#'
#' @description
#' The `parliament2017` data set contains the results of the 2017 Catalan
#' Parliament election aggregated by region.
#'
#' @format A data frame with 42 rows and 9 variables:
#' \describe{
#'   \item{com}{Region}
#'   \item{cs}{Votes for the Ciutadans party}
#'   \item{jxcat}{Votes for the Junts per Catalunya party}
#'   \item{erc}{Votes for the Esquerra Republicana de Catalunya party}
#'   \item{psc}{Votes for the Partit Socialista de Catalunya party}
#'   \item{catsp}{Votes for the Catalunya Sí que es Pot party}
#'   \item{cup}{Votes for the Candidatura d'Unitat Popular party}
#'   \item{pp}{Votes for the Partit Popular party}
#'   \item{other}{Votes for other parties}
#' }
#'
#' @source Idescat, statistics on Catalan Parliament elections.
#' @name parliament2017
"parliament2017"

#' Chemical compositions of Romano-British pottery
#'
#' @description
#' The `pottery` data set contains the chemical composition of 45 specimens of
#' Romano-British pottery. The measurements were obtained by atomic absorption
#' spectrophotometry and include nine oxides:
#' Al2O3, Fe2O3, MgO, CaO, Na2O, K2O, TiO2, MnO, and BaO.
#'
#' The specimens come from five different kiln sites.
#'
#' @name pottery
"pottery"

#' Arctic lake sediments at different depths
#'
#' @description
#' The `arctic_lake` data set records the three-part composition
#' \eqn{[sand, silt, clay]} of 39 sediment samples collected at different water
#' depths in an Arctic lake.
#'
#' @name arctic_lake
"arctic_lake"

#' Household budget patterns
#'
#' @description
#' In a sample survey of single persons living alone in rented accommodation,
#' twenty men and twenty women were randomly selected and asked to record their
#' expenditure over one month in the following four mutually exclusive and
#' exhaustive commodity groups:
#'
#' \itemize{
#'   \item `Hous`: housing, including fuel and light,
#'   \item `Food`: foodstuffs, including alcohol and tobacco,
#'   \item `Serv`: services, including transport and vehicles,
#'   \item `Other`: other goods, including clothing, footwear, and durable goods.
#' }
#'
#' @name household_budget
"household_budget"

#' Food consumption in European countries
#'
#' @description
#' The `alimentation` data set contains the percentage composition of food
#' consumption in 25 European countries during the 1980s. The food categories
#' are:
#'
#' \itemize{
#'   \item `RM`: red meat (pork, veal, beef),
#'   \item `WM`: white meat (chicken),
#'   \item `E`: eggs,
#'   \item `M`: milk,
#'   \item `F`: fish,
#'   \item `C`: cereals,
#'   \item `S`: starch (potatoes),
#'   \item `N`: nuts,
#'   \item `FV`: fruits and vegetables.
#' }
#'
#' The data set also contains categorical variables indicating whether the
#' country belongs to the North or South/Mediterranean group, and whether it is
#' an Eastern or Western European country.
#'
#' @name alimentation
"alimentation"

#' The MN blood system
#'
#' @description
#' In humans, the main blood group systems are the ABO system, the Rh system,
#' and the MN system. The MN blood system is related to proteins of the red
#' blood cell plasma membrane. Its inheritance pattern is autosomal with
#' codominance, meaning that the heterozygous phenotype is distinct from both
#' homozygous phenotypes.
#'
#' The three phenotypes are M, N, and MN. Their frequencies vary across
#' populations. Under the Hardy-Weinberg principle, allele and genotype
#' frequencies remain constant across generations in the absence of evolutionary
#' forces, implying that
#'
#' \eqn{\frac{x_{MM} x_{NN}}{x_{MN}^2} = \frac{1}{4}}{
#' (x_MM x_NN) / x_MN^2 = 1/4
#' }
#'
#' where \eqn{x_{MM}} and \eqn{x_{NN}} are the genotype frequencies of the
#' homozygotes and \eqn{x_{MN}} is the genotype frequency of heterozygotes.
#'
#' @name blood_mn
"blood_mn"

#' Physical activity and body mass index
#'
#' @description
#' The `bmi_activity` data set records the proportion of daily time spent in
#' sleep (`sleep`), sedentary behaviour (`sedent`), light physical activity
#' (`Lpa`), moderate physical activity (`Mpa`), and vigorous physical activity
#' (`Vpa`) for 393 children. The standardized body mass index (`zBMI`) of each
#' child is also included.
#'
#' This data set was used in the example of Dumuid et al. (2019) to examine the
#' expected differences in zBMI associated with reallocations of daily time
#' between sleep, sedentary behaviour, and physical activity. Because the
#' original data are confidential, `bmi_activity` contains simulated data that
#' mimic the main features of the original study.
#'
#' @references
#' Dumuid, D., Pedisic, Z., Stanford, T. E., Martín-Fernández, J. A., Hron, K.,
#' Maher, C., Lewis, L. K., & Olds, T. S. (2019).
#' \emph{The Compositional Isotemporal Substitution Model: a Method for
#' Estimating Changes in a Health Outcome for Reallocation of Time between
#' Sleep, Sedentary Behaviour, and Physical Activity}.
#' Statistical Methods in Medical Research, \strong{28}(3), 846--857.
#'
#' @name bmi_activity
"bmi_activity"

#' Employment distribution in EUROSTAT countries
#'
#' @description
#' According to the three-sector theory, employment shifts from the primary
#' sector (raw material extraction), to the secondary sector (industry, energy,
#' and construction), and then to the tertiary sector (services) as economies
#' develop. The `eurostat_employment` data set contains EUROSTAT data on
#' employment, aggregated for both sexes and all ages, distributed by economic
#' activity in 2008 for 29 EUROSTAT member countries.
#'
#' A related variable is the logarithm of gross domestic product per person in
#' EUR at current prices (`logGDP`). For exploratory purposes, it is also
#' categorised as a binary variable indicating values above or below the median
#' (`Binary GDP`).
#'
#' The employment composition has 11 parts:
#' \itemize{
#'   \item Primary sector
#'   \item Manufacturing
#'   \item Energy
#'   \item Construction
#'   \item Trade repair transport
#'   \item Hotels restaurants
#'   \item Financial intermediation
#'   \item Real estate
#'   \item Educ admin defense soc sec
#'   \item Health social work
#'   \item Other services
#' }
#'
#' @name eurostat_employment
"eurostat_employment"

#' Paleocological compositions
#'
#' @description
#' The `foraminiferals` data set (Aitchison, 1986) is a classical example of
#' paleocological compositional data. It contains the composition of four fossil
#' types (Neogloboquadrina atlantica, Neogloboquadrina pachyderma,
#' Globorotalia obesa, and Globigerinoides triloba) at 30 different depths.
#'
#' Because the data contain rounded zeros, zero-replacement techniques are
#' typically required before analysis. A natural goal is then to study the
#' association between fossil composition and depth.
#'
#' @name foraminiferals
"foraminiferals"

#' Household expenditures
#'
#' @description
#' The `house_expend` data set, obtained from Eurostat, records the composition
#' of mean household consumption expenditure across 12 expenditure categories in
#' 27 European Union countries. Some values are rounded zeros.
#'
#' In addition, the data set contains gross domestic product values for 2005
#' (`GDP05`) and 2014 (`GDP14`). A relevant analysis is the relationship
#' between expenditure compositions and GDP.
#'
#' @name house_expend
"house_expend"

#' Mammals' milk
#'
#' @description
#' The `mammals_milk` data set contains the percentages of five constituents of
#' the milk of 24 mammals:
#' \eqn{[W, P, F, L, A]}{[W, P, F, L, A]},
#' where `W` is water, `P` is protein, `F` is fat, `L` is lactose, and `A` is
#' ash.
#'
#' @name mammals_milk
"mammals_milk"

#' Milk composition study
#'
#' @description
#' In an attempt to improve the quality of cow milk, milk from thirty cows was
#' assessed before and after a controlled dietary and hormonal regime over eight
#' weeks. A control group of thirty cows kept under the usual regime was also
#' included.
#'
#' The `milk_cows` data set provides the complete before/after milk composition
#' data for the sixty cows, with the proportions of protein (`pr`), milk fat
#' (`mf`), carbohydrate (`ch`), calcium (`Ca`), sodium (`Na`), and potassium
#' (`K`).
#'
#' @name milk_cows
"milk_cows"

#' Concentration of minor elements in coal ashes
#'
#' @description
#' The `montana` data set contains 229 samples of the concentration (in ppm) of
#' five minor elements \eqn{[Cr, Cu, Hg, U, V]} in coal ashes from the Fort
#' Union formation (Montana, USA), in the Powder River Basin.
#'
#' The five measured elements form a fully observed subcomposition of a much
#' larger chemical composition. Since the data are given in parts per million and
#' all concentrations were measured, a residual component could in principle be
#' added to close the compositions to \eqn{10^6}.
#'
#' @name montana
"montana"

#' Calc-alkaline and tholeiitic volcanic rocks
#'
#' @description
#' The `petrafm` data set contains 100 classified volcanic rock samples from
#' Ontario (Canada). The three-part composition is
#'
#' \eqn{[A: Na_2O + K_2O;\ F: FeO + 0.8998\,Fe_2O_3;\ M: MgO]}{
#' [A: Na2O + K2O; F: FeO + 0.8998 Fe2O3; M: MgO]
#' }
#'
#' Rocks from the calc-alkaline magma series (25 samples) can be distinguished
#' from those of the tholeiitic magma series (75 samples) using an AFM diagram.
#'
#' @name petrafm
"petrafm"

#' Pollen composition in fossils
#'
#' @description
#' The `pollen` data set contains 30 fossil pollen samples from three different
#' locations (recorded in variable `group`). The measured composition is the
#' three-part composition \eqn{[pinus, abies, quercus]}.
#'
#' @name pollen
"pollen"

#' Serum proteins
#'
#' @description
#' The `serprot` data set records the percentages of four serum proteins from
#' blood samples of 30 patients. Fourteen patients have one disease and sixteen
#' have another.
#'
#' The four-part compositions are formed by
#' \eqn{[albumin, pre\text{-}albumin, globulin\ A, globulin\ B]}{
#' [albumin, pre-albumin, globulin A, globulin B]
#' }.
#'
#' @name serprot
"serprot"

#' A statistician's time budget
#'
#' @description
#' The `statistician_time` data set records the daily time budget of an academic
#' statistician across 20 working days. The six activities are teaching (`T`),
#' consultation (`C`), administration (`A`), research (`R`), other wakeful
#' activities (`O`), and sleep (`S`).
#'
#' These activities may also be grouped into work (`T`, `C`, `A`, `R`) and
#' leisure (`O`, `S`). The data allow investigation of the relationship between
#' detailed time-allocation patterns and the broader division between work and
#' leisure.
#'
#' @name statistician_time
"statistician_time"

#' Urban waste composition in Catalonia
#'
#' @description
#' The `waste` data set studies the relationship between waste composition and
#' floating population in Catalonia. The actual population of a municipality
#' combines census population and floating population (tourists, seasonal
#' visitors, temporary workers, and similar short-term residents), expressed as
#' equivalent full-time residents.
#'
#' The composition of urban solid waste is classified into five parts:
#' \itemize{
#'   \item `x1`: non-recyclable waste,
#'   \item `x2`: glass,
#'   \item `x3`: light containers,
#'   \item `x4`: paper and cardboard,
#'   \item `x5`: biodegradable waste.
#' }
#'
#' Waste generation and composition are influenced by floating population, which
#' makes waste composition a useful predictor of this difficult-to-measure
#' demographic quantity.
#'
#' @references
#' Coenders, G., Martín-Fernández, J. A., & Ferrer-Rosell, B. (2017).
#' \emph{When relative and absolute information matter: compositional predictor
#' with a total in generalized linear models}. Statistical Modelling,
#' \strong{17}(6), 494--512.
#'
#' @name waste
"waste"

#' Hotel posts in social media
#'
#' @description
#' The `weibo_hotels` data set compares the use of Weibo (the Chinese equivalent
#' of Facebook) in hospitality e-marketing between small and medium
#' establishments and larger hotel businesses in China.
#'
#' The 50 latest posts from the Weibo page of each hotel (\eqn{n = 10}) were
#' content-analysed and coded into a four-part composition:
#' \eqn{[facilities, food, events, promotions]}.
#' Hotels were also classified by size as large (`L`) or small (`S`).
#'
#' @name weibo_hotels
"weibo_hotels"

#' Chemical composition of volcanic rocks from Kilauea Iki
#'
#' @description
#' The `kilauea_iki` data set contains the chemical composition of volcanic
#' rocks sampled from the lava lake at Kilauea Iki (Hawaii). The data represent
#' major oxide concentrations in fractional form.
#'
#' @format A data frame with 17 observations and 11 variables:
#' \describe{
#'   \item{SiO2}{Silicon dioxide}
#'   \item{TiO2}{Titanium dioxide}
#'   \item{Al2O3}{Aluminium oxide}
#'   \item{Fe2O3}{Ferric oxide}
#'   \item{FeO}{Ferrous oxide}
#'   \item{MnO}{Manganese oxide}
#'   \item{MgO}{Magnesium oxide}
#'   \item{CaO}{Calcium oxide}
#'   \item{Na2O}{Sodium oxide}
#'   \item{K2O}{Potassium oxide}
#'   \item{P2O5}{Phosphorus pentoxide}
#' }
#'
#' @details
#' The variability in oxide concentrations is attributed to magnesian olivine
#' fractionation from a single magmatic mass, as suggested by Richter and Moore
#' (1966).
#'
#' @source
#' Richter, D. H., & Moore, J. G. (1966). \emph{Petrology of Kilauea Iki lava
#' lake, Hawaii}. Geological Survey Professional Paper 537-B.
#'
#' @name kilauea_iki
"kilauea_iki"
