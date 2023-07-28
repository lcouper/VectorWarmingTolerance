### Vector Thermal Safety Margin Estimation 

# This script is used to estimate vector thermal safety margins using ERA5 climate data as input
# Code overview: 
# 0. Load libraries
# 1. Pull in species occurrence data and subset to given latitudinal band
# 2. Calculate thermal vulnerability indices with and w/o behavior: 
# TSMs, total time in thermal danger, longest streak of days in thermal danger, 
# longest streak of hours in thermal danger
# Output data to csv

#### 0. Load libraries and data ####

setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")
library(mcera5)
library(ecmwfr)
library(NicheMapR)
library(dplyr)
library(ncdf4)
library(curl)
library(keyring)
library(abind)
library(lubridate)
library(tidync)
library(microclima) # to install: remotes::install_github("ilyamaclean/microclima")
library(NicheMapR) # to install: remotes::install_github("mrke/NicheMapR") or devtools::install_github('mrke/NicheMapR')


#### 1. Pull in occurrence data and subset to specific lat/lon #####
##### 1a. Aedes albopictus #####

aedes = read.csv("VectorOccurrence/aegypti_albopictus.csv")

# subset to just published, point records of Aedes albopicus 
aealbo = subset(aedes, VECTOR=="Aedes albopictus" & LOCATION_TYPE=="point" &
                  SOURCE_TYPE=="published")

# remove any duplicates (occurences from same lat, lon)
aealbo <- aealbo %>% dplyr::distinct(Latitude, Longitude, .keep_all = TRUE)

# latitidunal range = -37 to 49

# Subset to records within latitudinal band: 
aealbo = subset(aealbo, Latitude>=-38 & Latitude<=-35)

# sub-sample if hundreds of records
# aealbo = sample_n(aealbo, 25)


##### 1b. Aedes aegypti ######

# subset to just published, point records of Aedes aegypti or albopicus 
aeaegypti = subset(aedes, VECTOR=="Aedes aegypti" & LOCATION_TYPE=="point" &
                     SOURCE_TYPE=="published")

# remove any duplicates (occurences from same lat, lon)
aeaegypti<- aeaegypti %>% dplyr::distinct(Latitude, Longitude, .keep_all = TRUE)

# latitudinal range: -39 to -25

# Subset to records within a given latitudinal band: 
aeaegypti = subset(aeaegypti, Latitude>=-39 & Latitude<=-25)

# sub-sample if hundreds of records
# aeaegypti = sample_n(aeaegypti, 35)

##### 1c. Aedes camptorhynchus ####

AeCamp = read.delim("VectorOccurrence/Aedes camptorhynchus.csv")
aec = AeCamp[!is.na(AeCamp$decimalLatitude),]
aec = aec[!is.na(aec$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
aec = aec %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = -44 to -17
aec = subset(aec, decimalLatitude>=-25 & decimalLatitude<= -15)
# subset
aec = sample_n(aec, 30)


##### 1d. Aedes triseriatus ####

AeTri = read.csv("VectorOccurrence/AedesTriseriatus.csv")
aet = AeTri[!is.na(AeTri$decimalLatitude),]
aet = aet[!is.na(aet$decimalLongitude),]

# some dubious records in Cambodia/Loas region - filter these
aet = aet[aet$decimalLongitude < 90,]

# remove any duplicates (occurrences from same lat, lon)
aet = aet %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = 17 to 50
aet = subset(aet, decimalLatitude>=15 & decimalLatitude<= 25)
aet = sample_n(aet, 30) # sub-sample if hundreds of records

##### 1e. Aedes vexans ####

AeVex = read.csv("VectorOccurrence/AedesVexans_Full.csv")
aev = AeVex[!is.na(AeVex$decimalLatitude),]
aev = aev[!is.na(aev$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
aev = aev %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = -35 to 67
aev = subset(aev, decimalLatitude>=48 & decimalLatitude<=67)
aev = aev[aev$decimalLongitude >= -124,]

# To ensure that records outside of US & CA are also included
aev1 = aev[aev$countryCode == "US" | aev$countryCode == "CA",]
aev2 = aev[aev$countryCode != "US",]
aev2 = aev2[aev2$countryCode != "CA",]

aev = rbind.data.frame(sample_n(aev1, 5), sample_n(aev2, 5)) # sub-sample if hundreds of records

##### 1f. Anopheles gambiae ######

agambiae = read.csv("VectorOccurrence/Africa Vectors database_1898-2016.csv", header = T)
agambiae = agambiae[agambiae$An.gambiae.ss == "Y",]

# remove any duplicates (occurrences from same lat, lon)
agambi <- agambiae %>% dplyr::distinct(Lat, Long, .keep_all = TRUE)

agambi = agambi[!is.na(agambi$Lat),]
agambi = agambi[!is.na(agambi$Long),]

# range: -27 to 21

# Subset to records within latitudinal band: 
agambi = subset(agambi, Lat>=-15 & Lat<=-5)

# to ensure better spatial coverage
agambi <- agambi %>% group_by(Country) %>% sample_n(2)

#agambi = sample_n(agambi,25)


##### 1g. Anopheles stephensi ######

asteph = read.csv("VectorOccurrence/AnophelesStephensi_Additional.csv", header = T)

# remove any duplicates (occurrences from same lat, lon)
asteph <- asteph %>% dplyr::distinct(latitude, longitude, .keep_all = TRUE)

asteph = asteph[!is.na(asteph$latitude),]
asteph = asteph[!is.na(asteph$longitude),]

# lat range: 5 to 35

# Subset to records within latitudinal band: 
asteph = subset(asteph, latitude >=25 & latitude <=35)

asteph <- asteph %>% group_by(country) %>% sample_n(1)
#asteph = sample_n(asteph, 30)

##### 1h. Culex quinquefasciatus ######

cxq = read.csv("VectorOccurrence/Culex quinquefasciatus.csv")
cxq = cxq[!is.na(cxq$decimalLatitude),]
cxq = cxq[!is.na(cxq$decimalLongitude),]

# Latitudinal range: -44 to 47

# remove any duplicates (occurrences from same lat, lon)
cxq = cxq %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Subset to records within latitudinal band: 
cxq = subset(cxq, decimalLatitude>=35)  
cxq1 = subset(cxq, decimalLongitude <= -10)
cxq2 = subset(cxq, decimalLongitude > -5)

cxq = rbind.data.frame(sample_n(cxq1, 24), sample_n(cxq2,6))

##### 1i. Culex pipiens ######

cxp = read.csv("VectorOccurrence/CulexPipiens.csv")
cxp = cxp[!is.na(cxp$decimalLatitude),]
cxp = cxp[!is.na(cxp$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxp = cxp %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: -39 to 67
cxp = subset(cxp, decimalLatitude>=48 & decimalLatitude<=67)

# take random sample of 20-30 rows
#cxp = sample_n(cxp, 35)

##### 1j. Culex tarsalis ######

cxt = read.delim("VectorOccurrence/CulexTarsalis.csv")
cxt = cxt[!is.na(cxt$decimalLatitude),]
cxt = cxt[!is.na(cxt$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxt = cxt %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: 18 to 65
cxt = subset(cxt, decimalLatitude>=50 & decimalLatitude<= 67)
# take random sample of 20-30 rows
#cxt = sample_n(cxt, 10)

##### 1k. Culex annulirostris ######

cxa = read.csv("VectorOccurrence/CulexAnnulirostris.csv")
cxa = cxa[!is.na(cxa$decimalLatitude),]
cxa = cxa[!is.na(cxa$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxa = cxa %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: -39 to 16
cxa = subset(cxa, decimalLatitude>=5 & decimalLatitude<= 16)
# take random sample of 20-30 rows
#cxa = sample_n(cxa, 9)

##### 1l. Culex theileri ######

cxth = read.csv("VectorOccurrence/CulexTheileri.csv")
cxth = cxth[!is.na(cxth$decimalLatitude),]
cxth = cxth[!is.na(cxth$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxth = cxth %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: -35 to 43
cxth = subset(cxth, decimalLatitude>=35 & decimalLatitude<=45)
# take random sample of 20-30 rows
#cxth = sample_n(cxth, 9)


##### 1i. Aedes camptorhynchus ####

AeCamp = read.delim("VectorOccurrence/Aedes camptorhynchus.csv")
aec = AeCamp[!is.na(AeCamp$decimalLatitude),]
aec = aec[!is.na(aec$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
aec = aec %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = -44 to -17
aec = subset(aec, decimalLatitude>=-25 & decimalLatitude<= -15)
# subset
aec = sample_n(aec, 30)


##### 1j. Aedes triseriatus ####

AeTri = read.csv("VectorOccurrence/AedesTriseriatus.csv")
aet = AeTri[!is.na(AeTri$decimalLatitude),]
aet = aet[!is.na(aet$decimalLongitude),]

# some dubious records in Cambodia/Loas region - filter these
aet = aet[aet$decimalLongitude < 90,]

# remove any duplicates (occurrences from same lat, lon)
aet = aet %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = 17 to 50
aet = subset(aet, decimalLatitude>=15 & decimalLatitude<= 25)
aet = sample_n(aet, 30) # sub-sample if hundreds of records

##### 1k. Aedes vexans ####

AeVex = read.csv("VectorOccurrence/AedesVexans_Full.csv")
aev = AeVex[!is.na(AeVex$decimalLatitude),]
aev = aev[!is.na(aev$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
aev = aev %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# lat range = -35 to 67
aev = subset(aev, decimalLatitude>=48 & decimalLatitude<=67)
aev1 = aev[aev$countryCode == "US" | aev$countryCode == "CA",]
aev2 = aev[aev$countryCode != "US",]
aev2 = aev2[aev2$countryCode != "CA",]

aev = rbind.data.frame(sample_n(aev1, 5), sample_n(aev2, 5)) # sub-sample if hundreds of records

#### 2. Calculate thermal vulnerability for given species at each collection location #####

setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")

df1 = as.data.frame(matrix(nrow = nrow(asteph), ncol = 34))
colnames(df1) = c("Species.Pop", "HighestTrait", "lat", "lon", "maxtemp_point",
                  "maxtemp_lower", "maxtemp_upper",
                  "maxtemp_NoB_point", "maxtemp_NoB_lower", "maxtemp_NoB_upper",
                  "tolerance_point", "tolerance_lower", "tolerance_upper", 
                  "tolerance_NoB_point", "tolerance_NoB_lower", "tolerance_NoB_upper", 
                  "numhours_point", "numhours_lower", "numhours_upper", 
                  "numhours_NoB_point", "numhours_NoB_lower", "numhours_NoB_upper", 
                  "streak_point", "streak_lower", "streak_upper",
                  "streak_NoB_point", "streak_NoB_lower", "streak_NoB_upper",
                  "streak_days_point", "streak_days_lower", "streak_days_upper",
                  "streak_days_nob_point", "streak_days_nob_lower", "streak_days_nob_upper")

for (i in 1:nrow(asteph)){
 lat = asteph$latitude[i]
  lon = asteph$longitude[i]
  species = "Anopheles_stephensi"
  
  loc = c(lon, lat)
  dstart <- "01/01/2017"
  dfinish <- "31/12/2017"
  micro<-micro_era5(loc = loc, spatial = 'era5', Usrhyt = 1,
                    dstart = dstart, dfinish = dfinish,
                    minshade = 0, maxshade = 100)
  metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
  shadeout = as.data.frame(micro$shadmet) # above ground microclimatic conditions, full shade
  tzone<-paste("Etc/GMT+",0,sep="") # append dates
  dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
  met <- cbind(dates,metout)
  shade <- cbind(dates,shadeout)
  
  # Ectotherm model
  limitsdata = read.csv("VectorWarmingTolerance_SpeciesThermalLimits.csv")
  speciesdata = limitsdata[limitsdata$Species == species,]
  HighestTrait = unique(speciesdata$HighestTrait)
  traitdata = speciesdata[speciesdata$Trait == HighestTrait,]  
  
  CTmax_point = traitdata$CTmax_point
  CTmax_lower = traitdata$CTmax_lower
  CTmax_upper = traitdata$CTmax_upper
  
  Topt_point = traitdata$Topt_point
  Topt_lower = traitdata$Topt_lower
  Topt_upper = traitdata$Topt_upper
  
  CTmin_point = traitdata$CTmin_point
  CTmin_lower = traitdata$CTmin_lower
  CTmin_upper = traitdata$CTmin_upper
  
  MosqTempPoint = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                            T_pref = Topt_point, CT_max = CTmax_point, CT_min = CTmin_point, 
                            diurn = 0, nocturn = 1, crepus = 1, 
                            shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempPoint_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_point, CT_max = CTmax_point, CT_min = CTmin_point, 
                                diurn = 0, nocturn = 1, crepus = 1, 
                                shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempLower= ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                           T_pref = Topt_lower, CT_max = CTmax_lower, CT_min = CTmin_lower, 
                           diurn = 0, nocturn = 1, crepus = 1, 
                           shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempLower_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_lower, CT_max = CTmax_lower, CT_min = CTmin_lower, 
                                diurn = 0, nocturn = 1, crepus = 1, 
                                shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempUpper = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                            T_pref = Topt_upper, CT_max = CTmax_upper, CT_min = CTmin_upper, 
                            diurn = 0, nocturn = 1, crepus = 1, 
                            shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempUpper_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_upper, CT_max = CTmax_upper, CT_min = CTmin_upper, 
                                diurn = 0, nocturn = 1, crepus = 1, 
                                shade_seek = 0, burrow = 0, climb = 0)
  
  point = as.data.frame(MosqTempPoint$environ)
  lower = as.data.frame(MosqTempLower$environ)
  upper = as.data.frame(MosqTempUpper$environ)
  
  point_noB = as.data.frame(MosqTempPoint_NoB$environ)
  lower_noB = as.data.frame(MosqTempLower_NoB$environ)
  upper_noB = as.data.frame(MosqTempUpper_NoB$environ)
  
  ModeledTemps = cbind.data.frame(point, lower$TC, upper$TC,  point_noB$TC, lower_noB$TC, upper_noB$TC, met$TALOC, shade$TALOC)
  colnames(ModeledTemps)[c(5,29:35)] = c("Mosq_Temp_Point", "Mosq_Temp_Lower", "Mosq_Temp_Upper", 
                                         "Mosq_Temp_NoB_Point", "Mosq_Temp_NoB_Lower", "Mosq_Temp_NoB_Upper", 
                                         "Sun_Temp", "Shade_Temp")
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
  # Number hours where predicted body temp > CTmax
  numhours_point = sum(ModeledTemps$Mosq_Temp_Point > CTmax_point) 
  numhours_lower = sum(ModeledTemps$Mosq_Temp_Lower > CTmax_lower) 
  numhours_upper = sum(ModeledTemps$Mosq_Temp_Upper > CTmax_upper) 
  
  numhours_NoB_point = sum(ModeledTemps$Mosq_Temp_NoB_Point > CTmax_point) 
  numhours_NoB_lower = sum(ModeledTemps$Mosq_Temp_NoB_Lower > CTmax_lower) 
  numhours_NoB_upper = sum(ModeledTemps$Mosq_Temp_NoB_Upper > CTmax_upper) 
  
  # Number of consecutive days in thermal danger
  ModeledTemps$danger_point =  as.numeric((CTmax_point - ModeledTemps$Mosq_Temp_Point) <0)
  dangerdays_point = aggregate(ModeledTemps$danger_point, list(ModeledTemps$DAY), sum)$x
  runs_danger_point = rle(dangerdays_point > 0)
  streak_dangerdays_point = as.numeric(tapply(runs_danger_point$lengths, runs_danger_point$values, max))[2]
  ModeledTemps$danger_lower =  as.numeric((CTmax_lower - ModeledTemps$Mosq_Temp_Lower) <0)
  dangerdays_lower = aggregate(ModeledTemps$danger_lower, list(ModeledTemps$DAY), sum)$x
  runs_danger_lower = rle(dangerdays_lower > 0)
  streak_dangerdays_lower = as.numeric(tapply(runs_danger_lower$lengths, runs_danger_lower$values, max))[2]
  ModeledTemps$danger_upper =  as.numeric((CTmax_upper - ModeledTemps$Mosq_Temp_Upper) <0)
  dangerdays_upper = aggregate(ModeledTemps$danger_upper, list(ModeledTemps$DAY), sum)$x
  runs_danger_upper = rle(dangerdays_upper > 0)
  streak_dangerdays_upper = as.numeric(tapply(runs_danger_upper$lengths, runs_danger_upper$values, max))[2]
  
  # Number of consecutive days in thermal danger without behavior
  ModeledTemps$danger_nob_point =  as.numeric((CTmax_point - ModeledTemps$Mosq_Temp_NoB_Point) <0)
  dangerdays_nob_point = aggregate(ModeledTemps$danger_nob_point, list(ModeledTemps$DAY), sum)$x
  runs_danger_nob_point = rle(dangerdays_nob_point > 0)
  streak_dangerdays_nob_point = as.numeric(tapply(runs_danger_nob_point$lengths, runs_danger_nob_point$values, max))[2]
  ModeledTemps$danger_nob_lower =  as.numeric((CTmax_lower - ModeledTemps$Mosq_Temp_NoB_Lower) <0)
  dangerdays_nob_lower = aggregate(ModeledTemps$danger_nob_lower, list(ModeledTemps$DAY), sum)$x
  runs_danger_nob_lower = rle(dangerdays_nob_lower > 0)
  streak_dangerdays_nob_lower = as.numeric(tapply(runs_danger_nob_lower$lengths, runs_danger_nob_lower$values, max))[2]
  ModeledTemps$danger_nob_upper =  as.numeric((CTmax_upper - ModeledTemps$Mosq_Temp_NoB_Upper) <0)
  dangerdays_nob_upper = aggregate(ModeledTemps$danger_nob_upper, list(ModeledTemps$DAY), sum)$x
  runs_danger_nob_upper = rle(dangerdays_nob_upper > 0)
  streak_dangerdays_nob_upper = as.numeric(tapply(runs_danger_nob_upper$lengths, runs_danger_nob_upper$values, max))[2]
  
  
  # Length of longest streak above CTmax (in hours)
  runs_point = rle(ModeledTemps$Mosq_Temp_Point > CTmax_point)
  runs_lower = rle(ModeledTemps$Mosq_Temp_Lower > CTmax_lower)
  runs_upper = rle(ModeledTemps$Mosq_Temp_Upper > CTmax_upper)
  streak_point = as.numeric(tapply(runs_point$lengths, runs_point$values, max))[2]
  streak_lower = as.numeric(tapply(runs_lower$lengths, runs_lower$values, max))[2]
  streak_upper = as.numeric(tapply(runs_upper$lengths, runs_upper$values, max))[2]
  
  runs_NoB_point = rle(ModeledTemps$Mosq_Temp_NoB_Point > CTmax_point)
  runs_NoB_lower = rle(ModeledTemps$Mosq_Temp_NoB_Lower > CTmax_lower)
  runs_NoB_upper = rle(ModeledTemps$Mosq_Temp_NoB_Upper > CTmax_upper)
  streak_NoB_point = as.numeric(tapply(runs_NoB_point$lengths, runs_NoB_point$values, max))[2]
  streak_NoB_lower = as.numeric(tapply(runs_NoB_lower$lengths, runs_NoB_lower$values, max))[2]
  streak_NoB_upper = as.numeric(tapply(runs_NoB_upper$lengths, runs_NoB_upper$values, max))[2]
  
  df1[i,] = c(species, HighestTrait, lat, lon, maxtemp_point, maxtemp_lower, maxtemp_upper, 
              maxtemp_NoB_point, maxtemp_NoB_lower, maxtemp_NoB_upper,
              diff_point, diff_lower, diff_upper,
              diff_NoB_point, diff_NoB_lower, diff_NoB_upper,
              numhours_point, numhours_lower, numhours_upper,
              numhours_NoB_point, numhours_NoB_lower, numhours_NoB_upper,
              streak_point, streak_lower, streak_upper,
              streak_NoB_point, streak_NoB_lower, streak_NoB_upper,
              streak_dangerdays_point, streak_dangerdays_lower, streak_dangerdays_upper,
              streak_dangerdays_nob_point, streak_dangerdays_nob_lower, streak_dangerdays_nob_upper)
}

write.csv(df1, "~/Downloads/AnSteph_TSM_25to35.csv")

