#############################################################
########### CODE FOR VECTOR TSM ESTIMATION ##################
########### WITH DROUGHT MASK INCORPORATED ##################
########### WRITTEN BY LISA COUPER ##########################

# This script is used to estimate vector thermal safety margins using ERA5 climate data as input
# Code overview: 
# 0. Load libraries
# 1. Pull in species occurrence data and subset to given latitudinal band
# 2. Calculate thermal vulnerability indices with and w/o behavior: 
# TSMs, total time in thermal danger, longest streak of days in thermal danger, 
# longest streak of hours in thermal danger
# Output data to csv

#### 0. Load libraries and data ####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles")
library(mcera5)
library(ecmwfr)
library(dplyr)
library(ncdf4)
library(curl)
library(keyring)
library(abind)
library(data.table)
library(lubridate)
library(tidync)
library(maps)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(mapdata)
library(microclima) # to install: remotes::install_github("ilyamaclean/microclima")
library(NicheMapR) # to install: remotes::install_github("mrke/NicheMapR") or devtools::install_github('mrke/NicheMapR')


#### 1. Pull in occurrence data and subset to specific lat/lon #####
##### 1a. Aedes albopictus #####

aealbo = fread("Vector Occurrence Data/AeAlbopictus_SubSampledPoints.csv")
range(aealbo$Latitude)
aealbo = subset(aealbo, Latitude>=25 & Latitude<=35)

##### 1b. Aedes aegypti ######

aeaegypti = fread("Vector Occurrence Data/AeAegypti_SubSampledPoints.csv")
range(aeaegypti$Latitude)
aeaegypti = subset(aeaegypti, Latitude>=45) # & Latitude<=45)

##### 1c. Aedes triseriatus ####

AeTri = fread("Vector Occurrence Data/AeTri_SubSampledPoints.csv")
aet = subset(AeTri, Latitude>=15 & Latitude<=25)

##### 1d. Anopheles gambiae ######

agambi = fread("Vector Occurrence Data/AnGambi_SubSampledPoints.csv")
range(agambi$Latitude, na.rm = T)
agambi = subset(agambi, Latitude>=25 & Latitude<=35)

##### 1e. Anopheles stephensi ######

asteph = fread("Vector Occurrence Data/AnSteph_SubSampledPoints.csv")
range(asteph$Latitude)
asteph = subset(asteph, Latitude>=35 & Latitude<=45)

##### 1f. Culex quinquefasciatus ######

cxq = fread("Vector Occurrence Data/CxQuinque_SubSampledPoints.csv")
cxq = subset(cxq, Latitude>=-44 & Latitude<=-35)

##### 1g. Culex pipiens ######

cxp = fread("Vector Occurrence Data/CxPipiens_SubSampledPoints.csv")
cxp = subset(cxp, Latitude>=-44 & Latitude<=-35)

##### 1h. Culex tarsalis ######

cxt = fread("Vector Occurrence Data/CxTarsalis_SubSampledPoints.csv")
cxt = subset(cxt, Latitude>=5 & Latitude<=15)

##### 1i. Culex annulirostris ######

cxa = fread("Vector Occurrence Data/CxAnnul_SubSampledPoints.csv")
range(cxa$Latitude)
cxa = subset(cxa, Latitude>=15 & Latitude<=25)


#### 2. Calculate thermal vulnerability for given sprecies at each collection location #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles")

df = as.data.frame(matrix(nrow = nrow(cxq), ncol = 71))
colnames(df) = c("Species.Pop", "HighestTrait", "lat", "lon", "maxtemp_point",
                 "maxtemp_lower", "maxtemp_upper",
                 "maxtemp_NoB_point", "maxtemp_NoB_lower", "maxtemp_NoB_upper",
                 "tolerance_point", "tolerance_lower", "tolerance_upper", 
                 "tolerance_NoB_point", "tolerance_NoB_lower", "tolerance_NoB_upper", 
                 "tolerance_Tpref_point", "tolerance_Tpref_lower", "tolerance_Tpref_upper", 
                 "numhours_point", "numhours_lower", "numhours_upper", 
                 "numhours_NoB_point", "numhours_NoB_lower", "numhours_NoB_upper", 
                 "streak_point", "streak_lower", "streak_upper",
                 "streak_NoB_point", "streak_NoB_lower", "streak_NoB_upper",
                 "streak_days_point", "streak_days_lower", "streak_days_upper",
                 "streak_days_nob_point", "streak_days_nob_lower", "streak_days_nob_upper",
                 "AdultLifespan", 
                 "maxtemp_point2", "maxtemp_lower2", "maxtemp_upper2",
                 "maxtemp_NoB_point2", "maxtemp_NoB_lower2", "maxtemp_NoB_upper2",
                 "tolerance_point2", "tolerance_lower2", "tolerance_upper2", 
                 "tolerance_NoB_point2", "tolerance_NoB_lower2", "tolerance_NoB_upper2", 
                 "tolerance_Tpref_point2", "tolerance_Tpref_lower2", "tolerance_Tpref_upper2",
                 "numhours_point2", "numhours_lower2", "numhours_upper2", 
                 "numhours_NoB_point2", "numhours_NoB_lower2", "numhours_NoB_upper2",
                 "streak_point2", "streak_lower2", "streak_upper2",
                 "streak_NoB_point2", "streak_NoB_lower2", "streak_NoB_upper2",
                 "streak_days_point2", "streak_days_lower2", "streak_days_upper2",
                 "streak_days_nob_point2", "streak_days_nob_lower2", "streak_days_nob_upper2")

# Drought mask code
for (i in 1:nrow(cxq)){
  lat = cxq$Latitude[i] 
  lon = cxq$Longitude[i]
  species = "Culex_quinquefasciatus"
  
  loc = c(lon, lat)
  dstart <- "01/01/2017"
  dfinish <- "31/12/2017"
  
  # Note that era5 file must be in the working directoy and named "era5_2017.nc"
  micro<-micro_era5(loc = loc, spatial = 'era5', Usrhyt = 1,
                    dstart = dstart, dfinish = dfinish,
                    minshade = 0, maxshade = 100)
  metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, min shade
  shadeout = as.data.frame(micro$shadmet) # above ground microclimatic conditions, full shade
  tzone<-paste("Etc/GMT+",0,sep="") # append datesi 
  dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
  met <- cbind(dates,metout)
  shade <- cbind(dates,shadeout)
  
  # For drought mask: 
  # Identify days in which the prior 30 days each had soil moisture below 5%
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
    dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
    
  # We will remove days not in this 'daylist' later
    
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
  
  # For sensitivity analysis: consider a preferred body temp that is 5C lower than Topt
  Tpref_point = (traitdata$Topt_point - 3)
  Tpref_lower = (traitdata$Topt_lower - 3)
  Tpref_upper = (traitdata$Topt_upper - 3)
  
  MosqTempPoint = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                            T_pref = Topt_point, CT_max = CTmax_point, CT_min = CTmin_point, 
                            diurn = 0, nocturn = 1, crepus = 1, 
                            shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempPoint_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_point, CT_max = CTmax_point, CT_min = CTmin_point, 
                                diurn = 0, nocturn = 1, crepus = 1, 
                                shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefPoint = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Tpref_point, CT_max = CTmax_point, CT_min = CTmin_point, 
                                diurn = 0, nocturn = 1, crepus = 1,  
                                shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempLower= ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                           T_pref = Topt_lower, CT_max = CTmax_lower, CT_min = CTmin_lower, 
                           diurn = 0, nocturn = 1, crepus = 1,  
                           shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempLower_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_lower, CT_max = CTmax_lower, CT_min = CTmin_lower, 
                                diurn = 0, nocturn = 1, crepus = 1,  
                                shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefLower= ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                               T_pref = Tpref_lower, CT_max = CTmax_point, CT_min = CTmin_point, 
                               diurn = 0, nocturn = 1, crepus = 1,  
                               shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempUpper = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                            T_pref = Topt_upper, CT_max = CTmax_upper, CT_min = CTmin_upper, 
                            diurn = 0, nocturn = 1, crepus = 1,   
                            shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempUpper_NoB = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Topt_upper, CT_max = CTmax_upper, CT_min = CTmin_upper, 
                                diurn = 0, nocturn = 1, crepus = 1,  
                                shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefUpper = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Tpref_upper, CT_max = CTmax_point, CT_min = CTmin_point, 
                                diurn = 0, nocturn = 1, crepus = 1,  
                                shade_seek = 1, burrow = 0, climb = 0)
  
  point = as.data.frame(MosqTempPoint$environ)
  lower = as.data.frame(MosqTempLower$environ)
  upper = as.data.frame(MosqTempUpper$environ)
  
  point_noB = as.data.frame(MosqTempPoint_NoB$environ)
  lower_noB = as.data.frame(MosqTempLower_NoB$environ)
  upper_noB = as.data.frame(MosqTempUpper_NoB$environ)
  
  point_pref = as.data.frame(MosqTempPrefPoint$environ)
  lower_pref = as.data.frame(MosqTempPrefLower$environ)
  upper_pref = as.data.frame(MosqTempPrefUpper$environ)
  
  ModeledTemps = cbind.data.frame(point, lower$TC, upper$TC,  
                                  point_noB$TC, lower_noB$TC, upper_noB$TC,
                                  point_pref$TC, lower_pref$TC, upper_pref$TC,
                                  met$TALOC, shade$TALOC)
  colnames(ModeledTemps)[c(5,30:39)] = c("Mosq_Temp_Point", "Mosq_Temp_Lower", "Mosq_Temp_Upper", 
                                         "Mosq_Temp_NoB_Point", "Mosq_Temp_NoB_Lower", "Mosq_Temp_NoB_Upper", 
                                         "Mosq_Temp_Tpref_Point", "Mosq_Temp_Tpref_Lower", "Mosq_Temp_Tpref_Upper",
                                         "Sun_Temp", "Shade_Temp")
  
  # Remove all days not in 'daylist' (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  maxtemp_Pref_point = max(ModeledTemps_Adj$Mosq_Temp_Tpref_Point)
  maxtemp_Pref_lower = max(ModeledTemps_Adj$Mosq_Temp_Tpref_Lower)
  maxtemp_Pref_upper = max(ModeledTemps_Adj$Mosq_Temp_Tpref_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
  diff_Pref_point = CTmax_point - maxtemp_Pref_point
  diff_Pref_lower = CTmax_lower - maxtemp_Pref_lower
  diff_Pref_upper = CTmax_upper - maxtemp_Pref_upper
  
  # Number hours where predicted body temp > CTmax
  numhours_point = sum(ModeledTemps_Adj$Mosq_Temp_Point > CTmax_point) 
  numhours_lower = sum(ModeledTemps_Adj$Mosq_Temp_Lower > CTmax_lower) 
  numhours_upper = sum(ModeledTemps_Adj$Mosq_Temp_Upper > CTmax_upper) 
  
  numhours_NoB_point = sum(ModeledTemps_Adj$Mosq_Temp_NoB_Point > CTmax_point) 
  numhours_NoB_lower = sum(ModeledTemps_Adj$Mosq_Temp_NoB_Lower > CTmax_lower) 
  numhours_NoB_upper = sum(ModeledTemps_Adj$Mosq_Temp_NoB_Upper > CTmax_upper) 
  
  # Number of consecutive days in thermal danger
  ModeledTemps_Adj$danger_point =  as.numeric((CTmax_point - ModeledTemps_Adj$Mosq_Temp_Point) <0)
  dangerdays_point = aggregate(ModeledTemps_Adj$danger_point, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_point = rle(dangerdays_point > 0)
  streak_dangerdays_point = as.numeric(tapply(runs_danger_point$lengths, runs_danger_point$values, max))[2]
  ModeledTemps_Adj$danger_lower =  as.numeric((CTmax_lower - ModeledTemps_Adj$Mosq_Temp_Lower) <0)
  dangerdays_lower = aggregate(ModeledTemps_Adj$danger_lower, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_lower = rle(dangerdays_lower > 0)
  streak_dangerdays_lower = as.numeric(tapply(runs_danger_lower$lengths, runs_danger_lower$values, max))[2]
  ModeledTemps_Adj$danger_upper =  as.numeric((CTmax_upper - ModeledTemps_Adj$Mosq_Temp_Upper) <0)
  dangerdays_upper = aggregate(ModeledTemps_Adj$danger_upper, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_upper = rle(dangerdays_upper > 0)
  streak_dangerdays_upper = as.numeric(tapply(runs_danger_upper$lengths, runs_danger_upper$values, max))[2]
  
  # Number of consecutive days in thermal danger without behavior
  ModeledTemps_Adj$danger_nob_point =  as.numeric((CTmax_point - ModeledTemps_Adj$Mosq_Temp_NoB_Point) <0)
  dangerdays_nob_point = aggregate(ModeledTemps_Adj$danger_nob_point, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_nob_point = rle(dangerdays_nob_point > 0)
  streak_dangerdays_nob_point = as.numeric(tapply(runs_danger_nob_point$lengths, runs_danger_nob_point$values, max))[2]
  ModeledTemps_Adj$danger_nob_lower =  as.numeric((CTmax_lower - ModeledTemps_Adj$Mosq_Temp_NoB_Lower) <0)
  dangerdays_nob_lower = aggregate(ModeledTemps_Adj$danger_nob_lower, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_nob_lower = rle(dangerdays_nob_lower > 0)
  streak_dangerdays_nob_lower = as.numeric(tapply(runs_danger_nob_lower$lengths, runs_danger_nob_lower$values, max))[2]
  ModeledTemps_Adj$danger_nob_upper =  as.numeric((CTmax_upper - ModeledTemps_Adj$Mosq_Temp_NoB_Upper) <0)
  dangerdays_nob_upper = aggregate(ModeledTemps_Adj$danger_nob_upper, list(ModeledTemps_Adj$DAY), sum)$x
  runs_danger_nob_upper = rle(dangerdays_nob_upper > 0)
  streak_dangerdays_nob_upper = as.numeric(tapply(runs_danger_nob_upper$lengths, runs_danger_nob_upper$values, max))[2]
  
  # Length of longest streak above CTmax (in hours)
  runs_point = rle(ModeledTemps_Adj$Mosq_Temp_Point > CTmax_point)
  runs_lower = rle(ModeledTemps_Adj$Mosq_Temp_Lower > CTmax_lower)
  runs_upper = rle(ModeledTemps_Adj$Mosq_Temp_Upper > CTmax_upper)
  streak_point = as.numeric(tapply(runs_point$lengths, runs_point$values, max))[2]
  streak_lower = as.numeric(tapply(runs_lower$lengths, runs_lower$values, max))[2]
  streak_upper = as.numeric(tapply(runs_upper$lengths, runs_upper$values, max))[2]
  
  runs_NoB_point = rle(ModeledTemps_Adj$Mosq_Temp_NoB_Point > CTmax_point)
  runs_NoB_lower = rle(ModeledTemps_Adj$Mosq_Temp_NoB_Lower > CTmax_lower)
  runs_NoB_upper = rle(ModeledTemps_Adj$Mosq_Temp_NoB_Upper > CTmax_upper)
  streak_NoB_point = as.numeric(tapply(runs_NoB_point$lengths, runs_NoB_point$values, max))[2]
  streak_NoB_lower = as.numeric(tapply(runs_NoB_lower$lengths, runs_NoB_lower$values, max))[2]
  streak_NoB_upper = as.numeric(tapply(runs_NoB_upper$lengths, runs_NoB_upper$values, max))[2]
  
  # Repeat with adult life span
  adulttrait = speciesdata[speciesdata$Trait == "AdultLifespan",]  
  
  CTmax_point2 = adulttrait$CTmax_point
  CTmax_lower2 = adulttrait$CTmax_lower
  CTmax_upper2 = adulttrait$CTmax_upper
  
  Topt_point2 = adulttrait$Topt_point
  Topt_lower2 = adulttrait$Topt_lower
  Topt_upper2 = adulttrait$Topt_upper
  
  CTmin_point2 = adulttrait$CTmin_point
  CTmin_lower2 = adulttrait$CTmin_lower
  CTmin_upper2 = adulttrait$CTmin_upper
  
  # For sensitivity analysis: consider a preferred body temp that is 5C lower than Topt
  Tpref_point2 = (adulttrait$Topt_point - 3)
  Tpref_lower2 = (adulttrait$Topt_lower - 3)
  Tpref_upper2 = (adulttrait$Topt_upper - 3)
  
  MosqTempPoint2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                             T_pref = Topt_point2, CT_max = CTmax_point2, CT_min = CTmin_point2, 
                             diurn = 0, nocturn = 1, crepus = 1, 
                             shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempPoint_NoB2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                 T_pref = Topt_point2, CT_max = CTmax_point2, CT_min = CTmin_point2, 
                                 diurn = 0, nocturn = 1, crepus = 1, 
                                 shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefPoint2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                 T_pref = Tpref_point2, CT_max = CTmax_point2, CT_min = CTmin_point2, 
                                 diurn = 0, nocturn = 1, crepus = 1,  
                                 shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempLower2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                             T_pref = Topt_lower2, CT_max = CTmax_lower2, CT_min = CTmin_lower2, 
                             diurn = 0, nocturn = 1, crepus = 1,  
                             shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempLower_NoB2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                 T_pref = Topt_lower2, CT_max = CTmax_lower2, CT_min = CTmin_lower2, 
                                 diurn = 0, nocturn = 1, crepus = 1,  
                                 shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefLower2= ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                T_pref = Tpref_lower2, CT_max = CTmax_point2, CT_min = CTmin_point2, 
                                diurn = 0, nocturn = 1, crepus = 1,  
                                shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempUpper2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                             T_pref = Topt_upper2, CT_max = CTmax_upper2, CT_min = CTmin_upper2, 
                             diurn = 0, nocturn = 1, crepus = 1,   
                             shade_seek = 1, burrow = 0, climb = 0)
  
  MosqTempUpper_NoB2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                 T_pref = Topt_upper2, CT_max = CTmax_upper2, CT_min = CTmin_upper2, 
                                 diurn = 0, nocturn = 1, crepus = 1,  
                                 shade_seek = 0, burrow = 0, climb = 0)
  
  MosqTempPrefUpper2 = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                                 T_pref = Tpref_upper2, CT_max = CTmax_point2, CT_min = CTmin_point2, 
                                 diurn = 0, nocturn = 1, crepus = 1,  
                                 shade_seek = 1, burrow = 0, climb = 0)
  
  point2 = as.data.frame(MosqTempPoint2$environ)
  lower2 = as.data.frame(MosqTempLower2$environ)
  upper2 = as.data.frame(MosqTempUpper2$environ)
  
  point_noB2 = as.data.frame(MosqTempPoint_NoB2$environ)
  lower_noB2 = as.data.frame(MosqTempLower_NoB2$environ)
  upper_noB2 = as.data.frame(MosqTempUpper_NoB2$environ)
  
  point_pref2 = as.data.frame(MosqTempPrefPoint2$environ)
  lower_pref2 = as.data.frame(MosqTempPrefLower2$environ)
  upper_pref2 = as.data.frame(MosqTempPrefUpper2$environ)
  
  ModeledTemps2 = cbind.data.frame(point2, lower2$TC, upper2$TC,  
                                   point_noB2$TC, lower_noB2$TC, upper_noB2$TC,
                                   point_pref2$TC, lower_pref2$TC, upper_pref2$TC,
                                   met$TALOC, shade$TALOC)
  
  colnames(ModeledTemps2)[c(5,30:39)] = c("Mosq_Temp_Point2", "Mosq_Temp_Lower2", "Mosq_Temp_Upper2", 
                                          "Mosq_Temp_NoB_Point2", "Mosq_Temp_NoB_Lower2", "Mosq_Temp_NoB_Upper2", 
                                          "Mosq_Temp_Tpref_Point2", "Mosq_Temp_Tpref_Lower2", "Mosq_Temp_Tpref_Upper2",
                                          "Sun_Temp2", "Shade_Temp2")
  
  # Remove all days not in 'daylist' (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj2 = ModeledTemps2[ModeledTemps2$DOY %in% daylist,]
  
  
  # Maximum body temperature
  maxtemp_point2 = max(ModeledTemps_Adj2$Mosq_Temp_Point2)
  maxtemp_lower2 = max(ModeledTemps_Adj2$Mosq_Temp_Lower2)
  maxtemp_upper2 = max(ModeledTemps_Adj2$Mosq_Temp_Upper2)
  
  maxtemp_NoB_point2 = max(ModeledTemps_Adj2$Mosq_Temp_NoB_Point2)
  maxtemp_NoB_lower2 = max(ModeledTemps_Adj2$Mosq_Temp_NoB_Lower2)
  maxtemp_NoB_upper2 = max(ModeledTemps_Adj2$Mosq_Temp_NoB_Upper2)
  
  maxtemp_Pref_point2 = max(ModeledTemps_Adj2$Mosq_Temp_Tpref_Point2)
  maxtemp_Pref_lower2 = max(ModeledTemps_Adj2$Mosq_Temp_Tpref_Lower2)
  maxtemp_Pref_upper2 = max(ModeledTemps_Adj2$Mosq_Temp_Tpref_Upper2)
  
  # Thermal safety margin
  diff_point2 = CTmax_point2 - maxtemp_point2
  diff_lower2 = CTmax_lower2 - maxtemp_lower2
  diff_upper2 = CTmax_upper2 - maxtemp_upper2
  
  diff_NoB_point2 = CTmax_point2 - maxtemp_NoB_point2
  diff_NoB_lower2 = CTmax_lower2 - maxtemp_NoB_lower2
  diff_NoB_upper2 = CTmax_upper2 - maxtemp_NoB_upper2
  
  diff_Pref_point2 = CTmax_point2 - maxtemp_Pref_point2
  diff_Pref_lower2 = CTmax_lower2 - maxtemp_Pref_lower2
  diff_Pref_upper2 = CTmax_upper2 - maxtemp_Pref_upper2
  
  # Number hours where predicted body temp > CTmax
  numhours_point2 = sum(ModeledTemps_Adj2$Mosq_Temp_Point2 > CTmax_point2) 
  numhours_lower2 = sum(ModeledTemps_Adj2$Mosq_Temp_Lower2 > CTmax_lower2) 
  numhours_upper2 = sum(ModeledTemps_Adj2$Mosq_Temp_Upper2 > CTmax_upper2) 
  
  numhours_NoB_point2 = sum(ModeledTemps_Adj2$Mosq_Temp_NoB_Point2 > CTmax_point2) 
  numhours_NoB_lower2 = sum(ModeledTemps_Adj2$Mosq_Temp_NoB_Lower2 > CTmax_lower2) 
  numhours_NoB_upper2 = sum(ModeledTemps_Adj2$Mosq_Temp_NoB_Upper2 > CTmax_upper2) 
  
  # Number of consecutive days in thermal danger
  ModeledTemps_Adj2$danger_point2 =  as.numeric((CTmax_point2 - ModeledTemps_Adj2$Mosq_Temp_Point2) <0)
  dangerdays_point2 = aggregate(ModeledTemps_Adj2$danger_point2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_point2 = rle(dangerdays_point2 > 0)
  streak_dangerdays_point2 = as.numeric(tapply(runs_danger_point2$lengths, runs_danger_point2$values, max))[2]
  ModeledTemps_Adj2$danger_lower2 =  as.numeric((CTmax_lower2 - ModeledTemps_Adj2$Mosq_Temp_Lower2) <0)
  dangerdays_lower2 = aggregate(ModeledTemps_Adj2$danger_lower2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_lower2 = rle(dangerdays_lower2 > 0)
  streak_dangerdays_lower2 = as.numeric(tapply(runs_danger_lower2$lengths, runs_danger_lower2$values, max))[2]
  ModeledTemps_Adj2$danger_upper2 =  as.numeric((CTmax_upper2 - ModeledTemps_Adj2$Mosq_Temp_Upper2) <0)
  dangerdays_upper2 = aggregate(ModeledTemps_Adj2$danger_upper2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_upper2 = rle(dangerdays_upper2 > 0)
  streak_dangerdays_upper2 = as.numeric(tapply(runs_danger_upper2$lengths, runs_danger_upper2$values, max))[2]
  
  # Number of consecutive days in thermal danger without behavior
  ModeledTemps_Adj2$danger_nob_point2 =  as.numeric((CTmax_point2 - ModeledTemps_Adj2$Mosq_Temp_NoB_Point2) <0)
  dangerdays_nob_point2 = aggregate(ModeledTemps_Adj2$danger_nob_point2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_nob_point2 = rle(dangerdays_nob_point2 > 0)
  streak_dangerdays_nob_point2 = as.numeric(tapply(runs_danger_nob_point2$lengths, runs_danger_nob_point2$values, max))[2]
  ModeledTemps_Adj2$danger_nob_lower2 =  as.numeric((CTmax_lower2 - ModeledTemps_Adj2$Mosq_Temp_NoB_Lower2) <0)
  dangerdays_nob_lower2 = aggregate(ModeledTemps_Adj2$danger_nob_lower2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_nob_lower2 = rle(dangerdays_nob_lower2 > 0)
  streak_dangerdays_nob_lower2 = as.numeric(tapply(runs_danger_nob_lower2$lengths, runs_danger_nob_lower2$values, max))[2]
  ModeledTemps_Adj2$danger_nob_upper2 =  as.numeric((CTmax_upper2 - ModeledTemps_Adj2$Mosq_Temp_NoB_Upper2) <0)
  dangerdays_nob_upper2 = aggregate(ModeledTemps_Adj2$danger_nob_upper2, list(ModeledTemps_Adj2$DAY), sum)$x
  runs_danger_nob_upper2 = rle(dangerdays_nob_upper2 > 0)
  streak_dangerdays_nob_upper2 = as.numeric(tapply(runs_danger_nob_upper2$lengths, runs_danger_nob_upper2$values, max))[2]
  
  # Length of longest streak above CTmax (in hours)
  runs_point2 = rle(ModeledTemps_Adj2$Mosq_Temp_Point2 > CTmax_point2)
  runs_lower2 = rle(ModeledTemps_Adj2$Mosq_Temp_Lower2 > CTmax_lower2)
  runs_upper2 = rle(ModeledTemps_Adj2$Mosq_Temp_Upper2 > CTmax_upper2)
  streak_point2 = as.numeric(tapply(runs_point2$lengths, runs_point2$values, max))[2]
  streak_lower2 = as.numeric(tapply(runs_lower2$lengths, runs_lower2$values, max))[2]
  streak_upper2 = as.numeric(tapply(runs_upper2$lengths, runs_upper2$values, max))[2]
  
  runs_NoB_point2 = rle(ModeledTemps_Adj2$Mosq_Temp_NoB_Point2 > CTmax_point2)
  runs_NoB_lower2 = rle(ModeledTemps_Adj2$Mosq_Temp_NoB_Lower2 > CTmax_lower2)
  runs_NoB_upper2 = rle(ModeledTemps_Adj2$Mosq_Temp_NoB_Upper2 > CTmax_upper2)
  streak_NoB_point2 = as.numeric(tapply(runs_NoB_point2$lengths, runs_NoB_point2$values, max))[2]
  streak_NoB_lower2 = as.numeric(tapply(runs_NoB_lower2$lengths, runs_NoB_lower2$values, max))[2]
  streak_NoB_upper2 = as.numeric(tapply(runs_NoB_upper2$lengths, runs_NoB_upper2$values, max))[2]
  
  df[i,] = c(species, HighestTrait, lat, lon, maxtemp_point, maxtemp_lower, maxtemp_upper, 
             maxtemp_NoB_point, maxtemp_NoB_lower, maxtemp_NoB_upper,
             diff_point, diff_lower, diff_upper,
             diff_NoB_point, diff_NoB_lower, diff_NoB_upper,
             diff_Pref_point, diff_Pref_lower, diff_Pref_upper,
             numhours_point, numhours_lower, numhours_upper,
             numhours_NoB_point, numhours_NoB_lower, numhours_NoB_upper,
             streak_point, streak_lower, streak_upper,
             streak_NoB_point, streak_NoB_lower, streak_NoB_upper,
             streak_dangerdays_point, streak_dangerdays_lower, streak_dangerdays_upper,
             streak_dangerdays_nob_point, streak_dangerdays_nob_lower, streak_dangerdays_nob_upper,
             "AdultLifespan",
             maxtemp_point2, maxtemp_lower2, maxtemp_upper2,
             maxtemp_NoB_point2, maxtemp_NoB_lower2, maxtemp_NoB_upper2,
             diff_point2, diff_lower2, diff_upper2,
             diff_NoB_point2, diff_NoB_lower2, diff_NoB_upper2,
             diff_Pref_point2, diff_Pref_lower2, diff_Pref_upper2,
             numhours_point2, numhours_lower2, numhours_upper2,
             numhours_NoB_point2, numhours_NoB_lower2, numhours_NoB_upper2,
             streak_point2, streak_lower2, streak_upper2,
             streak_NoB_point2, streak_NoB_lower2, streak_NoB_upper2,
             streak_dangerdays_point2, streak_dangerdays_lower2, streak_dangerdays_upper2,
             streak_dangerdays_nob_point2, streak_dangerdays_nob_lower2, streak_dangerdays_nob_upper2)
}

write.csv(df, "~/Downloads/CxQuinque_TSM_-44to-35_WithDroughtMask.csv")

