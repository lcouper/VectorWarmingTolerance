#### Vector Thermal Safety TSM esimation with seasonality

# This script estimates TSMs and related indices with precipitation-driven seasonality incorporated
# Specifically if a 30-day period of low soil moisture (<5%) occurs,
# subsequent days are masked from calculations until soil moisture increases again
# (i.e., thermal extremes during this period are excluded)

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles")

library(dplyr)
library(mcera5)
library(ecmwfr)
library(NicheMapR)
library(ncdf4)
library(tidync)
library(microclima) 
library(NicheMapR) 


#### Estimate TSMs, excluding days where soil moisture <5% for past 30 days ####

# shown below for generic species info contained in 'data'
setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")

df1 = as.data.frame(matrix(nrow = nrow(data), ncol = 34))
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

for (i in 1:nrow(data)){
  lat = data$Lat[i]
  lon = data$Long[i]
  species = "Anopheles_gambiae"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")

# Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data

    for (j in 31:365)
    {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
    daylist = which(dayindex == FALSE)}

# Run NicheMapR as before. Later, we will exclude the dates identified above

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

# Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
    ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
    
    # Recalculate thermal vulnerability indicies
    
    # Maximum body temperature
    maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
    maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
    maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
    
    maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
    maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
    maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
    
    # Thermal safety margin
    diff_point = CTmax_point - maxtemp_point
    diff_lower = CTmax_lower - maxtemp_lower
    diff_upper = CTmax_upper - maxtemp_upper
    
    diff_NoB_point = CTmax_point - maxtemp_NoB_point
    diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
    diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
    
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

write.csv(df1, "~/Downloads/Species_TSM_WithSeasonality_Range.csv")



#### Anopheles stephensi #####

setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")
AnSteph = read.csv("Vector_TSM/AnSteph_TSM_WithoutSeasonality.csv")

# subset to given band, and run all points for that band
AnSteph = subset(AnSteph, latitude >=0 & latitude <=5)

df1 = as.data.frame(matrix(nrow = nrow(AnSteph), ncol = 34))
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

for (i in 1:nrow(AnSteph)){
  lat = AnSteph$latitude[i]
  lon = AnSteph$longitude[i]
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AnSteph_TSM_WithSeasonality_0to5.csv")







#### Anopheles gambiae #####
setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")
AnGambi = read.csv("Vector_TSM/AnGambiae_TSM_WithoutSeasonality.csv")

# subset to given band, and run all points for that band
AnGambi = subset(AnGambi, latitude >=15 & latitude <=25)

df1 = as.data.frame(matrix(nrow = nrow(AnGambi), ncol = 34))
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

for (i in 1:nrow(AnGambi)){
  lat = AnGambi$latitude[i]
  lon = AnGambi$longitude[i]
  species = "Anopheles_gambiae"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AnGambiae_TSM_WithSeasonality_25to35.csv")


#### Aedes aegypti #####

setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")
AeAg = read.csv("Vector_TSM/AeAegypti_TSM_WithoutSeasonality.csv")

# subset to given band, and run all points for that band
AeAg = subset(AeAg, latitude >=-45 & latitude <=-35)

df1 = as.data.frame(matrix(nrow = nrow(AeAg), ncol = 34))
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

for (i in 1:nrow(AeAg)){
  lat = AeAg$latitude[i]
  lon = AeAg$longitude[i]
  species = "Aedes_aegypti"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AeAegypti_TSM_WithSeasonality_-35to-45.csv")


#### Aedes albopictus ####

setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")
AeAl = read.csv("Vector_TSM/AeAlbopictus_TSM_WithoutSeasonality.csv")

# subset to given band, and run all points for that band
AeAl = subset(AeAl, latitude >=45 & latitude <=55)

df1 = as.data.frame(matrix(nrow = nrow(AeAl), ncol = 34))
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

for (i in 1:nrow(AeAl)){
  lat = AeAl$latitude[i]
  lon = AeAl$longitude[i]
  species = "Aedes_albopictus"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AeAlbo_TSM_WithSeasonality_45to55.csv")


#### Aedes camptorhynchus  ####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
AeCamp  = read.csv("Vector_TSM/AeCamp_TSM_WithoutSeasonality.csv")

# subset to given band, and run all points for that band
AeCamp = subset(AeCamp, latitude >=-45 & latitude <=-35)

df1 = as.data.frame(matrix(nrow = nrow(AeCamp), ncol = 34))
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

for (i in 1:nrow(AeCamp)){
  lat = AeCamp$latitude[i]
  lon = AeCamp$longitude[i]
  species = "Aedes_camptorhynchus"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AeCamp_TSM_WithSeasonality_-35to-45.csv")


#### Aedes triseriatus #####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
AeTri = read.csv("Vector_TSM/AeTri_TSM_WithoutSeasonality.csv")

# 25 to 47 lat
AeTri = subset(AeTri, latitude >=45 & latitude <=55)

df1 = as.data.frame(matrix(nrow = nrow(AeTri), ncol = 34))
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

for (i in 1:nrow(AeTri)){
  lat = AeTri$latitude[i]
  lon = AeTri$longitude[i]
  species = "Aedes_triseriatus"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AeTri_TSM_WithSeasonality_45to55.csv")


#### Aedes vexans #####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
AeVex = read.csv("Vector_TSM/AeVexans_TSM_WithoutSeasonality.csv")

# -35 to 56
AeVex = subset(AeVex, latitude >=45 & latitude <=56)

df1 = as.data.frame(matrix(nrow = nrow(AeVex), ncol = 34))
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

for (i in 6:nrow(AeVex)){
  lat = AeVex$latitude[i]
  lon = AeVex$longitude[i]
  species = "Aedes_vexans"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/AeVex_TSM_WithSeasonality_45to56.csv")


#### Culex quinquefasciatus #####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
CxQ= read.csv("Vector_TSM/CxQuinque_TSM_WithoutSeasonality.csv")

# - 43 to 46
CxQ = subset(CxQ, latitude >=35 & latitude <=46)

df1 = as.data.frame(matrix(nrow = nrow(CxQ), ncol = 34))
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

for (i in 1:nrow(CxQ)){
  lat = CxQ$latitude[i]
  lon = CxQ$longitude[i]
  species = "Culex_quinquefasciatus"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/CxQuinque_TSM_WithSeasonality_35to45.csv")




#### Culex pipiens #####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
CxP= read.csv("Vector_TSM/CxPipiens_TSM_WithoutSeasonailty.csv")

# - 38 to 58
CxP = subset(CxP, latitude >=25 & latitude <=35)

df1 = as.data.frame(matrix(nrow = nrow(CxP), ncol = 34))
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

for (i in 1:nrow(CxP)){
  lat = CxP$latitude[i]
  lon = CxP$longitude[i]
  species = "Culex_pipiens"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/CxPip_TSM_WithSeasonality_25to35_p2.csv")



#### Culex tarsalis ####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
CxT= read.csv("Vector_TSM/CxTarsalis_TSM_WithoutSesaonality.csv")

# lat range: 18 to 66
CxT = subset(CxT, latitude >=49.01 & latitude <=67)

df1 = as.data.frame(matrix(nrow = nrow(CxT), ncol = 34))
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

for (i in 1:nrow(CxT)){
  lat = CxT$latitude[i]
  lon = CxT$longitude[i]
  species = "Culex_tarsalis"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/CxTar_TSM_WithSeasonality_49to67.csv")

#### Culex annulirostris ####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
CxA= read.csv("Vector_TSM/CxAnnul_TSM_WithoutSeasonality.csv")

# lat range: -38 to 14
CxA = subset(CxA, latitude >=5 & latitude < 15)

df1 = as.data.frame(matrix(nrow = nrow(CxA), ncol = 34))
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

for (i in 1:nrow(CxA)){
  lat = CxA$latitude[i]
  lon = CxA$longitude[i]
  species = "Culex_annulirostris"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/CxAnnul_TSM_WithSeasonality_5to15.csv")



#### Culex theileri ####

setwd("/Users/lisacouper/Documents/Current Projects/WarmingTolerance/DataFiles")
CxTh = read.csv("Vector_TSM/CxTheileri_TSM_WithoutSeasonality.csv")

# lat range: -35 to 43
CxTh = subset(CxTh, latitude >=-35 & latitude < -25)

df1 = as.data.frame(matrix(nrow = nrow(CxTh), ncol = 34))
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

for (i in 1:nrow(CxTh)){
  lat = CxTh$latitude[i]
  lon = CxTh$longitude[i]
  species = "Culex_theileri"
  
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
  SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
  colnames(SoilMoisture) = c("DOY", "PCTWET")
  
  # Identify days in which the prior 30 days each had soil moisture below 5%
  dayindex = c(rep(FALSE, 30), rep(NA, 335)) # excluding first month since don't have pre-2017 data
  
  for (j in 31:365)
  {dayindex[j] = all(SoilMoisture$PCTWET[(j-31):(j-1)] < 5)
  daylist = which(dayindex == FALSE)}
  
  # Run NicheMapR as before. Later, we will exclude the dates identified above
  
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
  
  # Remove all days identified above (as these represent days when mosquitoes likely inactive due to aridity)
  ModeledTemps_Adj = ModeledTemps[ModeledTemps$DOY %in% daylist,]
  
  # Recalculate thermal vulnerability indicies
  
  # Maximum body temperature
  maxtemp_point = max(ModeledTemps_Adj$Mosq_Temp_Point)
  maxtemp_lower = max(ModeledTemps_Adj$Mosq_Temp_Lower)
  maxtemp_upper = max(ModeledTemps_Adj$Mosq_Temp_Upper)
  
  maxtemp_NoB_point = max(ModeledTemps_Adj$Mosq_Temp_NoB_Point)
  maxtemp_NoB_lower = max(ModeledTemps_Adj$Mosq_Temp_NoB_Lower)
  maxtemp_NoB_upper = max(ModeledTemps_Adj$Mosq_Temp_NoB_Upper)
  
  # Thermal safety margin
  diff_point = CTmax_point - maxtemp_point
  diff_lower = CTmax_lower - maxtemp_lower
  diff_upper = CTmax_upper - maxtemp_upper
  
  diff_NoB_point = CTmax_point - maxtemp_NoB_point
  diff_NoB_lower = CTmax_lower - maxtemp_NoB_lower
  diff_NoB_upper = CTmax_upper - maxtemp_NoB_upper
  
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

write.csv(df1, "~/Downloads/CxTheileri_TSM_WithSeasonality_-35to-25.csv")






#### Comparing Soil Moisture Estimate against total precipitatoin for validation #####

# Using Anopheles gambiae record from Niger: 
agambiae = read.csv("VectorOccurrence/Africa Vectors database_1898-2016.csv", header = T)
agambiae = agambiae[agambiae$An.gambiae.ss == "Y",]
agambi <- agambiae %>% dplyr::distinct(Lat, Long, .keep_all = TRUE)
agambi = agambi[!is.na(agambi$Lat),]
agambi = agambi[!is.na(agambi$Long),]
agambi = agambi[agambi$Country == "Niger",]


# Calculate total precipitation at this location, across the year

Climate = nc_open("era5_2017.nc")
Testlat = 16.992
Testlon = 7.979
lab = "tp"
lon <- ncvar_get(Climate, "longitude")
lat <- ncvar_get(Climate, "latitude")
ts.array <- ncvar_get(Climate, lab)

# Extract total precip values at this location across the year
# Gives index of closest lat/lon in the ERA5 data
closestlat = which(abs(lat-Testlat)==min(abs(lat-Testlat)))
closestlon = which(abs(lon-Testlon)==min(abs(lon-Testlon)))
# extract precip at that location
Precips = as.data.frame(ts.array[closestlon,closestlat,])
DOY = rep(1:364, each = 24)
Precips$DOY = DOY
colnames(Precips)[1] = "Precip"
PrecipByDay = aggregate(Precips$Precip ~ Precips$DOY, FUN = sum)
colnames(PrecipByDay) = c("DOY", "TP")
# Convert units from m to mm
PrecipByDay$TPmm = PrecipByDay$TP * 1000
dry2 = PrecipByDay[PrecipByDay$TPmm == 0,]


# Compare to NicheMapR soil moisture
species = "Anopheles_gambiae"
loc = c(Testlon, Testlat)
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

SoilMoisture = aggregate(metout$PCTWET ~ metout$DOY, FUN = mean)
colnames(SoilMoisture) = c("DOY", "PCTWET")
# Pull out days where soil moisture was <5
dry1 = SoilMoisture[SoilMoisture$PCTWET < 5,]


#### Plot daily total precip and soil moisture (Supp. Fig) #####
par(mar = c(5, 6, 4, 6) + 0.3)  # Leave space for z axis
plot(PrecipByDay$TPmm ~ Days, type = "l", lwd = 1.5, cex.axis = 1.5, cex.lab = 2,
     xlab = "Julian day", ylab = "total precip (mm)") # first plot
par(new = TRUE)
plot(SoilMoisture$PCTWET ~ Days, col = "red", type = "l", lwd = 1.5, lty = 2,
     axes = FALSE, bty = "n", xlab = "", ylab = "")
axis(side=4, at = pretty(range(SoilMoisture$PCTWET)), cex.axis = 1.5)
mtext("soil moisture (%)", side = 4, line = 3, cex = 2)
legend(300, 45, legend=c("total precip", "soil moisture"),
       col=c("black", "red"), lwd = 2, lty=1:2, cex=1.5)


