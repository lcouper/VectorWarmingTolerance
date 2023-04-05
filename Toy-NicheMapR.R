### Toy version: Calculating vector warming tolerance using era5 climate data ####

# from : https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html

#### 1. Load libraries and data ####
library(mcera5)
library(dplyr)
library(ecmwfr)
library(ncdf4)
library(curl)
library(keyring)
library(abind)
library(lubridate)
library(tidync)
library(microclima) # if you need to install this package: remotes::install_github("ilyamaclean/microclima")
library(NicheMapR) # if you need to install this package: remotes::install_github("mrke/NicheMapR")

# Set this to wherever you'd like to store the downloaded climate data files
setwd("~/Documents/Current_Projects/WarmingTolerance/DataFiles")

#### 2. Request ERA5 data (needed for microclimate model) ####

# Accessing era5 requires authentication key (below is Lisa's, but you can request your own here:
# https://cds.climate.copernicus.eu/user/register?destination=%2F%23!%2Fhome

uid <- "188426"
cds_api_key <- "4b2d9d42-30fe-4443-bef4-c8c9fc8819db"
ecmwfr::wf_set_key(user = uid, key = cds_api_key, service = "cds")

# set desired bounding coordinates for climate data (in WGS84 / EPSG:4326)
xmn <- 130
xmx <- 132
ymn <- -26
ymx <- -24

# set desired temporal extent
st_time <- lubridate::ymd("2020:01:01")
en_time <- lubridate::ymd("2020:12:31")

file_prefix <- "era5"
file_path <- getwd()

# build a request (covering the spatial and temporal extent defiend above)
req <- build_era5_request(xmin = xmn, xmax = xmx,
                          ymin = ymn, ymax = ymx,
                          start_time = st_time,
                          end_time = en_time,
                          outfile_name = file_prefix)
# climatevars = req[[1]][["variable"]] # To see what climate variables you get data on
# Submit request to access the ERA5 climate data
request_era5(request = req, uid = uid, out_path = file_path, overwrite = TRUE)

# With the now downloaded ERA5 data, run the microcliamte model 
# for a particular date range and point location (lat/long)

#### 3. Run the microclimate model for a given lat/lon ####

loc <- c(131, -25) # location in australia
dstart <- "01/01/2020" # first day to include in analysis
dfinish <- "31/12/2020" # last day to include in analysis

micro<-micro_era5(loc = loc, spatial = 'era5',
                    dstart = dstart, dfinish = dfinish,
                    minshade = 0, maxshade = 100) # specify a possible fully shaded and fully exposed environment

metout<-as.data.frame(micro$metout) # above ground microclimatic conditions, minimum shade
shadeout = as.data.frame(micro$shadmet) # above ground microclimatic conditions, full shade

tzone<-paste("Etc/GMT+",0,sep="") # append dates
dates<-seq(as.POSIXct(dstart, format="%d/%m/%Y",tz=tzone)-3600*12, as.POSIXct(dfinish, format="%d/%m/%Y",tz=tzone)+3600*11, by="hours")
met <- cbind(dates,metout)
shade <- cbind(dates,shadeout)

#### 4. Run the ectotherm model for a given species and trait

# ectotherm model takes many inputs including organism weight, thermal preference (here using Topt)
# thermal limits, capacity to seek shade, climb, burrow, daily activity (e.g., noctural, diurnal)
# Below, I've used the theraml limits for Aedes aegypti, trait: fecundity (from Mordecai et al. 2019)
CTmax = 34.4
Topt =  29.6
CTmin = 14.7
MosqTemp = ectotherm(Ww_g = 0.003, shape = 2, alpha_min = 0.85, alpha_max = 0.85,
                       T_pref = Topt, CT_max = CTmax, CT_min = CTmin, 
                       diurn = 1, nocturn = 0, crepus = 1, 
                       shade_seek = 1, burrow = 0, climb = 0)
  x = as.data.frame(MosqTemp$environ) # Pull oyt the output from the ectotherm model
  ModeledTemps = cbind.data.frame(MosqTemp$environ, met[,4:5], shade[,4]) # Bind it with output from microclimate model
  colnames(ModeledTemps)[c(5,29:31)] = c("Mosq_Temp", "Sun_Temp", "Ref_Temp", "Shade_Temp")
  
#### 5. Plot output of ectotherm and microclimate model #####

ModeledTemps$hour = 1:8784 # add column to track hour in year

# Plot modeled temperature of full sun, full shade, and mosquito body
plot(ModeledTemps$Sun_Temp ~ ModeledTemps$hour, col = alpha("red", 0.4), type = "l", cex = 0.5, 
     xlab = "Day/Hour of Year", ylab = "Temp (C)")
points(y = ModeledTemps$Shade_Temp, x = ModeledTemps$hour, col = alpha("blue", 0.6), type = "l", cex = 0.5)
points(y = ModeledTemps$Mosq_Temp, x = ModeledTemps$hour, pch =16, col = alpha("grey20", 0.6), type = "l", cex = 0.6)
abline(h = CTmax, lty = 3, lwd = 2)
legend("topleft", fill = c(alpha("red", 0.4), alpha("blue", 0.4), alpha("grey20", 0.6)), 
       legend= c("Full sun", "Full shade", "Body temp"))
