#############################################################
########### CODE FOR VECTOR TSM OCCURRENCE POINTS ###########
########### WRITTEN BY LISA COUPER ##########################

# This script is used to 
# a) Plot occurrence records used in the analysis
# b) Apply a gridded sampling approach to select occurrence records for analysis

#### 0a. Load libraries ####

# For GAMS
library(mgcv)
library(magrittr) 
library(scales)
library(elevatr)
library(tidyverse)
library(ggpmisc)
library(data.table)
library(maps)
library(mapdata)
library(maptools)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

#### 0b. Load colors and data frames #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles")

aeaegypti = fread("Vector_TSM/AeAegypti_TSM_DroughtMask_Combined_WithElevation.csv")
aealbo = fread("Vector_TSM/AeAlbo_TSM_DroughtMask_Combined_WithElevation.csv")
agambi = fread("Vector_TSM/AnGambiae_TSM_DroughtMask_Combined_WithElevation.csv")
asteph = fread("Vector_TSM/AnSteph_TSM_DroughtMask_Combined_WithElevation.csv")
cxa = fread("Vector_TSM/CxAnnul_TSM_DroughtMask_Combined_WithElevation.csv")
cxpip = fread("Vector_TSM/CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")
cxtar = fread("Vector_TSM/CxTarsalis_TSM_DroughtMask_Combined_WithElevation.csv")
cxq = fread("Vector_TSM/CxQuinque_TSM_DroughtMask_Combined_WithElevation.csv")


#aealbo = fread("Vector Occurrence Data/AeAlbopictus_SubSampledPoints.csv")
#aeaegypti = fread("Vector Occurrence Data/AeAegypti_SubSampledPoints.csv")
#agambi = fread("Vector Occurrence Data/AnGambi_SubSampledPoints.csv")
#asteph = fread("Vector Occurrence Data/AnSteph_SubSampledPoints.csv")
#cxq = fread("Vector Occurrence Data/CxQuinque_SubSampledPoints.csv")
#cxp = fread("Vector Occurrence Data/CxPipiens_SubSampledPoints.csv")
#cxt = fread("Vector Occurrence Data/CxTarsalis_SubSampledPoints.csv")
#cxa = fread("Vector Occurrence Data/CxAnnul_SubSampledPoints.csv")

SpeciesColors = c("#9b0000", "#ed2939",
                  "#cc5801", "#f4bb00",
                  "#7cae00", "#004B00",
                  "#52b2bf", "#016064", 
                  "#57a0d3", "#3944bc",
                  "#311465", "#dda0dd")
SpeciesList  = c("Aedes_aegypti", "Aedes_albopictus", 
                 "Anopheles_gambiae", "Anopheles_stephensi", 
                 "Culex_pipiens", "Culex_quinquefasciatus",
                 "Aedes_camptorhynchus", "Culex_annulirostris", 
                 "Aedes_triseriaatus", "Cx_tarsalis",
                 "Aedes_vexans", "Culex_theileri")

################[START OF ONLY RUN THIS PART ONCE] ############ 
##### Gridded sampling approach & obtain elevation data ######

###### Aedes albopictus ######
aedes = read.csv("Vector Occurrence Data/aegypti_albopictus.csv")

# subset to just published, point records of Aedes albopicus 
aealboOG = subset(aedes, VECTOR=="Aedes albopictus" & LOCATION_TYPE=="point" &
                    SOURCE_TYPE=="published")

# remove any duplicates (occurences from same lat, lon)
aealbo <- aealboOG %>% dplyr::distinct(Latitude, Longitude, .keep_all = TRUE)

# Ae. albopictus latitudinal range: -37 to 49

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

aealbo = subset(aealbo, Latitude>=45 & Latitude<=49)
aealbo$LongitudeF <- as.factor(as.integer(aealbo$Longitude))
aealbo <- aealbo %>% group_by(LongitudeF) %>% slice_sample(n=1)

# Identify elevation for each occurrence point
aealbo_sites <- st_as_sf(aealbo, coords = c("Longitude", "Latitude"), 
                         crs = 4326, agr = "constant")
aealbo <- get_elev_point(aealbo_sites, src = "aws")

# Separate out latitude and longitude
aealbo <- aealbo %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                                   Latitude = sf::st_coordinates(.)[,2])

# output this list
#fwrite(aealbo, "Vector Occurrence Data/AeAlbopictus_SubSampledPoints.csv")

###### Aedes aegypti ######

# subset to just published, point records of Aedes aegypti or albopicus 
aeaegypti = subset(aedes, VECTOR=="Aedes aegypti" & LOCATION_TYPE=="point" &
                     SOURCE_TYPE=="published")

# remove any duplicates (occurences from same lat, lon)
aeaegypti <- aeaegypti %>% dplyr::distinct(Latitude, Longitude, .keep_all = TRUE)

# Ae. aegypti total latitudinal range: -39 to 47

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

aeaegypti = subset(aeaegypti, Latitude>=45 & Latitude<=49)
aeaegypti$LongitudeF <- as.factor(as.integer(aeaegypti$Longitude))
aeaegypti <- aeaegypti %>% group_by(LongitudeF) %>% slice_sample(n=1)

# if still too many:
odds = seq(from = 1, to = nrow(aeaegypti), by = 2)
aeaegypti <- aeaegypti[odds,]

# Identify elevation for each occurrence point
aeaegypti_sites <- st_as_sf(aeaegypti, coords = c("Longitude", "Latitude"), 
                            crs = 4326, agr = "constant")
aeaegypti <- get_elev_point(aeaegypti_sites, src = "aws")

# Separate out latitude and longitude
aeaegypti <- aeaegypti %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                                         Latitude = sf::st_coordinates(.)[,2])

# output this list
#fwrite(aeaegypti, "Vector Occurrence Data/AeAegypti_SubSampledPoints.csv")


###### Anopheles gambiae #######

agambiae = read.csv("Vector Occurrence Data/Africa Vectors database_1898-2016.csv", header = T)
agambiae = agambiae[agambiae$An.gambiae.ss == "Y",]

# remove any duplicates (occurrences from same lat, lon)
agambi <- agambiae %>% dplyr::distinct(Lat, Long, .keep_all = TRUE)

agambi = agambi[!is.na(agambi$Lat),]
agambi = agambi[!is.na(agambi$Long),]

# range: -27 to 21

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

agambi = subset(agambi, Lat >=15 & Lat <= 25)
agambi$LongitudeF <- as.factor(as.integer(agambi$Long))
agambi <- agambi %>% group_by(LongitudeF) %>% slice_sample(n=1)

# Identify elevation for each occurrence point
agambi_sites <- st_as_sf(agambi, coords = c("Long", "Lat"), 
                         crs = 4326, agr = "constant")
agambi <- get_elev_point(agambi_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
agambi <- agambi %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                                   Latitude = sf::st_coordinates(.)[,2])

# output this list
#fwrite(agambi, "Vector Occurrence Data/AnGambi_SubSampledPoints.csv")


### Add-on data source of GBIF ###

agambi2 = read.csv("Vector Occurrence Data/Anopheles Gambiae.csv")
agambi2 = agambi2[!is.na(agambi2$decimalLatitude),]
agambi2 = agambi2[!is.na(agambi2$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
agambi2 = agambi2  %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# range: -30 to 20

# Apply gridded sampling approach: 
agambi2 = subset(agambi2, decimalLatitude >=15 & decimalLatitude <= 25)
agambi2$LongitudeF <- as.factor(as.integer(agambi2$decimalLongitude))
agambi2 <- agambi2 %>% group_by(LongitudeF) %>% slice_sample(n=1)

# if still too many:
odds = seq(from = 1, to = nrow(agambi2), by = 2)
agambi2 <- agambi2[odds,]


# Identify elevation for each occurrence point
agambi_sites <- st_as_sf(agambi2, coords = c("decimalLongitude", "decimalLatitude"), 
                         crs = 4326, agr = "constant")
agambi <- get_elev_point(agambi_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
agambi <- agambi %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                                   Latitude = sf::st_coordinates(.)[,2])

#fwrite(agambi, "Vector Occurrence Data/AnGambi_SubSampledPoints_2.csv")

# Note: manually combined output of these two data sources outside of R

###### Anopheles stephensi ######

asteph1 = read.csv("Vector Occurrence Data/Anopheles Stephensi.csv", header = T)
asteph2 = read.csv("Vector Occurrence Data/AnophelesStephensi_Additional.csv", header = T)

# remove any duplicates (occurrences from same lat, lon)
asteph1 <- asteph1 %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)
asteph2 <- asteph2 %>% dplyr::distinct(latitude, longitude, .keep_all = TRUE)

asteph1 = asteph1[!is.na(asteph1$decimalLatitude),]
asteph1 = asteph1[!is.na(asteph1$decimalLongitude),]
asteph2 = asteph2[!is.na(asteph2$latitude),]
asteph2 = asteph2[!is.na(asteph2$longitude),]

# lat range: 5 to 35

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

# Doing first for additional points (i.e for locations clearly missing from GBIF data)
asteph2 = subset(asteph2, latitude >=25 & latitude <= 35)
asteph2$LongitudeF <- as.factor(as.integer(asteph2$longitude))
asteph2 <- asteph2  %>% group_by(LongitudeF) %>% slice_sample(n=1)

# Repeat for GBIF data
asteph1 = subset(asteph1, decimalLatitude >=25 & decimalLatitude <= 35)
asteph1$LongitudeF <- as.factor(as.integer(asteph1$decimalLongitude))
asteph1 <- asteph1  %>% group_by(LongitudeF) %>% slice_sample(n=3)


# Identify elevation for each occurrence point
asteph_sites <- st_as_sf(asteph1, coords = c("decimalLongitude", "decimalLatitude"), 
                         crs = 4326, agr = "constant")
asteph <- get_elev_point(asteph_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
asteph <- asteph %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                                   Latitude = sf::st_coordinates(.)[,2])

# combined manually out of R

# output this list
#fwrite(asteph, "Vector Occurrence Data/AnSteph_SubSampledPoints.csv")


###### Culex quinquefasciatus ######

cxq = read.csv("Vector Occurrence Data/Culex quinquefasciatus.csv")
cxq = cxq[!is.na(cxq$decimalLatitude),]
cxq = cxq[!is.na(cxq$decimalLongitude),]

# Latitudinal range: -44 to 47

# remove any duplicates (occurrences from same lat, lon)
cxq = cxq %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

cxq = subset(cxq, decimalLatitude >=45 & decimalLatitude <= 47)
cxq$LongitudeF <- as.factor(as.integer(cxq$decimalLongitude))
cxq <- cxq %>% group_by(LongitudeF) %>% slice_sample(n=1)

# if still too many:
odds = seq(from = 1, to = nrow(cxq), by = 2)
cxq <- cxq[odds,]

# Identify elevation for each occurrence point
cxq_sites <- st_as_sf(cxq, coords = c("decimalLongitude", "decimalLatitude"), 
                      crs = 4326, agr = "constant")
cxq <- get_elev_point(cxq_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
cxq <- cxq %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                             Latitude = sf::st_coordinates(.)[,2])

# output this list
#fwrite(cxq, "Vector Occurrence Data/CxQuinque_SubSampledPoints.csv")

###### Culex pipiens ######

cxp = read.csv("Vector Occurrence Data/CulexPipiens.csv")
cxp = cxp[!is.na(cxp$decimalLatitude),]
cxp = cxp[!is.na(cxp$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxp = cxp %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: -39 to 67

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

cxp = subset(cxp, decimalLatitude >=65 & decimalLatitude <= 75)
cxp$LongitudeF <- as.factor(as.integer(cxp$decimalLongitude))
cxp <- cxp %>% group_by(LongitudeF) %>% slice_sample(n=1)

# if still too many:
odds = seq(from = 1, to = nrow(cxp), by = 2)
cxp <- cxp[odds,]

# Identify elevation for each occurrence point
cxp_sites <- st_as_sf(cxp, coords = c("decimalLongitude", "decimalLatitude"), 
                      crs = 4326, agr = "constant")
cxp <- get_elev_point(cxp_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
cxp <- cxp %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                             Latitude = sf::st_coordinates(.)[,2])


# output this list
#fwrite(cxp, "Vector Occurrence Data/CxPipiens_SubSampledPoints.csv")

###### Culex tarsalis ######

cxt = read.delim("Vector Occurrence Data/CulexTarsalis.csv")
cxt = cxt[!is.na(cxt$decimalLatitude),]
cxt = cxt[!is.na(cxt$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxt = cxt %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: 18 to 66

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

cxt = subset(cxt, decimalLatitude >=55 & decimalLatitude <= 66)
cxt$LongitudeF <- as.factor(as.integer(cxt$decimalLongitude))
cxt <- cxt %>% group_by(LongitudeF) %>% slice_sample(n=1)

# if still too many:
odds = seq(from = 1, to = nrow(cxt), by = 2)
cxt <- cxt[odds,]

# Identify elevation for each occurrence point
cxt_sites <- st_as_sf(cxt, coords = c("decimalLongitude", "decimalLatitude"), 
                      crs = 4326, agr = "constant")
cxt <- get_elev_point(cxt_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
cxt <- cxt %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                             Latitude = sf::st_coordinates(.)[,2])


# output this list
#fwrite(cxt, "Vector Occurrence Data/CxTarsalis_SubSampledPoints.csv")


###### Culex annulirostris ######

cxa = read.csv("Vector Occurrence Data/CulexAnnulirostris.csv")
cxa = cxa[!is.na(cxa$decimalLatitude),]
cxa = cxa[!is.na(cxa$decimalLongitude),]

# remove any duplicates (occurrences from same lat, lon)
cxa = cxa %>% dplyr::distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

# Latitudinal range: -39 to 16

# Apply a gridded sampling approach: 
# First, subset to a single 10 degree latitudinal band
# Then take random sample, capturing at least one from each degree longitude
# for which there is an occurrence record

cxa = subset(cxa, decimalLatitude>=5 & decimalLatitude<= 16)
cxa$LongitudeF <- as.factor(as.integer(cxa$decimalLongitude))
cxa <- cxa %>% group_by(LongitudeF) %>% slice_sample(n=2)

# if still too many:
odds = seq(from = 1, to = nrow(cxa), by = 2)
cxa <- cxa[odds,]

# Identify elevation for each occurrence point
cxa_sites <- st_as_sf(cxa, coords = c("decimalLongitude", "decimalLatitude"), 
                      crs = 4326, agr = "constant")
cxa <- get_elev_point(cxa_sites, src = "aws", overwrite = T)

# Separate out latitude and longitude
cxa <- cxa %>% dplyr::mutate(Longitude = sf::st_coordinates(.)[,1],
                             Latitude = sf::st_coordinates(.)[,2])


# output this list
#fwrite(cxa, "Vector Occurrence Data/CxAnnul_SubSampledPoints.csv")


############ [END OF: ONLY RUN THIS PART ONCE] ##############


#### 1. Plot Occurrence Points for each species #####

##### 1a. Aedes aegypti ######

#  Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aeaegypti_sites <- st_as_sf(aeaegypti, coords = c("lon", "lat"), 
                            crs = 4326, agr = "constant")
ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aeaegypti_sites, color = "#f46d43",  size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") 


##### 1b. Aedes albopictus ######

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aealbo_sites <- st_as_sf(aealbo, coords = c("lon", "lat"), 
                         crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() + labs(y = "") +
  geom_sf(data = aealbo_sites, 
          size = 3, shape = 16, col = "darkblue", alpha = 0.7) +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") 

##### 1f. Anopheles gambiae #####

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
agambi_sites <- st_as_sf(agambi, coords = c("lon", "lat"), 
                         crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = agambi_sites, 
          size = 3, shape = 16, col = "darkgreen", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-50, 70), ylim = c(-50, 30), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") 

##### 1g. Anopheles stephensi ######

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
asteph_sites <- st_as_sf(asteph, coords = c("lon", "lat"), 
                         crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = asteph_sites, 
          size = 3, shape = 16, col = "darkorange", alpha = 0.7) +
  coord_sf(xlim = c(10, 120), ylim = c(-20, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") 

##### 1h. Culex quinquefasciatus #####

# Plot occurrence records used in analysis
world <- ne_countries(scale = "medium", returnclass = "sf")
cxq_sites <- st_as_sf(cxq, coords = c("lon", "lat"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxq_sites, 
          size = 3, shape = 16, col = "purple", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") 

##### 1i. Culex pipiens ######

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxp_sites <- st_as_sf(cxpip, coords = c("lon", "lat"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxp_sites, 
          size = 3, shape = 16, col = "pink", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 80), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 60, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 70, lwd = 0.1, color = "black") 

##### 1j. Culex tarsalis #####

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxt_sites <- st_as_sf(cxtar, coords = c("lon", "lat"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxt_sites, 
          size = 3, shape = 16, col = "pink", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 80), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 60, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 70, lwd = 0.1, color = "black") 

##### 1l. Culex annulirostris #####

# Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxa_sites <- st_as_sf(cxa, coords = c("lon", "lat"), 
                      crs = 4326, agr = "constant")
ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxa_sites, 
          size = 3, shape = 16, col = "pink", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(80, 180), ylim = c(-50, 20), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 


#### 2. Plot Occurrence Points Colored by TSMs for each species #####
##### 2a. Aedes aegypti #####

# Note: aeaegypti_sites created in 1a

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aeaegypti_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
AeAegyptiColors = rgb(cols, maxColorValue=256)
world <- ne_countries(scale = "medium", returnclass = "sf")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aeaegypti_sites, color = AeAegyptiColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-50, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") 

##### 2b. Aedes albopictus #####

# note: aealbo_sites created in 1b
colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aealbo_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
AeAlboColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aealbo_sites, color = AeAlboColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-40, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") 

##### 2f. Anopheles gambiae #####

# note: agambi_sites created in 1f.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(agambi_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
AnGambiColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = agambi_sites, color = AnGambiColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-20, 60), ylim = c(-40, 30), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") 

##### 2g. Anopheles stephensi ######

# note: asteph_sites was created in 1g.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(asteph_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
AnStephColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = asteph_sites, color = AnStephColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(30, 105), ylim = c(-10, 40), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") 

##### 2h. Culex quinquefasciatus #####

# Note: cxq_sites was crated in 1h.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxq_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
CxQColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  theme(legend.position = c(.95, .95),
        legend.justification = c("right", "top"),
        legend.box.just = "right",
        legend.margin = margin(6, 6, 6, 6)) + 
  geom_sf(data = cxq_sites, color = CxQColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-170, 180), ylim = c(-50, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") 

##### 2i. Culex pipiens ######

# Note: cxp_sites created in 1i.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxp_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
CxPColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxp_sites, color = CxPColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-170, 180), ylim = c(-50, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 


##### 2j. Culex tarsalis #####

# Note: cxt_sites was created in 1j.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxt_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
CxTColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxt_sites, color = CxTColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-125, -75), ylim = c(15, 55), expand = FALSE) +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 

##### 2l. Culex annulirostris #####

# note: cxa_sites created in 1l.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxa_sites, (tolerance_point2 - min(tolerance_point2)) / diff(range(tolerance_point2)))
cols <- colorscale(colorbreaks)
CxAnnulColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = st_shift_longitude(cxa_sites), color = CxAnnulColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(110, 220), ylim = c(-50, 10), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black")

##### Create combined occurrence plots ##### 
###### Combined Aedes occurrence plots #####

# Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
aeaegypti_sites <- st_as_sf(aaegypti, coords = c("lon", "lat"), crs = 4326, agr = "constant")
aealbo_sites <- st_as_sf(aealbo, coords = c("lon", "lat"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aeaegypti_sites, size = 2, shape = 16, col = "#9b0000", alpha = 0.8) +
  geom_sf(data = aealbo_sites, size = 2, shape = 15, col = "#ed2939", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 13),axis.text.x = element_blank()) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") 

##### Combined Anopheles occurrence plots #####

# Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
agambi_sites <- st_as_sf(agambi, coords = c("lon", "lat"), crs = 4326, agr = "constant")
asteph_sites <- st_as_sf(asteph, coords = c("lon", "lat"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = agambi_sites,  size = 3, shape = 16, col = "#cc5801", alpha = 0.8) +
  geom_sf(data = asteph_sites, size = 3, shape = 15, col = "#f4bb00", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 13),axis.text.x = element_blank()) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-40, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 


##### Combined Culex occurrence plots ####

# Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
cxa_sites <- st_as_sf(cxa, coords = c("lon", "lat"), crs = 4326, agr = "constant")
cxp_sites <- st_as_sf(cxpip, coords = c("lon", "lat"), crs = 4326, agr = "constant")
cxq_sites <- st_as_sf(cxq, coords = c("lon", "lat"), crs = 4326, agr = "constant")
cxt_sites <- st_as_sf(cxtar, coords = c("lon", "lat"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxa_sites, size = 3, shape = 16, col = "#74add1", alpha = 0.8) +
  geom_sf(data = cxp_sites, size = 3, shape = 15, col = "#313695", alpha = 0.8) + 
  geom_sf(data = cxq_sites, size = 3, shape = 17, col = "#800080", alpha = 0.8) + 
  geom_sf(data = cxt_sites, size = 3, shape = 18, col = "#dda0dd", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 13),axis.text.x = element_blank()) + 
  coord_sf(xlim = c(-180, 180), ylim = c(-60, 70), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 60, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 70, lwd = 0.1, color = "black") 


##### All species occurrence plots ####

SpeciesColors = c("#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#309143", 
                  "#ab041b", "#ec3c30", 
                  "#abd9e9", "#74add1", "#313695", "#800080", "#dda0dd")
SpeciesList  = c("Aedes_aegypti", "Aedes_albopictus", "Aedes_camptorhynchus", "Aedes_triseriaatus", "Aedes_vexans",
                 "Anopheles_gambiae", "Anopheles_stephensi", 
                 "Culex_annulirostris", "Culex_pipiens", "Culex_quinquefasciatus", "Culex_tarsalis", "Culex_theileri")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aeaegypti_sites, size = 3.5, shape = 16, col = "#f46d43", alpha = 0.8) +
  geom_sf(data = aealbo_sites, size = 3.5, shape = 16, col = "#fdae61", alpha = 0.8) +
  geom_sf(data = aec_sites, size = 3.5, shape = 16, col = "#fed439ff", alpha = 0.8) +
  geom_sf(data = aet_sites, size = 3.5, shape = 16, col = "#7cae00", alpha = 0.8) +
  geom_sf(data = aev_sites, size = 3.5, shape = 16, col = "#309143", alpha = 0.8) +
  geom_sf(data = agambi_sites, size = 3.5, shape = 15, col = "#ab041b", alpha = 0.8) +
  geom_sf(data = asteph_sites, size = 3.5, shape = 15, col = "#ec3c30", alpha = 0.8) + 
  geom_sf(data = cxq_sites, size = 3.5, shape = 17, col = "#313695", alpha = 0.8) +
  geom_sf(data = cxp_sites, size = 3.5, shape = 17, col = "#74add1", alpha = 0.8) + 
  geom_sf(data = cxt_sites, size = 3.5, shape = 17, col = "#800080", alpha = 0.8) + 
  geom_sf(data = cxth_sites, size = 3.5, shape = 17, col = "#dda0dd", alpha = 0.8) + 
  labs(y = "") + 
  theme(axis.text.y = element_text(size = 16),axis.text.x = element_blank()) + 
  coord_sf(xlim = c(-100, 100), ylim = c(-90, 90), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.1, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.05, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.05, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.05, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.05, color = "black") +
  geom_hline(yintercept = -50, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 60, lwd = 0.05, color = "black") +
  geom_hline(yintercept = 70, lwd = 0.05, color = "black") 
