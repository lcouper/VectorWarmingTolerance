#### Vector thermal safety: Occurrence points ####

# This script is used to plot occurrence records used in the analysis

#### 0a. Load libraries ####

# For GAMS
library(mgcv)
library(scales)
library(elevatr)
library(tidyverse)
library(ggpmisc)
# For plotting
library(maps)
library(mapdata)
library(maptools)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)

#### 0b. Load colors and data frames #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files with seasonality
AeAegypti = read.csv("AeAegypti_TSM.csv")
AeAlbopictus = read.csv("AeAlbopictus_TSM.csv")
AeCamp = read.csv("AeCamp_TSM.csv")
AeTri = read.csv("AeTri_TSM.csv")
AeVexans = read.csv("AeVexans_TSM.csv")
AnGambiae = read.csv("AnGambiae_TSM.csv")
AnSteph = read.csv("AnSteph_TSM.csv")
CxQuinque = read.csv("CxQuinque_TSM.csv")  
CxPip = read.csv("CxPipiens_TSM.csv")
CxTar = read.csv("CxTarsalis_TSM.csv")
CxTh = read.csv("CxTheileri_TSM.csv")
CxAnnul = read.csv("CxAnnul_TSM.csv")

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


#### 1. Plot Occurrence Points for each species #####

##### 1a. Aedes aegypti ######

# remove any samples without elevation data
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aeaegypti_sites <- st_as_sf(AeAegypti, coords = c("longitude", "latitude"), 
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
##### 1b. Aedes albopictus #####

# remove any samples without elevation data
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aealbo_sites <- st_as_sf(AeAlbopictus, coords = c("longitude", "latitude"), 
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

##### 1c. Aedes camptorhynchus ######

# remove any samples without elevation data
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aec_sites <- st_as_sf(AeCamp, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aec_sites, 
          size = 3, shape = 16, col = "#fed439ff", alpha = 0.7) +
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
##### 1d. Aedes triseriatus #####

# remove any samples without elevation data
AeTri = AeTri[!is.na(AeTri$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aet_sites <- st_as_sf(AeTri, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aet_sites, 
          size = 3, shape = 16, col = "#7cae00", alpha = 0.7) +
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

##### 1e. Aedes vexans #####

# remove any samples without elevation data
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
aev_sites <- st_as_sf(AeVexans, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aev_sites, 
          size = 3, shape = 16, col = "#7cae00", alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-50, 60), expand = FALSE) +
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

##### 1f. Anopheles gambiae #####

# remove any samples without elevation data
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
agambi_sites <- st_as_sf(AnGambiae, coords = c("longitude", "latitude"), 
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

# remove any samples without elevation data
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
asteph_sites <- st_as_sf(AnSteph, coords = c("longitude", "latitude"), 
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

# remove any samples without elevation data
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]

# 1. Plot occurrence records used in analysis
world <- ne_countries(scale = "medium", returnclass = "sf")
cxq_sites <- st_as_sf(CxQuinque, coords = c("longitude", "latitude"), 
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

# remove any samples without elevation data
CxPip = CxPip[!is.na(CxPip$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxp_sites <- st_as_sf(CxPip, coords = c("longitude", "latitude"), 
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

# remove any samples without elevation data
CxTar = CxTar[!is.na(CxTar$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxt_sites <- st_as_sf(CxTar, coords = c("longitude", "latitude"), 
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

##### 1k. Culex theileri #####

CxTh = CxTh[!is.na(CxTh$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxth_sites <- st_as_sf(CxTh, coords = c("longitude", "latitude"), 
                       crs = 4326, agr = "constant")
ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxth_sites, 
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

CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]

# 1. Plot occurrence points used in analysis 
world <- ne_countries(scale = "medium", returnclass = "sf")
cxa_sites <- st_as_sf(CxAnnul, coords = c("longitude", "latitude"), 
                      crs = 4326, agr = "constant")
ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxa_sites, 
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
#### 2. Plot Occurrence Points Colored by TSMs for each species #####
##### 2a. Aedes aegypti #####

# Note: aeaegypti_sites created in 1a

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aeaegypti_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
AeAegyptiColors = rgb(cols, maxColorValue=256)

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
colorbreaks <- with(aealbo_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
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

##### 2c. Aedes camptorhynchus #####

# note: aec_sites created in 1c.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aec_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
AeCampColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aec_sites, color = AeCampColors , 
          size = 3, shape = 16, alpha = 0.7) +
  #  geom_sf(data = aealbo_sites, size = 3, shape = 16, col = "darkblue", alpha = 0.7) + 
  labs(y = "") +
  coord_sf(xlim = c(100, 160), ylim = c(-50,-10), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black") 

##### 2d. Aedes triseriatus ######

# Note: aet_sites was created in 1d
colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aet_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
AeTriColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aet_sites, color = AeTriColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-120, -60), ylim = c(20, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black")

##### 2e. Aedes vexans #####

# Note aev_sites created in 1e.
colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(aev_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
AeVexansColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aev_sites, color = AeVexansColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-180, 180), ylim = c(-30, 60), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 

##### 2f. Anopheles gambiae #####

# note: agambi_sites created in 1f.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(agambi_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
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
colorbreaks <- with(asteph_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
AnStephColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = asteph_sites, color = AnStephColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(20, 110), ylim = c(-10, 40), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") 

##### 2h. Culex quinquefasciatus #####

# Note: cxq_sites was crated in 1h.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxq_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
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
colorbreaks <- with(cxp_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
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
colorbreaks <- with(cxt_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
CxTColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxt_sites, color = CxTColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-130, -70), ylim = c(15, 55), expand = FALSE) +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 50, lwd = 0.1, color = "black") 

##### 2k. Culex theileri #####

# Note: cxth_sites created in 1k.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxth_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
CxThColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxth_sites, color = CxThColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(-30, 100), ylim = c(-40, 50), expand = FALSE) +
  geom_hline(yintercept = -10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 0, lwd = 0.2, color = "black") + 
  geom_hline(yintercept = 10, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -20, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -30, lwd = 0.1, color = "black") +
  geom_hline(yintercept = 40, lwd = 0.1, color = "black") +
  geom_hline(yintercept = -40, lwd = 0.1, color = "black")

##### 2l. Culex annulirostris #####

# note: cxa_sites created in 1l.

colorscale <- colorRamp(c("#b31818","#fdf4f4"))
colorbreaks <- with(cxa_sites, (tolerance_point - min(tolerance_point)) / diff(range(tolerance_point)))
cols <- colorscale(colorbreaks)
CxAnnulColors = rgb(cols, maxColorValue=256)

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxa_sites, color = CxAnnulColors , 
          size = 3, shape = 16, alpha = 0.7) +
  labs(y = "") +
  coord_sf(xlim = c(100, 180), ylim = c(-50, 20), expand = FALSE) +
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


Aedes = rbind.data.frame(AeAegypti, AeAlbo, AeCamp, AeTri, AeVex)
colors = c("#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#309143")
AedesColors = rep(colors, each = 80)

# 1. Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
aeaegypti_sites <- st_as_sf(AeAegypti, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
aealbo_sites <- st_as_sf(AeAlbo, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
aec_sites <- st_as_sf(AeCamp, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
aet_sites <- st_as_sf(AeTri, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
aev_sites <- st_as_sf(AeVex, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = aeaegypti_sites, size = 3, shape = 16, col = "#f46d43", alpha = 0.8) +
  geom_sf(data = aealbo_sites, size = 3, shape = 16, col = "#fdae61", alpha = 0.8) + 
  geom_sf(data = aec_sites, size = 3, shape = 16, col = "#fed439ff", alpha = 0.8) + 
  geom_sf(data = aet_sites, size = 3, shape = 16, col = "#7cae00", alpha = 0.8) + 
  geom_sf(data = aev_sites, size = 3, shape = 16, col = "#309143", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 16),axis.text.x = element_blank()) + 
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
  geom_hline(yintercept = -50, lwd = 0.1, color = "black") 

##### Combined Anopheles occurrence plots #####

Anopheles = rbind.data.frame(AnGam, AnSteph)
colorsAnopheles = c(rep("#ab041b", 80), rep("#ec3c30", 80))

# 1. Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
agambi_sites <- st_as_sf(AnGam, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
asteph_sites <- st_as_sf(AnSteph, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = agambi_sites,  size = 3, shape = 16, col = "#ab041b", alpha = 0.8) +
  geom_sf(data = asteph_sites, size = 3, shape = 16, col = "#ec3c30", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 16),axis.text.x = element_blank()) + 
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

Culex = rbind.data.frame(CxA, CxPip, CxQ, CxTar, CxTh)
colors = c("#abd9e9", "#74add1", "#313695", "#800080", "#dda0dd")
colorsCulex = rep(colors, each = 80)

# 1. Plot combined occurrence records
world <- ne_countries(scale = "medium", returnclass = "sf")
cxa_sites <- st_as_sf(CxA, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
cxp_sites <- st_as_sf(CxPip, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
cxq_sites <- st_as_sf(CxQ, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
cxt_sites <- st_as_sf(CxTar, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")
cxth_sites <- st_as_sf(CxTh, coords = c("longitude", "latitude"), crs = 4326, agr = "constant")

ggplot(data = world) +
  geom_sf() + theme_bw() +
  geom_sf(data = cxa_sites, size = 3, shape = 16, col = "#abd9e9", alpha = 0.8) +
  geom_sf(data = cxp_sites, size = 3, shape = 16, col = "#74add1", alpha = 0.8) + 
  geom_sf(data = cxq_sites, size = 3, shape = 16, col = "#313695", alpha = 0.8) + 
  geom_sf(data = cxt_sites, size = 3, shape = 16, col = "#800080", alpha = 0.8) + 
  geom_sf(data = cxth_sites, size = 3, shape = 16, col = "#dda0dd", alpha = 0.8) + 
  labs(y = "") +
  theme(axis.text.y = element_text(size = 16),axis.text.x = element_blank()) + 
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
