#### Vector Thermal Safety GAMs #####

# This script is used to model and plot results for each vector species including:
# 1. Run GAM of TSM ~ latitude
# 2. Run GAM of TSM ~ latitude (high v low elevation)
# 3. Run GAM of TSM ~ latitude (with and w/o behavior)
# 4. Run GAM of TSM ~ latitude (with and w/o drought mask)

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


#### 1. GAM of TSM across latitude for each species ####

##### 1a. Aedes aegypti #####

# First remove high elevation samples (>3000 ft)
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegypti = AeAegypti[AeAegypti$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAg = gam(tolerance_point ~ s(latitude, k = 8), data = AeAegypti, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti))

# Predict using GAM model
y = predict(gamAg, newdata = data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam 
ggplot() +  theme_minimal() +
  geom_point(aes(x = AeAegypti$latitude, y = AeAegypti$tolerance_point), size = 1.2, color = SpeciesColors[1]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[1], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[1], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[1], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes aegypti") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 1b. Aedes albopictus #####

# First remove any high elevation samples
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbopictus = AeAlbopictus[AeAlbopictus$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAb = gam(tolerance_point ~ s(latitude, k = 8), data = AeAlbopictus, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAlbopictus$latitude), to = max(AeAlbopictus$latitude), length.out = nrow(AeAlbopictus)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamAb, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AeAlbopictus, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[2]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[2], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[2], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[2], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes albopictus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 1c. Aedes camptorhynchus #####

# First remove any high elevation samples
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]
AeCamp = AeCamp[AeCamp$elevation_ft < 3000,]
gamAc = gam(tolerance_point ~ s(latitude, k = 8), data = AeCamp, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeCamp$latitude), to = max(AeCamp$latitude), length.out = nrow(AeCamp)) 

# Predict using GAM model
y = predict(gamAc, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AeCamp, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[7]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[7], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[7], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[7], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes camptorhynchus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1d. Aedes triseriatus #####

# First remove high elevation samples (>3000 ft)
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTri = AeTri[AeTri$elevation_ft < 3000,]
gamAt = gam(tolerance_point ~ s(latitude, k = 8), data = AeTri, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeTri$latitude), to = max(AeTri$latitude), length.out = nrow(AeTri)) 

# Predict using GAM model
y = predict(gamAt, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AeTri, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[9]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[9], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[9], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[9], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes triseriatus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1e. Aedes vexans #####

# First remove high elevation samples (>3000 ft)
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexans = AeVexans[AeVexans$elevation_ft < 3000,] # 16 high-elevation samples. 
gamAv = gam(tolerance_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamAv, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AeVexans, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[11]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[11], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[11], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[11], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes vexans") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1f. Anopheles gambiae #####

# First remove high elevation samples (>3000 ft)
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiae = AnGambiae[AnGambiae$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAm = gam(tolerance_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamAm, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AnGambiae, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[3]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[3], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[3], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[3], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles gambiae") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 1g. Anopheles stephensi #####

# First remove high elevation samples (>3000 ft)
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnSteph= AnSteph[AnSteph$elevation_ft < 3000,] # 9 high-elevation samples. 

# Run GAM estimating thermal safety margins across latitude
gamAs = gam(tolerance_point ~ s(latitude, k = 8), data = AnSteph, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph))

# Predict using GAM model
y = predict(gamAs, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnSteph, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[4]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[4], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[4], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[4], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles stephensi") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1h. Culex quinquefasciatus ####

# First remove high elevation samples (>3000 ft)
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinque = CxQuinque[CxQuinque$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCq = gam(tolerance_point ~ s(latitude, k = 8), data = CxQuinque, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque)) 

# Predict using GAM model
y = predict(gamCq, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(CxQuinque, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[6]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[6], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[6], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[6], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex quinquefasciatus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1i. Culex pipiens ####

# First separate any samples from high elevation
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPip = CxPip[CxPip$elevation_ft < 3000,] # 30 high-elevation samples 

# First remove high elevation samples (>3000 ft)
gamCp = gam(tolerance_point ~ s(latitude, k = 8), data = CxPip, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip)) 

# Predict using GAM model
y = predict(gamCp, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxPip, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[5]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[5], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[5], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[5], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex pipiens") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 1j. Culex tarsalis ####

# First remove high elevation samples (>3000 ft)
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTar = CxTar[CxTar$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCt = gam(tolerance_point ~ s(latitude, k = 8), data = CxTar, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamCt, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam 
ggplot(CxTar, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[10]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[10], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[10], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[10], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex tarsalis") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1k. Culex theileri ####

# First remove high elevation samples (>3000 ft)
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxTh = CxTh[CxTh$elevation_ft < 3000,] # 20 high-elevation samples. 

# Run GAM estimating thermal safety margins across latitude
gamCth = gam(tolerance_point ~ s(latitude, k = 8), data = CxTh, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh)) 

# Predict using GAM model
y = predict(gamCth, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(CxTh, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[12]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[12], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[12], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[12], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex theileri") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1l. Culex annulirostris ####

# First remove high elevation samples (>3000 ft)
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnul = CxAnnul[CxAnnul$elevation_ft < 3000,] # remove high-elevation sample

# Run GAM estimating thermal safety margins across latitude
gamCa = gam(tolerance_point ~ s(latitude, k = 8), data = CxAnnul, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxAnnul$latitude), to = max(CxAnnul$latitude), length.out = nrow(CxAnnul)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamCa, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(CxAnnul, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = SpeciesColors[8]) +
  geom_line(aes(x=newlats, y = y$fit), colour=SpeciesColors[8], lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha(SpeciesColors[8], 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha(SpeciesColors[8], 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Cx annulirostris") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




#### 2. GAM of TSM across latitude: high and low elevation records ####
##### 2a. Aedes aegypti #####

# Bring back previously removed high elevation samples
AeAegypti = read.csv("AeAegypti_TSM.csv")

# Classify records as "high" or "low" elevation
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegyptiHigh = AeAegypti[AeAegypti$elevation_ft >= 3000,] # for coloring points on plot below
AeAegypti$ElevClass = NA
AeAegypti$ElevClass[AeAegypti$elevation_ft >= 3000] <- "High"
AeAegypti$ElevClass[AeAegypti$elevation_ft < 3000] <- "Low"
AeAegypti$ElevClass = as.factor(AeAegypti$ElevClass)

# run gam, using a random intercept for elevational class
gamAg = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = AeAegypti, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamAg, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = AeAegypti$latitude, y = AeAegypti$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(AeAegyptiHigh$latitude, y = AeAegyptiHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:114]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:114]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:114]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[115:228]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[115:228]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[115:228]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes aegypti") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 






##### 2b. Aedes albopictus #####

# Bring back previously removed high elevation samples
AeAlbopictus = read.csv("AeAlbopictus_TSM.csv")

# Classify records as "high" or "low" elevation
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbopictusHigh = AeAlbopictus[AeAlbopictus$elevation_ft >= 3000,] # for coloring points on plot below
AeAlbopictus$ElevClass = NA
AeAlbopictus$ElevClass[AeAlbopictus$elevation_ft >= 3000] <- "High"
AeAlbopictus$ElevClass[AeAlbopictus$elevation_ft < 3000] <- "Low"
AeAlbopictus$ElevClass = as.factor(AeAlbopictus$ElevClass)

# run gam, using a random intercept for elevational class
gamAb = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = AeAlbopictus, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAlbopictus$latitude), to = max(AeAlbopictus$latitude), length.out = nrow(AeAlbopictus)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamAb, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = AeAlbopictus$latitude, y = AeAlbopictus$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(AeAlbopictusHigh$latitude, y = AeAlbopictusHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:114]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:114]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:114]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[115:228]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[115:228]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[115:228]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes albopictus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 











##### 2c. Aedes camptorhynchus (skip) #####

# Bring back previously removed high elevation samples
AeCamp = read.csv("AeCamp_TSM.csv")

# Classify records as "high" or "low" elevation
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]
AeCampHigh = AeCamp[AeCamp$elevation_ft >= 3000,] # for coloring points on plot below

## Only 2 high-elevation records: Not enough to read GAM



##### 2d. Aedes triseriatus (skip) #####

# Bring back previously removed high elevation samples
AeTri = read.csv("AeTri_TSM.csv")

# Classify records as "high" or "low" elevation
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTriHigh = AeTri[AeTri$elevation_ft >= 3000,] 

## 0 high-elevation records




##### 2e. Aedes vexans #####

# Bring back previously removed high elevation samples
AeVexans = read.csv("AeVexans_TSM.csv")

# Classify records as "high" or "low" elevation
AeVexans  = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexansHigh = AeVexans[AeVexans$elevation_ft >= 3000,] # for coloring points on plot below
AeVexans$ElevClass = NA
AeVexans$ElevClass[AeVexans$elevation_ft >= 3000] <- "High"
AeVexans$ElevClass[AeVexans$elevation_ft < 3000] <- "Low"
AeVexans$ElevClass = as.factor(AeVexans$ElevClass)

# run gam, using a random intercept for elevational class
gamAv = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = AeVexans, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamAv, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = AeVexans$latitude, y = AeVexans$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(AeVexansHigh$latitude, y = AeVexansHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:84]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:84]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:84]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[85:168]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[85:168]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[85:168]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes vexans") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 



##### 2f. Anopheles gambiae (no effect of elevation) #####

# Bring back previously removed high elevation samples
AnGambiae = read.csv("AnGambiae_TSM.csv")

# Classify records as "high" or "low" elevation
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiaeHigh = AnGambiae[AnGambiae$elevation_ft >= 3000,] # for coloring points on plot below
AnGambiae$ElevClass = NA
AnGambiae$ElevClass[AnGambiae$elevation_ft >= 3000] <- "High"
AnGambiae$ElevClass[AnGambiae$elevation_ft < 3000] <- "Low"
AnGambiae$ElevClass = as.factor(AnGambiae$ElevClass)

# run gam, using a random intercept for elevational class
gamAm = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = AnGambiae, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamAm, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = AnGambiae$latitude, y = AnGambiae$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(AnGambiaeHigh$latitude, y = AnGambiaeHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:91]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:91]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:91]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[92:182]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[92:182]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[92:182]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles gambiae") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 










##### 2g. Anopheles stephensi #####

# Bring back previously removed high elevation samples
AnSteph = read.csv("AnSteph_TSM.csv")

# Classify records as "high" or "low" elevation
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnStephHigh = AnSteph[AnSteph$elevation_ft >= 3000,] # for coloring points on plot below
AnSteph$ElevClass = NA
AnSteph$ElevClass[AnSteph$elevation_ft >= 3000] <- "High"
AnSteph$ElevClass[AnSteph$elevation_ft < 3000] <- "Low"
AnSteph$ElevClass = as.factor(AnSteph$ElevClass)

# run gam, using a random intercept for elevational class
gamAs = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = AnSteph, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamAs, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = AnSteph$latitude, y = AnSteph$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(AnStephHigh$latitude, y = AnStephHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:52]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:52]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:52]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[53:104]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[53:104]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[53:104]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles stephensi") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 



##### 2h. Culex quinquefasciatus #####

# Bring back previously removed high elevation samples
CxQuinque = read.csv("CxQuinque_TSM.csv")

# Classify records as "high" or "low" elevation
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinqueHigh = CxQuinque[CxQuinque$elevation_ft >= 3000,] # for coloring points on plot below
CxQuinque$ElevClass = NA
CxQuinque$ElevClass[CxQuinque$elevation_ft >= 3000] <- "High"
CxQuinque$ElevClass[CxQuinque$elevation_ft < 3000] <- "Low"
CxQuinque$ElevClass = as.factor(CxQuinque$ElevClass)

# run gam, using a random intercept for elevational class
gamCq = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = CxQuinque, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamCq, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = CxQuinque$latitude, y = CxQuinque$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(CxQuinqueHigh$latitude, y = CxQuinqueHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:135]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:135]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:135]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[136:270]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[136:270]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[136:270]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex quinquefasciatus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 


##### 2i. Culex pipiens (no effect of elevation) #####

# Bring back previously removed high elevation samples
CxPip = read.csv("CxPipiens_TSM.csv")

# Classify records as "high" or "low" elevation
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPipHigh = CxPip[CxPip$elevation_ft >= 3000,] # for coloring points on plot below
CxPip$ElevClass = NA
CxPip$ElevClass[CxPip$elevation_ft >= 3000] <- "High"
CxPip$ElevClass[CxPip$elevation_ft < 3000] <- "Low"
CxPip$ElevClass = as.factor(CxPip$ElevClass)

# run gam, using a random intercept for elevational class
gamCp = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = CxPip, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamCp, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = CxPip$latitude, y = CxPip$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(CxPipHigh$latitude, y = CxPipHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:91]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:91]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:91]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[92:182]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[92:182]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[92:182]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex pipiens") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 







##### 2j. Culex tarsalis #####

# Bring back previously removed high elevation samples
CxTar = read.csv("CxTarsalis_TSM.csv")

# Classify records as "high" or "low" elevation
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTarHigh = CxTar[CxTar$elevation_ft >= 3000,] # for coloring points on plot below
CxTar$ElevClass = NA
CxTar$ElevClass[CxTar$elevation_ft >= 3000] <- "High"
CxTar$ElevClass[CxTar$elevation_ft < 3000] <- "Low"
CxTar$ElevClass = as.factor(CxTar$ElevClass)

# run gam, using a random intercept for elevational class
gamCt = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = CxTar, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamCt, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = CxTar$latitude, y = CxTar$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(CxTarHigh$latitude, y = CxTarHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:62]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:62]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:62]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[63:124]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[63:124]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[63:124]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex tarsalis") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 









##### 2k. Culex theileri #####

# Bring back previously removed high elevation samples
CxTh = read.csv("CxTheileri_TSM.csv")

# Classify records as "high" or "low" elevation
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxThHigh = CxTh[CxTh$elevation_ft >= 3000,] # for coloring points on plot below
CxTh$ElevClass = NA
CxTh$ElevClass[CxTh$elevation_ft >= 3000] <- "High"
CxTh$ElevClass[CxTh$elevation_ft < 3000] <- "Low"
CxTh$ElevClass = as.factor(CxTh$ElevClass)

# run gam, using a random intercept for elevational class
gamCth = gam(tolerance_point ~ s(latitude, k = 8) + s(ElevClass, bs = "re"), data = CxTh, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh)/2) 
newelevclass = factor(c("Low", "High"))
newdata = expand.grid(newlats, newelevclass)
colnames(newdata) = c("latitude", "ElevClass")

# Predict using GAM model
y = predict(gamCth, newdata = newdata, se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam (showing low and high elevation fits if applicable)
ggplot() +  theme_minimal() +
  geom_point(aes(x = CxTh$latitude, y = CxTh$tolerance_point), size = 1.2, color = "#026633") +
  geom_point(aes(CxThHigh$latitude, y = CxThHigh$tolerance_point), col = "#846954") + 
  geom_line(aes(x=newlats, y = y$fit[1:41]), colour="#026633", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[1:41]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[1:41]), colour=alpha("#026633", 1.0), lty = 2, lwd = 1.1) +
  geom_line(aes(x=newlats, y = y$fit[42:82]), colour="#846954", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr[42:82]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr[42:82]), colour=alpha("#846954", 1.0), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex theileri") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0)) 








##### 2l. Culex annulirostris (skip) #####

# Bring back previously removed high elevation samples
CxAnnul = read.csv("CxAnnul_TSM.csv")

# Classify records as "high" or "low" elevation
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnulHigh = CxAnnul[CxAnnul$elevation_ft >= 3000,] # for coloring points on plot below

## Only 1 high-elevation record: Not enough to read GAM



#### 3. GAM of TSM across latitude: With and without behavior ####
##### 3a. Aedes aegypti #####

# First remove high elevation samples (>3000 ft)
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegypti = AeAegypti[AeAegypti$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAg = gam(tolerance_point ~ s(latitude, k = 8), data = AeAegypti, method = "REML")
gamAgNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeAegypti, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti)) 

# Predict using GAM model
y = predict(gamAg, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAgNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeAegypti, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes aegypti") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 3b. Aedes albopictus #####

# First remove high elevation samples (>3000 ft)
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbopictus  = AeAlbopictus[AeAlbopictus$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAb = gam(tolerance_point ~ s(latitude, k = 8), data = AeAlbopictus, method = "REML")
gamAbNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeAlbopictus, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAlbopictus$latitude), to = max(AeAlbopictus$latitude), length.out = nrow(AeAlbopictus)) 

# Predict using GAM model
y = predict(gamAb,data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAbNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeAlbopictus, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes albopictus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 3c. Aedes camptorhynchus #####

# First remove high elevation samples (>3000 ft)
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]
AeCamp = AeCamp[AeCamp$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAc = gam(tolerance_point ~ s(latitude, k = 8), data = AeCamp, method = "REML")
gamAcNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeCamp, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeCamp$latitude), to = max(AeCamp$latitude), length.out = nrow(AeCamp)) 

# Predict using GAM model
y = predict(gamAc, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAcNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeCamp, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes camptorhynchus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3d. Aedes triseriatus #####

# First remove high elevation samples (>3000 ft)
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTri = AeTri[AeTri$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAt = gam(tolerance_point ~ s(latitude, k = 8), data = AeTri, method = "REML")
gamAtNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeTri, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeTri$latitude), to = max(AeTri$latitude), length.out = nrow(AeTri)) 

# Predict using GAM model
y = predict(gamAt, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAtNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeTri, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes triseriatus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 3e. Aedes vexans #####

# First remove high elevation samples (>3000 ft)
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexans = AeVexans[AeVexans$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAv = gam(tolerance_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")
gamAvNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans)) 

# Predict using GAM model
y = predict(gamAv, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAvNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeVexans, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes vexans") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 3f. Anopheles gambiae #####

# First remove high elevation samples (>3000 ft)
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiae = AnGambiae[AnGambiae$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAm = gam(tolerance_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")
gamAmNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae)) 

# Predict using GAM model
y = predict(gamAm, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAmNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnGambiae, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles gambiae") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 3g. Anopheles stephensi #####

# First remove high elevation samples (>3000 ft)
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnSteph = AnSteph[AnSteph$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamAs = gam(tolerance_point ~ s(latitude, k = 8), data = AnSteph, method = "REML")
gamAsNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AnSteph, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph)) 

# Predict using GAM model
y = predict(gamAs, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamAsNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnSteph, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles stephensi") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 3h. Culex quinquefasciatus #####

# First remove high elevation samples (>3000 ft)
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinque = CxQuinque[CxQuinque$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCq = gam(tolerance_point ~ s(latitude, k = 8), data = CxQuinque, method = "REML")
gamCqNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxQuinque, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque)) 

# Predict using GAM model
y = predict(gamCq, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamCqNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxQuinque, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.1) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex quinquefasciatus") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3i. Culex pipiens #####

# First remove high elevation samples (>3000 ft)
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPip = CxPip[CxPip$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCp = gam(tolerance_point ~ s(latitude, k = 8), data = CxPip, method = "REML")
gamCpNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxPip, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip)) 

# Predict using GAM model
y = predict(gamCp, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamCpNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxPip, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex pipiens") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3j. Culex tarsalis #####

# First remove high elevation samples (>3000 ft)
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTar = CxTar[CxTar$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCt = gam(tolerance_point ~ s(latitude, k = 8), data = CxTar, method = "REML")
gamCtNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxTar, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar)) 

# Predict using GAM model
y = predict(gamCt, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamCtNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxTar, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex tarsalis") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3k. Culex theileri #####

# First remove high elevation samples (>3000 ft)
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxTh = CxTh[CxTh$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCth = gam(tolerance_point ~ s(latitude, k = 8), data = CxTh, method = "REML")
gamCthNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxTh, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh)) 

# Predict using GAM model
y = predict(gamCth, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamCthNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxTh, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex theileri") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 3l. Culex annulirostris #####

# First remove high elevation samples (>3000 ft)
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnul = CxAnnul[CxAnnul$elevation_ft < 3000,] 

# Run GAM estimating thermal safety margins across latitude
gamCa = gam(tolerance_point ~ s(latitude, k = 8), data = CxAnnul, method = "REML")
gamCaNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxAnnul, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxAnnul$latitude), to = max(CxAnnul$latitude), length.out = nrow(CxAnnul)) 

# Predict using GAM model
y = predict(gamCa, data.frame(latitude = newlats), se.fit = TRUE)
ynoB = predict(gamCaNoB, data.frame(latitude = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxAnnul, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_nob_point), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex annulirostris") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

#### 4. GAM of TSM across latitude: With and without drought mask ####
##### 4a. Aedes aegypti #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeAegyptiNS = read.csv("WithoutDroughtMask/AeAegypti_TSM_WithoutDroughtMask.csv")
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegypti = AeAegypti[AeAegypti$elevation_ft < 3000,] 
AeAegyptiNS = AeAegyptiNS[!is.na(AeAegyptiNS$elevation_ft),]
AeAegyptiNS = AeAegyptiNS[AeAegyptiNS$elevation_ft < 3000,] 

# with drought mask
gamAg = gam(tolerance_point ~ s(latitude, k = 8), data = AeAegypti, method = "REML")
newlats = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti)) 
newdata = data.frame(latitude = newlats)
y = predict(gamAg, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAgNS = gam(tolerance_point ~ s(latitude, k = 8), data = AeAegyptiNS, method = "REML")
yNS = predict(gamAgNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeAegypti, AeAegyptiNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes aegypti")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4b. Aedes albopictus  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeAlboNS = read.csv("WithoutDroughtMask/AeAlbopictus_TSM_WithoutDroughtMask.csv")
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbopictus = AeAlbopictus[AeAlbopictus$elevation_ft < 3000,] 
AeAlboNS = AeAlboNS[!is.na(AeAlboNS$elevation_ft),]
AeAlboNS = AeAlboNS[AeAlboNS$elevation_ft < 3000,] 

# with drought mask
gamAb = gam(tolerance_point ~ s(latitude, k = 8), data = AeAlbopictus, method = "REML")
newlats = seq(from = min(AeAlbopictus$latitude), to = max(AeAlbopictus$latitude), length.out = nrow(AeAlbopictus)) 
newdata = data.frame(latitude = newlats)
y = predict(gamAb, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAbNS = gam(tolerance_point ~ s(latitude, k = 8), data = AeAlboNS, method = "REML")
yNS = predict(gamAbNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeAlbopictus, AeAlboNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes albopictus")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4c. Aedes camptorhynchus #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeCampNS = read.csv("WithoutDroughtMask/AeCamp_TSM_WithoutDroughtMask.csv")
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]
AeCamp = AeCamp[AeCamp$elevation_ft < 3000,] 
AeCampNS = AeCampNS[!is.na(AeCampNS$elevation_ft),]
AeCampNS = AeCampNS[AeCampNS$elevation_ft < 3000,] 

# with drought mask
gamAc = gam(tolerance_point ~ s(latitude, k = 8), data = AeCamp, method = "REML")
newlats = seq(from = min(AeCamp$latitude), to = max(AeCamp$latitude), length.out = nrow(AeCamp))
newdata = data.frame(latitude = newlats)
y = predict(gamAc, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAcNS = gam(tolerance_point ~ s(latitude, k = 8), data = AeCampNS, method = "REML")
yNS = predict(gamAcNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeCamp, AeCampNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes camptorhynchus")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 4d. Aedes triseriatus  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeTriNS = read.csv("WithoutDroughtMask/AeTri_TSM_WithoutDroughtMask.csv")
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTri = AeTri[AeTri$elevation_ft < 3000,] 
AeTriNS = AeTriNS[!is.na(AeTriNS$elevation_ft),]
AeTriNS = AeTriNS[AeTriNS$elevation_ft < 3000,] 

# with drought mask
gamAt = gam(tolerance_point ~ s(latitude, k = 8), data = AeTri, method = "REML")
newlats = seq(from = min(AeTri$latitude), to = max(AeTri$latitude), length.out = nrow(AeTri))
newdata = data.frame(latitude = newlats)
y = predict(gamAt, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAtNS = gam(tolerance_point ~ s(latitude, k = 8), data = AeTriNS, method = "REML")
yNS = predict(gamAtNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeTri, AeTriNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes triseriatus")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4e. Aedes vexans  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeVexansNS = read.csv("WithoutDroughtMask/AeVexans_TSM_WithoutDroughtMask.csv")
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexans = AeVexans[AeVexans$elevation_ft < 3000,] 
AeVexansiNS = AeVexansNS[!is.na(AeVexansNS$elevation_ft),]
AeVexansNS = AeVexansNS[AeVexansNS$elevation_ft < 3000,] 

# with drought mask
gamAv = gam(tolerance_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")
gamAvNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")
newlats = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans))
newdata = data.frame(latitude = newlats)
y = predict(gamAv, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAvNS = gam(tolerance_point ~ s(latitude, k = 8), data = AeVexansNS, method = "REML")
yNS = predict(gamAvNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeVexans, AeVexansNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes vexans")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4f. Anopheles gambiae  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AnGambiaeNS = read.csv("WithoutDroughtMask/AnGambiae_TSM_WithoutDroughtMask.csv")
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiae = AnGambiae[AnGambiae$elevation_ft < 3000,] 
AnGambiaeNS = AnGambiaeNS[!is.na(AnGambiaeNS$elevation_ft),]
AnGambiaeNS = AnGambiaeNS[AnGambiaeNS$elevation_ft < 3000,]

# with drought mask
gamAm = gam(tolerance_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")
gamAmNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")
newlats = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae))
newdata = data.frame(latitude = newlats)
y = predict(gamAm, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAmNS = gam(tolerance_point ~ s(latitude, k = 8), data = AnGambiaeNS, method = "REML")
yNS = predict(gamAmNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AnGambiae, AnGambiaeNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Anopheles gambiae")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4g. Anopheles stephensi  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AnStephNS = read.csv("WithoutDroughtMask/AnSteph_TSM_WithoutDroughtMask.csv")
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnSteph = AnSteph[AnSteph$elevation_ft < 3000,]
AnStephNS = AnStephNS[!is.na(AnStephNS$elevation_ft),]
AnStephNS = AnStephNS[AnStephNS$elevation_ft < 3000,]

# with drought mask
gamAs = gam(tolerance_point ~ s(latitude, k = 8), data = AnSteph, method = "REML")
newlats = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph)) 
newelevs = seq(from = min(AnSteph$elevation_ft), to = max(AnSteph$elevation_ft), length.out = nrow(AnSteph))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamAs, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAsNS = gam(tolerance_point ~ s(latitude, k = 8), data = AnStephNS, method = "REML")
yNS = predict(gamAsNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AnSteph, AnStephNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Anopheles stephensi")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4h. Cx quinquefasciatus  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxQuinqueNS = read.csv("WithoutDroughtMask/CxQuinque_TSM_WithoutDroughtMask.csv")
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinque = CxQuinque[CxQuinque$elevation_ft < 3000,]
CxQuinqueNS = CxQuinqueNS[!is.na(CxQuinqueNS$elevation_ft),]
CxQuinqueNS = CxQuinqueNS[CxQuinqueNS$elevation_ft < 3000,]

# with drought mask
gamCq = gam(tolerance_point ~ s(latitude, k = 8), data = CxQuinque, method = "REML")
newlats = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque)) 
newelevs = seq(from = min(CxQuinque$elevation_ft), to = max(CxQuinque$elevation_ft), length.out = nrow(CxQuinque))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamCq, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCqNS = gam(tolerance_point ~ s(latitude, k = 8), data = CxQuinqueNS, method = "REML")
yNS = predict(gamCqNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxQuinque, CxQuinqueNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex quinquefasciatus")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4i. Cx pipiens  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxPipNS = read.csv("WithoutDroughtMask/CxPipiens_TSM_WithoutDroughtMask.csv")
CxPipNS = CxPipNS[!is.na(CxPipNS$elevation_ft),]
CxPipNS = CxPipNS[CxPipNS$elevation_ft < 3000,]
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPip = CxPip[CxPip$elevation_ft < 3000,]

# with drought mask
gamCp = gam(tolerance_point ~ s(latitude, k = 8), data = CxPip, method = "REML")
newlats = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip)) 
newelevs = seq(from = min(CxPip$elevation_ft), to = max(CxPip$elevation_ft), length.out = nrow(CxPip))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamCp, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCpNS = gam(tolerance_point ~ s(latitude, k = 8), data = CxPipNS, method = "REML")
yNS = predict(gamCpNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxPip, CxPipNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex pipiens")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4j. Cx tarsalis  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxTarNS = read.csv("WithoutDroughtMask/CxTarsalis_TSM_WithoutDroughtMask.csv")
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTar = CxTar[CxTar$elevation_ft < 3000,]
CxTarNS = CxTarNS[!is.na(CxTarNS$elevation_ft),]
CxTarNS = CxTarNS[CxTarNS$elevation_ft < 3000,]

# with drought mask
gamCt = gam(tolerance_point ~ s(latitude, k = 8), data = CxTar, method = "REML")
gamCtNoB = gam(tolerance_nob_point ~ s(latitude, k = 8), data = CxTar, method = "REML")
newlats = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar)) 
newelevs = seq(from = min(CxTar$elevation_ft), to = max(CxTar$elevation_ft), length.out = nrow(CxTar))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamCt, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCtNS = gam(tolerance_point ~ s(latitude, k = 8), data = CxTarNS, method = "REML")
yNS = predict(gamCtNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxTar, CxTarNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) +  
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex tarsalis")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4k. Cx theileri #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxThNS = read.csv("WithoutDroughtMask/CxTheileri_TSM_WithoutDroughtMask.csv")
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxTh = CxTh[CxTh$elevation_ft < 3000,]
CxThNS = CxThNS[!is.na(CxThNS$elevation_ft),]
CxThNS = CxThNS[CxThNS$elevation_ft < 3000,]

# with drought mask
gamCth = gam(tolerance_point ~ s(latitude, k = 8), data = CxTh, method = "REML")
newlats = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh)) 
newelevs = seq(from = min(CxTh$elevation_ft), to = max(CxTh$elevation_ft), length.out = nrow(CxTh))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamCth, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCthNS = gam(tolerance_point ~ s(latitude, k = 8), data = CxThNS, method = "REML")
yNS = predict(gamCthNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxTh, CxThNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex theileri")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 4l. Cx annulirostris  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxAnnulNS = read.csv("WithoutDroughtMask/CxAnnul_TSM_WithoutDroughtMask.csv")
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnul = CxAnnul[CxAnnul$elevation_ft < 3000,]
CxAnnulNS = CxAnnulNS[!is.na(CxAnnulNS$elevation_ft),]
CxAnnulNS = CxAnnulNS[CxAnnulNS$elevation_ft < 3000,]

# with drought mask
gamCa = gam(tolerance_point ~ s(latitude, k = 8), data = CxAnnul, method = "REML")
newlats = seq(from = min(CxAnnul$latitude), to = max(CxAnnul$latitude), length.out = nrow(CxAnnul)) 
newelevs = seq(from = min(CxAnnul$elevation_ft), to = max(CxAnnul$elevation_ft), length.out = nrow(CxAnnul))
newdata = data.frame(latitude = newlats, elevation_ft = newelevs)
y = predict(gamCa, data.frame(latitude = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCaNS = gam(tolerance_point ~ s(latitude, k = 8), data = CxAnnulNS, method = "REML")
yNS = predict(gamCaNS, data.frame(latitude = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxAnnul, CxAnnulNS, by = c("latitude", "longitude"))

ggplot(df, aes(x = latitude)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex annulirostris")  + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



#### 5. Create combined TSM plot (Figure 1) #####
##### 5a. Aedes aegypti and albopictus ######

newlats1 = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti)) # for projection onto a regular vector of latitudes
newelevs1 = seq(from = min(AeAegypti$elevation_ft), to = max(AeAegypti$elevation_ft), length.out = nrow(AeAegypti))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1)
y1 = predict(gamAg, data.frame(latitude = newlats1, elevation_ft = newelevs1), se.fit = TRUE) # Note: gamAm created in 3a
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(AeAlbopictus$latitude), to = max(AeAlbopictus$latitude), length.out = nrow(AeAlbopictus)) # for projection onto a regular vector of latitudes
newelevs2 = seq(from = min(AeAlbopictus$elevation_ft), to = max(AeAlbopictus$elevation_ft), length.out = nrow(AeAlbopictus))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamAb, data.frame(latitude = newlats2, elevation_ft = newelevs2), se.fit = TRUE) # Note: gamAb created in 3b
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[1], lwd = 1.3) + 
geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[1], 0.6), lty = 2, lwd = 0.7) + 
geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[1], 0.6), lty = 2, lwd = 0.7) +
geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[2], lwd = 1.3) + 
geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[2], 0.6), lwd = 0.7, lty = 2) + 
geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[2], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[1]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[2]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
  #ggtitle("Aedes aegypti & albopictus") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))
  
##### 5b. Anopheles gambiae and stephensi ######

newlats1 = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae)) # for projection onto a regular vector of latitudes
newelevs1 = seq(from = min(AnGambiae$elevation_ft), to = max(AnGambiae$elevation_ft), length.out = nrow(AnGambiae))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1) 
y1 = predict(gamAm, data.frame(latitude = newlats1), se.fit = TRUE) # Note: gamAm created in 3c
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph)) # for projection onto a regular vector of latitudes
newelevs2 = seq(from = min(AnSteph$elevation_ft), to = max(AnSteph$elevation_ft), length.out = nrow(AnSteph))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamAs, data.frame(latitude = newlats2), se.fit = TRUE)
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
  geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[3], lwd = 1.3) + 
  geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[3], 0.6), lty = 2, lwd = 0.7) + 
  geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[3], 0.6), lty = 2, lwd = 0.7) +
  geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[4], lwd = 1.3) + 
  geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[4], 0.6), lwd = 0.7, lty = 2) + 
  geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[4], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[3]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[4]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
  #ggtitle("Anopheles") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 5c. Culex pipiens and quinquefasciatus ######

newlats1 = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque)) # for projection onto a regular vector of latitudes
newelevs1 = seq(from = min(CxQuinque$elevation_ft), to = max(CxQuinque$elevation_ft), length.out = nrow(CxQuinque))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1) 
y1 = predict(gamCq, data.frame(latitude = newlats1), se.fit = TRUE) # Note: gamCq created in 3h
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip)) # for projection onto a regular vector of latitudes
newelevs2 = seq(from = min(CxPip$elevation_ft), to = max(CxPip$elevation_ft), length.out = nrow(CxPip))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamCp, data.frame(latitude = newlats2), se.fit = TRUE) # Note: gamCp created in 3i
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
  geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[6], lwd = 1.3) + 
  geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[6], 0.6), lty = 2, lwd = 0.7) + 
  geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[6], 0.6), lty = 2, lwd = 0.7) +
  geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[5], lwd = 1.3) + 
  geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[5], 0.6), lwd = 0.7, lty = 2) + 
  geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[5], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[6]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[5]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
 # ggtitle("Culex quinque & pipiens") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 5d. Ae. camptorhynchus and Cx annulirostris ######

newlats1 = seq(from = min(AeCamp$latitude), to = max(AeCamp$latitude), length.out = nrow(AeCamp)) 
newelevs1 = seq(from = min(AeCamp$elevation_ft), to = max(AeCamp$elevation_ft), length.out = nrow(AeCamp))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1) 
y1 = predict(gamAc, data.frame(latitude = newlats1), se.fit = TRUE) # Note: gamAc created in 3c
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(CxAnnul$latitude), to = max(CxAnnul$latitude), length.out = nrow(CxAnnul)) 
newelevs2 = seq(from = min(CxAnnul$elevation_ft), to = max(CxAnnul$elevation_ft), length.out = nrow(CxAnnul))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamCa, data.frame(latitude = newlats2), se.fit = TRUE) # Note: gamCa created in 3l
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
  geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[7], lwd = 1.3) + 
  geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[7], 0.6), lty = 2, lwd = 0.7) + 
  geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[7], 0.6), lty = 2, lwd = 0.7) +
  geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[8], lwd = 1.3) + 
  geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[8], 0.6), lwd = 0.7, lty = 2) + 
  geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[8], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[7]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[8]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
 # ggtitle("Ae camptorhyncus and Cx annulirostris") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 5e. Ae. triseriatus and Cx tarsalis ######

newlats1 = seq(from = min(AeTri$latitude), to = max(AeTri$latitude), length.out = nrow(AeTri)) 
newelevs1 = seq(from = min(AeTri$elevation_ft), to = max(AeTri$elevation_ft), length.out = nrow(AeTri))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1) 
y1 = predict(gamAt, data.frame(latitude = newlats1), se.fit = TRUE) # Note: gamAt created in 3d
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar)) 
newelevs2 = seq(from = min(CxTar$elevation_ft), to = max(CxTar$elevation_ft), length.out = nrow(CxTar))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamCt, data.frame(latitude = newlats2), se.fit = TRUE) # Note: gamCt created in 3j
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
  geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[9], lwd = 1.3) + 
  geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[9], 0.6), lty = 2, lwd = 0.7) + 
  geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[9], 0.6), lty = 2, lwd = 0.7) +
  geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[10], lwd = 1.3) + 
  geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[10], 0.6), lwd = 0.7, lty = 2) + 
  geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[10], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[9]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[10]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
 # ggtitle("Ae triseriatus and Cx tarsalis") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 5f. Ae. vexans and Cx theileri ######

newlats1 = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans)) 
newelevs1 = seq(from = min(AeVexans$elevation_ft), to = max(AeVexans$elevation_ft), length.out = nrow(AeVexans))
newdata1 = data.frame(latitude = newlats1, elevation_ft = newelevs1) 
y1 = predict(gamAv, data.frame(latitude = newlats1), se.fit = TRUE) # Note: gamAv created in 3e
upr1 <- y1$fit + (1.96 * y1$se.fit)
lwr1 <- y1$fit - (1.96 * y1$se.fit)

newlats2 = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh)) 
newelevs2 = seq(from = min(CxTh$elevation_ft), to = max(CxTh$elevation_ft), length.out = nrow(CxTh))
newdata2 = data.frame(latitude = newlats2, elevation_ft = newelevs2)
y2 = predict(gamCth, data.frame(latitude = newlats2), se.fit = TRUE) # Note: gamCth created in 3k
upr2 <- y2$fit + (1.96 * y2$se.fit)
lwr2 <- y2$fit - (1.96 * y2$se.fit)

ggplot() + theme_bw() +
  geom_line(aes(x=newlats1, y = y1$fit), colour=SpeciesColors[11], lwd = 1.3) + 
  geom_line(aes(x=newlats1, y = upr1), colour=alpha(SpeciesColors[11], 0.6), lty = 2, lwd = 0.7) + 
  geom_line(aes(x=newlats1, y = lwr1), colour=alpha(SpeciesColors[11], 0.6), lty = 2, lwd = 0.7) +
  geom_line(aes(x=newlats2, y = y2$fit), color = SpeciesColors[12], lwd = 1.3) + 
  geom_line(aes(x=newlats2, y = upr2), colour= alpha(SpeciesColors[12], 0.6), lwd = 0.7, lty = 2) + 
  geom_line(aes(x=newlats2, y = lwr2), colour= alpha(SpeciesColors[12], 0.6), lwd = 0.7, lty = 2) +
  geom_ribbon(aes(x = newlats1, ymin = lwr1, ymax = upr1), alpha = 0.1, fill = SpeciesColors[11]) +
  geom_ribbon(aes(x = newlats2, ymin = lwr2, ymax = upr2), alpha = 0.1, fill = SpeciesColors[12]) +
  geom_hline(yintercept=0, linetype="dashed", color = "black", lwd=0.6) + 
  labs(x = " ", y = "Thermal safety margin (\u00B0C)") + 
 #  ggtitle("Ae vexans and Cx theileri") + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18),
        legend.text=element_text(size=18), legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


#### All species maximum body temps across latitude plots ####
AllSp = rbind.data.frame(AeAegypti, AeAlbo, AeCamp, AeTri, AeVex,
                         AnGam, AnSteph, CxA, CxPip, CxQ, CxTar, CxTh)
plot(AllSp$maxtemp_point ~ AllSp$latitude)

SpeciesColors = c("#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#309143", 
                  "#ab041b", "#ec3c30", 
                  "#abd9e9", "#74add1", "#313695", "#800080", "#dda0dd")
SpColors = rep(SpeciesColors, each = 80)


# Max Temps
ggplot(AllSp, aes(x = latitude, y = maxtemp_point, color = "black")) +  theme_minimal() +   
  geom_smooth(method = "loess", span = 0.5, col = "black", alpha = 0.3) +
  geom_point(aes(colour=species.pop), size = 1.9, alpha = 0.5) + 
scale_colour_manual(values = SpeciesColors) + 
  scale_fill_manual(values = SpeciesColors) + 
  labs(x = "Latitude", y = "Maximum body temperature (\u00B0C)") + 
 scale_y_continuous(limits = c(28,48)) + 
  scale_x_continuous(limits = c(-40, 60), breaks = seq(-40, 60, 10), position = "bottom") + 
  ggtitle(" ") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(colour = "black", fill = NA)) 


# Elevation by latitude
ggplot(AllSp, aes(x = latitude, y = elevation_ft, color = "black")) +  theme_minimal() +   
  geom_smooth(method = "loess", span = 0.5, col = "black", alpha = 0.3) +
  geom_point(aes(colour=species.pop), size = 1.5, alpha = 0.5) + 
  scale_colour_manual(values = SpeciesColors) + 
  scale_fill_manual(values = SpeciesColors) + 
  labs(x = "Latitude", y = "Elevation (ft)") + 
  scale_y_continuous(limits = c(0,9000)) + 
  scale_x_continuous(limits = c(-40, 60), breaks = seq(-40, 60, 10), position = "bottom") + 
  ggtitle(" ") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(colour = "black", fill = NA)) 


  
