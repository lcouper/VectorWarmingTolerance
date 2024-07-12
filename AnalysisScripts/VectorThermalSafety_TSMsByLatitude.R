#############################################################
########### CODE FOR VECTOR TSM GAMS ########################
########### WRITTEN BY LISA COUPER ##########################

# This script is used to model and plot results for each vector species including:
# 1. Run GAM of TSM ~ latitude
# 2. Run GAM of TSM ~ latitude (high v low elevation)
# 3. Run GAM of TSM ~ latitude (with and w/o behavior)
# 4. Run GAM of TSM ~ latitude (with and w/o drought mask)
# 5. Max body temp across latitude 
# 6. Generate table of average TSM under varying specifications

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

#### 0b. Load data frames #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files with seasonality
AeAegypti = read.csv("AeAegypti_TSM_DroughtMask_Combined_WithElevation.csv")
AeAlbopictus = read.csv("AeAlbo_TSM_DroughtMask_Combined_WithElevation.csv")
AnGambiae = read.csv("AnGambiae_TSM_DroughtMask_Combined_WithElevation.csv")
AnSteph = read.csv("AnSteph_TSM_DroughtMask_Combined_WithElevation.csv")
CxQuinque = read.csv("CxQuinque_TSM_DroughtMask_Combined_WithElevation.csv")  
CxPip = read.csv("CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")
CxTar = read.csv("CxTarsalis_TSM_DroughtMask_Combined_WithElevation.csv")
CxAnnul = read.csv("CxAnnul_TSM_DroughtMask_Combined_WithElevation.csv")

#### 1. GAM of TSM across latitude for each species ####

##### 1a. Aedes aegypti #####

# First remove high elevation samples (>2,500 m)
range(AeAegypti$elevation) # none to remove
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation),]
AeAegypti = AeAegypti[AeAegypti$elevation < 2500,] 

# Run GAM estimating thermal safety margins across latitude
gamAg = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAegypti, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAegypti$lat), to = max(AeAegypti$lat), length.out = nrow(AeAegypti))

# Predict using GAM model
y = predict(gamAg, newdata = data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam 
ggplot() +  theme_minimal() +
  geom_point(aes(x = AeAegypti$lat, y = AeAegypti$tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes aegypti") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 1b. Aedes albopictus #####

# First remove any high elevation samples (>2,500 m)
range(AeAlbopictus$elevation)
AeAlbopictus = AeAlbopictus[!is.na(AeAlbopictus$elevation),] # none to remove
AeAlbopictus = AeAlbopictus[AeAlbopictus$elevation < 2500,] 

# Run GAM estimating thermal safety margins across latitude
gamAb = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAlbopictus, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAlbopictus$lat), to = max(AeAlbopictus$lat), length.out = nrow(AeAlbopictus)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamAb, data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AeAlbopictus, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes albopictus") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 1f. Anopheles gambiae #####

# First remove high elevation samples (>2,500 m)
range(AnGambiae$elevation) 
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation),] # none to remove
AnGambiae = AnGambiae[AnGambiae$elevation < 2500,] 

# Run GAM estimating thermal safety margins across latitude
gamAm = gam(tolerance_point2 ~ s(lat, k = 8), data = AnGambiae, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnGambiae$lat), to = max(AnGambiae$lat), length.out = nrow(AnGambiae)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamAm, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(AnGambiae, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles gambiae") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 1g. Anopheles stephensi #####

# First remove high elevation samples (>2500 m)
range(AnSteph$elevation)
AnSteph = AnSteph[!is.na(AnSteph$elevation),] # none to remove
AnSteph= AnSteph[AnSteph$elevation < 2500,] 

# Run GAM estimating thermal safety margins across latitude
gamAs = gam(tolerance_point2 ~ s(lat, k = 8), data = AnSteph, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnSteph$lat), to = max(AnSteph$lat), length.out = nrow(AnSteph))

# Predict using GAM model
y = predict(gamAs, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnSteph, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles stephensi") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1h. Culex quinquefasciatus ####

# First remove high elevation samples (>2500 m)
range(CxQuinque$elevation)
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation),] 
CxQuinque = CxQuinque[CxQuinque$elevation < 2500,] # removed 1 row

# Run GAM estimating thermal safety margins across latitude
gamCq = gam(tolerance_point2 ~ s(lat, k = 8), data = CxQuinque, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxQuinque$lat), to = max(CxQuinque$lat), length.out = nrow(CxQuinque)) 

# Predict using GAM model
y = predict(gamCq, data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(CxQuinque, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex quinquefasciatus") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 1i. Culex pipiens ####

# First separate any samples from high elevation (>2500 m)
range(CxPip$elevation)
CxPip = CxPip[!is.na(CxPip$elevation),]
CxPip = CxPip[CxPip$elevation < 2500,] # none to remove

# First remove high elevation samples (>2500 ft)
gamCp = gam(tolerance_point2 ~ s(lat, k = 8), data = CxPip, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxPip$lat), to = max(CxPip$lat), length.out = nrow(CxPip)) 

# Predict using GAM model
y = predict(gamCp, data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxPip, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex pipiens") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 1j. Culex tarsalis ####

# First remove high elevation samples (>2500 m)
range(CxTar$elevation)
CxTar = CxTar[!is.na(CxTar$elevation),]
CxTar = CxTar[CxTar$elevation < 2500,] # removes 2

# Run GAM estimating thermal safety margins across latitude
gamCt = gam(tolerance_point2 ~ s(lat, k = 8), data = CxTar, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTar$lat), to = max(CxTar$lat), length.out = nrow(CxTar)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamCt, data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam 
ggplot(CxTar, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex tarsalis") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 1l. Culex annulirostris ####

# First remove high elevation samples (>2500 m)
range(CxAnnul$elevation)
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation),]
CxAnnul = CxAnnul[CxAnnul$elevation < 2500,] # none to remove

# Run GAM estimating thermal safety margins across latitude
gamCa = gam(tolerance_point2 ~ s(lat, k = 8), data = CxAnnul, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxAnnul$lat), to = max(CxAnnul$lat), length.out = nrow(CxAnnul)) # for projection onto a regular vector of latitudes

# Predict using GAM model
y = predict(gamCa, data.frame(lat = newlats), se.fit = TRUE)

# Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# Plot gam
ggplot(CxAnnul, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = alpha("black", 0.6)) +
  geom_line(aes(x=newlats, y = y$fit), colour="black", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("black", 0.8), lty = 2, lwd = 1.1) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex annulirostris") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


#### 3. GAM of TSM across latitude: With and without behavior ####
##### 3a. Aedes aegypti #####

# no high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamAg = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAegypti, method = "REML")
gamAgNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = AeAegypti, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAegypti$lat), to = max(AeAegypti$lat), length.out = nrow(AeAegypti)) 

# Predict using GAM model
y = predict(gamAg, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamAgNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeAegypti, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes aegypti") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

##### 3b. Aedes albopictus #####

# no high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamAb = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAlbopictus, method = "REML")
gamAbNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = AeAlbopictus, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AeAlbopictus$lat), to = max(AeAlbopictus$lat), length.out = nrow(AeAlbo)) 

# Predict using GAM model
y = predict(gamAb,data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamAbNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AeAlbopictus, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Aedes albopictus") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 3f. Anopheles gambiae #####

# no high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamAm = gam(tolerance_point2 ~ s(lat, k = 8), data = AnGambiae, method = "REML")
gamAmNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = AnGambiae, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnGambiae$lat), to = max(AnGambiae$lat), length.out = nrow(AnGambiae)) 

# Predict using GAM model
y = predict(gamAm, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamAmNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnGambiae, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.3) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles gambiae") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 3g. Anopheles stephensi #####

# no high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamAs = gam(tolerance_point2 ~ s(lat, k = 8), data = AnSteph, method = "REML")
gamAsNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = AnSteph, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(AnSteph$lat), to = max(AnSteph$lat), length.out = nrow(AnSteph)) 

# Predict using GAM model
y = predict(gamAs, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamAsNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(AnSteph, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Anopheles stephensi") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))




##### 3h. Culex quinquefasciatus #####

# remove the 1 high elevation record
CxQuinque = CxQuinque[CxQuinque$elevation < 2500,] 

# Run GAM estimating thermal safety margins across latitude
gamCq = gam(tolerance_point2 ~ s(lat, k = 8), data = CxQuinque, method = "REML")
gamCqNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = CxQuinque, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxQuinque$lat), to = max(CxQuinque$lat), length.out = nrow(CxQuinque)) 

# Predict using GAM model
y = predict(gamCq, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamCqNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxQuinque, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.1) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex quinquefasciatus") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3i. Culex pipiens #####

# no high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamCp = gam(tolerance_point2 ~ s(lat, k = 8), data = CxPip, method = "REML")
gamCpNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = CxPip, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxPip$lat), to = max(CxPip$lat), length.out = nrow(CxPip)) 

# Predict using GAM model
y = predict(gamCp, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamCpNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxPip, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex pipiens") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3j. Culex tarsalis #####

# remove the 2 high elevation records
CxTar = CxTar[CxTar$elevation < 2500,] # removes 2

# Run GAM estimating thermal safety margins across latitude
gamCt = gam(tolerance_point2 ~ s(lat, k = 8), data = CxTar, method = "REML")
gamCtNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = CxTar, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxTar$lat), to = max(CxTar$lat), length.out = nrow(CxTar)) 

# Predict using GAM model
y = predict(gamCt, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamCtNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxTar, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex tarsalis") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 3l. Culex annulirostris #####

# No high elevation records to remove

# Run GAM estimating thermal safety margins across latitude
gamCa = gam(tolerance_point2 ~ s(lat, k = 8), data = CxAnnul, method = "REML")
gamCaNoB = gam(tolerance_NoB_point2 ~ s(lat, k = 8), data = CxAnnul, method = "REML")

# Create new data frame on which to predict TSM
newlats = seq(from = min(CxAnnul$lat), to = max(CxAnnul$lat), length.out = nrow(CxAnnul)) 

# Predict using GAM model
y = predict(gamCa, data.frame(lat = newlats), se.fit = TRUE)
ynoB = predict(gamCaNoB, data.frame(lat = newlats), se.fit = TRUE)

# Fourth, Create confidence intervals
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)
uprnoB <- ynoB$fit + (1.96 * ynoB$se.fit)
lwrnoB <- ynoB$fit - (1.96 * ynoB$se.fit)

# Fifth, Plot with and without behavioral thermoregulation
ggplot(CxAnnul, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2), size = 1.2, color = "#0D0887FF") +
  geom_line(aes(x=newlats, y = y$fit), colour="#0D0887FF", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#0D0887FF", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_NoB_point2), size = 1.2, color = "darkred") +
  geom_line(aes(x=newlats, y = ynoB$fit), colour=alpha("darkred", 1.1)) +
  geom_line(aes(x=newlats, y = uprnoB), colour=alpha("darkred", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrnoB), colour=alpha("darkred", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed", 
             color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C)
  ggtitle("Culex annulirostris") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

#### 4. GAM of TSM across latitude: With and without drought mask ####
##### 4a. Aedes aegypti #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AeAegyptiNS = read.csv("WithoutDroughtMask/AeAegypti_TSM_NoDroughtMask_Combined_WithElevation.csv")
AeAegyptiNS$lat == AeAegypti$lat # check they capture the same points

# with drought mask
gamAg = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAegypti, method = "REML")
newlats = seq(from = min(AeAegypti$lat), to = max(AeAegypti$lat), length.out = nrow(AeAegypti)) 
newdata = data.frame(lat = newlats)
y = predict(gamAg, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAgNS = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAegyptiNS, method = "REML")
yNS = predict(gamAgNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeAegypti, AeAegyptiNS, by = c("lat", "lon"))
# df <- df %>% distinct(lat, lon, .keep_all = TRUE) 

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes aegypti")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4b. Aedes albopictus  #####

# Bring in drought mask data set
AeAlboNS = read.csv("WithoutDroughtMask/AeAlbopictus_TSM_NoDroughtMask_Combined_WithElevation.csv")

# with drought mask
gamAb = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAlbopictus, method = "REML")
newlats = seq(from = min(AeAlbopictus$lat), to = max(AeAlbopictus$lat), length.out = nrow(AeAlbopictus)) 
newdata = data.frame(lat = newlats)
y = predict(gamAb, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAbNS = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAlboNS, method = "REML")
yNS = predict(gamAbNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AeAlbopictus, AeAlboNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Aedes albopictus")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4f. Anopheles gambiae  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AnGambiaeNS = read.csv("WithoutDroughtMask/AnGambiae_TSM_NoDroughtMask_Combined_WithElevation.csv")

# with drought mask
gamAm = gam(tolerance_point2 ~ s(lat, k = 8), data = AnGambiae, method = "REML")
newlats = seq(from = min(AnGambiae$lat), to = max(AnGambiae$lat), length.out = nrow(AnGambiae))
newdata = data.frame(lat = newlats)
y = predict(gamAm, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAmNS = gam(tolerance_point2 ~ s(lat, k = 8), data = AnGambiaeNS, method = "REML")
yNS = predict(gamAmNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AnGambiae, AnGambiaeNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Anopheles gambiae")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4g. Anopheles stephensi  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
AnStephNS = read.csv("WithoutDroughtMask/AnSteph_TSM_NoDroughtMask_Combined_WithElevation.csv")

# with drought mask
gamAs = gam(tolerance_point2 ~ s(lat, k = 8), data = AnSteph, method = "REML")
newlats = seq(from = min(AnSteph$lat), to = max(AnSteph$lat), length.out = nrow(AnSteph)) 
y = predict(gamAs, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamAsNS = gam(tolerance_point2 ~ s(lat, k = 8), data = AnStephNS, method = "REML")
yNS = predict(gamAsNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(AnSteph, AnStephNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Anopheles stephensi")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4h. Cx quinquefasciatus  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxQuinqueNS = read.csv("WithoutDroughtMask/CxQuinque_TSM_NoDroughtMask_Combined_WithElevation.csv")
CxQuinqueNS = CxQuinqueNS[CxQuinqueNS$elevation < 2500,] # removed 1 row
CxQuinque = CxQuinque[CxQuinque$elevation < 2500,] # removed 1 row

# with drought mask
gamCq = gam(tolerance_point2 ~ s(lat, k = 8), data = CxQuinque, method = "REML")
newlats = seq(from = min(CxQuinque$lat), to = max(CxQuinque$lat), length.out = nrow(CxQuinque)) 
y = predict(gamCq, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCqNS = gam(tolerance_point2 ~ s(lat, k = 8), data = CxQuinqueNS, method = "REML")
yNS = predict(gamCqNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxQuinque, CxQuinqueNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex quinquefasciatus")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4i. Cx pipiens  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxPipNS = read.csv("WithoutDroughtMask/CxPipiens_TSM_NoDroughtMask_Combined_WithElevation.csv")
CxPip = read.csv("CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")

# keep only CxPip also in CxPipNS and vice versa
CxPip <- CxPip %>% filter(lat %in% unique(CxPipNS$lat) & lon %in% unique(CxPipNS$lon))
CxPipNS <- CxPipNS %>% filter(lat %in% unique(CxPip$lat) & lon %in% unique(CxPip$lon))

# with drought mask
gamCp = gam(tolerance_point2 ~ s(lat, k = 8), data = CxPip, method = "REML")
newlats = seq(from = min(CxPip$lat), to = max(CxPip$lat), length.out = nrow(CxPip)) 
y = predict(gamCp, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCpNS = gam(tolerance_point2 ~ s(lat, k = 8), data = CxPipNS, method = "REML")
yNS = predict(gamCpNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxPip, CxPipNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex pipiens")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


##### 4j. Cx tarsalis  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxTarNS = read.csv("WithoutDroughtMask/CxTarsalis_TSM_NoDroughtMask_Combined_WithElevation.csv")
CxTarNS = CxTarNS[CxTarNS$elevation < 2500,] # removes 2
CxTar = CxTar[CxTar$elevation < 2500,] # removes 2

# keep only CxTar also in CxTarNS and vice versa
CxTar <- CxTar %>% filter(lat %in% unique(CxTarNS$lat) & lon %in% unique(CxTarNS$lon))
CxTarNS <- CxTarNS %>% filter(lat %in% unique(CxTar$lat) & lon %in% unique(CxTar$lon))

# with drought mask
gamCt = gam(tolerance_point2 ~ s(lat, k = 8), data = CxTar, method = "REML")
newlats = seq(from = min(CxTar$lat), to = max(CxTar$lat), length.out = nrow(CxTar)) 
y = predict(gamCt, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCtNS = gam(tolerance_point2 ~ s(lat, k = 8), data = CxTarNS, method = "REML")
yNS = predict(gamCtNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxTar, CxTarNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) +  
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex tarsalis")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### 4l. Cx annulirostris  #####

# Bring in drought mask data set. Remove high elevation samples in both this and original
CxAnnulNS = read.csv("WithoutDroughtMask/CxAnnul_TSM_NoDroughtMask_Combined_WithElevation.csv")
CxAnnulNS$lon == CxAnnul$lon

# with drought mask
gamCa = gam(tolerance_point2 ~ s(lat, k = 8), data = CxAnnul, method = "REML")
newlats = seq(from = min(CxAnnul$lat), to = max(CxAnnul$lat), length.out = nrow(CxAnnul)) 
y = predict(gamCa, data.frame(lat = newlats), se.fit = TRUE)
upr <- y$fit + (1.96 * y$se.fit)
lwr <- y$fit - (1.96 * y$se.fit)

# without drought mask
gamCaNS = gam(tolerance_point2 ~ s(lat, k = 8), data = CxAnnulNS, method = "REML")
yNS = predict(gamCaNS, data.frame(lat = newlats), se.fit = TRUE)
uprNS <- yNS$fit + (1.96 * yNS$se.fit)
lwrNS <- yNS$fit - (1.96 * yNS$se.fit)

# Plot with and without drought mask
df = merge(CxAnnul, CxAnnulNS, by = c("lat", "lon"))

ggplot(df, aes(x = lat)) +  theme_minimal() +
  geom_point(aes(y = tolerance_point2.x), size = 1.2, color = "#4b0076") +
  geom_line(aes(x=newlats, y = y$fit), colour="#4b0076", lwd = 1.2) + 
  geom_line(aes(x=newlats, y = upr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) + 
  geom_line(aes(x=newlats, y = lwr), colour=alpha("#4b0076", 1.0), lty = 2, lwd = 1.1) +
  geom_point(aes(y = tolerance_point2.y), size = 1.2, color = "#dda0dd") +
  geom_line(aes(x=newlats, y = yNS$fit), colour=alpha("#dda0dd", 1.2)) +
  geom_line(aes(x=newlats, y = uprNS), colour=alpha("#dda0dd", 1.0), lty = 2) + 
  geom_line(aes(x=newlats, y = lwrNS), colour=alpha("#dda0dd", 1.0), lty = 2) +
  geom_hline(yintercept=0, linetype="dashed",  color = "black", lwd=0.6) + 
  labs(x = " ", y = " ") + # Thermal safety margin (\u00B0C) +
  ggtitle("Culex annulirostris")  + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))



##### Plots: All species max body temps across latitude ####
AllSp = rbind.data.frame(AeAegypti, AeAlbopictus,
                         AnGambiae, AnSteph, CxAnnul, CxPip, CxQuinque, CxTar)
plot(AllSp$maxtemp_point ~ AllSp$latitude)

SpeciesColors = c("#ab041b", "#cc5801", "#f4bb00","#7cae00", "#309143",
                  "#74add1","#313695",  "#dda0dd")

# Max Temps
ggplot(AllSp, aes(x = lat, y = maxtemp_point2, color = "black")) +  theme_minimal() +   
  geom_smooth(method = "loess", span = 0.5, col = "black", alpha = 0.3) +
  geom_point(aes(colour=Species.Pop), size = 1.3, alpha = 0.5) + 
  scale_colour_manual(values = SpeciesColors) + 
  scale_fill_manual(values = SpeciesColors) + 
  labs(x = "Latitude", y = "Maximum body temperature (\u00B0C)") + 
  scale_y_continuous(limits = c(20,50)) + 
  scale_x_continuous(limits = c(-40, 60), breaks = seq(-40, 60, 10), position = "bottom") + 
  ggtitle(" ") + 
  theme(axis.text=element_text(size=20), 
        axis.title = element_text(size = 20),
        legend.position = "none",
        plot.title = element_text(hjust = 0.5, size = 24),
        panel.border = element_rect(colour = "black", fill = NA)) 


##### Supp Tables: Average TSM under different specifications #####

AverageTSMs = as.data.frame(matrix(nrow = 8, ncol = 5))
colnames(AverageTSMs) = c("species", "mean", "sd", "lowerCI", "upperCI")
AverageTSMs[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                    "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")
AverageTSMs[1,2:5] = c(mean(AeAegypti$tolerance_point2), sd(AeAegypti$tolerance_point2),
                       t.test(AeAegypti$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[2,2:5] = c(mean(AeAlbopictus$tolerance_point2), sd(AeAlbopictus$tolerance_point2),
                       t.test(AeAlbopictus$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[3,2:5] = c(mean(AnGambiae$tolerance_point2), sd(AnGambiae$tolerance_point2),
                       t.test(AnGambiae$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[4,2:5] = c(mean(AnSteph$tolerance_point2), sd(AnSteph$tolerance_point2),
                       t.test(AnSteph$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[5,2:5] = c(mean(CxAnnul$tolerance_point2), sd(CxAnnul$tolerance_point2),
                       t.test(CxAnnul$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[6,2:5] = c(mean(CxPip$tolerance_point2), sd(CxPip$tolerance_point2),
                       t.test(CxPip$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[7,2:5] = c(mean(CxQuinque$tolerance_point2), sd(CxQuinque$tolerance_point2),
                       t.test(CxQuinque$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMs[8,2:5] = c(mean(CxTar$tolerance_point2), sd(CxTar$tolerance_point2),
                       t.test(CxTar$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$tolerance_point2, conf.level = 0.95)$conf.int[2])

#fwrite(AverageTSMs, "SuppTable_AverageTSMs_MainModelSpecification.csv")
avgTSMs <- c(1.62, -3.61, 1.66, -1.48, -2.62, -0.34, -1.02, 1.30)

# Average TMS without behavior 
AverageTSMsNoB = as.data.frame(matrix(nrow = 8, ncol = 5))
colnames(AverageTSMsNoB) = c("species", "mean", "sd", "lowerCI", "upperCI")
AverageTSMsNoB[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                    "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageTSMsNoB[1,2:5] = c(mean(AeAegypti$tolerance_NoB_point2), sd(AeAegypti$tolerance_NoB_point2),
                       t.test(AeAegypti$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[2,2:5] = c(mean(AeAlbopictus$tolerance_NoB_point2), sd(AeAlbopictus$tolerance_NoB_point2),
                       t.test(AeAlbopictus$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[3,2:5] = c(mean(AnGambiae$tolerance_NoB_point2), sd(AnGambiae$tolerance_point2),
                       t.test(AnGambiae$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[4,2:5] = c(mean(AnSteph$tolerance_NoB_point2), sd(AnSteph$tolerance_point2),
                       t.test(AnSteph$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[5,2:5] = c(mean(CxAnnul$tolerance_NoB_point2), sd(CxAnnul$tolerance_point2),
                       t.test(CxAnnul$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[6,2:5] = c(mean(CxPip$tolerance_NoB_point2), sd(CxPip$tolerance_point2),
                       t.test(CxPip$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[7,2:5] = c(mean(CxQuinque$tolerance_NoB_point2), sd(CxQuinque$tolerance_point2),
                       t.test(CxQuinque$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNoB[8,2:5] = c(mean(CxTar$tolerance_NoB_point2), sd(CxTar$tolerance_point2),
                       t.test(CxTar$tolerance_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$tolerance_NoB_point2, conf.level = 0.95)$conf.int[2])

#fwrite(AverageTSMsNoB, "SuppTable_AverageTSMs_NoBehavior.csv")

avgTSMsNoB <- c(-1.34, -6.45, -2.04, -5.03, -5.81, -3.13, -3.93, -1.97)


# Average TMS without drought mask
AverageTSMsNS = as.data.frame(matrix(nrow = 8, ncol = 5))
colnames(AverageTSMsNS) = c("species", "mean", "sd", "lowerCI", "upperCI", "nax")
AverageTSMsNS[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                       "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageTSMsNS[1,2:5] = c(mean(AeAegyptiNS$tolerance_point2), sd(AeAegyptiNS$tolerance_point2),
                       t.test(AeAegyptiNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegyptiNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[2,2:5] = c(mean(AeAlboNS$tolerance_point2), sd(AeAlboNS$tolerance_point2),
                       t.test(AeAlboNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlboNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[3,2:5] = c(mean(AnGambiaeNS$tolerance_point2), sd(AnGambiaeNS$tolerance_point2),
                       t.test(AnGambiaeNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiaeNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[4,2:5] = c(mean(AnStephNS$tolerance_point2), sd(AnStephNS$tolerance_point2),
                       t.test(AnStephNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(AnStephNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[5,2:5] = c(mean(CxAnnulNS$tolerance_point2), sd(CxAnnulNS$tolerance_point2),
                       t.test(CxAnnulNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnulNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[6,2:5] = c(mean(CxPipNS$tolerance_point2), sd(CxPipNS$tolerance_point2),
                       t.test(CxPipNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxPipNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[7,2:5] = c(mean(CxQuinqueNS$tolerance_point2), sd(CxQuinqueNS$tolerance_point2),
                       t.test(CxQuinqueNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinqueNS$tolerance_point2, conf.level = 0.95)$conf.int[2])
AverageTSMsNS[8,2:5] = c(mean(CxTarNS$tolerance_point2), sd(CxTarNS$tolerance_point2),
                       t.test(CxTarNS$tolerance_point2, conf.level = 0.95)$conf.int[1], t.test(CxTarNS$tolerance_point2, conf.level = 0.95)$conf.int[2])

#fwrite(AverageTSMsNS, "SuppTable_AverageTSMs_NoDroughtMask.csv")

avgTSMsNS <- c(1.58, -3.65, 1.57, -1.68, -2.62, -0.91, -1.08, 1.26)


##### Supp Tables: Average streak (hours) in thermal danger under different specifications #####

AverageHours = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageHours) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageHours[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                     "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageHours[1,2:6] = c(mean(AeAegypti$streak_point2, na.rm = T), sd(AeAegypti$streak_point2, na.rm = T),
                      t.test(AeAegypti$streak_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$streak_point2, conf.level = 0.95)$conf.int[2],
                      max(AeAegypti$streak_point2, na.rm = T))
AverageHours[2,2:6] = c(mean(AeAlbopictus$streak_point2, na.rm = T), sd(AeAlbopictus$streak_point2, na.rm = T),
                      t.test(AeAlbopictus$streak_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$streak_point2, conf.level = 0.95)$conf.int[2],
                      max(AeAlbo$streak_point2, na.rm = T))
AverageHours[3,2:6] = c(mean(AnGambiae$streak_point2, na.rm = T), sd(AnGambiae$streak_point2,na.rm = T),
                     t.test(AnGambiae$streak_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$streak_point2, conf.level = 0.95)$conf.int[2],
                     max(AnGambiae$streak_point2, na.rm = T))
AverageHours[4,2:6] = c(mean(AnSteph$streak_point2, na.rm = T), sd(AnSteph$streak_point2,na.rm = T),
                       t.test(AnSteph$streak_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$streak_point2, conf.level = 0.95)$conf.int[2],
                       max(AnSteph$streak_point2, na.rm = T))
AverageHours[5,2:6] = c(mean(CxAnnul$streak_point2, na.rm = T), sd(CxAnnul$streak_point2,na.rm = T),
                       t.test(CxAnnul$streak_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$streak_point2, conf.level = 0.95)$conf.int[2],
                       max(CxAnnul$streak_point2, na.rm = T))
AverageHours[6,2:6] = c(mean(CxPip$streak_point2, na.rm = T), sd(CxPip$streak_point2,na.rm = T),
                     t.test(CxPip$streak_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$streak_point2, conf.level = 0.95)$conf.int[2],
                     max(CxPip$streak_point2, na.rm = T))
AverageHours[7,2:6] = c(mean(CxQuinque$streak_point2, na.rm = T), sd(CxQuinque$streak_point2,na.rm = T),
                        t.test(CxQuinque$streak_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$streak_point2, conf.level = 0.95)$conf.int[2],
                        max(CxQuinque$streak_point2, na.rm = T))
AverageHours[8,2:6] = c(mean(CxTar$streak_point2, na.rm = T), sd(CxTar$streak_point2,na.rm = T),
                      t.test(CxTar$streak_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$streak_point2, conf.level = 0.95)$conf.int[2],
                      max(CxTar$streak_point2, na.rm = T))
                        
#fwrite(AverageHours, "SuppTable_AverageHoursThermalDanger_MainModelSpecification.csv")

avgHours <- c(6.66, 9.36, 7.26, 8.56, 9.90, 7.37, 8.99, 5.94)

## Average hours without behavioral thermoregulation
AverageHoursNoB = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageHoursNoB) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageHoursNoB[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                     "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageHoursNoB[1,2:6] = c(mean(AeAegypti$streak_NoB_point2, na.rm = T), sd(AeAegypti$streak_NoB_point2, na.rm = T),
                                t.test(AeAegypti$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                                max(AeAegypti$streak_NoB_point2, na.rm = T))
AverageHoursNoB[2,2:6] = c(mean(AeAlbopictus$streak_NoB_point2, na.rm = T), sd(AeAlbopictus$streak_NoB_point2, na.rm = T),
                           t.test(AeAlbopictus$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                            max(AeAlbo$streak_NoB_point2, na.rm = T))
AverageHoursNoB[3,2:6] = c(mean(AnGambiae$streak_NoB_point2, na.rm = T), sd(AnGambiae$streak_NoB_point2,na.rm = T),
                           t.test(AnGambiae$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                            max(AnGambiae$streak_NoB_point2, na.rm = T))
 AverageHoursNoB[4,2:6] = c(mean(AnSteph$streak_NoB_point2, na.rm = T), sd(AnSteph$streak_NoB_point2,na.rm = T),
                             t.test(AnSteph$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                             max(AnSteph$streak_NoB_point2, na.rm = T))
 AverageHoursNoB[5,2:6] = c(mean(CxAnnul$streak_NoB_point2, na.rm = T), sd(CxAnnul$streak_NoB_point2,na.rm = T),
                           t.test(CxAnnul$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                           max(CxAnnul$streak_NoB_point2, na.rm = T))
AverageHoursNoB[6,2:6] = c(mean(CxPip$streak_NoB_point2, na.rm = T), sd(CxPip$streak_NoB_point2,na.rm = T),
                         t.test(CxPip$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                       max(CxPip$streak_NoB_point2, na.rm = T))
AverageHoursNoB[7,2:6] = c(mean(CxQuinque$streak_NoB_point2, na.rm = T), sd(CxQuinque$streak_NoB_point2,na.rm = T),
                        t.test(CxQuinque$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                      max(CxQuinque$streak_NoB_point2, na.rm = T))
AverageHoursNoB[8,2:6] = c(mean(CxTar$streak_NoB_point2, na.rm = T), sd(CxTar$streak_NoB_point2,na.rm = T),
                          t.test(CxTar$streak_NoB_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$streak_NoB_point2, conf.level = 0.95)$conf.int[2],
                           max(CxTar$streak_NoB_point2, na.rm = T))
                          
# fwrite(AverageHoursNoB, "SuppTable_AverageHoursThermalDanger_WithoutBehavior.csv")
avgHoursNoB <- c(6.00, 10, 6.96, 8.19, 8.93, 7.88, 8.66, 6.64)

AverageHoursNS = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageHoursNS) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageHoursNS[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                        "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")
AverageHoursNS[1,2:6] = c(mean(AeAegyptiNS$streak_point2, na.rm =T), sd(AeAegyptiNS$streak_point2, na.rm =T),
                          t.test(AeAegyptiNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(AeAegyptiNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(AeAegyptiNS$streak_point2, na.rm = T))
AverageHoursNS[2,2:6] = c(mean(AeAlboNS$streak_point2, na.rm =T), sd(AeAlboNS$streak_point2, na.rm =T),
                          t.test(AeAlboNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(AeAlboNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(AeAlboNS$streak_point2, na.rm = T))
AverageHoursNS[3,2:6] = c(mean(AnGambiaeNS$streak_point2, na.rm =T), sd(AnGambiaeNS$streak_point2, na.rm =T),
                          t.test(AnGambiaeNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(AnGambiaeNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(AnGambiaeNS$streak_point2, na.rm = T))
AverageHoursNS[4,2:6] = c(mean(AnStephNS$streak_point2, na.rm =T), sd(AnStephNS$streak_point2, na.rm =T),
                          t.test(AnStephNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(AnStephNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(AnStephNS$streak_point2, na.rm = T))
AverageHoursNS[5,2:6] = c(mean(CxAnnulNS$streak_point2, na.rm =T), sd(CxAnnulNS$streak_point2, na.rm =T),
                          t.test(CxAnnulNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(CxAnnulNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(CxAnnulNS$streak_point2, na.rm = T))
AverageHoursNS[6,2:6] = c(mean(CxPipNS$streak_point2, na.rm =T), sd(CxPipNS$streak_point2, na.rm =T),
                          t.test(CxPipNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(CxPipNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(CxPip$streak_point2, na.rm = T))
AverageHoursNS[7,2:6] = c(mean(CxQuinqueNS$streak_point2, na.rm =T), sd(CxQuinqueNS$streak_point2, na.rm =T),
                          t.test(CxQuinqueNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(CxQuinqueNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(CxQuinqueNS$streak_point2, na.rm = T))
AverageHoursNS[8,2:6] = c(mean(CxTarNS$streak_point2, na.rm =T), sd(CxTarNS$streak_point2, na.rm =T),
                          t.test(CxTarNS$streak_point2,  conf.level = 0.95)$conf.int[1], t.test(CxTarNS$streak_point2,  conf.level = 0.95)$conf.int[2],
                          max(CxTarNS$streak_point2, na.rm = T))
#fwrite(AverageHoursNS, "SuppTable_AverageHoursThermalDanger_WithoutDroughtMask.csv")

avgHoursNS <- c(6.82, 9.42, 7.52, 9.01, 9.92, 7.67, 9.44, 6.13)

##### Supp Tables: Average streak (days) in thermal danger under different specifications #####

AverageDays = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageDays) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageDays[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                    "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageDays[1,2:6] = c(mean(AeAegypti$streak_days_point2, na.rm = T), sd(AeAegypti$streak_days_point2, na.rm = T),
                       t.test(AeAegypti$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(AeAegypti$streak_days_point2, na.rm = T))
AverageDays[2,2:6] = c(mean(AeAlbopictus$streak_days_point2, na.rm = T), sd(AeAlbopictus$streak_days_point2, na.rm = T),
                       t.test(AeAlbopictus$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(AeAlbo$streak_days_point2, na.rm = T))
AverageDays[3,2:6] = c(mean(AnGambiae$streak_days_point2, na.rm = T), sd(AnGambiae$streak_days_point2,na.rm = T),
                       t.test(AnGambiae$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(AnGambiae$streak_days_point2, na.rm = T))
AverageDays[4,2:6] = c(mean(AnSteph$streak_days_point2, na.rm = T), sd(AnSteph$streak_days_point2,na.rm = T),
                       t.test(AnSteph$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(AnSteph$streak_days_point2, na.rm = T))
AverageDays[5,2:6] = c(mean(CxAnnul$streak_days_point2, na.rm = T), sd(CxAnnul$streak_days_point2,na.rm = T),
                       t.test(CxAnnul$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(CxAnnul$streak_days_point2, na.rm = T))
AverageDays[6,2:6] = c(mean(CxPip$streak_days_point2, na.rm = T), sd(CxPip$streak_days_point2,na.rm = T),
                       t.test(CxPip$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(CxPip$streak_days_point2, na.rm = T))
AverageDays[7,2:6] = c(mean(CxQuinque$streak_days_point2, na.rm = T), sd(CxQuinque$streak_days_point2,na.rm = T),
                       t.test(CxQuinque$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(CxQuinque$streak_days_point2, na.rm = T))
AverageDays[8,2:6] = c(mean(CxTar$streak_days_point2, na.rm = T), sd(CxTar$streak_days_point2,na.rm = T),
                       t.test(CxTar$streak_days_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$streak_days_point2, conf.level = 0.95)$conf.int[2],
                       max(CxTar$streak_days_point2, na.rm = T))
#fwrite(AverageDays, "SuppTable_AverageDaysThermalDanger_MainModelSpecification.csv")

avgdays <- c(9.7, 21.93, 11.47, 17.53, 16.45, 7.73, 11.25, 6.48)

# Average days streak Without behavior thermoregulation

AverageDaysNoB = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageDaysNoB) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageDaysNoB[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                       "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")

AverageDaysNoB[1,2:6] = c(mean(AeAegypti$streak_days_nob_point2, na.rm = T), sd(AeAegypti$streak_days_nob_point2, na.rm = T),
                          t.test(AeAegypti$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(AeAegypti$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(AeAegypti$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[2,2:6] = c(mean(AeAlbopictus$streak_days_nob_point2, na.rm = T), sd(AeAlbopictus$streak_days_nob_point2, na.rm = T),
                          t.test(AeAlbopictus$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(AeAlbopictus$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(AeAlbo$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[3,2:6] = c(mean(AnGambiae$streak_days_nob_point2, na.rm = T), sd(AnGambiae$streak_days_nob_point2,na.rm = T),
                          t.test(AnGambiae$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(AnGambiae$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(AnGambiae$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[4,2:6] = c(mean(AnSteph$streak_days_nob_point2, na.rm = T), sd(AnSteph$streak_days_nob_point2,na.rm = T),
                          t.test(AnSteph$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(AnSteph$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(AnSteph$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[5,2:6] = c(mean(CxAnnul$streak_days_nob_point2, na.rm = T), sd(CxAnnul$streak_days_nob_point2,na.rm = T),
                          t.test(CxAnnul$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(CxAnnul$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(CxAnnul$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[6,2:6] = c(mean(CxPip$streak_days_nob_point2, na.rm = T), sd(CxPip$streak_days_nob_point2,na.rm = T),
                          t.test(CxPip$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(CxPip$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(CxPip$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[7,2:6] = c(mean(CxQuinque$streak_days_nob_point2, na.rm = T), sd(CxQuinque$streak_days_nob_point2,na.rm = T),
                          t.test(CxQuinque$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(CxQuinque$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(CxQuinque$streak_days_nob_point2, na.rm = T))
AverageDaysNoB[8,2:6] = c(mean(CxTar$streak_days_nob_point2, na.rm = T), sd(CxTar$streak_days_nob_point2,na.rm = T),
                          t.test(CxTar$streak_days_nob_point2, conf.level = 0.95)$conf.int[1], t.test(CxTar$streak_days_nob_point2, conf.level = 0.95)$conf.int[2],
                          max(CxTar$streak_days_nob_point2, na.rm = T))

#fwrite(AverageDaysNoB, "SuppTable_AverageDaysThermalDanger_WithoutBehavior.csv")

avgdaysNoB <- c(12, 31.56, 17.72, 29.19, 23.78, 11.91, 16.15, 7.85)

# Average streaks (days) in thermal danger without drought mask
AverageDaysNS = as.data.frame(matrix(nrow = 8, ncol = 6))
colnames(AverageDaysNS) = c("species", "mean", "sd", "lowerCI", "upperCI", "max")
AverageDaysNS[,1] = c("Ae. aegypti", "Ae. albopictus",  "An. gambiae", "An. stephensi", 
                      "Cx. annul", "Cx. pipiens", "Cx. quinque", "Cx. tarsalis")
AverageDaysNS[1,2:6] = c(mean(AeAegyptiNS$streak_days_point2, na.rm =T), sd(AeAegyptiNS$streak_days_point2, na.rm =T),
                         t.test(AeAegyptiNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(AeAegyptiNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(AeAegyptiNS$streak_days_point2, na.rm = T))
AverageDaysNS[2,2:6] = c(mean(AeAlboNS$streak_days_point2, na.rm =T), sd(AeAlboNS$streak_days_point2, na.rm =T),
                         t.test(AeAlboNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(AeAlboNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(AeAlboNS$streak_days_point2, na.rm = T))
AverageDaysNS[3,2:6] = c(mean(AnGambiaeNS$streak_days_point2, na.rm =T), sd(AnGambiaeNS$streak_days_point2, na.rm =T),
                         t.test(AnGambiaeNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(AnGambiaeNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(AnGambiaeNS$streak_days_point2, na.rm = T))
AverageDaysNS[4,2:6] = c(mean(AnStephNS$streak_days_point2, na.rm =T), sd(AnStephNS$streak_days_point2, na.rm =T),
                         t.test(AnStephNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(AnStephNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(AnStephNS$streak_days_point2, na.rm = T))
AverageDaysNS[5,2:6] = c(mean(CxAnnulNS$streak_days_point2, na.rm =T), sd(CxAnnulNS$streak_days_point2, na.rm =T),
                         t.test(CxAnnulNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(CxAnnulNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(CxAnnulNS$streak_days_point2, na.rm = T))
AverageDaysNS[6,2:6] = c(mean(CxPipNS$streak_days_point2, na.rm =T), sd(CxPipNS$streak_days_point2, na.rm =T),
                         t.test(CxPipNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(CxPipNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(CxPip$streak_days_point2, na.rm = T))
AverageDaysNS[7,2:6] = c(mean(CxQuinqueNS$streak_days_point2, na.rm =T), sd(CxQuinqueNS$streak_days_point2, na.rm =T),
                         t.test(CxQuinqueNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(CxQuinqueNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(CxQuinqueNS$streak_days_point2, na.rm = T))
AverageDaysNS[8,2:6] = c(mean(CxTarNS$streak_days_point2, na.rm =T), sd(CxTarNS$streak_days_point2, na.rm =T),
                         t.test(CxTarNS$streak_days_point2,  conf.level = 0.95)$conf.int[1], t.test(CxTarNS$streak_days_point2,  conf.level = 0.95)$conf.int[2],
                         max(CxTarNS$streak_days_point2, na.rm = T))

#fwrite(AverageDaysNS, "SuppTable_AverageDaysThermalDanger_WithoutDroughtMask.csv")

avgdaysns <- c(11.58, 22.21, 14.03, 28.21, 18.96, 18.25, 14.13, 9.03)
