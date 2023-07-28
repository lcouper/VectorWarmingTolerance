#### Vector Thermal Safety Finding Peaks & Valleys in TSMs 

# This script is used to detects peaks and valleys in fitted GAMS of TSMs across latitude 
# Code overview:
# 1. Run GAM to model relationship between TSM and latitude 
# 2. Predict TSM across latitude. Create confidence intervals
# 3. Estimate locations of peaks and valleys and uncertainty
# 4. Visually inspect locations of peaks/valleys
# 5. Identify location of peaks in fitted GAM and estimate uncertainty (for 1000 samples)
# 6. Identify location of valleys in fitted GAM and estimate uncertainty (for 1000 samples)

#### Load libraries and pull in data ####
library(mgcv)
library(ggplot2)
library(Rmisc)

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files
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

SpeciesColors = c("#f46d43", "#fdae61", "#fed439ff", "#7cae00", "#309143", 
                  "#ab041b", "#ec3c30", 
                  "#abd9e9", "#74add1", "#313695", "#800080", "#dda0dd")
SpeciesList  = c("Aedes_aegypti", "Aedes_albopictus", "Aedes_camptorhynchus", 
                 "Aedes_triseriaatus", "Aedes_vexans",
                 "Anopheles_gambiae", "Anopheles_stephensi", 
                 "Culex_annulirostris", "Culex_pipiens", 
                 "Culex_quinquefasciatus", "Culex_tarsalis", "Culex_theileri")

##### Functions for finding peaks and valleys #####

# code adapted from Pinskey et al. 2019 
findPeaksGAM_CI <- function(mod, newdat, n){
  library(mgcv)
  set.seed(1)
  Xp <- predict(mod, data.frame(latitude = newdat),type="lpmatrix") # prediction decomposed into the linear predictors
  br <- rmvn(n,coef(mod),mod$Vp) ## n replicate parameter vectors sampled with uncertainty
  res <- matrix(NA, ncol=n, nrow=length(newdat)) # to hold the predictions
  for(i in 1:n) res[,i] <- Xp %*% br[i,] ## replicate predictions
  peaks <- data.frame(npeaks=rep(NA, n), peaklocs=I(vector('list', n)), y=I(vector('list', n))) # hold the number and location of peaks
  for(i in 1:n){
    temp <- ggpmisc:::find_peaks(res[,i], span = 3)
    index = which(temp == TRUE)
    peaks$npeaks[i] <- length(index)
    peaks$peaklocs[[i]] <- newdat[index]
    peaks$y[[i]] <- res[,i]
  }
  return(peaks)
}


# For GAMs with elevation as well
findPeaksGAM_CI_elev <- function(mod, newdat, n){
  library(mgcv)
  set.seed(1)
  Xp <- predict(mod,newdat,type="lpmatrix") # prediction decomposed into the linear predictors
  br <- rmvn(n,coef(mod),mod$Vp) ## n replicate parameter vectors sampled with uncertainty
  res <- matrix(NA, ncol=n, nrow=nrow(newdat)) # to hold the predictions
  for(i in 1:n) res[,i] <- Xp %*% br[i,] ## replicate predictions
  peaks <- data.frame(npeaks=rep(NA, n), peaklocs=I(vector('list', n)), y=I(vector('list', n))) # hold the number and location of peaks
  for(i in 1:n){
    temp <- ggpmisc:::find_peaks(res[,i], span = 3)
    index = which(temp == TRUE)
    peaks$npeaks[i] <- length(index)
    peaks$peaklocs[[i]] <- newdat[index,1]
    peaks$y[[i]] <- res[,i]
  }
  return(peaks)
}

# to find valleys: multiple the thermal safety margins for each species 
# (column labaled 'tolerance_point' by -1 and use find_peaks function above)


##### Ae. aegypti ######

# Subset to just records <3000 ft elevation
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegypti = AeAegypti[AeAegypti$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAg = gam(tolerance_point ~ s(latitude, k = 8), data = AeAegypti, method = "REML")
newlats = seq(from = min(AeAegypti$latitude), to = max(AeAegypti$latitude), length.out = nrow(AeAegypti))
y = predict(gamAg, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AgPeaks = findPeaksGAM_CI(gamAg, newlats, 1000)
peaklocations = unlist(as.vector(AgPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AgPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -8 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -15 & peaklocations < 10]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
AeAegypti$neg_tolerance_point = -AeAegypti$tolerance_point
gamAgNeg =  gam(neg_tolerance_point ~ s(latitude), data = AeAegypti, method = "REML")
AgValleys = findPeaksGAM_CI(gamAgNeg, newlats, 1000)
valleylocations = unlist(as.vector(AgValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AgValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. Located around -30 and 25 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -35 & valleylocations < -20]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 19 & valleylocations < 40]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAg, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAg),gamAg$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### Ae. albopictus ######

# Subset to just records <3000 ft elevation
AeAlbo = AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbo = AeAlbopictus[AeAlbopictus$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAb = gam(tolerance_point ~ s(latitude, k = 8), data = AeAlbo, method = "REML")
newlats = seq(from = min(AeAlbo$latitude), to = max(AeAlbo$latitude), length.out = nrow(AeAlbo))
y = predict(gamAb, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AbPeaks = findPeaksGAM_CI(gamAb, newlats, 1000)
peaklocations = unlist(as.vector(AbPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AbPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -10 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -20 & peaklocations < 0]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
AeAlbo$neg_tolerance_point = -AeAlbo$tolerance_point
gamAbNeg =  gam(neg_tolerance_point ~ s(latitude), data = AeAlbo, method = "REML")
AbValleys = findPeaksGAM_CI(gamAbNeg, newlats, 1000)
valleylocations = unlist(as.vector(AbValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AbValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valley. Located around 25 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > 15 & valleylocations < 40]
round(mean(v1), 2) 
round(CI(v1), 2) 


# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAb, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAb),gamAb$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,520),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) # histogram of peak locations
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### Ae. camptorhynchus ######

# Subset to just records <3000 ft elevation
AeCamp = AeCamp[!is.na(AeCamp$elevation_ft),]
AeCamp = AeCamp[AeCamp$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAc = gam(tolerance_point ~ s(latitude, k = 8), data = AeCamp, method = "REML")
newlats = seq(from = min(AeCamp$latitude), to = max(AeCamp$latitude), length.out = nrow(AeCamp))
y = predict(gamAc, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AcPeaks = findPeaksGAM_CI(gamAc, newlats, 1000)
peaklocations = unlist(as.vector(AcPeaks$peaklocs))

# No peaks detected within range (i.e. range edges are peaks)

# Find valleys
AeCamp$neg_tolerance_point = -AeCamp$tolerance_point
gamAcNeg =  gam(neg_tolerance_point ~ s(latitude), data = AeCamp, method = "REML")
AcValleys = findPeaksGAM_CI(gamAcNeg, newlats, 1000)
valleylocations = unlist(as.vector(AcValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AcValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valleys. Located around -30 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > -35 & valleylocations < -25]
round(mean(v1), 2) 
round(CI(v1), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAc, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAc),gamAc$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,400),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Ae. triseriatus ######

# Subset to just records <3000 ft elevation
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTri = AeTri[AeTri$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAt = gam(tolerance_point ~ s(latitude, k = 8), data = AeTri, method = "REML")
newlats = seq(from = min(AeTri$latitude), to = max(AeTri$latitude), length.out = nrow(AeTri))
y = predict(gamAt, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AtPeaks = findPeaksGAM_CI(gamAt, newlats, 1000)
peaklocations = unlist(as.vector(AtPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AtPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# No peaks detected within range (i.e. range edges are peaks)

# Find valleys
AeTri$neg_tolerance_point = -AeTri$tolerance_point
gamAtNeg =  gam(neg_tolerance_point ~ s(latitude), data = AeTri, method = "REML")
AtValleys = findPeaksGAM_CI(gamAtNeg, newlats, 1000)
valleylocations = unlist(as.vector(AtValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AtValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valley. Located around 32 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > 29 & valleylocations < 37]
round(mean(v1), 2) 
round(CI(v1), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAc, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAc),gamAc$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,600),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))










##### Ae. vexans ######

# Subset to just records <3000 ft elevation
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexans = AeVexans[AeVexans$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAv = gam(tolerance_point ~ s(latitude, k = 8), data = AeVexans, method = "REML")
newlats = seq(from = min(AeVexans$latitude), to = max(AeVexans$latitude), length.out = nrow(AeVexans))
y = predict(gamAv, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AvPeaks = findPeaksGAM_CI(gamAv, newlats, 1000)
peaklocations = unlist(as.vector(AvPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AvPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -5 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -15 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
AeVexans$neg_tolerance_point = -AeVexans$tolerance_point
gamAvNeg =  gam(neg_tolerance_point ~ s(latitude), data = AeVexans, method = "REML")
AvValleys = findPeaksGAM_CI(gamAvNeg, newlats, 1000)
valleylocations = unlist(as.vector(AvValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AvValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valley. Located around 30 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > 17 & valleylocations < 40]
round(mean(v1), 2) 
round(CI(v1), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAv, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAv),gamAv$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,400),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### An. gambiae ######

# Subset to just records <3000 ft elevation
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiae = AnGambiae[AnGambiae$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAm = gam(tolerance_point ~ s(latitude, k = 8), data = AnGambiae, method = "REML")
newlats = seq(from = min(AnGambiae$latitude), to = max(AnGambiae$latitude), length.out = nrow(AnGambiae))
y = predict(gamAm, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AmPeaks = findPeaksGAM_CI(gamAm, newlats, 1000)
peaklocations = unlist(as.vector(AmPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AmPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around 0 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -5 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
AnGambiae$neg_tolerance_point = -AnGambiae$tolerance_point
gamAmNeg =  gam(neg_tolerance_point ~ s(latitude), data = AnGambiae, method = "REML")
AmValleys = findPeaksGAM_CI(gamAmNeg, newlats, 1000)
valleylocations = unlist(as.vector(AmValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AmValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valley. Located around -20 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > -30 & valleylocations < -10]
round(mean(v1), 2) 
round(CI(v1), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAm, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAm),gamAm$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,500),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### An. stephensi ######

# Subset to just records <3000 ft elevation
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnSteph = AnSteph[AnSteph$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAs = gam(tolerance_point ~ s(latitude, k = 8), data = AnSteph, method = "REML")
newlats = seq(from = min(AnSteph$latitude), to = max(AnSteph$latitude), length.out = nrow(AnSteph))
y = predict(gamAs, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AsPeaks = findPeaksGAM_CI(gamAs, newlats, 1000)
peaklocations = unlist(as.vector(AsPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AsPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# No peaks detected within range (i.e. decreasing across range)

# Find valleys
AnSteph$neg_tolerance_point = -AnSteph$tolerance_point
gamAsNeg =  gam(neg_tolerance_point ~ s(latitude), data = AnSteph, method = "REML")
AsValleys = findPeaksGAM_CI(gamAsNeg, newlats, 1000)
valleylocations = unlist(as.vector(AsValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AsValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# No valleys detected within range (i.e. decreasing across range)

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAs, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAs),gamAs$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,100),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))

##### Cx. quinque ######

# Subset to just records <3000 ft elevation
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinque = CxQuinque[CxQuinque$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCq = gam(tolerance_point ~ s(latitude, k = 8), data = CxQuinque, method = "REML")
newlats = seq(from = min(CxQuinque$latitude), to = max(CxQuinque$latitude), length.out = nrow(CxQuinque))
y = predict(gamCq, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CqPeaks = findPeaksGAM_CI(gamCq, newlats, 1000)
peaklocations = unlist(as.vector(CqPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CqPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -5 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -15 & peaklocations < 0]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
CxQuinque$neg_tolerance_point = -CxQuinque$tolerance_point
gamCqNeg =  gam(neg_tolerance_point ~ s(latitude), data = CxQuinque, method = "REML")
CqValleys = findPeaksGAM_CI(gamCqNeg, newlats, 1000)
valleylocations = unlist(as.vector(CqValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CqValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. Located around -30 and 30 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -40 & valleylocations < -20]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 22 & valleylocations < 40]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCq, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCq),gamCq$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Cx. pipiens #####

# Subset to just records <3000 ft elevation
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPip = CxPip[CxPip$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCp = gam(tolerance_point ~ s(latitude, k = 8), data = CxPip, method = "REML")
newlats = seq(from = min(CxPip$latitude), to = max(CxPip$latitude), length.out = nrow(CxPip))
y = predict(gamCp, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CpPeaks = findPeaksGAM_CI(gamCp, newlats, 1000)
peaklocations = unlist(as.vector(CpPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CpPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -2 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -10 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
CxPip$neg_tolerance_point = -CxPip$tolerance_point
gamCpNeg =  gam(neg_tolerance_point ~ s(latitude), data = CxPip, method = "REML")
CpValleys = findPeaksGAM_CI(gamCpNeg, newlats, 1000)
valleylocations = unlist(as.vector(CpValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CpValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 3 valleys. Located around -35, 20, and 45 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -40 & valleylocations < -20]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 19 & valleylocations < 30]
round(mean(v2), 2) 
round(CI(v2), 2)
v3 = valleylocations[valleylocations > 30 & valleylocations < 50]
round(mean(v3), 2) 
round(CI(v3), 2)

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCp, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCp),gamCp$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Cx. tarsalis #####

# Subset to just records <3000 ft elevation
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTar = CxTar[CxTar$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCt = gam(tolerance_point ~ s(latitude, k = 8), data = CxTar, method = "REML")
newlats = seq(from = min(CxTar$latitude), to = max(CxTar$latitude), length.out = nrow(CxTar))
y = predict(gamCt, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CtPeaks = findPeaksGAM_CI(gamCt, newlats, 1000)
peaklocations = unlist(as.vector(CtPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CtPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

## No peaks detected. Increasing across range

# Find valleys
CxTar$neg_tolerance_point = -CxTar$tolerance_point
gamCtNeg =  gam(neg_tolerance_point ~ s(latitude), data = CxTar, method = "REML")
CtValleys = findPeaksGAM_CI(gamCtNeg, newlats, 1000)
valleylocations = unlist(as.vector(CtValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CtValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valley. Around 33 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > 30 & valleylocations < 37]
round(mean(v1), 2) 
round(CI(v1), 2) 


# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCt, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCt),gamCt$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,300),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))







##### Cx. theileri #####

# Subset to just records <3000 ft elevation
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxTh = CxTh[CxTh$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCth = gam(tolerance_point ~ s(latitude, k = 8), data = CxTh, method = "REML")
newlats = seq(from = min(CxTh$latitude), to = max(CxTh$latitude), length.out = nrow(CxTh))
y = predict(gamCth, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CthPeaks = findPeaksGAM_CI(gamCth, newlats, 1000)
peaklocations = unlist(as.vector(CthPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CthPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# no peaks detected 

# Find valleys
CxTh$neg_tolerance_point = -CxTh$tolerance_point
gamCthNeg =  gam(neg_tolerance_point ~ s(latitude), data = CxTh, method = "REML")
CthValleys = findPeaksGAM_CI(gamCthNeg, newlats, 1000)
valleylocations = unlist(as.vector(CthValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CthValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# no valleys detected

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCth, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCth),gamCth$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,300),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Cx. annulirostris ######

# Subset to just records <3000 ft elevation
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnul = CxAnnul[CxAnnul$elevation_ft < 3000,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCa = gam(tolerance_point ~ s(latitude, k = 8), data = CxAnnul, method = "REML")
newlats = seq(from = min(CxAnnul$latitude), to = max(CxAnnul$latitude), length.out = nrow(CxAnnul))
y = predict(gamCa, data.frame(latitude = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CaPeaks = findPeaksGAM_CI(gamCa, newlats, 1000)
peaklocations = unlist(as.vector(CaPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CaPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 2 peaks. Located around -20 and -2 degrees latitude
# Identify mean and CI for these peaks
p1 = peaklocations[peaklocations > -35 & peaklocations < -15]
round(mean(p1), 2) 
round(CI(p1), 2) 
p2 = peaklocations[peaklocations > -8 & peaklocations < 8]
round(mean(p2), 2) 
round(CI(p2), 2) 

# Find valleys
CxAnnul$neg_tolerance_point = -CxAnnul$tolerance_point
gamCaNeg =  gam(neg_tolerance_point ~ s(latitude), data = CxAnnul, method = "REML")
CaValleys = findPeaksGAM_CI(gamCaNeg, newlats, 1000)
valleylocations = unlist(as.vector(CaValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CaValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. Located around -30 and-12 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -40 & valleylocations < -25]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > -22 & valleylocations < -8]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCa, data.frame(latitude = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCa),gamCa$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "Latitude", ylab = "TSM (\u00B0C)", cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000),
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))

#### Multi-line plot of peaks and valleys #####

# Pull in dataframe created from code above
pv = read.csv("~/Documents/Current Projects/WarmingTolerance/DataFiles/VectorThermalSafety_PeaksValleys.csv", header = T)
pv$Species= factor(pv$Species, levels = rev(unique(pv$Species)))

colors = vector(mode = "character", length = nrow(pv))
colors[which(pv$Parameter == "LatRange")] <- alpha("black", 0.3)
colors[which(pv$Parameter == "Peak")] <- "#0D0887FF"
colors[which(pv$Parameter == "Valley")] <- "darkred"

ggplot(pv, aes(x=Species, y=Mean)) + 
 scale_y_continuous("Latitude", seq(from = -40, to = 70, by = 10)) + 
  scale_x_discrete(position = "top") +
  geom_point(stat="identity",  size = 3, col = colors,
             position=position_dodge(width = 1)) +
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95), size = 1, col = colors, 
                position=position_dodge(.9), width = 0) + 
  ggtitle("Peaks and valleys in thermal safety") + 
  theme_minimal() +  coord_flip()  + 
  labs(x = " ", y = " ") + 
  theme(axis.text.y = element_text(size=14, face = "italic"), 
        axis.text.x = element_text(size = 16),
        axis.title = element_text(size = 16),
        legend.text=element_text(size=15),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 24),
        legend.position= "none",
        panel.border = element_rect(colour = "black", fill = NA))
