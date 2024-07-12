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
AeAegypti = read.csv("AeAegypti_TSM_DroughtMask_Combined_WithElevation.csv")
AeAlbopictus = read.csv("AeAlbo_TSM_DroughtMask_Combined_WithElevation.csv")
AnGambiae = read.csv("AnGambiae_TSM_DroughtMask_Combined_WithElevation.csv")
AnSteph = read.csv("AnSteph_TSM_DroughtMask_Combined_WithElevation.csv")
CxQuinque = read.csv("CxQuinque_TSM_DroughtMask_Combined_WithElevation.csv")  
CxPip = read.csv("CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")
CxTar = read.csv("CxTarsalis_TSM_DroughtMask_Combined_WithElevation.csv")
CxAnnul = read.csv("CxAnnul_TSM_DroughtMask_Combined_WithElevation.csv")

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
  Xp <- predict(mod, data.frame(lat = newdat),type="lpmatrix") # prediction decomposed into the linear predictors
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

# to find valleys: multiply the thermal safety margins for each species 
# (column labeled 'tolerance_point2' by -1 and use find_peaks function above)


##### Ae. aegypti ######

# Subset to just records <2500 m elevation
AeAegypti = AeAegypti[AeAegypti$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAg = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAegypti, method = "REML")
newlats = seq(from = min(AeAegypti$lat), to = max(AeAegypti$lat), length.out = nrow(AeAegypti))
y = predict(gamAg, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AgPeaks = findPeaksGAM_CI(gamAg, newlats, 1000)
peaklocations = unlist(as.vector(AgPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AgPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 2 peaks. 1 located around -8 degrees latitude, 1 located around 40
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -15 & peaklocations < 10]
round(mean(p1), 2) 
round(CI(p1), 2) 
p2 = peaklocations[peaklocations > 30 & peaklocations < 50]
round(mean(p2), 2) 
round(CI(p2), 2) 

# Find valleys
AeAegypti$neg_tolerance_point2 = -AeAegypti$tolerance_point2
gamAgNeg =  gam(neg_tolerance_point2 ~ s(lat), data = AeAegypti, method = "REML")
AgValleys = findPeaksGAM_CI(gamAgNeg, newlats, 1000)
valleylocations = unlist(as.vector(AgValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AgValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. Located around -25 and 25 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -35 & valleylocations < -15]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 19 & valleylocations < 40]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAg, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAg),gamAg$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "", ylab = "",
     #ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000), 
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### Ae. albopictus ######

# Subset to just records <2500 m elevation
AeAlbo = AeAlbopictus[AeAlbopictus$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAb = gam(tolerance_point2 ~ s(lat, k = 8), data = AeAlbo, method = "REML")
newlats = seq(from = min(AeAlbo$lat), to = max(AeAlbo$lat), length.out = nrow(AeAlbo))
y = predict(gamAb, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AbPeaks = findPeaksGAM_CI(gamAb, newlats, 1000)
peaklocations = unlist(as.vector(AbPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AbPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 2 peaks. 1 located around -10 degrees latitude, another around 40 degrees
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -20 & peaklocations < 0]
round(mean(p1), 2) 
round(CI(p1), 2) 
p2 = peaklocations[peaklocations > 30 & peaklocations < 50]
round(mean(p2), 2) 
round(CI(p2), 2) 

# Find valleys
AeAlbo$neg_tolerance_point2 = -AeAlbo$tolerance_point2
gamAbNeg =  gam(neg_tolerance_point2 ~ s(lat), data = AeAlbo, method = "REML")
AbValleys = findPeaksGAM_CI(gamAbNeg, newlats, 1000)
valleylocations = unlist(as.vector(AbValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AbValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. 1 located around -25, 1 around 25 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > 15 & valleylocations < 40]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > -35 & valleylocations < -15]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAb, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAb),gamAb$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "", ylab = "",
     #ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,520), main = "", yaxt = 'n',
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) # histogram of peak locations
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))

##### An. gambiae ######

# Subset to just records <3000 ft elevation
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation),]
AnGambiae = AnGambiae[AnGambiae$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAm = gam(tolerance_point2 ~ s(lat, k = 8), data = AnGambiae, method = "REML")
newlats = seq(from = min(AnGambiae$lat), to = max(AnGambiae$lat), length.out = nrow(AnGambiae))
y = predict(gamAm, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AmPeaks = findPeaksGAM_CI(gamAm, newlats, 1000)
peaklocations = unlist(as.vector(AmPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AmPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Located around -5 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -10 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
AnGambiae$neg_tolerance_point2 = -AnGambiae$tolerance_point2
gamAmNeg =  gam(neg_tolerance_point2 ~ s(lat), data = AnGambiae, method = "REML")
AmValleys = findPeaksGAM_CI(gamAmNeg, newlats, 1000)
valleylocations = unlist(as.vector(AmValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AmValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. One around -25, another around 20 degrees latitude
# Identify mean and CI for this valley
v1 = valleylocations[valleylocations > -30 & valleylocations < -15]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 10 & valleylocations < 25]
round(mean(v2), 2) 
round(CI(v2), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamAm, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamAm),gamAm$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "", ylab = "",
     #ylab = "TSM (\u00B0C)",
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,500), yaxt = 'n', main = "",
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))



##### An. stephensi ######

AnSteph = AnSteph[!is.na(AnSteph$elevation),]
AnSteph = AnSteph[AnSteph$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamAs = gam(tolerance_point2 ~ s(lat, k = 8), data = AnSteph, method = "REML")
newlats = seq(from = min(AnSteph$lat), to = max(AnSteph$lat), length.out = nrow(AnSteph))
y = predict(gamAs, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
AsPeaks = findPeaksGAM_CI(gamAs, newlats, 1000)
peaklocations = unlist(as.vector(AsPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(AsPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# No peaks detected within range (i.e. decreasing across range)

# Find valleys
AnSteph$neg_tolerance_point2 = -AnSteph$tolerance_point2
gamAsNeg =  gam(neg_tolerance_point2 ~ s(lat), data = AnSteph, method = "REML")
AsValleys = findPeaksGAM_CI(gamAsNeg, newlats, 1000)
valleylocations = unlist(as.vector(AsValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(AsValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# No valleys detected within range (i.e. decreasing monotonically across range)

# No plot for this species

##### Cx. quinque ######

# Subset to just records 2,500 m
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation),]
CxQuinque = CxQuinque[CxQuinque$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCq = gam(tolerance_point2 ~ s(lat, k = 8), data = CxQuinque, method = "REML")
newlats = seq(from = min(CxQuinque$lat), to = max(CxQuinque$lat), length.out = nrow(CxQuinque))
y = predict(gamCq, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CqPeaks = findPeaksGAM_CI(gamCq, newlats, 1000)
peaklocations = unlist(as.vector(CqPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CqPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. Around 0
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -10 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
CxQuinque$neg_tolerance_point2 = -CxQuinque$tolerance_point2
gamCqNeg =  gam(neg_tolerance_point2 ~ s(lat), data = CxQuinque, method = "REML")
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
Xp <- predict(gamCq, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCq),gamCq$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "", ylab = "",
   #  ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000), main = "",
     yaxt = "n", xlab = "",
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Cx. pipiens #####

# Subset to just records <2,500 m elevation
CxPip = CxPip[!is.na(CxPip$elevation),]
CxPip = CxPip[CxPip$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCp = gam(tolerance_point2 ~ s(lat, k = 8), data = CxPip, method = "REML")
newlats = seq(from = min(CxPip$lat), to = max(CxPip$lat), length.out = nrow(CxPip))
y = predict(gamCp, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CpPeaks = findPeaksGAM_CI(gamCp, newlats, 1000)
peaklocations = unlist(as.vector(CpPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CpPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peak. One around -2 degrees latitude
# Identify mean and CI for this peak
p1 = peaklocations[peaklocations > -15 & peaklocations < 5]
round(mean(p1), 2) 
round(CI(p1), 2) 


# Find valleys
CxPip$neg_tolerance_point2 = -CxPip$tolerance_point2
gamCpNeg =  gam(neg_tolerance_point2 ~ s(lat), data = CxPip, method = "REML")
CpValleys = findPeaksGAM_CI(gamCpNeg, newlats, 1000)
valleylocations = unlist(as.vector(CpValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CpValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 2 valleys. Located around -35 and 30 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -40 & valleylocations < -20]
round(mean(v1), 2) 
round(CI(v1), 2) 
v2 = valleylocations[valleylocations > 19 & valleylocations < 40]
round(mean(v2), 2) 
round(CI(v2), 2)

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCp, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCp),gamCp$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), main = "", 
     xlab = "", ylab = "",
    # ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000), yaxt = "n",
     main = "", xlab = "",
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))


##### Cx. tarsalis #####

# Subset to just records <2500 m elevation
CxTar = CxTar[!is.na(CxTar$elevation),]
CxTar = CxTar[CxTar$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCt = gam(tolerance_point2 ~ s(lat, k = 8), data = CxTar, method = "REML")
newlats = seq(from = min(CxTar$lat), to = max(CxTar$lat), length.out = nrow(CxTar))
y = predict(gamCt, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CtPeaks = findPeaksGAM_CI(gamCt, newlats, 1000)
peaklocations = unlist(as.vector(CtPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CtPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

## 2 peaks detected - one around 20, one around 45
p1 = peaklocations[peaklocations > 20 & peaklocations < 25]
round(mean(p1), 2) 
round(CI(p1), 2) 
p2 = peaklocations[peaklocations > 42 & peaklocations < 50]
round(mean(p2), 2) 
round(CI(p2), 2) 

# Find valleys
CxTar$neg_tolerance_point2 = -CxTar$tolerance_point2
gamCtNeg =  gam(neg_tolerance_point2 ~ s(lat), data = CxTar, method = "REML")
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
Xp <- predict(gamCt, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCt),gamCt$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), 
     xlab = "", ylab = "",
     #ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,300), yaxt = 'n',
     xlab = "", ylab = "", main = "",
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))







##### Cx. annulirostris ######

# Subset to just records <2,500 m elevation
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation),]
CxAnnul = CxAnnul[CxAnnul$elevation < 2500,]

# Run GAM of TSM using occurence record data and use to predict across latitude
gamCa = gam(tolerance_point2 ~ s(lat, k = 8), data = CxAnnul, method = "REML")
newlats = seq(from = min(CxAnnul$lat), to = max(CxAnnul$lat), length.out = nrow(CxAnnul))
y = predict(gamCa, data.frame(lat = newlats), se.fit = TRUE)

# Identify locations of peaks and valleys
CaPeaks = findPeaksGAM_CI(gamCa, newlats, 1000)
peaklocations = unlist(as.vector(CaPeaks$peaklocs))

# Number of peaks detected
numpeaks = as.integer(unlist(as.vector(CaPeaks$npeaks)))
hist(numpeaks) # Histogram of number of peaks
hist(peaklocations) # histogram of peak locations

# 1 peaks. Located around -10 degrees latitude
# Identify mean and CI for these peaks
p1 = peaklocations[peaklocations > -15 & peaklocations < 0]
round(mean(p1), 2) 
round(CI(p1), 2) 

# Find valleys
CxAnnul$neg_tolerance_point2 = -CxAnnul$tolerance_point2
gamCaNeg =  gam(neg_tolerance_point2 ~ s(lat), data = CxAnnul, method = "REML")
CaValleys = findPeaksGAM_CI(gamCaNeg, newlats, 1000)
valleylocations = unlist(as.vector(CaValleys$peaklocs))
numvalleys = as.integer(unlist(as.vector(CaValleys$npeaks))) # Number of valleys detected
hist(numvalleys) # Histogram of number of valleys
hist(valleylocations) # Histogram of valleys

# 1 valleys. Located around -30 degrees latitude
# Identify mean and CI for these valleys
v1 = valleylocations[valleylocations > -40 & valleylocations < -25]
round(mean(v1), 2) 
round(CI(v1), 2) 

# Plot simulated GAMs (only 100)
par(mar= c(3.5, 3.5, 1, 1), mgp = c(2.2,1,0))
set.seed(1)
Xp <- predict(gamCa, data.frame(lat = newlats),type="lpmatrix") 
br <- rmvn(100,coef(gamCa),gamCa$Vp) 
res <- matrix(NA, ncol=100, nrow=length(newlats)) 
for(i in 1:100) {res[,i] <- Xp %*% br[i,]}
res2 <- list()   # convert to list for easier plotting
for(i in 1:ncol(res)) {            
  res2[[i]] <- res[ , i]}

plot(NA, xlim = range(newlats), ylim = range(res), main = "",
     xlab = "", ylab = "",
     #ylab = "TSM (\u00B0C)", 
     cex.lab = 1.8, cex.axis = 1.8)
sapply(seq_along(res2), function(z) lines(res2[[z]] ~ newlats, col = "gray80"))
lines(x=newlats, y = y$fit, col = "black", lwd = 1.2)

# To show histogram above this plot 
hist(peaklocations, xlim = range(newlats), ylim = c(0,1000), yaxt = 'n', 
     main = "", xlab = "",
     cex.lab = 1.3, cex.axis = 1.3, col = alpha("#276DC2", 0.6)) 
hist(valleylocations, add = T, col = alpha("#F20000", 0.6))

#### Multi-line plot of peaks and valleys (figure 2 of manuscript) #####

# Pull in dataframe created from code above
pv = read.csv("~/Documents/Current Projects/WarmingTolerance/DataFiles/VectorThermalSafety_PeaksValleys.csv", header = T)
pv$Species= factor(pv$Species, levels = rev(unique(pv$Species)))

colors = vector(mode = "character", length = nrow(pv))
colors[which(pv$Parameter == "LatRange")] <- alpha("black", 0.0)
colors[which(pv$Parameter == "Peak1")] <- "#113486"
colors[which(pv$Parameter == "Peak2")] <- "#113486"
colors[which(pv$Parameter == "Peak3")] <- "#113486"
colors[which(pv$Parameter == "Valley1")] <- "firebrick"
colors[which(pv$Parameter == "Valley2")] <- "firebrick"

ggplot(pv, aes(x=Species, y=Mean)) + 
 scale_y_continuous("Latitude", seq(from = -40, to = 70, by = 10)) + 
  scale_x_discrete(position = "top") +
  geom_point(stat="identity",  size = 3, col = colors,
             position=position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin=Lower95, ymax=Upper95), size = 1, col = colors, 
                position=position_dodge(.9), width = 0.9) + 
  ggtitle("Maxima and minima in thermal safety") + 
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
