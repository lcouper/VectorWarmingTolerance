##### Vector thermal safety Data pre-processing ####

# This script is used to:
# Rescale TSMs from 0 - 1 for each species
# Classify the region associated with each occurrence point
# Plot relative TSMs by region

##### Load data frames #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files with seasonality
AeAegypti = read.csv("AeAegypti_TSM_WithSeasonality.csv")
AeAlbopictus = read.csv("AeAlbopictus_TSM_WithSeasonality.csv")
AeCamp = read.csv("AeCamp_TSM_WithSeasonality.csv")
AeTri = read.csv("AeTri_TSM_WithSeasonality.csv")
AeVexans = read.csv("AeVexans_TSM_WithSeasonality.csv")
AnGambiae = read.csv("AnGambiae_TSM_WithSeasonality.csv")
AnSteph = read.csv("AnSteph_TSM_WithSeasonality.csv")
CxQuinque = read.csv("CxQuinque_TSM_WithSeasonality.csv")  
CxPip = read.csv("CxPipiens_TSM_WithSeasonality.csv")
CxTar = read.csv("CxTarsalis_TSM_WithSeasonality.csv")
CxTh = read.csv("CxTheileri_TSM_WithSeasonality.csv")
CxAnnul = read.csv("CxAnnul_TSM_WithSeasonality.csv")

# Remove any samples without elevation data or from high elevation (>3000 ft)
AeAegypti = AeAegypti[!is.na(AeAegypti$elevation_ft),]
AeAegypti = AeAegypti[AeAegypti$elevation_ft < 3000,] 
AeAlbopictus= AeAlbopictus[!is.na(AeAlbopictus$elevation_ft),]
AeAlbopictus = AeAlbopictus[AeAlbopictus$elevation_ft < 3000,] 
AeCamp= AeCamp[!is.na(AeCamp$elevation_ft),]
AeCamp = AeCamp[AeCamp$elevation_ft < 3000,] 
AeTri = AeTri[!is.na(AeTri$elevation_ft),]
AeTri = AeTri[AeTri$elevation_ft < 3000,] 
AeVexans = AeVexans[!is.na(AeVexans$elevation_ft),]
AeVexans = AeVexans[AeVexans$elevation_ft < 3000,] 
AnGambiae = AnGambiae[!is.na(AnGambiae$elevation_ft),]
AnGambiae = AnGambiae[AnGambiae$elevation_ft < 3000,] 
AnSteph = AnSteph[!is.na(AnSteph$elevation_ft),]
AnSteph = AnSteph[AnSteph$elevation_ft < 3000,] 
CxQuinque = CxQuinque[!is.na(CxQuinque$elevation_ft),]
CxQuinque = CxQuinque[CxQuinque$elevation_ft < 3000,] 
CxPip = CxPip[!is.na(CxPip$elevation_ft),]
CxPip = CxPip[CxPip$elevation_ft < 3000,] 
CxTar = CxTar[!is.na(CxTar$elevation_ft),]
CxTar = CxTar[CxTar$elevation_ft < 3000,] 
CxTh = CxTh[!is.na(CxTh$elevation_ft),]
CxTh = CxTh[CxTh$elevation_ft < 3000,] 
CxAnnul = CxAnnul[!is.na(CxAnnul$elevation_ft),]
CxAnnul = CxAnnul[CxAnnul$elevation_ft < 3000,] 

#### Calculate relative TSMs for each species #####

# Scale each species TSM from 0 to 1 so they are more directly comparable
library(scales)

AeAegypti$ScaledTSM = rescale(AeAegypti$tolerance_point)
AeAlbopictus$ScaledTSM = rescale(AeAlbopictus$tolerance_point)
AeCamp$ScaledTSM = rescale(AeCamp$tolerance_point)
AeTri$ScaledTSM = rescale(AeTri$tolerance_point)
AeVexans$ScaledTSM = rescale(AeVexans$tolerance_point)
AnGambiae$ScaledTSM = rescale(AnGambiae$tolerance_point)
AnSteph$ScaledTSM = rescale(AnSteph$tolerance_point)
CxQuinque$ScaledTSM = rescale(CxQuinque$tolerance_point)
CxPip$ScaledTSM = rescale(CxPip$tolerance_point)
CxTar$ScaledTSM = rescale(CxTar$tolerance_point)
CxTh$ScaledTSM = rescale(CxTh$tolerance_point)
CxAnnul$ScaledTSM = rescale(CxAnnul$tolerance_point)

# Bind all species into single dataframe
Mosqs = rbind.data.frame(AeAegypti, AeAlbopictus, AeCamp, AeTri, AeVexans,
        AnGambiae, AnSteph, CxQuinque, CxPip, CxTar, CxTh, CxAnnul)

#### TSMs by region (tropic, subtropic, temperate) ####
# classify lat/longs into different regions
Region = vector(length = nrow(Mosqs))
for (i in 1:length(Region)){
     if (abs(Mosqs$latitude[i]) <= 23.5)
     {Region[i] = "Tropical"}
  else if (abs(Mosqs$latitude[i]) > 23.5 & abs(Mosqs$latitude[i]) < 35.0)
  {Region[i] = "Subtropical"}
  else if (abs(Mosqs$latitude[i]) >= 35.0 & abs(Mosqs$latitude[i]) < 66.5)
  {Region[i] = "Temperate"}
}

Mosqs$Region = factor(Region, levels = c("Tropical", "Subtropical", "Temperate"))
Mosqs$BinnedLat = cut(Mosqs$latitude, breaks = c(-45, -35, -23.5, -13, 0, 13, 23.5, 35, 45, 67))

boxplot(data = Mosqs, ScaledTSM ~ Region, 
        xlab = " ", ylab = "Relative Thermal Safety Margin",
        col = c("#F52549", "#FFD54D", "#99BE1B"),
        cex.lab=1.6, cex.axis=1.6, frame = F)

boxplot(data = Mosqs, ScaledTSM ~ BinnedLat, 
        xlab = " ", ylab = "Relative Thermal Safety Margin",
       col = c("#99BE1B", "#FFD54D", 
               "#F52549", "#F52549", "#F52549", "#F52549",
               "#FFD54D","#99BE1B", "#99BE1B"),
        cex.lab=1.6, cex.axis=1.6, frame = F, boxwex = 0.6)

##### TSMs by Biome #####

# Note: assignment of occurrence points to biomes was done in QGIS
# See https://github.com/lcouper/VectorThermalSafety/blob/main/Biogeography/README.md
# for description of process and shapefile

biomes = read.csv("~/Documents/Current Projects/WarmingTolerance/DataFiles/Biogeography/VectorTSM_Biome_Classifications.csv", header = T)
# merge with TSM data in Mosqs dataframe
MosqsBiomes = merge(Mosqs, biomes, by = c("species.pop", "latitude", "longitude"))
MosqsBiomes$biome = as.factor(MosqsBiomes$biome)

# remove records with missing biomes      
MosqsBiomes = MosqsBiomes[!is.na(MosqsBiomes$biome),]

# remove records with likely mis-classificiont (i.e., "Mangroves")
# or very rare (tundra, taiga, montane grasslands, tropical & subtropical coniferous forests)
MosqsBiomes = MosqsBiomes[MosqsBiomes$biome != "Mangroves",]
MosqsBiomes = MosqsBiomes[MosqsBiomes$biome != "Tundra",]
MosqsBiomes = MosqsBiomes[MosqsBiomes$biome != "Boreal Forests/Taiga",]
MosqsBiomes = MosqsBiomes[MosqsBiomes$biome != "Boreal Forests/Taiga",]


# Reclassify biome name (quite long) with number (using WWF biome #'s)
MosqsBiomes$biome = recode(MosqsBiomes$biome, 
               "Tropical & Subtropical Moist Broadleaf Forests" = '1',
               "Tropical & Subtropical Dry Broadleaf Forests" = '2',
                  "Tropical & Subtropical Coniferous Forests" = '3',
               "Temperate Broadleaf & Mixed Forests" = '4',
               "Temperate Conifer Forests" = '5',
                   "Boreal Forests/Taiga" = '6',
               "Tropical & Subtropical Grasslands, Savannas & Shrublands" = '7',
               "Temperate Grasslands, Savannas & Shrublands" = '8',
               "Flooded Grasslands & Savannas" = '9',
                  "Montane Grasslands & Shrublands" = '10',
                  "Tundra" = '11',
               "Mediterranean Forests, Woodlands & Scrub" = '12',
               "Deserts & Xeric Shrublands" = '13',
                  "Mangroves" = '14')

# drop unused levels
MosqsBiomes$biome = droplevels(MosqsBiomes$biome)
# re-order factor levels
MosqsBiomes$biome = factor(MosqsBiomes$biome, levels = c(1,2,4,5,12, 7:9, 13))

par(mar = c(4.5, 4.5, 1, 1))
boxplot(data = MosqsBiomes, ScaledTSM ~ biome,
        xlab = "Biome", ylab = "Relative Thermal Safety Margin",
        col = c("#113900", "#2e7921",
                "#326d76","#00aeac",
                "#cc6e1f",
                "#6c4d02","#9d8613","#cdc904",
                "#eba34d"),
        cex.lab=1.6, cex.axis=1.6, frame = F,
        at = c(0, 1, 3, 4, 6, 8,9,10, 12))


