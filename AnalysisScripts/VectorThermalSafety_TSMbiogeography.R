##### Vector thermal safety marings biogeographical analysis ####

# This script is used to:
# Rescale TSMs from 0 - 1 for each species
# Classify the region associated with each occurrence point
# Plot relative TSMs by region

##### Load packages ######

library(ggplot2)
library(maptools)
library(terra)
library(sf)
library(dplyr)
library(rnaturalearth)
library(gghalves)
library(scales)
library(data.table)

##### Load data frames #####

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files with drought mask applied
AeAegypti = read.csv("AeAegypti_TSM_DroughtMask_Combined_WithElevation.csv")
AeAlbo = read.csv("AeAlbo_TSM_DroughtMask_Combined_WithElevation.csv")
AnGambiae = read.csv("AnGambiae_TSM_DroughtMask_Combined_WithElevation.csv")
AnSteph = read.csv("AnSteph_TSM_DroughtMask_Combined_WithElevation.csv")
CxQuinque = read.csv("CxQuinque_TSM_DroughtMask_Combined_WithElevation.csv")  
CxPip = read.csv("CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")
CxTar = read.csv("CxTarsalis_TSM_DroughtMask_Combined_WithElevation.csv")
CxAnnul = read.csv("CxAnnul_TSM_DroughtMask_Combined_WithElevation.csv")
#Mosqs = read.csv("AllSpecies_TSM_DroughtMask_Combined_WithElevation.csv")[,-1] # this data frame also created below

# Remove any high elevation samples (>2500 m)
AeAegypti = AeAegypti[AeAegypti$elevation < 2500,] 
AeAlbo= AeAlbo[AeAlbo$elevation < 2500,] 
AnGambiae = AnGambiae[AnGambiae$elevation < 2500,] 
AnSteph = AnSteph[AnSteph$elevation < 2500,] 
CxQuinque = CxQuinque[CxQuinque$elevation < 2500,] 
CxPip = CxPip[CxPip$elevation < 2500,] 
CxTar = CxTar[CxTar$elevation < 2500,] 
CxAnnul = CxAnnul[CxAnnul$elevation < 2500,] 

###### Calculate relative TSMs for each species #####

# Scale each species TSM from 0 to 1 so they are more directly comparable
AeAegypti$ScaledTSM = scales::rescale(AeAegypti$tolerance_point2)
AeAlbo$ScaledTSM = scales::rescale(AeAlbo$tolerance_point2)
AnGambiae$ScaledTSM = scales::rescale(AnGambiae$tolerance_point2)
AnSteph$ScaledTSM = scales::rescale(AnSteph$tolerance_point2)
CxQuinque$ScaledTSM = scales::rescale(CxQuinque$tolerance_point2)
CxPip$ScaledTSM = scales::rescale(CxPip$tolerance_point2)
CxTar$ScaledTSM = scales::rescale(CxTar$tolerance_point2)
CxAnnul$ScaledTSM = scales::rescale(CxAnnul$tolerance_point2)

# Bind all species into single dataframe
Mosqs = rbind.data.frame(AeAegypti, AeAlbo, 
        AnGambiae, AnSteph, CxQuinque, CxPip, CxTar,CxAnnul)

### Add regional designation for all species  (tropic, subtropic, temperate) 

# Note this part already done if loading in the 'Mosq' dataframe
# classify lat/longs into different regions
Region = vector(length = nrow(Mosqs))
for (i in 1:length(Region)){
  if (abs(Mosqs$lat[i]) <= 23.5)
  {Region[i] = "Tropical"}
  else if (abs(Mosqs$lat[i]) > 23.5 & abs(Mosqs$lat[i]) < 35.0)
  {Region[i] = "Subtropical"}
  else if (abs(Mosqs$lat[i]) >= 35.0 & abs(Mosqs$lat[i]) < 66.5)
  {Region[i] = "Temperate"}
}

Mosqs$Region = factor(Region, levels = c("Tropical", "Subtropical", "Temperate"))

#fwrite(Mosqs, "AllSpecies_TSM_Scaled.csv")

##### TSMs by grid cell (can skip to here) #####

Mosqs = fread("AllSpecies_TSM_Scaled.csv")

# Code below creates a grid across the globe, subset by locations with
# mosquito occurrence. Then calculates the mean relative TSM within each grid cell

# Convert mosquito occurrence points to sf object
MosqsSf = st_as_sf(Mosqs[,c(1,3,4,73)], coords = c("lon", "lat"), crs = 4326)

# Create grid. Units are in degrees lat/long
grid = st_make_grid(MosqsSf, cellsize = c(4, 4)) 
# Keep only grid cells with vector occurrences
index = which(lengths(st_intersects(grid, MosqsSf)) >0) 
fnet = grid[index] 
# Convert to sf object and add grid ID as a column
sffnet = st_as_sf(fnet) 
sffnet_id = mutate(sffnet, id = row_number()) 
# Combine with TSM data
griddata = st_join(MosqsSf, sffnet_id, left = TRUE) 

meandata = griddata %>%
  group_by(id) %>%
  summarize(mean_TSM = mean(ScaledTSM)) # Average TSM by grid cell
GeoData = cbind.data.frame(sffnet, meandata$mean_TSM)
colnames(GeoData) = c("geometry", "ScaledTSM")
GeoDataSF = st_as_sf(GeoData, sf_column_name = "geometry")

# create base map
world <-  rnaturalearth::ne_countries(returnclass = "sf")

ggplot() + 
  geom_sf(data = world, col = "white") + 
  geom_sf(data = GeoDataSF, aes(fill = ScaledTSM)) + theme_bw() + 
  coord_sf(ylim = c(-60,80)) + 
  theme( panel.border = element_blank()) + 
  scale_fill_viridis_c(option = "magma", direction = -1)

##### TSMs by region (tropic, subtropic, temperate) ####

# Bin occurrence records by latitude
Mosqs$BinnedLat = cut(Mosqs$lat, breaks = c(-45, -35, -23.5, -13, 0, 13, 23.5, 35, 45, 67))

# Plot relative TSM by region
boxplot(data = Mosqs, ScaledTSM ~ Region, 
        xlab = " ", ylab = "Relative Thermal Safety Margin",
        col = c("#F52549", "#FFD54D", "#99BE1B"),
        cex.lab=1.6, cex.axis=1.6, frame = F)

# Plot relative TSM by binned latitude
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
MosqsBiomes = merge(Mosqs, biomes[,c("Species.Pop","lat", "lon", "BIOME")], by = c("Species.Pop", "lat", "lon"))
MosqsBiomes$BIOME = as.factor(MosqsBiomes$BIOME)

# remove records with likely mis-classificiont (i.e., "Mangroves" (14), and "Lake" (98))
# or very rare (tundra, taiga, montane grasslands, tropical & subtropical coniferous forests)
MosqsBiomes = MosqsBiomes[MosqsBiomes$BIOME != "14",] # Mangroves
MosqsBiomes = MosqsBiomes[MosqsBiomes$BIOME != "98",] # Lake

BiomeColors = c("#113900", "#2e7921", "#326d76","#00aeac", "#cc6e1f",
                "#6c4d02","#9d8613","#cdc904","#eba34d")

# drop unused levels
MosqsBiomes$BIOME = droplevels(MosqsBiomes$BIOME)
MosqsBiomes <- MosqsBiomes[-which(is.na(MosqsBiomes$BIOME)),]
# re-order factor levels
MosqsBiomes$BIOME = factor(MosqsBiomes$BIOME, levels = c(1,2,4,5,12, 7:9, 13))

par(mar = c(4.5, 4.5, 1, 1))  
ggplot(MosqsBiomes, aes(x = BIOME, y = ScaledTSM)) + ylim(0, 1.0) +
  geom_violin(trim = FALSE, aes(fill = factor(BIOME)), draw_quantiles = 0.5) + 
  scale_fill_manual(values = BiomeColors) + 
  geom_jitter(height = 0, width = 0.1, col = alpha("black", 0.1)) +
  theme_bw() + 
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14),
        panel.border = element_blank()) + 
  xlab("Biome") + ylab("Relative Thermal Safety Margin")

# How many records in each biome?
NumPerBiome <- MosqsBiomes %>% 
  group_by(BIOME) %>% 
  summarise(n = n())

##### TSMs by Global Burden of Disease regions (not using) #####

# Global burden of disease classifications listed here: 
# https://www.healthdata.org/sites/default/files/files/Projects/GBD/GBDRegions_countries.pdf

CentralAsia = c("Armenia", "Azerbaijan", "Georgia", "Kazakhstan", "Kyrgyzstan",
                "Mongolia", "Tajikistan", "Turkmenistan", "Uzbekistan")
CentralEurope = c("Albania", "Bosnia and Herzegovina", "Croatia", "Czech Republic", 
                  "Hungary", "Macdeonia", "Montenegro", "Poland", "Romania", 
                  "Serbia", "Slovakia", "Slovenia")
EasternEurope = c("Belarus", "Estonia", "Latvia", "Lithuania", "Moldova", 
                  "Russia", "Ukraine")
Australasia = c("Australia", "New Zealand")
HighIncomeAP = c("Brunei", "Japan", "Singapore", "South Korea")
HighIncomeNA = c("Canada", "United States of America")
SouthernLatAm = c("Argentina", "Chile", "Uruguay")
WesternEurope = c("Andorra", "Austria", "Belgium", "Cyprus", "Denmark", 
                  "Finland", "France", "Germany", "Greece", "Greeland", 
                  "Iceland", "Ireland", "Isreal", "Italy", "Luxembourg",
                  "Malta", "Netherlands", "Norway", "Portugal", "Spain", 
                  "Sweden", "Switzerland", "United Kingdom")
AndeanLatAm = c("Bolivia", "Ecuador", "Peru")
Caribbean = c("Antigua and Barbua", "Bahamas", "Barbados", "Belize", "Bermuda",
              "Cuba", "Dominica", "Dominican Republic", "Grenada", "Guyana", 
              "Haiti", "Jamaica", "Puerto Rico", "Saint Lucia", 
              "Saint Vincent and the Grenadines", "Suriname", "Trinidad and Tobago")
CentralLatAm = c("Colombia", "Costa Rica", "El Salvador", "Guatemala", "Honduras",
                 "Mexico", "Nicaragua", "Panama", "Venezuela")
TropicalLatAm = c("Brazil", "Paraguay")
NorthAfricaMidEast = c("Afghanistan", "Algeria", "Bahrain", "Egypti", "Iran",
                       "Iraq", "Jordan", "Kuwait", "Lebanon", "Libya", "Morocco",
                       "Palestine", "Oman", "Qatar", "Saudi Arabia", "Sudan", 
                       "Syria", "Tunisia", "Turkey",  
                       "United Arab Emirates", "Yemen")
SouthAsia = c("Bangladesh", "Bhutan", "India", "Nepal", "Pakistan")
CentralSubSahara = c("Angola", "Central African Republic", "Congo", "Gabon",
                     "Democratic Republic of the Congo", "Equatorial Guinea")
EasternSubSahara = c("Burundi", "Comoros", "Djibouti", "Eritrea","Zambia",
                     "Ethiopia", "Kenya", "Madagascar", "Malawi", "Mozambique",
                     "Rwanda", "Somalia", "South Sudan", "Tanzania", "Uganda")
SouthernSubSahara = c("Botswana", "Lesotho", "Namibia", "South Africa",
                      "Swaziland", "Zimbabwe")
WesternSubSahara = c("Benin", "Burkina Faso", "Cameroon", "Cape Verde",
                     "Chad", "Cote d'Ivoire", "The Gambia","Ghana", "Guinea", 
                     "Guinea-Bissau", "Liberia", "Mali", "Mauritania", 
                     "Niger", "Nigeria", "Sao Tome and Principe", "Senegal", 
                     "Sierra Leon", "Togo")
EastAsia = c("China", "North Korea", "Taiwan")
SEAsia = c("Cambodia", "Indonesia", "Laos", "Malaysia", "Maldives",
           "Mauritius", "Myanmar", "Phillipines", "Seychelles", "Sri Lanka",
           "Thailand", "Timor-Leste", "Vietnam")
Oceania = c("American Samoa", "Micronesia", "Fiji", "Guam", "Kiribati", 
            "Marshall Islands", "Papua New Guinea", "Samoa", "Solomon Islands",
            "Tonga", "Vanuatu")

# Pull in dataset where each occurrence record has been assigned to country
# Note this dataset was generated in QGIS
countries = read.csv("~/Downloads/AllVectorSpecies_CountryLabel.csv", header = T)

# Assign countries to Global Burden of Disease Regions
GBD = vector(length = nrow(countries))
for (i in 1:nrow(countries)){
  if (countries$name[i] %in% CentralAsia)
  {GBD[i] = "CentralAsia"}
  else if (countries$name[i] %in% CentralEurope)
  {GBD[i] = "CentralEurope"}
  else if (countries$name[i] %in% EasternEurope)
  {GBD[i] = "EasternEurope"}
  else if (countries$name[i] %in% Australasia)
  {GBD[i] = "Australasia"}
  else if (countries$name[i] %in% HighIncomeAP)
  {GBD[i] = "HighIncomeAP"}
  else if (countries$name[i] %in% HighIncomeNA)
  {GBD[i] = "HighIncomeNA"}
  else if (countries$name[i] %in% SouthernLatAm)
  {GBD[i] = "SouthernLatAm"}
  else if (countries$name[i] %in% WesternEurope)
  {GBD[i] = "WesternEurope"}
  else if (countries$name[i] %in% AndeanLatAm)
  {GBD[i] = "AndeanLatAm"}
  else if (countries$name[i] %in% Caribbean)
  {GBD[i] = "Caribbean"}
  else if (countries$name[i] %in% CentralLatAm)
  {GBD[i] = "CentralLatAm"}
  else if (countries$name[i] %in% TropicalLatAm)
  {GBD[i] = "TropicalLatAm"}
  else if (countries$name[i] %in% NorthAfricaMidEast)
  {GBD[i] = "NAfricaMidEast"}
  else if (countries$name[i] %in% SouthAsia)
  {GBD[i] = "SouthAsia"}
  else if (countries$name[i] %in% CentralSubSahara)
  {GBD[i] = "CenSubSahara"}
  else if (countries$name[i] %in% EasternSubSahara)
  {GBD[i] = "EastSubSahara"}
  else if (countries$name[i] %in% SouthernSubSahara)
  {GBD[i] = "SouthSubSahara"}
  else if (countries$name[i] %in% WesternSubSahara)
  {GBD[i] = "WestSubSahara"}
  else if (countries$name[i] %in% EastAsia)
  {GBD[i] = "EastAsia"}
  else if (countries$name[i] %in% SEAsia)
  {GBD[i] = "SEAsia"}
  else if (countries$name[i] %in% Oceania)
  {GBD[i] = "Oceania"}
  else {GBD[i] = "NotListed"}
}

countries$GBD = GBD

# Merge with relative TSM data
MosqsGBD = merge(Mosqs, countries, by = c("latitude", "longitude"))
MosqsGBD$GBD = as.factor(MosqsGBD$GBD)

# Order by mean scaled TSM
Ordering = MosqsGBD %>% group_by(GBD) %>% summarize(meanTSM = mean(ScaledTSM))
Ordering = Ordering[order(Ordering$meanTSM),]

MosqsGBD$GBD = factor(MosqsGBD$GBD, levels = unique(Ordering$GBD))

par(mar = c(10,5,1,1))
boxplot(data = MosqsGBD, ScaledTSM ~ GBD, 
        xlab = " ", ylab = "Relative Thermal Safety Margin",
        #  col = c("#F52549", "#FFD54D", "#99BE1B"),
        cex.lab=1.3, cex.axis=1.5, frame = F, las = 2)

