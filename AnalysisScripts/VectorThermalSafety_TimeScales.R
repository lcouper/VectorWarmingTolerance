#### Vector Thermal Safety Time Scales of Thermal Danger #####

# This script is used to model and plot results relating to 
# the time scales of thermal danger (# consecutive hours or day)

#### 0. Load libraries and data frames #####

library(ggplot2)

setwd("~/Documents/Current Projects/WarmingTolerance/DataFiles/Vector_TSM")

# Pull in species data files (with seasonality)
AeAegypti = read.csv("AeAegypti_TSM_DroughtMask_Combined_WithElevation.csv")
AeAlbo = read.csv("AeAlbo_TSM_DroughtMask_Combined_WithElevation.csv")
AnGambiae = read.csv("AnGambiae_TSM_DroughtMask_Combined_WithElevation.csv")
AnSteph = read.csv("AnSteph_TSM_DroughtMask_Combined_WithElevation.csv")
CxQuinque = read.csv("CxQuinque_TSM_DroughtMask_Combined_WithElevation.csv")  
CxPip = read.csv("CxPipiens_TSM_DroughtMask_Combined_WithElevation.csv")
CxTar = read.csv("CxTarsalis_TSM_DroughtMask_Combined_WithElevation.csv")
CxAnnul = read.csv("CxAnnul_TSM_DroughtMask_Combined_WithElevation.csv")
 
Mosqs = fread("AllSpecies_TSM_Scaled.csv") # high elevation observations already removed
SpeciesList  = c("Aedes_aegypti", "Aedes_albopictus", 
                 "Anopheles_gambiae", "Anopheles_stephensi", 
                 "Culex_annulirostris",  "Culex_pipiens", "Culex_quinquefasciatus", 
                 "Culex_tarsalis")
Species.Colors = c("Aedes_aegypti" = "#ab041b", "Aedes_albopictus" = "#cc5801", 
                   "Anopheles_gambiae" = "#f4bb00",  "Anopheles_stephensi" = "#7cae00", 
                   "Culex_annulirostris" = "#309143", "Culex_pipiens" = "#74add1", 
                   "Culex_quinquefasciatus" = "#313695","Culex_tarsalis" = "#dda0dd")

Mosqs$Species.Pop = factor(Mosqs$Species.Pop, levels = SpeciesList)
Mosqs$Region = factor(Mosqs$Region, levels = c("Tropical", "Subtropical", "Temperate"))


#### Consecutive hours in thermal danger ####

Mosqs$BinnedStreakHours = Mosqs$streak_point2
Mosqs$BinnedStreakHours[Mosqs$streak_point2 > 20] <- 21 # to group outliers

###### Color by species ######

# Convert NAs to 0s
Mosqs$streak_point2[is.na(Mosqs$streak_point2)] <-0
Mosqs$streak_NoB_point2[is.na(Mosqs$streak_NoB_point2)] <-0

ggplot(data = Mosqs, aes(x = BinnedStreakHours, fill = factor(Species.Pop, levels = SpeciesList))) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 1) + theme_minimal() + 
  scale_fill_manual(values = Species.Colors,
                    name = "Species", 
                    labels = SpeciesList) + 
  xlim(c(0,22)) + # note setting this x-axis limit removes 4 outliers
  labs(x = "Consecutive hours", y = " ") + 
  ggtitle(" ") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18, face="italic"),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

# Boxplot
ggplot(data = Mosqs, aes(x = Species.Pop, y =streak_point2, fill = Species.Pop)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,25)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive hours") + 
  scale_fill_manual(values=Species.Colors) + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


###### Color by region ######

ggplot(data = Mosqs, aes(x = BinnedStreakHours, fill = Region)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 1) + theme_minimal() + 
  scale_fill_manual(values = c("#F52549", "#FFD54D", "#99BE1B"),
                    name = " ") + 
  labs(x = "Consecutive hours", y = " ") + 
 # xlim(c(0,25)) + # note setting this x-axis limit removes 4 outliers
  ggtitle(" ") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

aggregate(Mosqs$streak_point2 ~ Mosqs$Region, FUN = "mean")
aggregate(Mosqs$streak_point2 ~ Mosqs$Region, FUN = "median")

# Boxplot
ggplot(data = Mosqs, aes(x = Region, y =streak_point2, fill = Region)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,25)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive hours") + 
  scale_fill_manual(values=c("#F52549", "#FFD54D", "#99BE1B")) + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

###### Color by biome ####

biomes = read.csv("~/Documents/Current Projects/WarmingTolerance/DataFiles/Biogeography/VectorTSM_Biome_Classifications.csv", header = T)
# merge with TSM data in Mosqs dataframe
Mosqs$streak_point2[is.na(Mosqs$streak_point2)] <-0 # Make sure this is true before merging
MosqsBiomes = merge(Mosqs, biomes[,c("Species.Pop","lat", "lon", "BIOME")], by = c("Species.Pop", "lat", "lon"))
MosqsBiomes$BIOME = as.factor(MosqsBiomes$BIOME)

# remove records with likely mis-classificiont (i.e., "Mangroves" (14), and "Lake" (98))
# or very rare (tundra, taiga, montane grasslands, tropical & subtropical coniferous forests)
MosqsBiomes = MosqsBiomes[MosqsBiomes$BIOME != "14",] # Mangroves
MosqsBiomes = MosqsBiomes[MosqsBiomes$BIOME != "98",] # Lake

BiomeColors = c("#113900", "#2e7921", "#326d76","#00aeac", "#cc6e1f",
                "#6c4d02","#9d8613","#cdc904","#eba34d")

# drop unused levels & re-order
#MosqsBiomes$BIOME = droplevels(MosqsBiomes$BIOME)
#MosqsBiomes <- MosqsBiomes[-which(is.na(MosqsBiomes$BIOME)),]
MosqsBiomes$BIOME = factor(MosqsBiomes$BIOME, levels = c(1,2,4,5,12, 7:9, 13))

aggregate(MosqsBiomes$streak_point2 ~ MosqsBiomes$BIOME, FUN = "mean")
aggregate(MosqsBiomes$streak_point2 ~ MosqsBiomes$BIOME, FUN = "median")

# Boxplot
ggplot(data = MosqsBiomes, aes(x = BIOME, y =streak_point2, fill = BIOME)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,25)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive hours") + xlab("Biome") +
  scale_fill_manual(values=BiomeColors) + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

#### Consecutive days in thermal danger ####

Mosqs$BinnedStreakDays = Mosqs$streak_days_point2
Mosqs$BinnedStreakDays[Mosqs$streak_days_point2 > 40] <- 41 # to group outliers

###### Color by species #####

# Convert NAs to 0s
Mosqs$streak_days_point2[is.na(Mosqs$streak_days_point2)] <-0

ggplot(data = Mosqs, aes(x = BinnedStreakDays, fill = Species.Pop)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 1) + theme_minimal() +  
  scale_fill_manual(values = Species.Colors, name = "Species", 
                    labels = SpeciesList) + 
  labs(x = "Consecutive days", y = " ") + 
  ggtitle(" NOTE: 41 here is for any value >40 ! ") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18, face = "italic"),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


# boxplot
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "mean")
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "median")

# Boxplot
ggplot(data = Mosqs, aes(x = Species.Pop, y =streak_days_point2, fill = Species.Pop)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,100)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive hours") + 
  scale_fill_manual(values=Species.Colors) + 
  theme(axis.text.y=element_text(size=18),
        axis.text.x = element_text(angle = 60, vjust = 1, hjust=1, size = 12),
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


###### Color by region ######

ggplot(data = Mosqs, aes(x = BinnedStreakDays, fill = Region)) + 
  geom_histogram(aes(y = after_stat(count / sum(count))), binwidth = 1) + theme_minimal() + 
  labs(x = "Consecutive days", y = " ") + 
  scale_fill_manual(values = c("#F52549", "#FFD54D", "#99BE1B"),
                    name = " ") + 
 # ggtitle(" NOTE: 41 here is for any value >40 ! ") + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


aggregate(Mosqs$streak_days_point2 ~ Mosqs$Region, FUN = "mean")
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Region, FUN = "sd")


# Convert NAs to 0s
Mosqs$streak_days_point2[is.na(Mosqs$streak_days_point2)] <-0
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Region, FUN = "mean")
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Region, FUN = "median")

# Boxplot
ggplot(data = Mosqs, aes(x = Region, y =streak_days_point2, fill = Region)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,100)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive days") + 
  scale_fill_manual(values=c("#F52549", "#FFD54D", "#99BE1B")) + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))


###### Color by biome ######

# note uses MosqsBiomes created above
aggregate(MosqsBiomes$streak_days_point2 ~ MosqsBiomes$BIOME, FUN = "mean")
aggregate(MosqsBiomes$streak_days_point2 ~ MosqsBiomes$BIOME, FUN = "median")

# Boxplot
ggplot(data = MosqsBiomes, aes(x = BIOME, y =streak_days_point2, fill = BIOME)) + 
  geom_boxplot(outlier.size=0.4) + ylim(c(0,100)) + theme_minimal() +
  geom_jitter(shape=16, position=position_jitter(0.2), cex = 0.4) + 
  ylab("Consecutive days") + xlab("Biome") +
  scale_fill_manual(values=BiomeColors) + 
  theme(axis.text=element_text(size=18), 
        axis.title = element_text(size = 18),
        legend.text=element_text(size=18),
        legend.title = element_text(size=16),
        plot.title = element_text(hjust = 0.5, size = 20, face="italic"),
        panel.border = element_rect(colour = "black", fill = NA), 
        plot.margin = margin(20, 10, 0, 0))

#### For Supplemental Tables ######

Mosqs$streak_days_nob_point2[is.na(Mosqs$streak_days_nob_point2)] <-0 # convert NAs to 0s

Mosqs2 <- fread("WithoutDroughtMask/AllSpecies_TSM_Scaled_NoDroughtMask.csv")
Mosqs2$streak_days_point2[is.na(Mosqs2$streak_days_point2)] <-0 # convert NAs to 0s
Mosqs2$streak_point2[is.na(Mosqs2$streak_point2)] <-0 # convert NAs to 0s


# Average length of longest streak (hours) in thermal danger: main model specification
aggregate(Mosqs$streak_point2 ~ Mosqs$Species.Pop, FUN = "mean")
aggregate(Mosqs$streak_point2 ~ Mosqs$Species.Pop, FUN = "sd")
aggregate(Mosqs$streak_point2 ~ Mosqs$Species.Pop, FUN = "max")
aggregate(Mosqs$streak_point2 ~ Mosqs$Species.Pop, FUN = "median")
AcrossSp <- aggregate(Mosqs$streak_point2 ~ Mosqs$Species.Pop, FUN = "mean")[,2]
mean(AcrossSp)
sd(AcrossSp)

# Average length of longest streak (hours) in thermal danger: no behavior
aggregate(Mosqs$streak_NoB_point2 ~ Mosqs$Species.Pop, FUN = "mean")
aggregate(Mosqs$streak_NoB_point2 ~ Mosqs$Species.Pop, FUN = "sd")
aggregate(Mosqs$streak_NoB_point2 ~ Mosqs$Species.Pop, FUN = "max")
AcrossSp <- aggregate(Mosqs$streak_NoB_point2 ~ Mosqs$Species.Pop, FUN = "mean")[,2]

# Average length of longest streak (hours) in thermal danger: no drought mask #####
aggregate(Mosqs2$streak_point2 ~ Mosqs2$Species.Pop, FUN = "mean")
aggregate(Mosqs2$streak_point2 ~ Mosqs2$Species.Pop, FUN = "sd")
aggregate(Mosqs2$streak_point2 ~ Mosqs2$Species.Pop, FUN = "max")
AcrossSp <- aggregate(Mosqs2$streak_point2 ~ Mosqs2$Species.Pop, FUN = "mean")[,2]
mean(AcrossSp)
sd(AcrossSp)

# Average length of longest streak (days) in thermal danger: main model specification
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "mean")
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "sd")
aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "max")
AcrossSp <- aggregate(Mosqs$streak_days_point2 ~ Mosqs$Species.Pop, FUN = "mean")[,2]
mean(AcrossSp)
sd(AcrossSp)

# Average length of longest streak (days) in thermal danger: no behavior
aggregate(Mosqs$streak_days_nob_point2 ~ Mosqs$Species.Pop, FUN = "mean")
aggregate(Mosqs$streak_days_nob_point2 ~ Mosqs$Species.Pop, FUN = "sd")
aggregate(Mosqs$streak_days_nob_point2 ~ Mosqs$Species.Pop, FUN = "max")
AcrossSp <- aggregate(Mosqs$streak_days_nob_point2 ~ Mosqs$Species.Pop, FUN = "mean")[,2]
mean(AcrossSp)
sd(AcrossSp)

# Average length of longest streak (days) in thermal danger: no drought mask
aggregate(Mosqs2$streak_days_point2 ~ Mosqs2$Species.Pop, FUN = "mean")
aggregate(Mosqs2$streak_days_point2 ~ Mosqs2$Species.Pop, FUN = "sd")
aggregate(Mosqs2$streak_days_point2 ~ Mosqs2$Species.Pop, FUN = "max")
AcrossSp <- aggregate(Mosqs2$streak_days_point2 ~ Mosqs2$Species.Pop, FUN = "mean")[,2]
mean(AcrossSp)
sd(AcrossSp)
