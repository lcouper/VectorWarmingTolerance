# Data sources used in analysis of vector thermal safety

### Vector occurrence data

- GBIF (all species
- Kraemer et al 2015 (*Aedes aegypti* and *Ae. albopictus*)
- Korba et al. 2015, Moradi-Asl et al. 2019, Al Ahmad et al. 2011, Biteye et al. 2019 (*Ae. vexans*)
- Carlson et al. 2023 (*Anopheles gambiae*)
- Sinka et al. 2020, Ochomo et al. 2023 (*An. stephensi*)


### Vector thermal limits

- Villena et al. 2022 (*An. gambiae* and *An. stephensi*)
- Mordecai et al. 2019 (*Aedes aegypti, albopictus, vexans, camptorhynchus, triseriatus; Culex annulirostris, pipiens, quinquefasciatus, tarsalis, theileri*)


## Climate data

- ERA5 climate re-analysis data, accessed here: https://cds.climate.copernicus.eu
- Resolution: 1 hour,  .25°
- Coverage: 1950 to 5 days from present    
- Parameter values used and shared across all species: 
```
# For microclimate model:
# Usrhyt = 1 # to estimate microcliamte conditions at 1m (roughly average height of activity)
# minshade = 0 # 0% canopy cover (full sun)
# maxshade = 100 # 100% canopy cover (full shade)

# For ectotherm model:
# Ww_g = 0.003 # for now, using average of 3 mg as mosquito weight for all species
# shape = 2 keeping at 2 (ellipsis) for all species
# alpha_min : 0.85 (default, solar absorptivity)
# alpha_max : 0.85 (default, solar absorptivity)
# shade_seek = 1 # Running mdoels under shade-seeking and no-shade seeking for each species
# burrow = 0 (not allowing burrowing for any species)
# climb = 0 (not allowing climbing for any species)
# minshade = 0 (minimum available shade. 0 = full sun)
# maxshade = 100 (full shade)
T_pref, CT_max, and CT_min, diurnal, nocturnal, crepuscular are set to be species specific
```
