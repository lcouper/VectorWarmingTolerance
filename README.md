# VectorWarmingTolerance
Repository for code and data sources used to analyze vector warming tolerance

Project outline doc:
https://docs.google.com/document/d/1tI3qwjM38633R3KNWEdQIeQdMpTHPs_BhFMrloRZyVI/edit 

Questions:
- in ectotherm model, what to set for 'T-pref' - the optimal body temp? Topt?
- what life history trait(s) to use? Only have adult lifespan for a couple species. Could use a variety of adult traits? (e.g. fecundity, biting rate, vector competence)
- to look up: are any mosquito species nocturnal? currently specifiying all as diurnal and crepuscular, but not nocturnal

Parameter values used and shared across all species
```
# Ww_g = 0.3 # for now, using average of 0.3 g as mosquito weight for all species
# shape = 2 keeping at 2 (ellipsis) for all species
# alpha_min : 0.85 (default, solar absorptivity)
# alpha_max : 0.85 (default, solar absorptivity)
# diurn = 1 using for all species
# nocturn = 0 using for all species
# crepus = 1 using for all species
# shade_seek = 1 # allowing shade-seeking for all species
# burrow = 0 (not allowing burrowing for any species)
# climb = 0 (not allowing climbing for any species)
# minshade = 0 (minimum available shade. 0 = full sun)
# maxshade = 100 (full shade)
T_pref, CT_max, and CT_min are set to be species specific
```

Source of Aedes aegypti and Aedes albopictus occurrences:
Kraemer et al. 2017: https://datadryad.org/stash/dataset/doi:10.5061/dryad.47v3c
