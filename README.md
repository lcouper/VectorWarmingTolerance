# VectorWarmingTolerance
Repository for code and data sources used to analyze vector warming tolerance

Project outline doc:
https://docs.google.com/document/d/1tI3qwjM38633R3KNWEdQIeQdMpTHPs_BhFMrloRZyVI/edit 

Manuscript doc:
https://docs.google.com/document/d/1xCS3hHwLBeA_4W-JbToOJdZhmBASHR2Y1lQmUz2uQUo/edit

Supplemental Files:
'Supp Table 1': DataFiles > VectorWarmingTolerance_SpeciesThermalLimits.csv
'Supp Table 2': DataFiles > VectorOccurrence > VectorWarmingTolerance_OccurrenceDataSources.csv
Models run so far:
- Aedes aegypti, -5 to 5 latitude.   
Saved as WarmingTolerance > DataFiles > AeAegyptiTSM_neg5to5.csv

TO DO: 
-find Anopheles occurrence records
-adjust models for Anopheles so that nocturnal activity allowed!

## Sources 

### Vector occurrence data
Aedes aegypti and Aedes albopictus:
Kraemer et al. 2017: https://datadryad.org/stash/dataset/doi:10.5061/dryad.47v3c
Anopheles global: 
Sinka et al. 2012: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3349467/
Anopheles (Africa only): 
Kyalo et al. 2017:  https://dataverse.harvard.edu/dataset.xhtml?persistentId=doi:10.7910/DVN/NQ6CUN


### Vector thermal limits
Anopheles gambiae and Anopheles stephensi: 
Villena et al. 2022: https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3685&file=ecy3685-sup-0001-AppendixS1.pdf

Aedes aegypti, albopictus, vexans, camptorhynchus, notoscriptus, triseriatus; 
Culex annulirostris, pipiens, quinquefasciatus, tarsalis, theileri; 
Culiseta melanura: 
Mordecai et al. 2019: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.13335&file=ele13335-sup-0001-SupInfo.pdf

### Climate data 
https://www.climdex.org/access/
using era5 data (available at a temporal resolution of 1 hr, a spatial resolution of 0.25° and span from 1950 to 5 days from present).  
info about using micro_era5: https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html


Parameter values used and shared across all species
```
# Ww_g = 0.003 # for now, using average of 0.3 g as mosquito weight for all species
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




