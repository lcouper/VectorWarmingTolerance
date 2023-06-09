# Vector Warming Tolerance
Repository for code and data sources used to analyze vector warming tolerance

Manuscript doc:
https://docs.google.com/document/d/1xCS3hHwLBeA_4W-JbToOJdZhmBASHR2Y1lQmUz2uQUo/edit (this is the most up-to-date version)


## Data Sources 

#### Vector occurrence data
(in addition to GBIF)

- Aedes aegypti and Aedes albopictus:
Kraemer et al. 2017: https://datadryad.org/stash/dataset/doi:10.5061/dryad.47v3c

- Aedes Vexans: 
compiled here: https://onlinelibrary.wiley.com/doi/full/10.1111/tbed.14404?saml_referrer
Algeria: https://orbi.uliege.be/bitstream/2268/196649/1/AEB%20288-294.pdf
Iran: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6885136/
Saudi Arabia: https://academic.oup.com/jme/article/48/4/717/898766?login=false
Senegeal: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0215194

Anopheles stephensi:
Kenya: https://www.researchsquare.com/article/rs-2498485/v2

- Culex species (not used, but available):
Cleaned and filtered occurrence records (from GBIF) and other sources: 
https://academic.oup.com/jme/advance-article/doi/10.1093/jme/tjad027/7123791

#### Vector thermal limits

Anopheles gambiae and Anopheles stephensi: 
Villena et al. 2022: https://esajournals.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1002%2Fecy.3685&file=ecy3685-sup-0001-AppendixS1.pdf

Aedes aegypti, albopictus, vexans, camptorhynchus, triseriatus; 
Culex annulirostris, pipiens, quinquefasciatus, tarsalis, theileri; 
Mordecai et al. 2019: https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fele.13335&file=ele13335-sup-0001-SupInfo.pdf

#### Climate data 
ERA5 climate re-analysis data, accessed here: https://cds.climate.copernicus.eu
Available at a temporal resolution of 1 hr, a spatial resolution of 0.25° and span from 1950 to 5 days from present
Info about using era5 in NicheMapR microclimate models: https://rdrr.io/github/mrke/NicheMapR/man/micro_era5.html 
and https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13877
Downloaded data in 10 degree latitudinal bands at a time:
Submit request for 'ERA hourly data on single levels from 1940 to present'. Then select the following variables needed for the microclimate model: 
- 10m u-component of wind
- 10m v-component of wind
- land-sea mask
- total precipitation
- cloud cover
- 2m temperature
- 2m dewpoint temperature
- surface pressure
- mean surface downward longwave radiation flux
- mean surface net longwave radiation flux
- Total sky direct solar radiation at surface
- Surface solar radiation downwards, clear sky

Then select all days and months within 2017. Input the desired lat/long bounds and request data.
Given the large file size, requests can take ~4 hours to fulfill (and ~30 minutes to download)
Files need to be named "era5_2017.nc" on local computer to be run properly


Parameter values used and shared across all species
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




