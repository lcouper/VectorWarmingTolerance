# Assignment of vector occurrence records to biomes

We used QGIS to assign vector occurrence records to the WWF 14 biomes (as described in Olson et al. 2001) based on their latitude and longitude.
First, we added the "AllSpecies_LatLongs" as a vector (point) layer, and "WWF_Biomes.shp" as a vector (polygon layer). Note the WWF_Biomes shape file was obtained here: https://www.worldwildlife.org/publications/terrestrial-ecoregions-of-the-world and was not edited prior to upload in QGIS.  
Next, we usied the vector geoprocessing tool 'Intersection' to calculate the intersection of the point occurrence records to the polygon biome designations. This assigned a biome to each occurrence record, creating a new layer, exported as a .csv, and uploaded here for reference as "VectorTSM_Biome_Classifications.csv"

<p align="center">
  <img width="800"
    src="https://github.com/user-attachments/assets/3619c195-0957-49f5-8687-1f01e56cea2d">
  </p>    
<p align="center"> 
