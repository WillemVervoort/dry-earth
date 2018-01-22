#!/bin/bash -f

# This script downloads and creates a DEM with the right extent at the highest possible and at 0.05 degree resolution. It is based on Jarvis, A., Reuter, H. I., Nelson, A., & Guevara, E. (2008). Hole-filled SRTM for the globe Version 4. available from the CGIAR-CSI SRTM 90m Database (http://srtm. csi. cgiar. org).

echo 'Start Pre-processing step 1: download DEM'

# Setup info for the CGIAR STRM server
DIRECTORY='ftp://srtm.csi.cgiar.org/SRTM_V41/SRTM_Data_GeoTiff/'
TILES=(srtm_65_19 srtm_66_19 srtm_65_18 srtm_66_18)

for TILE in ${TILES[@]}; do
  echo ${DIRECTORY}${TILE}
  curl ${DIRECTORY}${TILE}.zip -o ./${TILE}.zip
  unzip -o ${TILE}.zip ${TILE}.tif
  rm ${TILE}.zip
done

# Create the dem at the highest 3 arc-second, and at the 0.05 degree resolution (155x145 cells), write output as geotiff and also as netcdf.
# NoData Value=-32767
# netCDF variable name due to gdal: 'Band1'
gdal_merge.py srtm* -ul_lr 140 -25 147.75 -32.25 -o merged_dem_highest_res.tif
gdal_merge.py merged_dem_highest_res.tif -of netCDF -o merged_dem_highest_res.nc # Not really a merge, but more a transformation of the .tif. In other cases an error occurs (halve the file is missing)
gdal_merge.py srtm* -ps 0.05 0.05 -ul_lr 140 -25 147.75 -32.25 -o merged_dem_p05deg.tif
gdal_merge.py srtm* -ps 0.05 0.05 -ul_lr 140 -25 147.75 -32.25 -of netCDF -o merged_dem_p05deg.nc

rm srtm*

