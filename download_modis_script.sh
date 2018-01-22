#!/bin/bash -f

# This script needs to be provided with a file of MODIS url's, a list that can be constructed in the earth-data poral of NASA. Currently the MOD13C2 product is chosen.
# Check https://lpdaac.usgs.gov/sites/default/files/public/product_documentation/mod13_user_guide.pdf for the MOD13C2 meta-data. 
# Note: 2008-06-01 and 2008-08-01 were missed during one run, url was found but download did not happen. curl error?

FILE="$1"

echo 'Start Pre-processing step 1: download MODIS'

# Preparation for logging in and downloading the data, create files inside folder of script.
touch .urs_netrc
echo 'machine urs.earthdata.nasa.gov login chiemvs password Earthdata123' > .urs_netrc
chmod 0600 .urs_netrc

touch .urs_cookies

# These control the timestamps, begin and end should match the provided urls, otherwise a warning will be given.
YEARB="2000"
MONTHB="02"
YEARE="2017"
MONTHE="09"
DAY="01"

# Loop through time and for each year/month combination: find the url, download, extract only EVI and the quality variable, convert from hdf to netCDF, manipulate timestamp, and save. Temporal stacking will happen later, after the loop.
YEAR=${YEARB}

while [ ${YEAR} -le ${YEARE} ]; do
  # subroutine to assign correct months to each year, based on whether it is the first (feb to dec), an intermediate (jan to dec) or the last year (jan to sep)
  if [ ${YEAR} == ${YEARB} ]; then
    MONTHS=($(seq $MONTHB 1 12))
  else 
    if [ ${YEAR} == ${YEARE} ]; then
      MONTHS=($(seq 1 1 $MONTHE))
    else
      MONTHS=($(seq 1 1 12))
    fi
  fi
  echo "${YEAR} 'with months' ${MONTHS[*]}"
  
  # Now we want to execute something for each month in that year.
  for MONTH in ${MONTHS[@]}; do
    if [ ${MONTH} -lt 10 ]; then MONTH=0${MONTH}; fi # Small correction for extra zero in front
    # Loop through the provided file with urls and find the corresponding line.
    while read LINE; do
      if echo ${LINE} | grep -q "${YEAR}.${MONTH}.${DAY}"; then
        URL=${LINE}
      fi
    done <${FILE}
    
    # If the url is found:
    if [ ! -z ${URL} ]; then
      # Download, can take a while (global files, many variables)
      curl -L -c .urs_cookies -b .usr_cookies --netrc-file .urs_netrc ${URL} -o MOD13C2_${YEAR}${MONTH}.hdf
      # Extract two important variables: EVI: range [-2000,10000], fill -3000; Quality bitfield: [0,65534], fill 65535
      gdal_translate -ot Int16 -of netCDF HDF4_EOS:EOS_GRID:"MOD13C2_${YEAR}${MONTH}.hdf":'MOD_Grid_monthly_CMG_VI:CMG 0.05 Deg Monthly EVI' MOD13C2_EVI_${YEAR}${MONTH}.nc
      gdal_translate -ot Int32 -of netCDF HDF4_EOS:EOS_GRID:"MOD13C2_${YEAR}${MONTH}.hdf":'MOD_Grid_monthly_CMG_VI:CMG 0.05 Deg Monthly VI Quality' MOD13C2_QUAL_${YEAR}${MONTH}.nc
      cdo -O sellonlatbox,140,147.75,-25,-32.25 -setname,evi MOD13C2_EVI_${YEAR}${MONTH}.nc MOD13C2_evi_${YEAR}${MONTH}.nc
      cdo -O sellonlatbox,140,147.75,-25,-32.25 -setname,quality MOD13C2_QUAL_${YEAR}${MONTH}.nc MOD13C2_qual_${YEAR}${MONTH}.nc
      # Timestamp is set to first of that month, for later temporal merging. Global files are removed
      cdo -O merge -setdate,${YEAR}-${MONTH}-${DAY},00:00:00 MOD13C2_evi_${YEAR}${MONTH}.nc -setdate,${YEAR}-${MONTH}-${DAY},00:00:00 MOD13C2_qual_${YEAR}${MONTH}.nc MOD_${YEAR}${MONTH}.nc
      find . -name "MOD13C2_*" | xargs rm
    else
      echo "'url for' ${YEAR} ${MONTH} 'could not be found'"
      exit 1
    fi
    
    URL= # Reset the parameter to empty, so a warning can be given if nothing is found in the next step.
  done
  YEAR=$[${YEAR}+1]
done

# Merge all the files to a temporal stack, which can be read into R.
# Currently the variable 'evi' remains unscaled. The needed scaling is 0.0001 so the range becomes its original [-0.2, 1]
# Integer variable 'quality' needs to be converted to 16 bit binary sequence to read, check metadata for the meaning of each bit.
# There are 145 rows, 155 columns and 210 timesteps, for Analysis 1
cdo -O mergetime MOD_*.nc stack_p05deg.nc

# This file is also averaged to a 0.25degree version, for comparison with model output in Analysis 2. 29 rows, 31 columns.
cdo -O gridboxmean,5,5 stack_p05deg.nc stack_p25deg.nc

# Clean-up.
find . -name "MOD_*" | xargs rm
