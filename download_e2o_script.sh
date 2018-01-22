#!/bin/bash -f

echo 'Start Pre-processing step 1: download E2O'

# Setup info for the E2O ftp server
touch .e2o_netrc
echo 'machine wci.earth2observe.eu login e2o_uu password ipa9aeRae7ahC4ae' > .e2o_netrc
chmod 0600 .e2o_netrc
DIRECTORY='ftp://wci.earth2observe.eu/data/primary/public'
# Setup the chosen variables and models
INSTS=(anu polytechfr ecmwf uu) # leads to models: W3RA, ORCHIDEE, H-TESSEL, PCR-GLOBWB
CODES=(anu cnrs ecmwf univu) # Extra alias inside the project (needed for the correct url)
VARIABLES=(PotEvap Precip TVeg TotMoist)

# For the chosen models this script downloads the desired variables (monthly average), crops spatially and temporally. It sets the timestamp to the first day of the month, similar to download_modis_script.sh
# NOTE:
# Two of the combinations are non-existent on the server: anu-Precip and polytechfr-Precip, this produces warnings.
# Multiple datasets: anu-TVeg, polytechfr-TVeg and polytechfr-PotEvap, are stored on the portal unaccording to the sign-conventions! This is later corrected in the dataset_building.R script.
for i in `seq 0 1 $[${#INSTS[@]}-1]`; do # Index (0:3) to loop over INSTS and CODES at once, other option: for ((i=0; i < ${#INSTS[@]}; i++)); do
  echo 'download variables for' ${INSTS[i]}
  for VARIABLE in ${VARIABLES[@]}; do
    echo ${VARIABLE}
    curl --netrc-file .e2o_netrc ${DIRECTORY}/${INSTS[i]}/wrr2/e2o_${CODES[i]}_wrr2_glob15_mon_${VARIABLE}_1980-2014.nc -o ./temp.nc
    # CDO is used for the many operations: spatial and temporal cropping and the manipulation of the timestamp. However, for anu and uu an extra inversion of latitudes is needed to make the reading into R equal (done in dataset_building.R). The if statement below leads thus only to ecmwf and polytechfr:
    if [ ${i} -eq 1 ] || [ ${i} -eq 2 ]; then
      cdo -O setreftime,1901-01-01,00:00:00 -settime,00:00:00 -setday,1 -settunits,days -seldate,2000-02-01,2014-12-31 -sellonlatbox,140,147.75,-25,-32.25 ./temp.nc ./e2o_${INSTS[i]}_${VARIABLE}_2000-2014_p25deg.nc
      # Special treatment of polytechfr, their z-dimension is named 'time_counter' instead of 'time'
      if [ ${i} -eq 1 ]; then
        ncrename -d time_counter,time ./e2o_${INSTS[i]}_${VARIABLE}_2000-2014_p25deg.nc -O ./e2o_${INSTS[i]}_${VARIABLE}_2000-2014_p25deg.nc
      fi
    else # this is for anu and uu
      cdo -O invertlat -settime,00:00:00 -setday,1 -seldate,2000-02-01,2014-12-31 -sellonlatbox,140,147.75,-25,-32.25 ./temp.nc ./e2o_${INSTS[i]}_${VARIABLE}_2000-2014_p25deg.nc
    fi
    rm ./temp.nc
  done
done

