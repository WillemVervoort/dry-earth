#!/bin/bash -f

# This is the meta script for Chiem's Australia research.
# Place this in a desired folder with these other scripts:
# download_modis_script.sh, download_e2o_script.sh, download_dem_script.sh, adapt_presup_script.sh, dataset_building.R, analysis1.R, analysis2.R
# this document:
# modis_urls.txt 
# and these presupplied data: 
# All_gw_data.rds Precip_2000-2014_p05deg.nc PotEvap_2000-2014_p05deg.nc UtrechtSiloData_ncdf4.nc
# Then the rest of the downloading, processing and analyses will unfold there.

# Create a dedicated DATA folder, as working directory for the pre-processing, and move presupplied data there
echo 'Create DATA directory in' ${PWD}
mkdir -p ./DATA/
mv All_gw_data.rds Precip_2000-2014_p05deg.nc PotEvap_2000-2014_p05deg.nc UtrechtSiloData_ncdf4.nc ./DATA/

# Data-Preprocessing step 1
# Download the data, crop spatial extent, average if needed, clip to temporal extent and create netCDF stacks.
cd ./DATA/
echo 'Start Pre-processing step 1'
# MODIS MOD13C2 contains monthly values at 0.05 deg, they need also to be averaged to the models 0.25 degree format.
# This script downloads the products one by one, does conversion from HDF to netcdf and stacks. The global files are large, so it can take a while:
bash ../download_modis_script.sh ../modis_urls.txt
# The eartH2Observe portal is accessed to download: transpiration, potential-evapotranspiration and Total Soil Moisture of the four models. Precipitation (which is the MSWEP product) is only stored for PCR-GLOBWB and W3RA. Therefore, in step 2, the presupplied highres set (which is the MSWEP product) will be transformed to the model's 0.25 deg fields as a substitute.
bash ../download_e2o_script.sh
# Creates DEM at highest possible resolution and at the target 0.05deg resolution of Analysis 1.
bash ../download_dem_script.sh

# Data-Preprocessing step 2
# Adaptations to presupplied data, when neccesary
echo 'Start Pre-processing step 2'
# This script transforms the highres Precipitation set to 0.25 degree as the substitute for PCR-GLOBWB and W3RA precipitation. It also does some corrections on the netCDF with SILO observations.
bash ../adapt_presup_script.sh
# Nothing is done to the groundwater data, which already are presupplied in merged and filtered form because the original datasets are large and not easily downloaded. For transparency however, and to complement the report, the dedicated script for this filtering can be found in ./other_material/prepare_gw_data.R


# Data-Preprocessing step 3
echo 'Start Pre-processing step 3'
# All the netcdfs are converted to tidy datasets to be used in R, a seperate one for each Analysis. It is also here that aridity variables are created from presuplied highres values or from the e2o model output.
# Dataset 1 contains the highres aridity index of presupplied data (Prec/PotEvap), EVI and EVI-quality. All at 0.05 degree, which gives 22475 cells that are indexed.
# Dataset 2 contains the modeled variables, including an aridity variable (Prec - PotEvap) as each model estimates it's own PotEvap. Some extra corrections take place because the sign-conventions are not strictly followed on the dataportal. The averaged EVI is also present, just like the SILO observed Precip-PotEvap. All at 0.25 degree, which gives 899 cells that are indexed.
# Dataset 3 is a tidy DEM dataset at highest resolution, for deriving Elevation at each groundwater observation point and also has an aggregated version of p05deg for possible kriging to the effective EVI scale.
Rscript ../dataset_building.R

# Data-Preprocessing step 4
echo 'Start Pre-processing step 4'
# Here a table is written that at multiple resolutions links the indexed cellnr's (899 at 0.25 degree and 22475 at 0.05 degree) to certain regions of interest. It gives an option to split the later analyses spatially and single out the known hydrological mechanisms of these areas.
Rscript ../spatial_classification.R

# Create a dedicated RESULTS folder, as working directory for the analyses and in which figures and tables can be written. 
echo 'Create RESULTS directory in' ${PWD}
mkdir -p ../RESULTS/
cd ../RESULTS/

# Analysis 1
echo 'Start Analysis 1'
# At this point the temporal dimension of the data gets modified. The full EVI time series of the tify datasets get aggregated to summary metrics for each cell. The groundwater values are added, as these are also temporal summaries. There is a quality control.
Rscript ../analysis1.R

 
# Analysis 2
# Comparison 2 EVI with E2O models at native model resolution (0.25 degree)
Rscript ../analysis2.R

# Clean-up, 


