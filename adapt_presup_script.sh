#!/bin/bash -f

# The highres MSWEP precipitation set (presupplied and constructed at Deltares) recieves an additional treatment. It is averaged to p25deg for use in the Precipitation - PotEvap variable, as Anu and polytechfr do not supply Precip on the E2O portal.
cdo -O gridboxmean,5,5 Precip_2000-2014_p05deg.nc Precip_2000-2014_p25deg.nc
cp Precip_2000-2014_p25deg.nc e2o_anu_Precip_2000-2014_p25deg.nc
cp Precip_2000-2014_p25deg.nc e2o_polytechfr_Precip_2000-2014_p25deg.nc

# SILO data drills, prepared by Willem
# Correct a typo in the time units for correct reading in cdo. 
ncatted -O -a units,time,m,c,"days since 1999-12-31" UtrechtSiloData_ncdf4.nc
# Correct the units attributes to kg m-2 s-1, temporal averaging and set timestamp to first of the month.
cdo -O setday,1 -monmean  UtrechtSiloData_ncdf4.nc silo_temp_2000-2015.nc
ncatted -O -a units,rain,m,c,'kg m-2 s-1' silo_temp_2000-2015.nc
ncatted -O -a units,evap,m,c,'kg m-2 s-1' silo_temp_2000-2015.nc
ncatted -O -a units,MaET,m,c,'kg m-2 s-1' silo_temp_2000-2015.nc

# The latitude midpoints of these cells are not correct. They should be shifted by 0.1 degree
cdo griddes silo_temp_2000-2015.nc > temp_grid
sed -i 's/-25.025/-25.125/' temp_grid
cdo setgrid,temp_grid silo_temp_2000-2015.nc silo_shifted_2000-2015.nc

# Select the variables, do the multiplication for unit conversion according E2O sign conventions (positive downward)
cdo mulc,0.000011574 -selvar,rain silo_shifted_2000-2015.nc silo_Precip_2000-2015.nc
cdo mulc,-0.000011574 -selvar,evap silo_shifted_2000-2015.nc silo_PotEvap_2000-2015.nc
cdo mulc,-0.000011574 -selvar,MaET silo_shifted_2000-2015.nc silo_MaET_2000-2015.nc # This is Morton's actual ET estimate


rm silo_temp_2000-2015.nc silo_shifted_2000-2015.nc temp_grid
