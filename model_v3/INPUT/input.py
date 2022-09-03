## import python packages
import xarray as xr
import netCDF4 as nc4
from pathlib import Path


####################################
## import essential AMCLIM modules
####################################
from CONFIG.config import *        ## animal data file names specified in CONFIG.config.py
from MODULES.FUNC import *

## directing input file directory
# infile_path = '/home/s1576984/scratch/working_directory/AMCLIM/input_files/'
# infile_path = '/gws/nopw/j04/macaque/JJz/jjz_virtual_env/working_directory/new_AMCLIM/AMCLIM/input_files/' 
animal_data_path = 'animal_data/'
met_data_path = 'met_data/met_data/'
crop_data_path = 'crop_data/'
soil_data_path ='soil_data/'
## Open animal data files and meteorological files
# animal_file = xr.open_dataset(infile_path+animal_data_path+animal_file_name)

## open meteorology/parameter/variable files
## meteorological data in 2018: 1) temperature (2m air; ground), 2) relative humidity, 3) 10m wind speed, 4) evaporation from soil,
##     5) soil moisture data, 6) percentage of saturation soil moisture 7) sensible heat flux (J/m^2/s)(not used) 8) rainfall (kg/m^2) 
##     9) aerogynamic and boundary layer resistance (s/m)
temp_file = xr.open_dataset(infile_path+met_data_path+'ERA5-temp2m-2018d-05.nc')
groundtemp_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1temp-2018d-05.nc')
groundtemp_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2temp-2018d-05.nc')
rhum_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-rhum2m-2018d-05.nc')
wind_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-wind10m-2018d-05.nc')
evap_file = xr.open_dataset(infile_path+met_data_path+'ERA5-evapfromsoil-2018d-05.nc')
# soilmoist_file = xr.open_dataset(infile_path+met_data_path+'SOILMOISTURE-L3S-SSMV-COMBINED-DAILY-2018-360x720.nc')
soilmoist_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1moist-2018d-05.nc')
soilmoist_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2moist-2018d-05.nc')
# soilsm_file = xr.open_dataset(infile_path+met_data_path+'SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2018-360x720.nc')
rain_file = xr.open_dataset(infile_path+met_data_path+'ERA5-totalcolumn_rainwater-2018d-05.nc')
ratm_file = xr.open_dataset(infile_path+met_data_path+'I2000clm50_RAM_output.clm2.h1.0018-01-01-00000-05.nc')
runoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-srfrunoff-2018d-05.nc')
subrunoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-subrunoff-2018d-05.nc')

#################################
## crop and fertilizer data
#################################
fertfilename = 'IFA_fert_info_21st_century.nc'
cropfileformat = '_cropping_info.nc'
crop_filledcalendar = 'ALL_CROPS_netCDF_0.5deg_filled_dir/'
crop_unfilledcalendar = 'ALL_CROPS_netCDF_0.5deg_unfilled_dir/'
crop_filledcalendarformat = '.crop.calendar.fill.nc'
crop_unfilledcalendarformat = '.crop.calendar.nc'
manure_appcalendar ='manure_app_calendar.nc'

##################################
## soil data
##################################
soilpHfile = 'topsoil_pH_H2O.nc'
soilCECfile = 'topsoil_CEC_clay.nc'
soilclayfile = 'topsoil_clay.nc'
soilsandfile = 'topsoil_sand.nc'
soilsiltfile = 'topsoil_silt.nc'
soilorgCfile = 'topsoil_OC.nc'
soilbdfile = 'topsoil_BD.nc'
