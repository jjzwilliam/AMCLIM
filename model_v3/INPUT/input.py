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
animal_data_path = 'animal_data/'
met_data_path = 'met_data/'
crop_data_path = 'crop_data/'
soil_data_path ='soil_data/'

## open meteorology/parameter/variable files
## meteorological data in 2018: 
##  1) 2m air temperature (K), 
##  2) relative humidity (%), 
##  3) 10m wind speed (m/s), 
##  4) evaporation from soil (m; accumulated), 
##  5-6) soil temperature data (2 levels, K), 
##  7-8) soil moisture data (2 levels, m3/m3),
##  9) rainfall (kg/m^2), 
##  10) surface runoff (m; accumulated) ,
##  11) subsurface runoff (m; accumulated), 
##  12) aerogynamic and boundary layer resistance (s/m)
if CONFIG_machine == "Stream":
    temp_file = xr.open_dataset(infile_path+met_data_path+'ERA5-temp2m-'+str(sim_year)+'d-05.nc')
    rhum_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-rhum2m-'+str(sim_year)+'d-05.nc')
    wind_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-wind10m-'+str(sim_year)+'d-05.nc')
    evap_file = xr.open_dataset(infile_path+met_data_path+'ERA5-evapfromsoil-'+str(sim_year)+'d-05.nc')
    groundtemp_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1temp-'+str(sim_year)+'d-05.nc')
    groundtemp_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2temp-'+str(sim_year)+'d-05.nc')
    soilmoist_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1moist-'+str(sim_year)+'d-05.nc')
    soilmoist_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2moist-'+str(sim_year)+'d-05.nc')
    rain_file = xr.open_dataset(infile_path+met_data_path+'ERA5-totalcolumn_rainwater-'+str(sim_year)+'d-05.nc')
    runoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-srfrunoff-'+str(sim_year)+'d-05.nc')
    subrunoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-subrunoff-'+str(sim_year)+'d-05.nc')
    ratm_file = xr.open_dataset(infile_path+met_data_path+'I2000clm50_RAM_output.clm2.h1.'+str(sim_year-2000)+'-01-01-00000-05.nc')
else:
    temp_file = xr.open_dataset(infile_path+met_data_path+'ERA5-temp2m-'+str(sim_year)+'h-05.nc')
    rhum_file = xr.open_dataset(infile_path+met_data_path+'ERA5-rhum2m-'+str(sim_year)+'h-05.nc')
    wind_file = xr.open_dataset(infile_path+met_data_path+'ERA5-wind10m-'+str(sim_year)+'h-05.nc')
    evap_file = xr.open_dataset(infile_path+met_data_path+'ERA5-evapfromsoil-'+str(sim_year)+'h-05.nc')
    groundtemp_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1temp-'+str(sim_year)+'h-05.nc')
    groundtemp_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2temp-'+str(sim_year)+'h-05.nc')
    soilmoist_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1moist-'+str(sim_year)+'h-05.nc')
    soilmoist_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2moist-'+str(sim_year)+'h-05.nc')
    rain_file = xr.open_dataset(infile_path+met_data_path+'ERA5-totalcolumn_rainwater-'+str(sim_year)+'h-05.nc')
    runoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-surfrunoff-'+str(sim_year)+'h-05.nc')
    subrunoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-subrunoff-'+str(sim_year)+'h-05.nc')
    ratm_file = xr.open_dataset(infile_path+met_data_path+'I2000clm50-RATM-'+str(sim_year)+'h-05.nc')

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
