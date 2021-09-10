## import python packages
import xarray as xr
import netCDF4 as nc4
from pathlib import Path


####################################
## import essential AMCLIM modules
####################################
from CONFIG.config import *        ## animal data file names specified in CONFIG.config.py

## directing input file directory
file_path = '/home/s1576984/scratch/working_directory/AMCLIM/input_files/' 
animal_data_path = 'animal_data/'
met_data_path = 'met_data/'

## Open animal data files and meteorological files
animal_file = xr.open_dataset(file_path+animal_data_path+animal_file_name)

## open meteorology files
## meteorological data in 2018: 1) temperature (2m air), 2) relative humidity, 3) 10m wind speed, 4) evaporation from soil,
##     5) soil moisture data, 6) percentage of saturation soil moisture 7) sensible heat flux
temp_file = xr.open_dataset(file_path+met_data_path+'ERA5_2018d_meant2m05.nc')
rhum_file = xr.open_dataset(file_path+met_data_path+'AgERA5_2018d_meanRH2m05.nc')
wind_file = xr.open_dataset(file_path+met_data_path+'AgERA5_2018d_10mwind05.nc')
evap_file = xr.open_dataset(file_path+met_data_path+'ERA5_2018d_evapfromsoil_dailytotal.nc')
soilmoist_file = xr.open_dataset(file_path+met_data_path+'SOILMOISTURE-L3S-SSMV-COMBINED-DAILY-2018-360x720.nc')
soilsm_file = xr.open_dataset(file_path+met_data_path+'SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2018-360x720.nc')
sshf_file = xr.open_dataset(file_path+met_data_path+'ERA5_2018d_sshf05.nc')

temp_data = temp_file['t2m'] - 273.15
rhum_data = rhum_file['Relative_Humidity_2m_06h']
wind_data = wind_file['Wind_Speed_10m_Mean']
evap_data = evap_file['evabs']*(-1000)
soilmoist_data = soilmoist_file['sm']
persm_data = soilsm_file['sm']
sshf_data = sshf_file['sshf']

#temp_file = xr.open_dataset(file_path+met_data_path+'Regridded_airT_2010.nc')
#rhum_file = xr.open_dataset(file_path+met_data_path+'Regridded_rhum_2010.nc')
#temp_data = temp_file['air'] - 273.15
#rhum_data = rhum_file['rhum']
#wind_data = np.zeros(mtrx)
#evap_data = np.zeros(mtrx)

## open MMS data
try:
    MMS_file = xr.open_dataset(file_path+animal_data_path+MMS_file_name)
    #MMS_data = MMS_file[MMS_type]
except:
    pass
