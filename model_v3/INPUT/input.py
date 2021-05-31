## import python packages
import numpy as np
import xarray as xr
import netCDF4 as nc4
from pathlib import Path
import os

####################################
## import essential AMCLIM modules
####################################
from CONFIG.config import *        ## animal data file names specified in CONFIG.config.py

## directing input file directory
file_path = '/Users/gilgameshj/Desktop/PhD/project/Chapter1/model_development/global model/AMCLIM/input_files/' 
animal_data_path = 'animal_data/'
met_data_path = 'met_data/'

## Open animal data files and meteorological files
animal_file_name = 'Pig_FAO_Gleam.nc'
animal_file = xr.open_dataset(file_path+animal_data_path+animal_file_name)

## open meteorology files
## test by using 2010 met data
temp_file = xr.open_dataset(file_path+met_data_path+'Regridded_airT_2010.nc')
rhum_file = xr.open_dataset(file_path+met_data_path+'Regridded_rhum_2010.nc')
# wind_file = xr.open_dataset()
# evap_file = xr.open_dataset()
temp_data = temp_file['air'] - 273.15
rhum_data = rhum_file['rhum']
wind_data = np.zeros(mtrx)
evap_data = np.zeros(mtrx)

## open MMS data
try:
    MMS_file = xr.open_dataset(file_path+animal_data_path+MMS_file_name)
    #MMS_data = MMS_file[MMS_type]
except:
    pass
