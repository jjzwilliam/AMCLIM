## import python packages
import numpy as np
import xarray as xr
import netCDF4 as nc4
from pathlib import Path
import os

## Open animal data files and meteorological files
file_path = '/home/s1576984/scratch/working_directory/AMCLIM/INPUT/'
animal_file_name = 'Pig_FAO_Gleam.nc'

animal_data = xr.open_dataset(file_path+animal_file_name)
excretN_info = animal_data['Excreted_N'][0]
animal_head = animal_data['Animal_head'][0]
try:
    animal_weight = animal_data['Animal_weight'][0]
except:
    pass
#print(np.nanmedian(animal_weight.values[np.where(animal_weight!=0)]))
animal_weight[np.where((animal_head!=0)&(animal_weight==0))] = np.nanmedian(animal_weight.values[np.where(animal_weight!=0)])
animal_Nrate = animal_data['N_rate'][0]
#print(np.nanmedian(animal_Nrate.values[np.where(animal_Nrate!=0)]))
animal_Nrate[np.where((animal_head!=0)&(animal_Nrate==0))] = np.nanmedian(animal_Nrate.values[np.where(animal_Nrate!=0)])
## pig density 1: 120 kg/m^2
animal_density = 120.0 
massgrid = animal_head*animal_weight
farming_area = animal_head*animal_weight/animal_density
## total animal mass per grid
# massgrid = 1000*excretN_info/(animal_Nrate*365)
# farming_area = massgrid/animal_density
## pig density 2: 1 head/m^2
# animal_density = 1.0
# farming_area = animal_head/animal_density

## open meteorology files
## test by using 2010 met data
temp_file = xr.open_dataset(file_path+'Regridded_airT_2010.nc')
rhum_file = xr.open_dataset(file_path+'Regridded_rhum_2010.nc')
# wind_file = xr.open_dataset()
temp_data = temp_file['air'] - 273.15
rhum_data = rhum_file['rhum']
wind_data = 0.0

## test code:
#print(np.nansum(farming_area)/1e9)
