import numpy as np
import sys
import os
from pathlib import Path

## specify CONFIG_machine: STREAM, JASMIN, ARCHER
CONFIG_machine = 'STREAM'

CONFIG_inputpathdict = {
        "STREAM": '/home/s1576984/scratch/working_directory/AMCLIM/input_files/',
        "JASMIN": '/gws/nopw/j04/macaque/JJz/jjz_virtual_env/working_directory/new_AMCLIM/AMCLIM/input_files/',
        "ARCHER": '/work/n02/n02/jjz/working_dir/AMCLIM/input_files/',
}

CONFIG_outputpathdict = {
        "STREAM": '/exports/csce/datastore/geos/users/s1576984/test_transfer/output_ncfiles/',
        "JASMIN": '/gws/nopw/j04/macaque/JJz/jjz_virtual_env/working_directory/new_AMCLIM/AMCLIM/outputs/output_ncfiles/',
        "ARCHER": '/work/n02/n02/jjz/working_dir/AMCLIM/outputs/output_ncfiles/',
}

infile_path = CONFIG_inputpathdict[CONFIG_machine]
output_path = CONFIG_outputpathdict[CONFIG_machine]

sim_year = 2018
## daily simulation or hourly simulation
nhours = 24
if (sim_year%4)==0: Days = 366
else: Days = 365
Hours = nhours * Days
CONFIG_time = Days
timestep = 1

## resolution; e.g. 0.5 degree
CONFIG_dlat = 0.5
CONFIG_dlon = 0.5

## lat/lon dimensions
CONFIG_lats = int(180.0/CONFIG_dlat)
CONFIG_lons = int(360.0/CONFIG_dlon)

## array dimensions;
CONFIG_mtrx = (CONFIG_time,CONFIG_lats,CONFIG_lons)
CONFIG_mtrx2 = (CONFIG_time*2,CONFIG_lats,CONFIG_lons)

## define months info
Months_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
        'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
Months_days = [31,28,31,30,31,30,31,31,30,31,30,31]
Months_idx = [0,31,59,90,120,151,181,212,243,273,304,334,365]
# Months_idx = [1,32,60,91,121,152,182,213,244,274,305,335,366]

## specify livestock type and production sytem
## livestock_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
# livestock = 'BEEF_CATTLE'
# livestock = 'POULTRY'

## level index: 
## PIG: industrial-0; intermediate-1, backyard-2
## BEEF: grassland-0; mixed-1
## DAIRY: grassland-0; mixed-1
## OTHER CATTLE: grassland-0; mixed-1
## SHEEP: grassland-0; mixed-1
## POULTRY: broiler-0, layer-1, backyard-2
# lvl_idx = 1
CONFIG_production_system_dict = {
        'PIG':['industrial','intermediate','backyard'],
        'BEEF_CATTLE':['grassland','mixed',None],
        'DAIRY_CATTLE':['grassland','mixed',None],
        'OTHER_CATTLE':['grassland','mixed',None],
        'SHEEP':['grassland','mixed',None],
        'POULTRY':['broiler','layer','backyard']
        }
## manure N to land production systems
CONFIG_livestock_prodsyst_Nappland_lict = {
        'PIG':['industrial','intermediate','backyard'],
        'BEEF_CATTLE':['mixed'],
        'DAIRY_CATTLE':['mixed'],
        'OTHER_CATTLE':['mixed'],
        'SHEEP':['mixed'],
        'POULTRY':['broiler','layer','backyard']
        }
# production_system = CONFIG_production_system_dict[livestock][lvl_idx]
## housing_system: 1. insulated building with pit (or without pit) 2. open/naturally ventilated barn 3. poultry houses
## housing_system_list = ['slat/pit house','barn','poultry_house','bck_poultry']
CONFIG_housing_system_dict = {
        'PIG':['insulated','naturally ventilated','naturally ventilated'],
        'BEEF_CATTLE':[None,'naturally ventilated',None],
        'DAIRY_CATTLE':[None,'naturally ventilated',None],
        'OTHER_CATTLE':[None,'naturally ventilated',None],
        'SHEEP':[None,'naturally ventilated',None],
        'POULTRY':['insulated','insulated','naturally ventilated']
        }
# housing_system = housing_system_dict[livestock][lvl_idx]

CONFIG_animal_file_dict = {
        'PIG':'piginfo_faogleam.nc',
        'BEEF_CATTLE':'beefinfo_faogleam.nc',
        'DAIRY_CATTLE':'dairyinfo_faogleam.nc',
        'OTHER_CATTLE':'otherdairyinfo_faogleam.nc',
        'SHEEP':'sheepinfo_faogleam.nc',
        'POULTRY':'chickeninfo_faogleam.nc'
        }
# animal_file_name = animal_file_dict[livestock]      ## input files should be put in AMCLIM/INPUT/

CONFIG_MMS_file_dict = {
        'PIG':'pigmms_faogleam.nc',
        'BEEF_CATTLE':'beefmms_faogleam.nc',
        'DAIRY_CATTLE':'dairymms_faogleam.nc',
        'OTHER_CATTLE':'otherdairymms_faogleam.nc',
        'SHEEP':'sheepmms_faogleam.nc',
        'POULTRY':'chickenmms_faogleam.nc'
        }