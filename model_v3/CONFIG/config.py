import numpy as np
import sys
import os
from pathlib import Path

## specify machine: STREAM, JASMIN, ARCHER
machine = ''

inputpathdict = {
        "STREAM": '/home/s1576984/scratch/working_directory/AMCLIM/input_files/',
        "JASMIN": '/gws/nopw/j04/macaque/JJz/jjz_virtual_env/working_directory/new_AMCLIM/AMCLIM/input_files/',
        "ARCHER": '/work/n02/n02/jjz/working_dir/AMCLIM/input_files/',
}

outputpathdict = {
        "STREAM": '/exports/csce/datastore/geos/users/s1576984/test_transfer/output_ncfiles/',
        "JASMIN": '/gws/nopw/j04/macaque/JJz/jjz_virtual_env/working_directory/new_AMCLIM/AMCLIM/outputs/output_ncfiles/',
        "ARCHER": '/work/n02/n02/jjz/working_dir/AMCLIM/outputs/output_ncfiles/',
}

infile_path = inputpathdict[machine]
output_path = outputpathdict[machine]

sim_year = 2018
## daily simulation or hourly simulation
nhours = 24
if (sim_year%4)==0: Days = 366
else: Days = 365
Hours = nhours * Days
time = Days
timestep = 1

## resolution; e.g. 0.5 degree
dlat = 0.5
dlon = 0.5

## lat/lon dimensions
lats = int(180.0/dlat)
lons = int(360.0/dlon)

## levels: this refers to how many types of practices
## currently two: housing and MMS (manure management system)
levels = 2

## array dimensions;
#mtrx = [levels,time,lats,lons]
mtrx = (time+1,lats,lons)
mtrx2 = (time*2+1,lats,lons)

## define months info
Months_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
        'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
Months_days = [31,28,31,30,31,30,31,31,30,31,30,31]
# Months_idx = [0,31,59,90,120,151,181,212,243,273,304,334,365]
Months_idx = [1,32,60,91,121,152,182,213,244,274,305,335,366]

## specify livestock type and production sytem
## livestock_list = ['CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
livestock = 'PIG'
animal_file_dict = {
        'PIG':'Pig_FAO_Gleam_new.nc',
        'CATTLE':'',
        'POULTRY':''
        }
animal_file_name = animal_file_dict[livestock]      ## input files should be put in AMCLIM/INPUT/
## level index: 
## PIG: industrial-0; intermediate-1, backyard-2
lvl_idx = 2
production_system_dict = {
        'PIG':['industrial','intermediate','backyard'],
        'CATTLE':[],
        'POULTRY':['broiler','layer','backyard']
        }
production_system = production_system_dict[livestock][lvl_idx]
## housing_system: 1. insulated building with pit (or without pit) 2. open/naturally ventilated barn 3. poultry houses
## housing_system_list = ['slat/pit house','barn','poultry_house','bck_poultry']
housing_system_dict = {
        'PIG':['slat/pit house','barn','barn'],
        'CATTLE':[],
        'POULTRY':['poultry_house','poultry_house','bck_poultry']
        }
housing_system = housing_system_dict[livestock][lvl_idx]
## MMS_type: ind, med, bck
MMS_type_dict = {
        'PIG':['ind','med','bck'],
        'CATTLE':[],
        'POULTRY':['bro','lay','bck']
        }
MMS_type = MMS_type_dict[livestock][lvl_idx]
# MMS_file_dict = {}
MMS_file_name ='ind_MMS.nc'                           ## input files should be put in AMCLIM/INPUT/