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
met_data_path = 'met_data/hourly_data'
crop_data_path = 'crop_data/'
soil_data_path ='soil_data/'
## Open animal data files and meteorological files
animal_file = xr.open_dataset(infile_path+animal_data_path+animal_file_name)

## open meteorology/parameter/variable files
## meteorological data in 2018: 1) temperature (2m air; ground), 2) relative humidity, 3) 10m wind speed, 4) evaporation from soil,
##     5) soil moisture data, 6) percentage of saturation soil moisture 7) sensible heat flux (J/m^2/s)(not used) 8) rainfall (kg/m^2) 
##     9) aerogynamic and boundary layer resistance (s/m)
# temp_file = xr.open_dataset(infile_path+met_data_path+'ERA5-temp2m-2018d-05.nc')
groundtemp_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1temp-2018h-05.nc')
groundtemp_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2temp-2018h-05.nc')
# rhum_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-rhum2m-2018d-05.nc')
# wind_file = xr.open_dataset(infile_path+met_data_path+'AgERA5-wind10m-2018d-05.nc')
# evap_file = xr.open_dataset(infile_path+met_data_path+'ERA5-evapfromsoil-2018d-05.nc')
# soilmoist_file = xr.open_dataset(infile_path+met_data_path+'SOILMOISTURE-L3S-SSMV-COMBINED-DAILY-2018-360x720.nc')
soilmoist_filelvl1 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl1moist-2018h-05.nc')
soilmoist_filelvl2 = xr.open_dataset(infile_path+met_data_path+'ERA5-soillvl2moist-2018h-05.nc')
# soilsm_file = xr.open_dataset(infile_path+met_data_path+'SOILMOISTURE-L3S-SSMS-ACTIVE-DAILY-2018-360x720.nc')
# rain_file = xr.open_dataset(infile_path+met_data_path+'ERA5-totalcolumn_rainwater-2018d-05.nc')
ratm_file = xr.open_dataset(infile_path+met_data_path+'I2000clm50-RATM-2018h-05.nc')
runoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-surfrunoff-2018h-05.nc')
subrunoff_file = xr.open_dataset(infile_path+met_data_path+'ERA5-subrunoff-2018h-05.nc')

# temp_data = temp_file['t2m'] - 273.15  ## degC
# groundtemp_datalvl1 = groundtemp_filelvl1['stl1'] - 273.15  ## degC
# groundtemp_datalvl2 = groundtemp_filelvl2['stl2'] - 273.15  ## degC
# rhum_data = rhum_file['Relative_Humidity_2m_06h']  ## per cent
# wind_data = wind_file['Wind_Speed_10m_Mean']  ## m/s
# evap_data = evap_file['evabs']*(-1e6)  ## g/day
# # soilmoist_data = soilmoist_file['sm']  ## m3/m3
# soilmoist_datalvl1 = soilmoist_filelvl1['swvl1']  ## m3/m3
# soilmoist_datalvl2 = soilmoist_filelvl2['swvl2']  ## m3/m3
# persm_data = soilsm_file['sm']  ## per cent
# # sshf_data = sshf_file['sshf']/(24*3600)  ## J/m2/s
# rain_data = rain_file['tcrw']*1000  ## g/m2
# ram1_data = ratm_file['RAM1']  ## s/m
# rb1_data = ratm_file['RB1']  ## s/m
# runoff_data = runoff_file['sro']  ## m/day
# subrunoff_data = subrunoff_file['ssro']  ## m/day

##################################
## fill land input data
##################################
# groundtemp_datalvl1 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=groundtemp_datalvl1)  ## degC
# groundtemp_datalvl2 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=groundtemp_datalvl2)  ## degC
# rhum_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=rhum_data)  ## per cent
# wind_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=wind_data)  ## m/s
# evap_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=evap_data) ## g/day
# soilmoist_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_data)  ## m3/m3
# # soilmoist_datalvl1 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_datalvl1)  ## m3/m3
# # soilmoist_datalvl2 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_datalvl2)  ## m3/m3
# persm_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=persm_data)  ## per cent
# ram1_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=ram1_data)  ## s/m
# rb1_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=rb1_data)  ## s/m
# runoff_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=runoff_data)  ## m/day
# subrunoff_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=subrunoff_data)  ## m/day

###############################################
## insert an extra time slice to met fields
###############################################
# temp_data = insert_time_slice(temp_data.values)
# groundtemp_datalvl1 = insert_time_slice(groundtemp_datalvl1)
# groundtemp_datalvl2 = insert_time_slice(groundtemp_datalvl2)
# rhum_data = insert_time_slice(rhum_data)
# wind_data = insert_time_slice(wind_data)
# evap_data = insert_time_slice(evap_data)
# # soilmoist_data = insert_time_slice(soilmoist_data)
# soilmoist_datalvl1 = insert_time_slice(soilmoist_datalvl1)
# soilmoist_datalvl2 = insert_time_slice(soilmoist_datalvl2)
# persm_data = insert_time_slice(persm_data)
# rain_data = insert_time_slice(rain_data)
# ram1_data = insert_time_slice(ram1_data)
# rb1_data = insert_time_slice(rb1_data)
# runoff_data = insert_time_slice(runoff_data)
# subrunoff_data = insert_time_slice(subrunoff_data)

#temp_file = xr.open_dataset(infile_path+met_data_path+'Regridded_airT_2010.nc')
#rhum_file = xr.open_dataset(infile_path+met_data_path+'Regridded_rhum_2010.nc')
#temp_data = temp_file['air'] - 273.15
#rhum_data = rhum_file['rhum']
#wind_data = np.zeros(mtrx)
#evap_data = np.zeros(mtrx)

## open MMS data
try:
    MMS_file = xr.open_dataset(infile_path+animal_data_path+MMS_file_name)
    #MMS_data = MMS_file[MMS_type]
except:
    pass

#################################
## crop and fertilizer data
#################################
fertfilename = 'IFA_fert_info_21st_century.nc'
cropfileformat = '_cropping_info.nc'
crop_filledcalendar = 'ALL_CROPS_netCDF_0.5deg_filled_dir/'
crop_unfilledcalendar = 'ALL_CROPS_netCDF_0.5deg_unfilled_dir/'
crop_filledcalendarformat = '.crop.calendar.fill.nc'
crop_unfilledcalendarformat = '.crop.calendar.nc'

##################################
## soil data
##################################
soilpHfile = 'topsoil_pH_H2O.nc'
soilCECfile = 'topsoil_CEC_clay.nc'
soilclayfile = 'topsoil_clay.nc'
soilbdfile = 'topsoil_BD.nc'
