import numpy as np
import xarray as xr

def open_ds(filename):
    with xr.open_dataset(filename) as ds:
        ds.load()
    return ds

def xr_to_np(arrayin):
    try:
        arrayout = arrayin.values
    except:
        pass
        arrayout = arrayin
    return arrayout

## function: fill nan values (very simple method)
def field_var_fill(sd_template,input_field):
    ## duplicate input field
    output_field = np.array(input_field)
    ## calculated filled values, 1) by lat, 2) by lon, 3) global average
#     lat_fill = np.nanmean(input_field,axis=2)
#     lon_fill = np.nanmean(input_field,axis=1)
    time_fill = np.nanmean(input_field,axis=(1,2))
    days = input_field.shape[0]
    for dd in np.arange(days):
        output_field[dd][(sd_template!=0)&np.isnan(input_field[dd])] = time_fill[dd]        
        output_field[dd][(sd_template!=0)&(input_field[dd]==0)] = time_fill[dd]    
    return output_field

## function: insert an extra slice (time) to array
def insert_time_slice(input_array):
    shape = input_array.shape
    new_time_dim = int(shape[0])+1
    new_array = np.zeros([new_time_dim,shape[1],shape[2]])
    new_array[1:,:,:] = input_array
    new_array[0,:,:] = input_array[-1,:,:]
    return new_array

## function: calculation total Gg of each pathways, and rounding to 3 decimal numbers
def sum_totalGg(input):
    total = np.round(np.nansum(input)/1e9,decimals=3)
    return total
