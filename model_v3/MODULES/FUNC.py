import numpy as np
import xarray as xr

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
    
    for dd in np.arange(365):
        output_field[dd][(sd_template!=0)&np.isnan(input_field[dd])] = time_fill[dd]        
        output_field[dd][(sd_template!=0)&(input_field[dd]==0)] = time_fill[dd]    
    return output_field