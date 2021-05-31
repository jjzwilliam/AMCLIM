import numpy as np
import xarray as xr

def xr_to_np(arrayin):
    try:
        arrayout = arrayin.values
    except:
        pass
        arrayout = arrayin
    return arrayout
