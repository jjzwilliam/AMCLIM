import numpy as np
import numpy.ma as ma
import math
import xarray as xr
import netCDF4 as nc4
import matplotlib.pyplot as plt
import xesmf as xe

## Earth's radius: 6370 km
Radius = 6371

## function: calculating grid area
def area_lat(dlat,dlon):
    nlat = int(180.0/dlat)
    nlon = int(360.0/dlon)
    area_latidx = np.zeros([nlat,nlon])
    ## define lat index
    lat_idx = 90 - dlat*np.arange(nlat)
    R = Radius * 1000
    ## calculating the surface area of same lat index
    for i in np.arange(nlat):
        area_latidx[i] = (np.pi/180*R**2)*abs(np.sin(np.deg2rad(lat_idx[i]))
                                              -np.sin(np.deg2rad(lat_idx[i]-dlat)))*abs(dlon)
    return area_latidx

## function: calculating grid area given a particular coordinate 
def grid_area(lat1,lat2,lon1,lon2):
    R = Radius * 1000
    ## calculating the surface area of same lat index
    area = (np.pi/180*R**2)*abs(np.sin(np.deg2rad(lat1)-np.sin(np.deg2rad(lat2)))*abs(lon1-lon2))
    return area

## function: simple upscaling
def simple_upscaling(data_0,size_lat,size_lon,plat,plon):
    size_lat = int(size_lat)
    size_lon = int(size_lon)
    rows = int(plat)
    columns = int(plon)
    data = np.zeros([plat,plon])
    data_0 = np.array(data_0)
    for row in np.arange(0,rows):
        for column in np.arange(0,columns):
            data[row,column] =np.nansum(np.nansum(data_0[row*size_lat:(row+1)*size_lat,
                                            column*size_lon:(column+1)*size_lon],
                                                axis=1),axis=0)
    return data

## function: simple weighted upscaling
def simple_weighted_upscaling(data_0,weight_0,size_lat,size_lon,plat,plon):
    size_lat = int(size_lat)
    size_lon = int(size_lon)
    rows = int(plat)
    columns = int(plon)
    data = np.zeros([plat,plon])
    weight = np.zeros([plat,plon])
    data_0 = np.array(data_0)
    weight_0 = np.array(weight_0)
    for row in np.arange(0,rows):
        for column in np.arange(0,columns):
            weight[row,column] = np.nansum(np.nansum(weight_0[row*size_lat:(row+1)*size_lat,
                                            column*size_lon:(column+1)*size_lon],
                                                axis=1),axis=0)
            data[row,column] = np.nansum(np.nansum(data_0[row*size_lat:(row+1)*size_lat,
                                            column*size_lon:(column+1)*size_lon]*weight_0[row*size_lat:(row+1)*size_lat,
                                            column*size_lon:(column+1)*size_lon],
                                                axis=1),axis=0)
    data = data/weight
    return data

## function: regridding using xesmf (Jiawei Zhuang's regridding package)
def xesmf_regrid(ds_in,dr_in,regrid_method,dlat,dlon):
    ds_out = xr.Dataset({'lat': (['lat'], np.arange(-90, 90, dlat)),
                     'lon': (['lon'], np.arange(-180, 180, dlon)),
                    }
                   )
    regridder = xe.Regridder(ds_in, ds_out, regrid_method)
    dr_out = regridder(dr_in)
    regridder.clean_weight_file()
    return dr_out
