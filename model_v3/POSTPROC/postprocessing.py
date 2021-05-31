import numpy as np
import xarray as xr
import matplotlib.pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as feature
from matplotlib import colors
from CONFIG.config import *

## data visualization: plotting monthly emissions and Pv rate; housing mainly
def plot_monthlysim(emiss_data,pv_data,simmonth):
    fig, axes = plt.subplots(2, 1, figsize=[20,20], subplot_kw={'projection':ccrs.PlateCarree()})
    lvls1 = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
    im1 = (emiss_data/1e9).plot(ax=axes[0], norm=colors.LogNorm(), vmin = 1e-10, vmax = 1, 
                                   cmap ='RdYlGn_r', 
                                   add_colorbar = False)

    cax1 = fig.add_axes([axes[0].get_position().x1+0.02,axes[0].get_position().y0,0.015,axes[0].get_position().height])
    cb1 = plt.colorbar(im1, cax=cax1)
    cb1.set_label(label='\n Gg N yr$\mathregular{^{-1}}$', size=20)
    cb1.ax.tick_params(labelsize='large')
    cb1.set_ticks(lvls1)
    cb1.set_ticklabels(['10$\mathregular{^{-10}}$','10$\mathregular{^{-9}}$','10$\mathregular{^{-8}}$',
                        '10$\mathregular{^{-7}}$','10$\mathregular{^{-6}}$','10$\mathregular{^{-5}}$',
                        '10$\mathregular{^{-4}}$','10$\mathregular{^{-3}}$','10$\mathregular{^{-2}}$',
                        '10$\mathregular{^{-1}}$','1'])
    cb1.ax.tick_params(labelsize=15)

    # lvls2 = [0,3,6,9,12,15,18,21,24,30,36,42,48,56,64,72]
    lvls2 = [0,8,16,24,32,40,48,56,64,72]
    im2= pv_data.plot(ax=axes[1], norm = colors.BoundaryNorm(lvls2,len(lvls2)),
                 cmap='coolwarm',vmin=0, vmax = 72, 
                     add_colorbar = False)

    cax2 = fig.add_axes([axes[1].get_position().x1+0.02,axes[1].get_position().y0,0.015,axes[1].get_position().height])
    cb2 = plt.colorbar(im2, cax= cax2)
    cb2.set_label(label='\n Pv (%)', size=20)
    cb2.ax.tick_params(labelsize='large')
    cb2.set_ticks(lvls2)

    cb2.ax.tick_params(labelsize=15)     

    for ax in axes:
        ax.gridlines(linestyle='--')
        ax.add_feature(feature.COASTLINE, linewidth=2)
        ax.add_feature(feature.BORDERS, linewidth=1)

    axes[0].set_title('AMCLIM-Urea '+Months_name[simmonth-1]+' simulation\n',size=20)
    plt.show()
    return

## data visualization: field plottings
def plot_fields(field_data,plot_title):
    
    fig, axes = plt.subplots(1, 1, figsize=[20,10], subplot_kw={'projection':ccrs.PlateCarree()})
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.coastlines()
    ax.add_feature(feature.BORDERS)
    ax.gridlines(linestyle='--')
    
    maxval = np.nanmax(field_data)
    print(maxval)
    if maxval <=1000:
        lvls = np.arange(0,np.round(maxval+10,0),5)
        im1 = (field_data).plot(norm = colors.BoundaryNorm(lvls,len(lvls)), vmin = 0, vmax = lvls[-1], 
                                       cmap ='rainbow', 
                                       cbar_kwargs={'shrink':0.9, 'label': ''})
        ax.set_title(str(plot_title),fontsize = 20)
    elif maxval>1000:
        lvls = [1e-10,1e-9,1e-8,1e-7,1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1]
        im1 = (field_data/1e9).plot(norm = colors.LogNorm(), vmin = 1e-10, vmax = 1, 
                                       cmap ='RdYlGn_r', 
                                       cbar_kwargs={'shrink':0.9, 'label': ''})

        ax.set_title(str(plot_title),fontsize = 20)
    plt.show()
    return

## writing netcdf file
## 
def write_netcdf_file(file_info,file_data):
    ## define dimensions
    ## file_shape should be a list: e.g. mtrx
    file_shape = file_info['shape']
    time = int(file_shape[0])
    lats = int(file_shape[1])
    lons = int(file_shape[2])
    
    dlat = 180.0/lats
    dlon = 360.0/lons

    lat_coords = 90.0 - np.arange(lats) * dlat 
    lon_coords = -180.0 + np.arange(lons) * dlon

    meta_data = xr.Dataset({'time': time, 'lat': lat, 'lon':lon})
    meta_data = meta_data.assign_coords(time = np.arange(1,time+1), lat = lat_coords, lon = lon_coords)
    meta_data.attrs = {'Title': str(file_info['title'])}

    for ii in np.arange(len(file_info['variable'])):
        if len(file_data[ii].shape) == 2:
            meta_data[file_info['variable'][ii]] = (('lat','lon'), file_data[ii])
            meta_data[file_info['variable'][ii]].attrs['unit'] = file_info['unit'][ii]
            meta_data[file_info['variable'][ii]].attrs['long_name'] =  file_info['long_name'][ii]
        elif len(file_data[ii].shape) == 3:
            meta_data[file_info['variable'][ii]] = (('time','lat','lon'), file_data[ii])
            meta_data[file_info['variable'][ii]].attrs['unit'] = file_info['unit'][ii]
            meta_data[file_info['variable'][ii]].attrs['long_name'] =  file_info['long_name'][ii]

    meta_data.to_netcdf(str(file_info['file_name'])+'.nc')
    print(str(file_info['file_name'])+' file complete!')
    return
