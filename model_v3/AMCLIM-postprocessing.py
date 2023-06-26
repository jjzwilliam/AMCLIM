import numpy as np
import xarray as xr
import time
import sys
import pandas as pd

import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.LAND as LAND


yearofstudy = int(sys.argv[1])

crop_list = ['barley','cassava','cotton','groundnut','maize','millet','potato','rapeseed',
            'rice','rye','sorghum','soybean','sugarbeet','sugarcane','sunflower','wheat']
ruminants_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','SHEEP','GOAT']

## define months info
Months_name = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul',
        'Aug', 'Sep', 'Oct', 'Nov', 'Dec']
Months_days = [31,28,31,30,31,30,31,31,30,31,30,31]
Months_idx = [0,31,58,89,119,150,180,211,242,272,303,333,364]
# Months_idx = [1,32,60,91,121,152,182,213,244,274,305,335,366]

ferts = ['ammonium','urea']
techs = ['surf','disk','injection']

refyear = 2000
nyears = 20
nlats = 360
nlons = 720

lats = 90 - 0.5*np.arange(nlats)
lons = -180 + 0.5*np.arange(nlons)
years = np.arange(2000,2020).astype(int)

sim_year = CONFIG.sim_year
Days = CONFIG.Days
nlat = CONFIG.CONFIG_lats
nlon = CONFIG.CONFIG_lons

ntime = Days
yearidx = str(yearofstudy)+'-01-01'
times = pd.date_range(yearidx,periods=ntime)

techs = ["surf","disk","injection"]
ntechs = 3
techlvls = np.arange(ntechs)
ncrops = 16
croplvls = np.arange(ncrops)
ferts = ["ammonium","urea"]
nferts = 2
fertlvls = np.arange(nferts)

nmms = 6
napp = 4
mmslvl = np.arange(nmms)
applvl = np.arange(napp)


# path = "/exports/csce/datastore/geos/users/s1576984/test_transfer/"
path = CONFIG.output_path + "requisite_files/"
ncfilepath = CONFIG.output_path + "hourly_sims/"
infile_path = CONFIG.infile_path
outpath = ncfilepath + "sectoral_emissions/"

def tech_matrix():
    ## tech used in countries determined by income level by WB data
    tech_algorithim = [[0.70,0.20,0.10],
                    [0.8,0.15,0.05],
                    [0.95,0.05,0.0],
                    [1.0,0.0,0.0]]

    ## geographical distribution of tech used in fertilizer application
    global_tech_used = np.zeros((3,nlat,nlon))

    ## read WB country income level
    wbincomelvlds = xr.open_dataset(path+'AMCLIM_fert_info/WB_21st_century_income_classification.nc')
    wb_incomelvl = wbincomelvlds.income_lvl.sel(year=yearofstudy)

    for ii in np.arange(3):
        for jj in np.arange(1,5):
            global_tech_used[ii][wb_incomelvl.values==jj] = tech_algorithim[jj-1][ii]
    return global_tech_used

def tech_matrix():
    ## tech used in countries determined by income level by WB data
    tech_algorithim = [[0.70,0.20,0.10],
                    [0.8,0.15,0.05],
                    [0.95,0.05,0.0],
                    [1.0,0.0,0.0]]

    ## geographical distribution of tech used in fertilizer application
    global_tech_used = np.zeros((3,nlat,nlon))

    ## read WB country income level
    wbincomelvlds = xr.open_dataset(path+'AMCLIM_fert_info/WB_21st_century_income_classification.nc')
    wb_incomelvl = wbincomelvlds.income_lvl.sel(year=yearofstudy)

    for ii in np.arange(3):
        for jj in np.arange(1,5):
            global_tech_used[ii][wb_incomelvl.values==jj] = tech_algorithim[jj-1][ii]
    return global_tech_used

def chemfert_output():
    ## N pathways
    NH3emiss = np.zeros((nferts, ntechs, ncrops, Days, nlat, nlon))
    TANrunoff = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    nitrif = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    NH4uptake = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    TANleaching = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    TANdiffaq = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    TANdiffgas = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    NO3uptake = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    NO3runoff = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    NO3diff = np.zeros((nferts, ntechs, ncrops, nlat, nlon))
    NO3leaching = np.zeros((nferts, ntechs, ncrops, nlat, nlon))

    ## read data from sim outputs
    for ii in fertlvls:
        fert = ferts[ii]
        print("fert",ii,fert)
        for jj in techlvls:
            tech = techs[jj]
            print("tech",jj,tech)
            for kk in croplvls:
                crop = crop_list[kk]
                print("crop",kk,crop)
                
                ## open nc files
                simdf = xr.open_dataset(ncfilepath+'chemfert/base.'+crop+'.'+fert+'.'+tech+'.'+str(yearofstudy)+'.pp.nc')
                ## read data
                NH3emiss[ii][jj][kk] = simdf.NH3emiss.fillna(0)
                TANrunoff[ii][jj][kk] = simdf.TANwashoff.fillna(0)
                nitrif[ii][jj][kk] = simdf.NH4nitrif.fillna(0)
                NH4uptake[ii][jj][kk] = simdf.NH4uptake.fillna(0)
                TANleaching[ii][jj][kk] = simdf.TANleaching.fillna(0)
                TANdiffaq[ii][jj][kk] = simdf.TANdiffaq.fillna(0)
                TANdiffgas[ii][jj][kk] = simdf.NH3diffgas.fillna(0)
                NO3uptake[ii][jj][kk] = simdf.NO3uptake.fillna(0)
                NO3runoff[ii][jj][kk] = simdf.NO3washoff.fillna(0)
                NO3diff[ii][jj][kk] = simdf.NO3diff.fillna(0)
                NO3leaching[ii][jj][kk] = simdf.NO3leaching.fillna(0)
                print("read over.")

    ## N app data
    Napp = np.zeros((nferts+1, ncrops, nlat, nlon))
    for kk in croplvls:
        crop = crop_list[kk]
        print("crop",kk,crop)
        cropinfo = xr.open_dataset(path+'AMCLIM_fert_info/chemfertinfo.'+crop+'.21st.nc')
        Napp[0][kk] = np.nan_to_num(cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.ammN_area.sel(year=yearofstudy))
        Napp[1][kk] = np.nan_to_num(cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.ureaN_area.sel(year=yearofstudy))
        Napp[2][kk] = np.nan_to_num(cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.nitN_area.sel(year=yearofstudy))

    Napp_ds = xr.Dataset(
        data_vars=dict(
            amm_Napp=(['croplvl','lat','lon'],Napp[0]/1e3),
            urea_Napp=(['croplvl','lat','lon'],Napp[1]/1e3),
            nit_Napp=(['croplvl','lat','lon'],Napp[2]/1e3),
        ),
        coords = dict(
            croplvl=(["croplvl"], croplvls),
            lon=(["lon"], lons),
            lat=(["lat"], lats),
                    ),
        attrs=dict(
            description="AMCLIM: chemical fertilizer application in " +str(yearofstudy),
            info = "Fertilizer level: 0-ammonium, 1-urea. 3-nitrate (Napps are listed seperately) ",
            units= "gN per grid (emission); kgN per grid (application)",
        ),
    )

    ## component N pathways
    compNH3emiss = np.zeros((nferts, ncrops, Days, nlat, nlon))
    compTANrunoff = np.zeros((nferts, ncrops, nlat, nlon))
    compnitrif = np.zeros((nferts, ncrops, nlat, nlon))
    compNH4uptake = np.zeros((nferts, ncrops, nlat, nlon))
    compTANleaching = np.zeros((nferts, ncrops, nlat, nlon))
    compTANdiffaq = np.zeros((nferts, ncrops, nlat, nlon))
    compTANdiffgas = np.zeros((nferts, ncrops, nlat, nlon))
    compNO3uptake = np.zeros((nferts, ncrops, nlat, nlon))
    compNO3runoff = np.zeros((nferts, ncrops, nlat, nlon))
    compNO3diff = np.zeros((nferts, ncrops, nlat, nlon))
    compNO3leaching = np.zeros((nferts, ncrops, nlat, nlon))

    global_tech_used = tech_matrix()

    compNH3emiss = global_tech_used[0]*NH3emiss[:,0]+\
                    global_tech_used[1]*NH3emiss[:,1]+\
                        global_tech_used[2]*NH3emiss[:,2]
    compTANrunoff = global_tech_used[0]*TANrunoff[:,0]+\
                        global_tech_used[1]*TANrunoff[:,1]+\
                            global_tech_used[2]*TANrunoff[:,2]
    compnitrif = global_tech_used[0]*nitrif[:,0]+\
                        global_tech_used[1]*nitrif[:,1]+\
                            global_tech_used[2]*nitrif[:,2]
    compNH4uptake = global_tech_used[0]*NH4uptake[:,0]+\
                        global_tech_used[1]*NH4uptake[:,1]+\
                            global_tech_used[2]*NH4uptake[:,2]
    compTANleaching = global_tech_used[0]*TANleaching[:,0]+\
                        global_tech_used[1]*TANleaching[:,1]+\
                            global_tech_used[2]*TANleaching[:,2]
    compTANdiffaq = global_tech_used[0]*TANdiffaq[:,0]+\
                        global_tech_used[1]*TANdiffaq[:,1]+\
                            global_tech_used[2]*TANdiffaq[:,2]
    compTANdiffgas = global_tech_used[0]*TANdiffgas[:,0]+\
                        global_tech_used[1]*TANdiffgas[:,1]+\
                            global_tech_used[2]*TANdiffgas[:,2]
    compNO3uptake = global_tech_used[0]*NO3uptake[:,0]+\
                        global_tech_used[1]*NO3uptake[:,1]+\
                            global_tech_used[2]*NO3uptake[:,2]
    compNO3runoff = global_tech_used[0]*NO3runoff[:,0]+\
                        global_tech_used[1]*NO3runoff[:,1]+\
                            global_tech_used[2]*NO3runoff[:,2]
    compNO3diff = global_tech_used[0]*NO3diff[:,0]+\
                        global_tech_used[1]*NO3diff[:,1]+\
                            global_tech_used[2]*NO3diff[:,2]
    compNO3leaching = global_tech_used[0]*NO3leaching[:,0]+\
                        global_tech_used[1]*NO3leaching[:,1]+\
                            global_tech_used[2]*NO3leaching[:,2]

    for ii in techlvls:
        NH3emiss[:,ii] = global_tech_used[ii]*NH3emiss[:,ii]
        
        TANrunoff[:,ii] = global_tech_used[ii]*TANrunoff[:,ii]
        
        nitrif[:,ii] = global_tech_used[ii]*nitrif[:,ii]
        
        NH4uptake[:,ii] = global_tech_used[ii]*NH4uptake[:,ii]
        
        TANleaching[:,ii] = global_tech_used[ii]*TANleaching[:,ii]
        
        TANdiffaq[:,ii] = global_tech_used[ii]*TANdiffaq[:,ii]
        
        TANdiffgas[:,ii] = global_tech_used[ii]*TANdiffgas[:,ii]
        
        NO3uptake[:,ii] = global_tech_used[ii]*NO3uptake[:,ii]
        
        NO3runoff[:,ii] = global_tech_used[ii]*NO3runoff[:,ii]
        
        NO3diff[:,ii] = global_tech_used[ii]*NO3diff[:,ii]
        
        NO3leaching[:,ii] = global_tech_used[ii]*NO3leaching[:,ii]
        
    alltech_simds = xr.Dataset(
        data_vars=dict(
            NH3emiss=(['fertlvl','techlvl','croplvl','time','lat','lon'],NH3emiss),
            TANwashoff=(['fertlvl','techlvl','croplvl','lat','lon'],TANrunoff),
            TANleaching=(['fertlvl','techlvl','croplvl','lat','lon'],TANleaching),
            TANdiffaq=(['fertlvl','techlvl','croplvl','lat','lon'],TANdiffaq),
            NH3diffgas=(['fertlvl','techlvl','croplvl','lat','lon'],TANdiffgas),
            NH4nitrif=(['fertlvl','techlvl','croplvl','lat','lon'],nitrif),
            NH4uptake=(['fertlvl','techlvl','croplvl','lat','lon'],NH4uptake),
            NO3uptake=(['fertlvl','techlvl','croplvl','lat','lon'],NO3uptake),
            NO3washoff=(['fertlvl','techlvl','croplvl','lat','lon'],NO3runoff),
            NO3leaching=(['fertlvl','techlvl','croplvl','lat','lon'],NO3leaching),
            NO3diff=(['fertlvl','techlvl','croplvl','lat','lon'],NO3diff),
        ),
        coords = dict(
            fertlvl=(["fertlvl"], fertlvls),
            techlvl=(["techlvl"], techlvls),
            croplvl=(["croplvl"], croplvls),
            time=(["time"], times),
            lon=(["lon"], lons),
            lat=(["lat"], lats),
                    ),
        attrs=dict(
            description="AMCLIM: Annaul NH3 emissions from chemical fertilizer application in " +str(yearofstudy),
            info = "NH3 emission from ammonium and urea fertilizer application.\
                    Fertilizer level: 0-ammonium, 1-urea. 3-nitrate (Napps are listed seperately) ",
            units= "gN per grid (emission); kgN per grid (application)",
        ),
    )

    ## save timeseries ncfile
    outds = xr.Dataset(
        data_vars=dict(
            NH3emiss=(['fertlvl','time','lat','lon'],alltech_simds.NH3emiss.sum(dim=['techlvl','croplvl']).values),
            amm_Napp=(['lat','lon'],Napp_ds.amm_Napp.sum(dim='croplvl').values),
            urea_Napp=(['lat','lon'],Napp_ds.urea_Napp.sum(dim='croplvl').values),
            nit_Napp=(['lat','lon'],Napp_ds.nit_Napp.sum(dim='croplvl').values),
        ),
        coords = dict(
            fertlvl=(["fertlvl"], fertlvls),
            time=(["time"], times),
            lon=(["lon"], lons),
            lat=(["lat"], lats),
                    ),
        attrs=dict(
            description="AMCLIM: NH3 emissions from chemical fertilizer application in " +str(yearofstudy),
            info = "NH3 emission from ammonium and urea fertilizer application.\
                    Fertilizer level: 0-ammonium, 1-urea. 3-nitrate (Napps are listed seperately) ",
            units= "gN per grid (emission); kgN per grid (application)",
        ),
    )

    outds.NH3emiss.attrs["unit"] = 'gN/day'
    outds.NH3emiss.attrs["long name"] = 'NH3 emission from fertilizer application'
    outds.amm_Napp.attrs["unit"] = 'kgN/year'
    outds.amm_Napp.attrs["long name"] = 'Chemical fertilizer (ammonium) application'  
    outds.urea_Napp.attrs["unit"] = 'kgN/year'
    outds.urea_Napp.attrs["long name"] = 'Chemical fertilizer (urea) application'  
    outds.nit_Napp.attrs["unit"] = 'kgN/year'
    outds.nit_Napp.attrs["long name"] = 'Chemical fertilizer (nitrate) application'

    outfilename = 'base.chemfert.NH3emission.timeseries.'+str(yearofstudy)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

    ## save annual data ncfile
    outds = xr.Dataset(
        data_vars=dict(
            amm_Napp=(['croplvl','lat','lon'],Napp[0]/1e3),
            urea_Napp=(['croplvl','lat','lon'],Napp[1]/1e3),
            nit_Napp=(['croplvl','lat','lon'],Napp[2]/1e3),
            NH3emiss=(['fertlvl','techlvl','croplvl','lat','lon'],alltech_simds.NH3emiss.sum(dim='time').values),
            TANwashoff=(['fertlvl','techlvl','croplvl','lat','lon'],TANrunoff),
            TANleaching=(['fertlvl','techlvl','croplvl','lat','lon'],TANleaching),
            TANdiffaq=(['fertlvl','techlvl','croplvl','lat','lon'],TANdiffaq),
            NH3diffgas=(['fertlvl','techlvl','croplvl','lat','lon'],TANdiffgas),
            NH4nitrif=(['fertlvl','techlvl','croplvl','lat','lon'],nitrif),
            NH4uptake=(['fertlvl','techlvl','croplvl','lat','lon'],NH4uptake),
            NO3uptake=(['fertlvl','techlvl','croplvl','lat','lon'],NO3uptake),
            NO3washoff=(['fertlvl','techlvl','croplvl','lat','lon'],NO3runoff),
            NO3leaching=(['fertlvl','techlvl','croplvl','lat','lon'],NO3leaching),
            NO3diff=(['fertlvl','techlvl','croplvl','lat','lon'],NO3diff),
        ),
        coords = dict(
            fertlvl=(["fertlvl"], fertlvls),
            techlvl=(["techlvl"], techlvls),
            croplvl=(["croplvl"], croplvls),
            time=(["time"], times),
            lon=(["lon"], lons),
            lat=(["lat"], lats),
                    ),
        attrs=dict(
            description="AMCLIM: Annaul NH3 emissions from chemical fertilizer application in " +str(yearofstudy),
            info = "NH3 emission from ammonium and urea fertilizer application.\
                    Fertilizer level: 0-ammonium, 1-urea. 3-nitrate (Napps are listed seperately) ",
            units= "gN per grid (emission); kgN per grid (application)",
        ),
    )

    outds.amm_Napp.attrs["unit"] = 'kgN/year'
    outds.amm_Napp.attrs["long name"] = 'Chemical fertilizer (ammonium) application'  
    outds.urea_Napp.attrs["unit"] = 'kgN/year'
    outds.urea_Napp.attrs["long name"] = 'Chemical fertilizer (urea) application'  
    outds.nit_Napp.attrs["unit"] = 'kgN/year'
    outds.nit_Napp.attrs["long name"] = 'Chemical fertilizer (nitrate) application'
    outds.NH3emiss.attrs["unit"] = 'gN/day'
    outds.NH3emiss.attrs["long name"] = 'NH3 emission from fertilizer application'
    outds.TANwashoff.attrs["unit"] = 'gN/day'
    outds.TANwashoff.attrs["long name"] = 'TAN washoff from fertilizer application'
    outds.TANleaching.attrs["unit"] = 'gN/day'
    outds.TANleaching.attrs["long name"] = 'TAN leaching from fertilizer application'
    outds.TANdiffaq.attrs["unit"] = 'gN/day'
    outds.TANdiffaq.attrs["long name"] = 'TAN diffusion to deep soil from fertilizer application'
    outds.NH3diffgas.attrs["unit"] = 'gN/day'
    outds.NH3diffgas.attrs["long name"] = 'NH3 diffusion to deep soil from fertilizer application'
    outds.NH4nitrif.attrs["unit"] = 'gN/day'
    outds.NH4nitrif.attrs["long name"] = 'Nitrification'
    outds.NH4uptake.attrs["unit"] = 'gN/day'
    outds.NH4uptake.attrs["long name"] = 'Uptake of NH4+ by crops'
    outds.NO3uptake.attrs["unit"] = 'gN/day'
    outds.NO3uptake.attrs["long name"] = 'Uptake of NO3- by crops'
    outds.NO3washoff.attrs["unit"] = 'gN/day'
    outds.NO3washoff.attrs["long name"] = 'Nitrate washoff from fertilizer application'
    outds.NO3leaching.attrs["unit"] = 'gN/day'
    outds.NO3leaching.attrs["long name"] = 'Nitrate leaching from fertilizer application'
    outds.NO3diff.attrs["unit"] = 'gN/day'
    outds.NO3diff.attrs["long name"] = 'Nitrate diffusion to deep soil from fertilizer application'
        
    outfilename = 'base.chemfert.NH3emission.Npathways.annual.'+str(yearofstudy)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

def livestock_output():

    livestock = "PIG"
    print("=======================================")
    print("Processing: ",livestock)

    ## production systems and levels
    prodsysts = CONFIG.CONFIG_production_system_dict[livestock]
    prodsyst_lvl = len(prodsysts)
    ## three types of NH3 emissions: housing, MMS, and land spreading
    housing_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    mms_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    landspreading_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    housing_compNH3 = np.zeros((prodsyst_lvl,nlat,nlon))
    housing_comptotalN = np.zeros((prodsyst_lvl,nlat,nlon))
    mms_compNH3 = np.zeros((prodsyst_lvl,6,nlat,nlon))
    mms_comptotalN = np.zeros((prodsyst_lvl,6,nlat,nlon))
    landspreading_compNH3 = np.zeros((prodsyst_lvl,4,nlat,nlon))
    landspreading_comptotalN = np.zeros((prodsyst_lvl,4,nlat,nlon))

    housing_types = ["slat-pit_house","barn","barn"]
    for prodsyst_idx in np.arange(prodsyst_lvl):
        prodsyst = prodsysts[prodsyst_idx]
        print("=======================================")
        print("Production system: ",prodsyst)
        ## HOUSING NH3 emissions
        housetype = housing_types[prodsyst_idx]
        if prodsyst == "industrial":
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                    str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss_slat + housingds.NH3emiss_pit
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission (slat-pit): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (slat-pit): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'_insitu.'+\
                    str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss_slat + housingds.NH3emiss_pit
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values + housing_compNH3[prodsyst_idx][:]
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret + housing_comptotalN[prodsyst_idx][:]
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values + housing_NH3emission[prodsyst_idx][:]
            
            print("Housing emission (slat-pit,in-situ): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (slat-pit,in-situ): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            print("Housing emission (total slat-pit): ",livestock,prodsyst,np.nansum(housing_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Housing total N (total slat-pit): ",livestock,prodsyst,np.nansum(housing_comptotalN[prodsyst_idx])/1e9," GgN.")

        else:
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                        str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values
            housing_emission = housingds.NH3emiss
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values 
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission: ",livestock,prodsyst,np.nansum(housing_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
        
        ## MMS and emissions
        mmsidx = 0
        print("[MMS]")
        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats: 
                try:
                    print("MMS files: ",mmsidx,mmsphase,mmscat)
                    mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                            str(mmsphase)+'.'+str(yearofstudy)+'.nc'
                    mmsds = xr.open_dataset(mmsfile)
                    mmsemiss = mmsds.NH3emiss
                    mms_NH3emission[prodsyst_idx][:] = mmsemiss.values + mms_NH3emission[prodsyst_idx][:]
                    mms_compNH3[prodsyst_idx][mmsidx][:] = mmsemiss.sum(dim="time").values
                    mms_comptotalN[prodsyst_idx][mmsidx][:] = mmsds.NtotalMMS.values
                    print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
                    print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
                    if prodsyst == "industrial":
                        if mmsphase == "solid":
                            if mmscat == "MMS_open":
                                ## include insit MMS
                                mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                                        str(mmsphase)+'.insitu.'+str(yearofstudy)+'.nc'
                                mmsds = xr.open_dataset(mmsfile)
                                mmsemiss = mmsds.NH3emiss
                                mms_NH3emission[prodsyst_idx][:] = mmsemiss.values + mms_NH3emission[prodsyst_idx][:]
                                mms_compNH3[prodsyst_idx][mmsidx][:] = mmsemiss.sum(dim="time").values + mms_compNH3[prodsyst_idx][mmsidx][:]
                                mms_comptotalN[prodsyst_idx][mmsidx] = np.nan_to_num(mmsds.NtotalMMS.values) + mms_comptotalN[prodsyst_idx][mmsidx]
                                print("MMS emission (insitu MMS): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mmsemiss.sum(dim="time").values)/1e9," GgN.")
                                print("MMS N (insitu MMS): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mmsds.NtotalMMS.values)/1e9," GgN.")
                                print("MMS emission (updated): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
                                print("MMS N (updated): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
                    mmsidx = mmsidx + 1
                except:
                    pass
                    print(livestock,"no "+mmsphase+' '+mmscat+' managemet.')
                    mmsidx = mmsidx + 1
        print("MMS total emission: ",livestock,prodsyst,np.nansum(mms_NH3emission[prodsyst_idx])/1e9," GgN.")
        print("MMS total N: ",livestock,prodsyst,np.nansum(mms_comptotalN[prodsyst_idx])/1e9," GgN.")
        
        ## Land spreading emissions
        mmsidx = 0
        print("[Manure land spreading]")
        for mmstype in ["MMS_cover.liquid","MMS_indoor.liquid","MMS_indoor.solid","MMS_open.liquid"]:
            print("MMS app, ",mmstype)
            manureappfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+'.'+str(mmstype)+\
                        '.surf.'+str(yearofstudy)+'.nc'
            manureappds = xr.open_dataset(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[prodsyst_idx][:] = manureappemiss.values + landspreading_NH3emission[prodsyst_idx][:]
            landspreading_compNH3[prodsyst_idx][mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(yearofstudy)+'.manureapp.nc'
            manureappNds = xr.open_dataset(manureappNfile)
            mmstotalappN = manureappNds.spring_TAN+manureappNds.spring_UAN+\
                        manureappNds.spring_availN+manureappNds.spring_resistN+manureappNds.spring_unavailN+\
                        manureappNds.winter_TAN+manureappNds.winter_UAN+\
                        manureappNds.winter_availN+manureappNds.winter_resistN+manureappNds.winter_unavailN
            landspreading_comptotalN[prodsyst_idx][mmsidx][:] = mmstotalappN.values
            print("Manure spreading emission: ",livestock,prodsyst,mmstype,np.nansum(landspreading_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
            print("Manure spreading N: ",livestock,prodsyst,mmstype,np.nansum(landspreading_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
            mmsidx = mmsidx + 1
        print("Manure spreading total emission: ",livestock,prodsyst,np.nansum(landspreading_NH3emission[prodsyst_idx])/1e9," GgN.")
        print("Manure spreading total applied N: ",livestock,prodsyst,np.nansum(landspreading_comptotalN[prodsyst_idx])/1e9," GgN.")
        
    prodsyst_lvl = np.arange(3)    
    outds = xr.Dataset(
            data_vars=dict(
                NH3emiss_housing=(['level','time','lat','lon'],housing_NH3emission),
                NH3emiss_MMS=(['level','time','lat','lon'],mms_NH3emission),
                NH3emiss_manureapp=(['level','time','lat','lon'],landspreading_NH3emission),
                        ),
            coords = dict(
                level=(["level"], prodsyst_lvl),
                time=(["time"], times),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                info = "NH3 emission from housing, manure management, manure application, and grazing.\
                        Level: 0-industrial, 1-intermediate, 2-backyard.",
                units="gN per grid; m2 per grid",
            ),
        )
    outds.NH3emiss_housing.attrs["unit"] = 'gN/day'
    outds.NH3emiss_housing.attrs["long name"] = 'NH3 emission from housing'
    outds.NH3emiss_MMS.attrs["unit"] = 'gN/day'
    outds.NH3emiss_MMS.attrs["long name"] = 'NH3 emission from MMS'
    outds.NH3emiss_manureapp.attrs["unit"] = 'gN/day'
    outds.NH3emiss_manureapp.attrs["long name"] = 'NH3 emission from manure applicaton'

    outfilename = str(livestock)+'.NH3emission.'+str(yearofstudy)+'.timeseries.nc'
    outds.to_netcdf(outpath+outfilename)
    print("timeseries ncfile saved.")

    outds = xr.Dataset(
            data_vars=dict(
                housing_NH3=(['level','lat','lon'],housing_compNH3),
                housing_N=(['level','lat','lon'],housing_comptotalN),
                mms_NH3=(['level','mmslevel','lat','lon'],mms_compNH3),
                mms_N=(['level','mmslevel','lat','lon'],mms_comptotalN),
                manureapp_NH3=(['level','applevel','lat','lon'],landspreading_compNH3),
                manureapp_N=(['level','applevel','lat','lon'],landspreading_comptotalN),
                        ),
            coords = dict(
                level=(["level"], prodsyst_lvl),
                mmslevel=(["mmslevel"], mmslvl),
                applevel=(["applevel"], applvl),
                time=(["time"], times),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                info = "NH3 emission and N from housing, manure management, manure application, and grazing.\
                        Level: 0-industrial, 1-intermediate, 2-backyard.\
                        MMS level:0-liquid MMS_indoor, 1-liquid MMS_open, 2-liquid MMS_cover, 3-liquid MMS_lagoon,4-solid MMS_indoor, 5-solid MMS_open.\
                        Manure app level: 0-MMS_cover liquid, 1-MMS_indoor liquid, 2-MMS_indoor solid, 3-MMS_open liquid.",
                units="gN per grid; m2 per grid",
            ),
        )
    outds.housing_NH3.attrs["unit"] = 'gN/day'
    outds.housing_NH3.attrs["long name"] = 'NH3 emission from housing'
    outds.housing_N.attrs["unit"] = 'gN/day'
    outds.housing_N.attrs["long name"] = 'Excreted N from housing'
    outds.mms_NH3.attrs["unit"] = 'gN/day'
    outds.mms_NH3.attrs["long name"] = 'NH3 emission from MMS'
    outds.mms_N.attrs["unit"] = 'gN/day'
    outds.mms_N.attrs["long name"] = 'N from MMS'
    outds.manureapp_NH3.attrs["unit"] = 'gN/day'
    outds.manureapp_NH3.attrs["long name"] = 'NH3 emission from manure applicaton'
    outds.manureapp_N.attrs["unit"] = 'gN/day'
    outds.manureapp_N.attrs["long name"] = 'Applied manure N for land application'

    outfilename = str(livestock)+'.NH3emission.annual.'+str(yearofstudy)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

    ## Poultry
    livestock = "POULTRY"
    print("=======================================")
    print("Processing: ",livestock)

    ## production systems and levels
    prodsysts = CONFIG.CONFIG_production_system_dict[livestock]
    prodsyst_lvl = len(prodsysts)
    ## three types of NH3 emissions: housing, MMS, and land spreading
    housing_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    mms_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    landspreading_NH3emission = np.zeros((prodsyst_lvl,Days,nlat,nlon))
    housing_compNH3 = np.zeros((prodsyst_lvl,nlat,nlon))
    housing_comptotalN = np.zeros((prodsyst_lvl,nlat,nlon))
    mms_compNH3 = np.zeros((prodsyst_lvl,6,nlat,nlon))
    mms_comptotalN = np.zeros((prodsyst_lvl,6,nlat,nlon))
    landspreading_compNH3 = np.zeros((prodsyst_lvl,4,nlat,nlon))
    landspreading_comptotalN = np.zeros((prodsyst_lvl,4,nlat,nlon))

    housing_types = ["poultry_house","poultry_house","poultry_house"]
    for prodsyst_idx in np.arange(prodsyst_lvl):
        prodsyst = prodsysts[prodsyst_idx]
        print("=======================================")
        print("Production system: ",prodsyst)
        
        ## HOUSING NH3 emissions
        housetype = housing_types[prodsyst_idx]
        if prodsyst != "backyard":
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                    str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission (no litter): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (no litter): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'_litter.'+\
                    str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values + housing_compNH3[prodsyst_idx][:]
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret + housing_comptotalN[prodsyst_idx][:]
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values + housing_NH3emission[prodsyst_idx][:]
            
            print("Housing emission (with litter): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (with litter): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            print("Housing emission (total poultry house): ",livestock,prodsyst,np.nansum(housing_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Housing total N (total poultry house): ",livestock,prodsyst,np.nansum(housing_comptotalN[prodsyst_idx])/1e9," GgN.")

        else:
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                        str(yearofstudy)+'.nc'
            housingds = xr.open_dataset(housingfile)
            housing_Nexcret = housingds.Nexcret.values
            housing_emission = housingds.NH3emiss
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret 
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission: ",livestock,prodsyst,np.nansum(housing_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
        
        ## MMS and emissions
        mmsidx = 0
        print("[MMS]")
        if prodsyst_idx < 2:
            for mmsphase in ["liquid","solid"]:
                if mmsphase == "liquid":
                    mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
                elif mmsphase == "solid":
                    mmscats = ['MMS_indoor',"MMS_open"]
                for mmscat in mmscats:
                    try:
                        print("MMS files: ",mmsidx,mmsphase,mmscat)
                        mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                                str(mmsphase)+'.'+str(yearofstudy)+'.nc'
                        mmsds = xr.open_dataset(mmsfile)
                        mmsemiss = mmsds.NH3emiss
                        mms_NH3emission[prodsyst_idx][:] = mmsemiss.values + mms_NH3emission[prodsyst_idx][:]
                        mms_compNH3[prodsyst_idx][mmsidx][:] = mmsemiss.sum(dim="time").values
                        mms_comptotalN[prodsyst_idx][mmsidx][:] = mmsds.NtotalMMS.values
                        print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
                        print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
                        if prodsyst != "backyard":
                            if mmsphase == "solid":
                                if mmscat == "MMS_open":
                                    ## include insit MMS
                                    mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                                            str(mmsphase)+'.insitu.'+str(yearofstudy)+'.nc'
                                    mmsds = xr.open_dataset(mmsfile)
                                    mmsemiss = mmsds.NH3emiss
                                    mms_NH3emission[prodsyst_idx][:] = np.nan_to_num(mmsemiss.values) + np.nan_to_num(mms_NH3emission[prodsyst_idx][:])
                                    mms_compNH3[prodsyst_idx][mmsidx][:] = mmsemiss.sum(dim="time").values + mms_compNH3[prodsyst_idx][mmsidx][:]
                                    mms_comptotalN[prodsyst_idx][mmsidx][:] = np.nan_to_num(mmsds.NtotalMMS.values) + np.nan_to_num(mms_comptotalN[prodsyst_idx][mmsidx][:])
                                    print("MMS emission (insitu MMS): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mmsemiss.sum(dim="time").values)/1e9," GgN.")
                                    print("MMS N (insitu MMS): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mmsds.NtotalMMS.values)/1e9," GgN.")
                                    print("MMS emission (updated): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
                                    print("MMS N (updated): ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
                        mmsidx = mmsidx + 1
                    except:
                        pass
                        print(livestock,"no "+mmsphase+' '+mmscat+' managemet.')
                        mmsidx = mmsidx + 1
            print("MMS total emission: ",livestock,prodsyst,np.nansum(mms_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret[prodsyst_idx])/1e9," GgN.")
        else:
            mmsphase = "solid.insitu"
            mmscat = "MMS_open"
            mmsidx = 5
            print("MMS files: ",mmsidx,mmsphase,mmscat)
            mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                    str(mmsphase)+'.'+str(yearofstudy)+'.nc'
            mmsds = xr.open_dataset(mmsfile)
            mmsemiss = mmsds.NH3emiss
            mms_NH3emission[prodsyst_idx][:] = mmsemiss.values + mms_NH3emission[prodsyst_idx][:]
            mms_compNH3[prodsyst_idx][mmsidx][:] = mmsemiss.sum(dim="time").values
            mms_comptotalN[prodsyst_idx][mmsidx][:] = mmsds.NtotalMMS.values
            print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
            print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")

        
        ## Land spreading emissions
        if prodsyst_idx < 2:
            mmsidx = 0
            print("[Manure land spreading]")
            for mmstype in ["MMS_cover.liquid","MMS_indoor.liquid","MMS_indoor.solid","MMS_open.liquid"]:
                print("MMS app, ",mmstype)
                manureappfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+'.'+str(mmstype)+\
                            '.surf.'+str(yearofstudy)+'.nc'
                manureappds = xr.open_dataset(manureappfile)
                manureappemiss = manureappds.NH3emiss
                landspreading_NH3emission[prodsyst_idx][:] = manureappemiss.values + landspreading_NH3emission[prodsyst_idx][:]
                landspreading_compNH3[prodsyst_idx][mmsidx][:] = manureappemiss.sum(dim="time").values

                manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                            str(yearofstudy)+'.manureapp.nc'
                manureappNds = xr.open_dataset(manureappNfile)
                mmstotalappN = manureappNds.spring_TAN+manureappNds.spring_UAN+\
                            manureappNds.spring_availN+manureappNds.spring_resistN+manureappNds.spring_unavailN+\
                            manureappNds.winter_TAN+manureappNds.winter_UAN+\
                            manureappNds.winter_availN+manureappNds.winter_resistN+manureappNds.winter_unavailN
                landspreading_comptotalN[prodsyst_idx][mmsidx][:] = mmstotalappN.values
                print("Manure spreading emission: ",livestock,prodsyst,mmstype,np.nansum(landspreading_compNH3[prodsyst_idx][mmsidx])/1e9," GgN.")
                print("Manure spreading N: ",livestock,prodsyst,mmstype,np.nansum(landspreading_comptotalN[prodsyst_idx][mmsidx])/1e9," GgN.")
                mmsidx = mmsidx + 1
            print("Manure spreading total emission: ",livestock,prodsyst,np.nansum(landspreading_NH3emission[prodsyst_idx])/1e9," GgN.")
            print("Manure spreading total applied N: ",livestock,prodsyst,np.nansum(landspreading_comptotalN[prodsyst_idx])/1e9," GgN.")
            
    prodsyst_lvl = np.arange(3)    
    outds = xr.Dataset(
            data_vars=dict(
                NH3emiss_housing=(['level','time','lat','lon'],housing_NH3emission),
                NH3emiss_MMS=(['level','time','lat','lon'],mms_NH3emission),
                NH3emiss_manureapp=(['level','time','lat','lon'],landspreading_NH3emission),
                        ),
            coords = dict(
                level=(["level"], prodsyst_lvl),
                time=(["time"], times),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                info = "NH3 emission from housing, manure management, manure application, and grazing.\
                        Level: 0-broiler, 1-layer, 2-backyard.",
                units="gN per grid; m2 per grid",
            ),
        )

    outds.NH3emiss_housing.attrs["unit"] = 'gN/day'
    outds.NH3emiss_housing.attrs["long name"] = 'NH3 emission from housing'
    outds.NH3emiss_MMS.attrs["unit"] = 'gN/day'
    outds.NH3emiss_MMS.attrs["long name"] = 'NH3 emission from MMS'
    outds.NH3emiss_manureapp.attrs["unit"] = 'gN/day'
    outds.NH3emiss_manureapp.attrs["long name"] = 'NH3 emission from manure applicaton'
    outfilename = str(livestock)+'.NH3emission.'+str(yearofstudy)+'.timeseries.nc'
    outds.to_netcdf(outpath+outfilename)
    print("timeseries ncfile saved.")

    outds = xr.Dataset(
            data_vars=dict(
                housing_NH3=(['level','lat','lon'],housing_compNH3),
                housing_N=(['level','lat','lon'],housing_comptotalN),
                mms_NH3=(['level','mmslevel','lat','lon'],mms_compNH3),
                mms_N=(['level','mmslevel','lat','lon'],mms_comptotalN),
                manureapp_NH3=(['level','applevel','lat','lon'],landspreading_compNH3),
                manureapp_N=(['level','applevel','lat','lon'],landspreading_comptotalN),
                        ),
            coords = dict(
                level=(["level"], prodsyst_lvl),
                mmslevel=(["mmslevel"], mmslvl),
                applevel=(["applevel"], applvl),
                time=(["time"], times),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                info = "NH3 emission and N from housing, manure management, manure application, and grazing.\
                        Level: 0-broiler, 1-layer, 2-backyard.\
                        MMS level:0-liquid MMS_indoor, 1-liquid MMS_open, 2-liquid MMS_cover, 3-liquid MMS_lagoon,4-solid MMS_indoor, 5-solid MMS_open.\
                        Manure app level: 0-MMS_cover liquid, 1-MMS_indoor liquid, 2-MMS_indoor solid, 3-MMS_open liquid.",
                units="gN per grid; m2 per grid",
            ),
        )

    outds.housing_NH3.attrs["unit"] = 'gN/day'
    outds.housing_NH3.attrs["long name"] = 'NH3 emission from housing'
    outds.housing_N.attrs["unit"] = 'gN/day'
    outds.housing_N.attrs["long name"] = 'Excreted N from housing'
    outds.mms_NH3.attrs["unit"] = 'gN/day'
    outds.mms_NH3.attrs["long name"] = 'NH3 emission from MMS'
    outds.mms_N.attrs["unit"] = 'gN/day'
    outds.mms_N.attrs["long name"] = 'N from MMS'
    outds.manureapp_NH3.attrs["unit"] = 'gN/day'
    outds.manureapp_NH3.attrs["long name"] = 'NH3 emission from manure applicaton'
    outds.manureapp_N.attrs["unit"] = 'gN/day'
    outds.manureapp_N.attrs["long name"] = 'Applied manure N for land application'
    outfilename = str(livestock)+'.NH3emission.annual.'+str(yearofstudy)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

    for livestock in ruminants_list:
        
        infile_path = CONFIG.infile_path
        animal_data_path = 'animal_data/'
        animal_file_name = CONFIG.CONFIG_animal_file_dict[livestock]
        
        print("=======================================")
        print("Processing: ",livestock)
        ## production systems and levels
        prodsysts = CONFIG.CONFIG_production_system_dict[livestock]
        prodsyst_lvl = len(prodsysts)

        ## three/four types of NH3 emissions: housing, MMS, land spreading, grazing
        housing_NH3emission = np.zeros((Days,nlat,nlon))
        mms_NH3emission = np.zeros((Days,nlat,nlon))
        landspreading_NH3emission = np.zeros((Days,nlat,nlon))
        allyeargrazing_NH3emission = np.zeros((Days,nlat,nlon))
        seasonalgrazing_NH3emission = np.zeros((Days,nlat,nlon))

        ## Mixed production system: NH3 from housing, MMS, land spreading and seasonal grazing
        prodsyst_idx = CONFIG.CONFIG_production_system_dict[livestock].index("mixed")
        prodsyst = prodsysts[prodsyst_idx]
        print("Production system: ",prodsyst)

        ## HOUSING NH3 emissions
        print("[Housing]")
        housetype = "barn"
        housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                    str(yearofstudy)+'.nc'
        housingds = xr.open_dataset(housingfile)
        animalds = xr.open_dataset(infile_path+animal_data_path+animal_file_name)
        livestockNfactords = xr.open_dataset(infile_path+animal_data_path+'livestockN_factor_updated.nc')
        yearly_factor = livestockNfactords.yearly_factor.sel(year=yearofstudy)
        housing_Nexcret = animalds['Excreted_N'][prodsyst_idx].values*1e3*yearly_factor - housingds.grazing_N.sum(dim="time").values
        housing_emission = housingds.NH3emiss
        housing_NH3emission[:] = housing_emission.values
        print("Housing emission: ",livestock,prodsyst,np.nansum(housing_NH3emission)/1e9," GgN.")
        print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
        
        ## MMS and emissions
        mms_compNH3 = np.zeros((6,nlat,nlon))
        mms_comptotalN = np.zeros((6,nlat,nlon))
        mmsidx = 0
        print("[MMS]")
        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats:
                try:
                    print("MMS files: ",mmsidx,mmsphase,mmscat)
                    mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                            str(mmsphase)+'.'+str(yearofstudy)+'.nc'
                    mmsds = xr.open_dataset(mmsfile)
                    mmsemiss = mmsds.NH3emiss
                    mms_NH3emission[:] = mmsemiss.values + mms_NH3emission[:]
                    mms_compNH3[mmsidx][:] = mmsemiss.sum(dim="time").values
                    mms_comptotalN[mmsidx][:] = mmsds.NtotalMMS.values
                    print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[mmsidx])/1e9," GgN.")
                    print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[mmsidx])/1e9," GgN.")
                    mmsidx = mmsidx + 1
                except:
                    pass
                    print(livestock,"no "+mmsphase+' '+mmscat+' managemet.')
                    mmsidx = mmsidx + 1
        print("MMS total emission: ",livestock,prodsyst,np.nansum(mms_NH3emission)/1e9," GgN.")
        print("MMS total N: ",livestock,prodsyst,np.nansum(mms_comptotalN)/1e9," GgN.")

        ## Land spreading emissions
        landspreading_compNH3 = np.zeros((4,nlat,nlon))
        landspreading_comptotalN = np.zeros((4,nlat,nlon))
        mmsidx = 0
        print("[Manure land spreading]")
        for mmstype in ["MMS_cover.liquid","MMS_indoor.liquid","MMS_indoor.solid","MMS_open.liquid"]:
            print("MMS app, ",mmstype)
            manureappfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+'.'+str(mmstype)+\
                        '.surf.'+str(yearofstudy)+'.nc'
            manureappds = xr.open_dataset(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[:] = manureappemiss.values + landspreading_NH3emission[:]
            landspreading_compNH3[mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(yearofstudy)+'.manureapp.nc'
            manureappNds = xr.open_dataset(manureappNfile)
            mmstotalappN = manureappNds.spring_TAN+manureappNds.spring_UAN+\
                        manureappNds.spring_availN+manureappNds.spring_resistN+manureappNds.spring_unavailN+\
                        manureappNds.winter_TAN+manureappNds.winter_UAN+\
                        manureappNds.winter_availN+manureappNds.winter_resistN+manureappNds.winter_unavailN
            landspreading_comptotalN[mmsidx][:] = mmstotalappN.values
            print("Manure spreading emission: ",livestock,prodsyst,mmstype,np.nansum(landspreading_compNH3[mmsidx])/1e9," GgN.")
            print("Manure spreading N: ",livestock,prodsyst,mmstype,np.nansum(landspreading_comptotalN[mmsidx])/1e9," GgN.")
            mmsidx = mmsidx + 1
        print("Manure spreading total emission: ",livestock,prodsyst,np.nansum(landspreading_NH3emission)/1e9," GgN.")
        print("Manure spreading total applied N: ",livestock,prodsyst,np.nansum(landspreading_comptotalN)/1e9," GgN.")
        
        ## Grazing emission
        allyeargrazing_compNH3 = np.zeros((nlat,nlon))
        allyeargrazinging_comptotalN = np.zeros((nlat,nlon))
        seasonalgrazing_compNH3 = np.zeros((nlat,nlon))
        seasonalgrazing_comptotalN = np.zeros((nlat,nlon))

        # grazing_compNH3 = np.zeros((2,nlat,nlon))
        # grazinging_comptotalN = np.zeros((2,nlat,nlon))
        print("[Grazing]")
        for prodsyst_idx in np.arange(2):
            prodsyst = prodsysts[prodsyst_idx]
            print("Production system: ",prodsyst)
            grazingfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+\
                        '.grazing.'+str(yearofstudy)+'.nc'
            grazingds = xr.open_dataset(grazingfile)
            grazingemiss = grazingds.NH3emiss_urinepatch + \
                            grazingds.NH3emiss_dungpat + \
                            grazingds.NH3emiss_mixed
            if prodsyst == "grassland":
                allyeargrazing_NH3emission[:] = grazingemiss.values
                allyeargrazing_compNH3[:] = grazingemiss.sum(dim="time")
                grasslandds = xr.open_dataset(infile_path+animal_data_path+animal_file_name)
                livestockNfactords = xr.open_dataset(infile_path+animal_data_path+'livestockN_factor_updated.nc')
                yearly_factor = livestockNfactords.yearly_factor.sel(year=yearofstudy)
                allyeargrazinging_comptotalN[:] = grasslandds['Excreted_N'][prodsyst_idx].values*1e3*yearly_factor
                print("Year-round grazing emission: ",livestock,prodsyst,np.nansum(allyeargrazing_NH3emission)/1e9," GgN.")
                print("Year-round grazing total N: ",livestock,prodsyst,np.nansum(allyeargrazinging_comptotalN)/1e9," GgN.")
            elif prodsyst == "mixed":
                seasonalgrazing_NH3emission[:] = grazingemiss.values
                seasonalgrazing_compNH3[:] = grazingemiss.sum(dim="time").values
                seasonalgrazing_comptotalN[:] = housingds.grazing_N.sum(dim="time").values
                print("Seasonal grazing emission: ",livestock,prodsyst,np.nansum(seasonalgrazing_NH3emission)/1e9," GgN.")
                print("Seasonal grazing total N: ",livestock,prodsyst,np.nansum(seasonalgrazing_comptotalN)/1e9," GgN.")
        
        
        outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss_housing=(['time','lat','lon'],housing_NH3emission),
                    NH3emiss_MMS=(['time','lat','lon'],mms_NH3emission),
                    NH3emiss_manureapp=(['time','lat','lon'],landspreading_NH3emission),
                    NH3emiss_seasonalgrazing=(['time','lat','lon'],seasonalgrazing_NH3emission),
                    NH3emiss_allyeargrazing=(['time','lat','lon'],allyeargrazing_NH3emission),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission from housing, manure management, manure application, and grazing",
                    units="gN per grid; m2 per grid",
                ),
            )
        
        outds.NH3emiss_housing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_housing.attrs["long name"] = 'NH3 emission from housing'
        outds.NH3emiss_MMS.attrs["unit"] = 'gN/day'
        outds.NH3emiss_MMS.attrs["long name"] = 'NH3 emission from MMS'
        outds.NH3emiss_manureapp.attrs["unit"] = 'gN/day'
        outds.NH3emiss_manureapp.attrs["long name"] = 'NH3 emission from manure applicaton'
        outds.NH3emiss_seasonalgrazing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_seasonalgrazing.attrs["long name"] = 'NH3 emission from seasonal grazing'
        outds.NH3emiss_allyeargrazing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_allyeargrazing.attrs["long name"] = 'NH3 emission from year-round grazing'

        outfilename = str(livestock)+'.NH3emission.'+str(yearofstudy)+'.timeseries.nc'
        outds.to_netcdf(outpath+outfilename)
        print("timeseries ncfile saved.")

        outds = xr.Dataset(
                data_vars=dict(
                    housing_NH3=(['lat','lon'],housing_emission.sum(dim="time").values),
                    housing_N=(['lat','lon'],housing_Nexcret.values),
                    mms_NH3=(['mmslevel','lat','lon'],mms_compNH3),
                    mms_N=(['mmslevel','lat','lon'],mms_comptotalN),
                    manureapp_NH3=(['applevel','lat','lon'],landspreading_compNH3),
                    manureapp_N=(['applevel','lat','lon'],landspreading_comptotalN),
                    seasonalgrazing_NH3=(['lat','lon'],seasonalgrazing_compNH3),
                    seasonalgrazing_N=(['lat','lon'],seasonalgrazing_comptotalN),
                    allyeargrazing_NH3=(['lat','lon'],allyeargrazing_compNH3),
                    allyeargrazing_N=(['lat','lon'],allyeargrazinging_comptotalN),
                            ),
                coords = dict(
                    mmslevel=(["mmslevel"], mmslvl),
                    applevel=(["applevel"], applvl),
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission and N from housing, manure management, manure application, and grazing. MMS level:\
                            0-liquid MMS_indoor, 1-liquid MMS_open, 2-liquid MMS_cover, 3-liquid MMS_lagoon,\
                            4-solid MMS_indoor, 5-solid MMS_open.  Manure app level:\
                            0-MMS_cover liquid, 1-MMS_indoor liquid, 2-MMS_indoor solid, 3-MMS_open liquid.",
                    units="gN per grid; m2 per grid",
                ),
            )

        outds.housing_NH3.attrs["unit"] = 'gN/day'
        outds.housing_NH3.attrs["long name"] = 'NH3 emission from housing'
        outds.housing_N.attrs["unit"] = 'gN/day'
        outds.housing_N.attrs["long name"] = 'Excreted N from housing'
        outds.mms_NH3.attrs["unit"] = 'gN/day'
        outds.mms_NH3.attrs["long name"] = 'NH3 emission from MMS'
        outds.mms_N.attrs["unit"] = 'gN/day'
        outds.mms_N.attrs["long name"] = 'N from MMS'
        outds.manureapp_NH3.attrs["unit"] = 'gN/day'
        outds.manureapp_NH3.attrs["long name"] = 'NH3 emission from manure applicaton'
        outds.manureapp_N.attrs["unit"] = 'gN/day'
        outds.manureapp_N.attrs["long name"] = 'Applied manure N for land application'
        outds.seasonalgrazing_NH3.attrs["unit"] = 'gN/day'
        outds.seasonalgrazing_NH3.attrs["long name"] = 'NH3 emission from seasonal grazing'
        outds.seasonalgrazing_N.attrs["unit"] = 'gN/day'
        outds.seasonalgrazing_N.attrs["long name"] = 'Excreted N during seasonal grazing'
        outds.allyeargrazing_NH3.attrs["unit"] = 'gN/day'
        outds.allyeargrazing_NH3.attrs["long name"] = 'NH3 emission from year-round grazing'
        outds.allyeargrazing_N.attrs["unit"] = 'gN/day'
        outds.allyeargrazing_N.attrs["long name"] = 'Excreted N during year-round grazing'

        outfilename = str(livestock)+'.NH3emission.annual.'+str(yearofstudy)+'.nc'
        outds.to_netcdf(outpath+outfilename)
        print("ncfile saved.")

    for livestock in ['BUFFALO_BEEF','BUFFALO_DAIRY']:
        infile_path = CONFIG.infile_path
        animal_data_path = 'animal_data/'
        animal_file_name = CONFIG.CONFIG_animal_file_dict[livestock]

        print("=======================================")
        print("Processing: ",livestock)

        ## production systems and levels
        prodsysts = CONFIG.CONFIG_production_system_dict[livestock]
        prodsyst_lvl = len(prodsysts)

        ## three/four types of NH3 emissions: housing, MMS, land spreading, grazing
        housing_NH3emission = np.zeros((Days,nlat,nlon))
        mms_NH3emission = np.zeros((Days,nlat,nlon))
        landspreading_NH3emission = np.zeros((Days,nlat,nlon))
        allyeargrazing_NH3emission = np.zeros((Days,nlat,nlon))
        seasonalgrazing_NH3emission = np.zeros((Days,nlat,nlon))


        ## Mixed production system: NH3 from housing, MMS, land spreading and seasonal grazing
        prodsyst_idx = CONFIG.CONFIG_production_system_dict[livestock].index("mixed")
        prodsyst = prodsysts[prodsyst_idx]
        print("Production system: ",prodsyst)

        ## HOUSING NH3 emissions
        print("[Housing]")
        housetype = "barn"
        
        housingfile = ncfilepath+'BUFFALO/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                    str(yearofstudy)+'.nc'

        housingds = xr.open_dataset(housingfile)
        animalds = xr.open_dataset(infile_path+animal_data_path+animal_file_name)
        livestockNfactords = xr.open_dataset(infile_path+animal_data_path+'livestockN_factor_updated.nc')
        yearly_factor = livestockNfactords.yearly_factor.sel(year=yearofstudy)
        housing_Nexcret = animalds['Excreted_N'][prodsyst_idx].values*1e3*yearly_factor - housingds.grazing_N.sum(dim="time").values
        housing_emission = housingds.NH3emiss
        housing_NH3emission[:] = housing_emission.values
        print("Housing emission: ",livestock,prodsyst,np.nansum(housing_NH3emission)/1e9," GgN.")
        print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
        
        ## MMS and emissions
        mms_compNH3 = np.zeros((6,nlat,nlon))
        mms_comptotalN = np.zeros((6,nlat,nlon))
        mmsidx = 0
        print("[MMS]")
        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats:
                try:
                    print("MMS files: ",mmsidx,mmsphase,mmscat)
                    mmsfile = ncfilepath+'BUFFALO/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                            str(mmsphase)+'.'+str(yearofstudy)+'.nc'
                    mmsds = xr.open_dataset(mmsfile)
                    mmsemiss = mmsds.NH3emiss
                    mms_NH3emission[:] = mmsemiss.values + mms_NH3emission[:]
                    mms_compNH3[mmsidx][:] = mmsemiss.sum(dim="time").values
                    mms_comptotalN[mmsidx][:] = mmsds.NtotalMMS.values
                    print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[mmsidx])/1e9," GgN.")
                    print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[mmsidx])/1e9," GgN.")
                    mmsidx = mmsidx + 1
                except:
                    pass
                    print(livestock,"no "+mmsphase+' '+mmscat+' managemet.')
                    mmsidx = mmsidx + 1
        print("MMS total emission: ",livestock,prodsyst,np.nansum(mms_NH3emission)/1e9," GgN.")
        print("MMS total N: ",livestock,prodsyst,np.nansum(mms_comptotalN)/1e9," GgN.")

        ## Land spreading emissions
        landspreading_compNH3 = np.zeros((4,nlat,nlon))
        landspreading_comptotalN = np.zeros((4,nlat,nlon))
        mmsidx = 0
        print("[Manure land spreading]")
        for mmstype in ["MMS_indoor.liquid","MMS_indoor.solid","MMS_open.liquid"]:
            print("MMS app, ",mmstype)
            manureappfile = ncfilepath+'BUFFALO/base.'+str(prodsyst)+'.'+str(livestock)+'.'+str(mmstype)+\
                        '.surf.'+str(yearofstudy)+'.nc'
            manureappds = xr.open_dataset(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[:] = manureappemiss.values + landspreading_NH3emission[:]
            landspreading_compNH3[mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+'BUFFALO/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(yearofstudy)+'.manureapp.nc'
            manureappNds = xr.open_dataset(manureappNfile)
            mmstotalappN = manureappNds.spring_TAN+manureappNds.spring_UAN+\
                        manureappNds.spring_availN+manureappNds.spring_resistN+manureappNds.spring_unavailN+\
                        manureappNds.winter_TAN+manureappNds.winter_UAN+\
                        manureappNds.winter_availN+manureappNds.winter_resistN+manureappNds.winter_unavailN
            landspreading_comptotalN[mmsidx][:] = mmstotalappN.values
            print("Manure spreading emission: ",livestock,prodsyst,mmstype,np.nansum(landspreading_compNH3[mmsidx])/1e9," GgN.")
            print("Manure spreading N: ",livestock,prodsyst,mmstype,np.nansum(landspreading_comptotalN[mmsidx])/1e9," GgN.")
            mmsidx = mmsidx + 1
        print("Manure spreading total emission: ",livestock,prodsyst,np.nansum(landspreading_NH3emission)/1e9," GgN.")
        print("Manure spreading total applied N: ",livestock,prodsyst,np.nansum(landspreading_comptotalN)/1e9," GgN.")
        
        ## Grazing emission
        allyeargrazing_compNH3 = np.zeros((nlat,nlon))
        allyeargrazinging_comptotalN = np.zeros((nlat,nlon))
        seasonalgrazing_compNH3 = np.zeros((nlat,nlon))
        seasonalgrazing_comptotalN = np.zeros((nlat,nlon))

        # grazing_compNH3 = np.zeros((2,nlat,nlon))
        # grazinging_comptotalN = np.zeros((2,nlat,nlon))
        print("[Grazing]")
        for prodsyst_idx in np.arange(2):
            prodsyst = prodsysts[prodsyst_idx]
            print("Production system: ",prodsyst)
            grazingfile = ncfilepath+'BUFFALO/base.'+str(prodsyst)+'.'+str(livestock)+\
                        '.grazing.'+str(yearofstudy)+'.nc'
            grazingds = xr.open_dataset(grazingfile)
            grazingemiss = grazingds.NH3emiss_urinepatch + \
                            grazingds.NH3emiss_dungpat + \
                            grazingds.NH3emiss_mixed
            if prodsyst == "grassland":
                allyeargrazing_NH3emission[:] = grazingemiss.values
                allyeargrazing_compNH3[:] = grazingemiss.sum(dim="time")
                grasslandds = xr.open_dataset(infile_path+animal_data_path+animal_file_name)
                livestockNfactords = xr.open_dataset(infile_path+animal_data_path+'livestockN_factor_updated.nc')
                yearly_factor = livestockNfactords.yearly_factor.sel(year=yearofstudy)
                allyeargrazinging_comptotalN[:] = grasslandds['Excreted_N'][prodsyst_idx].values*1e3*yearly_factor
                print("Year-round grazing emission: ",livestock,prodsyst,np.nansum(allyeargrazing_NH3emission)/1e9," GgN.")
                print("Year-round grazing total N: ",livestock,prodsyst,np.nansum(allyeargrazinging_comptotalN)/1e9," GgN.")
            elif prodsyst == "mixed":
                seasonalgrazing_NH3emission[:] = grazingemiss.values
                seasonalgrazing_compNH3[:] = grazingemiss.sum(dim="time").values
                seasonalgrazing_comptotalN[:] = housingds.grazing_N.sum(dim="time").values
                print("Seasonal grazing emission: ",livestock,prodsyst,np.nansum(seasonalgrazing_NH3emission)/1e9," GgN.")
                print("Seasonal grazing total N: ",livestock,prodsyst,np.nansum(seasonalgrazing_comptotalN)/1e9," GgN.")
        
        
        outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss_housing=(['time','lat','lon'],housing_NH3emission),
                    NH3emiss_MMS=(['time','lat','lon'],mms_NH3emission),
                    NH3emiss_manureapp=(['time','lat','lon'],landspreading_NH3emission),
                    NH3emiss_seasonalgrazing=(['time','lat','lon'],seasonalgrazing_NH3emission),
                    NH3emiss_allyeargrazing=(['time','lat','lon'],allyeargrazing_NH3emission),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission from housing, manure management, manure application, and grazing",
                    units="gN per grid; m2 per grid",
                ),
            )
        
        outds.NH3emiss_housing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_housing.attrs["long name"] = 'NH3 emission from housing'
        outds.NH3emiss_MMS.attrs["unit"] = 'gN/day'
        outds.NH3emiss_MMS.attrs["long name"] = 'NH3 emission from MMS'
        outds.NH3emiss_manureapp.attrs["unit"] = 'gN/day'
        outds.NH3emiss_manureapp.attrs["long name"] = 'NH3 emission from manure applicaton'
        outds.NH3emiss_seasonalgrazing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_seasonalgrazing.attrs["long name"] = 'NH3 emission from seasonal grazing'
        outds.NH3emiss_allyeargrazing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_allyeargrazing.attrs["long name"] = 'NH3 emission from year-round grazing'

        outfilename = str(livestock)+'.NH3emission.'+str(yearofstudy)+'.timeseries.nc'
        outds.to_netcdf(outpath+outfilename)
        print("timeseries ncfile saved.")

        outds = xr.Dataset(
                data_vars=dict(
                    housing_NH3=(['lat','lon'],housing_emission.sum(dim="time").values),
                    housing_N=(['lat','lon'],housing_Nexcret.values),
                    mms_NH3=(['mmslevel','lat','lon'],mms_compNH3),
                    mms_N=(['mmslevel','lat','lon'],mms_comptotalN),
                    manureapp_NH3=(['applevel','lat','lon'],landspreading_compNH3),
                    manureapp_N=(['applevel','lat','lon'],landspreading_comptotalN),
                    seasonalgrazing_NH3=(['lat','lon'],seasonalgrazing_compNH3),
                    seasonalgrazing_N=(['lat','lon'],seasonalgrazing_comptotalN),
                    allyeargrazing_NH3=(['lat','lon'],allyeargrazing_compNH3),
                    allyeargrazing_N=(['lat','lon'],allyeargrazinging_comptotalN),
                            ),
                coords = dict(
                    mmslevel=(["mmslevel"], mmslvl),
                    applevel=(["applevel"], applvl),
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission and N from housing, manure management, manure application, and grazing. MMS level:\
                            0-liquid MMS_indoor, 1-liquid MMS_open, 2-liquid MMS_cover, 3-liquid MMS_lagoon,\
                            4-solid MMS_indoor, 5-solid MMS_open.  Manure app level:\
                            0-MMS_cover liquid, 1-MMS_indoor liquid, 2-MMS_indoor solid, 3-MMS_open liquid.",
                    units="gN per grid; m2 per grid",
                ),
            )

        outds.housing_NH3.attrs["unit"] = 'gN/day'
        outds.housing_NH3.attrs["long name"] = 'NH3 emission from housing'
        outds.housing_N.attrs["unit"] = 'gN/day'
        outds.housing_N.attrs["long name"] = 'Excreted N from housing'
        outds.mms_NH3.attrs["unit"] = 'gN/day'
        outds.mms_NH3.attrs["long name"] = 'NH3 emission from MMS'
        outds.mms_N.attrs["unit"] = 'gN/day'
        outds.mms_N.attrs["long name"] = 'N from MMS'
        outds.manureapp_NH3.attrs["unit"] = 'gN/day'
        outds.manureapp_NH3.attrs["long name"] = 'NH3 emission from manure applicaton'
        outds.manureapp_N.attrs["unit"] = 'gN/day'
        outds.manureapp_N.attrs["long name"] = 'Applied manure N for land application'
        outds.seasonalgrazing_NH3.attrs["unit"] = 'gN/day'
        outds.seasonalgrazing_NH3.attrs["long name"] = 'NH3 emission from seasonal grazing'
        outds.seasonalgrazing_N.attrs["unit"] = 'gN/day'
        outds.seasonalgrazing_N.attrs["long name"] = 'Excreted N during seasonal grazing'
        outds.allyeargrazing_NH3.attrs["unit"] = 'gN/day'
        outds.allyeargrazing_NH3.attrs["long name"] = 'NH3 emission from year-round grazing'
        outds.allyeargrazing_N.attrs["unit"] = 'gN/day'
        outds.allyeargrazing_N.attrs["long name"] = 'Excreted N during year-round grazing'

        outfilename = str(livestock)+'.NH3emission.annual.'+str(yearofstudy)+'.nc'
        outds.to_netcdf(outpath+outfilename)
        print("ncfile saved.")


    for livestock in ["FEEDLOT_CATTLE"]:
        print("=======================================")
        print("Processing: ",livestock)

        ## production systems and levels
        prodsysts = CONFIG.CONFIG_production_system_dict[livestock]
        prodsyst_lvl = len(prodsysts)

        ## three/four types of NH3 emissions: housing, MMS, land spreading, grazing
        housing_NH3emission = np.zeros((Days,nlat,nlon))
        mms_NH3emission = np.zeros((Days,nlat,nlon))
        landspreading_NH3emission = np.zeros((Days,nlat,nlon))
        allyeargrazing_NH3emission = np.zeros((Days,nlat,nlon))
        seasonalgrazing_NH3emission = np.zeros((Days,nlat,nlon))


        ## Mixed production system: NH3 from housing, MMS, land spreading and seasonal grazing
        prodsyst_idx = CONFIG.CONFIG_production_system_dict[livestock].index("feedlot")
        prodsyst = prodsysts[prodsyst_idx]
        print("Production system: ",prodsyst)

        ## HOUSING NH3 emissions
        print("[Housing]")
        housetype = "barn"
        housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'.'+\
                    str(yearofstudy)+'.nc'
        housingds = xr.open_dataset(housingfile)
        housing_Nexcret = housingds.Nexcret.values
        housing_emission = housingds.NH3emiss
        housing_NH3emission[:] = housing_emission.values
        print("Housing emission: ",livestock,prodsyst,np.nansum(housing_NH3emission)/1e9," GgN.")
        print("Housing total N: ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
        
        ## MMS and emissions
        mms_compNH3 = np.zeros((6,nlat,nlon))
        mms_comptotalN = np.zeros((6,nlat,nlon))
        mmsidx = 0
        print("[MMS]")
        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats:
                try:
                    print("MMS files: ",mmsidx,mmsphase,mmscat)
                    mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                            str(mmsphase)+'.'+str(yearofstudy)+'.nc'
                    mmsds = xr.open_dataset(mmsfile)
                    mmsemiss = mmsds.NH3emiss
                    mms_NH3emission[:] = mmsemiss.values + mms_NH3emission[:]
                    mms_compNH3[mmsidx][:] = mmsemiss.sum(dim="time").values
                    mms_comptotalN[mmsidx][:] = mmsds.NtotalMMS.values
                    print("MMS emission: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_compNH3[mmsidx])/1e9," GgN.")
                    print("MMS N: ",livestock,prodsyst,mmscat,mmsphase,np.nansum(mms_comptotalN[mmsidx])/1e9," GgN.")
                    mmsidx = mmsidx + 1
                except:
                    pass
                    print(livestock,"no "+mmsphase+' '+mmscat+' managemet.')
                    mmsidx = mmsidx + 1
        print("MMS total emission: ",livestock,prodsyst,np.nansum(mms_NH3emission)/1e9," GgN.")
        print("MMS total N: ",livestock,prodsyst,np.nansum(mms_comptotalN)/1e9," GgN.")

        ## Land spreading emissions
        landspreading_compNH3 = np.zeros((4,nlat,nlon))
        landspreading_comptotalN = np.zeros((4,nlat,nlon))
        mmsidx = 0
        print("[Manure land spreading]")
        for mmstype in ["MMS_cover.liquid","MMS_indoor.liquid","MMS_indoor.solid","MMS_open.liquid"]:
            print("MMS app, ",mmstype)
            manureappfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+'.'+str(mmstype)+\
                        '.surf.'+str(yearofstudy)+'.nc'
            manureappds = xr.open_dataset(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[:] = manureappemiss.values + landspreading_NH3emission[:]
            landspreading_compNH3[mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(yearofstudy)+'.manureapp.nc'
            manureappNds = xr.open_dataset(manureappNfile)
            mmstotalappN = manureappNds.spring_TAN+manureappNds.spring_UAN+\
                        manureappNds.spring_availN+manureappNds.spring_resistN+manureappNds.spring_unavailN+\
                        manureappNds.winter_TAN+manureappNds.winter_UAN+\
                        manureappNds.winter_availN+manureappNds.winter_resistN+manureappNds.winter_unavailN
            landspreading_comptotalN[mmsidx][:] = mmstotalappN.values
            print("Manure spreading emission: ",livestock,prodsyst,mmstype,np.nansum(landspreading_compNH3[mmsidx])/1e9," GgN.")
            print("Manure spreading N: ",livestock,prodsyst,mmstype,np.nansum(landspreading_comptotalN[mmsidx])/1e9," GgN.")
            mmsidx = mmsidx + 1
        print("Manure spreading total emission: ",livestock,prodsyst,np.nansum(landspreading_NH3emission)/1e9," GgN.")
        print("Manure spreading total applied N: ",livestock,prodsyst,np.nansum(landspreading_comptotalN)/1e9," GgN.")
        
        
        
        outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss_housing=(['time','lat','lon'],housing_NH3emission),
                    NH3emiss_MMS=(['time','lat','lon'],mms_NH3emission),
                    NH3emiss_manureapp=(['time','lat','lon'],landspreading_NH3emission),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission from housing, manure management, manure application, and grazing",
                    units="gN per grid; m2 per grid",
                ),
            )
        
        outds.NH3emiss_housing.attrs["unit"] = 'gN/day'
        outds.NH3emiss_housing.attrs["long name"] = 'NH3 emission from housing'
        outds.NH3emiss_MMS.attrs["unit"] = 'gN/day'
        outds.NH3emiss_MMS.attrs["long name"] = 'NH3 emission from MMS'
        outds.NH3emiss_manureapp.attrs["unit"] = 'gN/day'
        outds.NH3emiss_manureapp.attrs["long name"] = 'NH3 emission from manure applicaton'

        outfilename = str(livestock)+'.NH3emission.'+str(yearofstudy)+'.timeseries.nc'
        outds.to_netcdf(outpath+outfilename)
        print("timeseries ncfile saved.")

        outds = xr.Dataset(
                data_vars=dict(
                    housing_NH3=(['lat','lon'],housing_emission.sum(dim="time").values),
                    housing_N=(['lat','lon'],housing_Nexcret),
                    mms_NH3=(['mmslevel','lat','lon'],mms_compNH3),
                    mms_N=(['mmslevel','lat','lon'],mms_comptotalN),
                    manureapp_NH3=(['applevel','lat','lon'],landspreading_compNH3),
                    manureapp_N=(['applevel','lat','lon'],landspreading_comptotalN),
                            ),
                coords = dict(
                    mmslevel=(["mmslevel"], mmslvl),
                    applevel=(["applevel"], applvl),
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(yearofstudy),
                    info = "NH3 emission and N from housing, manure management, manure application, and grazing. MMS level:\
                            0-liquid MMS_indoor, 1-liquid MMS_open, 2-liquid MMS_cover, 3-liquid MMS_lagoon,\
                            4-solid MMS_indoor, 5-solid MMS_open.  Manure app level:\
                            0-MMS_cover liquid, 1-MMS_indoor liquid, 2-MMS_indoor solid, 3-MMS_open liquid.",
                    units="gN per grid; m2 per grid",
                ),
            )

        outds.housing_NH3.attrs["unit"] = 'gN/day'
        outds.housing_NH3.attrs["long name"] = 'NH3 emission from housing'
        outds.housing_N.attrs["unit"] = 'gN/day'
        outds.housing_N.attrs["long name"] = 'Excreted N from housing'
        outds.mms_NH3.attrs["unit"] = 'gN/day'
        outds.mms_NH3.attrs["long name"] = 'NH3 emission from MMS'
        outds.mms_N.attrs["unit"] = 'gN/day'
        outds.mms_N.attrs["long name"] = 'N from MMS'
        outds.manureapp_NH3.attrs["unit"] = 'gN/day'
        outds.manureapp_NH3.attrs["long name"] = 'NH3 emission from manure applicaton'
        outds.manureapp_N.attrs["unit"] = 'gN/day'
        outds.manureapp_N.attrs["long name"] = 'Applied manure N for land application'

        outfilename = str(livestock)+'.NH3emission.annual.'+str(yearofstudy)+'.nc'
        outds.to_netcdf(outpath+outfilename)
        print("ncfile saved.")

if __name__ == "__main__":
    livestock_output()
    chemfert_output()