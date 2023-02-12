import numpy as np
import xarray as xr
import pandas as pd
import time
import sys

import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.PARAMETERS as PARA

sim_year = CONFIG.sim_year
yearofstudy = CONFIG.sim_year
Days = CONFIG.Days
nlat = CONFIG.CONFIG_lats
nlon = CONFIG.CONFIG_lons
ncfilepath = CONFIG.output_path + "hourly_sims/"
infile_path = CONFIG.infile_path
outpath = ncfilepath + "sectoral_emissions/"

ntime = Days
lats = 90 - 0.5*np.arange(nlat)
lons = -180 + 0.5*np.arange(nlon)
yearidx = str(sim_year)+'-01-01'
times = pd.date_range(yearidx,periods=ntime)
mmslvl = np.arange(6)
applvl = np.arange(4)

livestock_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','FEEDLOT_CATTLE','BUFFALO_BEEF','BUFFALO_DAIRY',
                    'SHEEP','GOAT',
                    'PIG',
                    'POULTRY']

ruminants_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','BUFFALO_BEEF','BUFFALO_DAIRY','SHEEP','GOAT']

CONFIG_production_system_dict = {
        'PIG':['industrial','intermediate','backyard'],
        'BEEF_CATTLE':['grassland','mixed',None],
        'DAIRY_CATTLE':['grassland','mixed',None],
        'OTHER_CATTLE':['grassland','mixed',None],
        'FEEDLOT_CATTLE':['feedlot',None,None],
        'BUFFALO_BEEF':['grassland','mixed',None],
        'BUFFALO_DAIRY':['grassland','mixed',None],
        'SHEEP':['grassland','mixed',None],
        'GOAT':['grassland','mixed',None],
        'POULTRY':['broiler','layer','backyard']
        }

housing_types = {
        'PIG':[["slat-pit_house_in_situ","slat-pit_house"],["barn"],["barn"]],
        'BEEF_CATTLE':[[None],["barn"]],
        'DAIRY_CATTLE':[[None],["barn"]],
        'OTHER_CATTLE':[[None],["barn"]],
        'BUFFALO_BEEF':[[None],["barn"]],
        'BUFFALO_DAIRY':[[None],["barn"]],
        'SHEEP':[[None],["barn"]],
        'GOAT':[[None],["barn"]],
        'POULTRY':[["poultry_house_litter","poultry_house"],
                    ["poultry_house_litter","poultry_house"],["poultry_house"]],
        }

####################################
## Source sector ncfiles: livestock
####################################
def livestock_output():
    print("###################################################")
    print("## livestock sector synthesis; ncfiles")
    print("###################################################")

    ## RUMINANTS: Housing, MMS, Land application, Grazing
    for livestock in ruminants_list:

        print("=======================================")
        print("Processing: ",livestock)

        infile_path = CONFIG.infile_path
        animal_data_path = 'animal_data/'
        animal_file_name = CONFIG.CONFIG_animal_file_dict[livestock]

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
                    str(sim_year)+'.nc'
        housingds = FUNC.open_ds(housingfile)
        housing_Nexcret = housingds.Nexcret.values - housingds.grazing_N.sum(dim="time").values
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
                            str(mmsphase)+'.'+str(sim_year)+'.nc'
                    mmsds = FUNC.open_ds(mmsfile)
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
                        '.surf.'+str(sim_year)+'.nc'
            manureappds = FUNC.open_ds(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[:] = manureappemiss.values + landspreading_NH3emission[:]
            landspreading_compNH3[mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(sim_year)+'.manureapp.nc'
            manureappNds = FUNC.open_ds(manureappNfile)
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

        print("[Grazing]")
        for prodsyst_idx in np.arange(2):
            prodsyst = prodsysts[prodsyst_idx]
            print("Production system: ",prodsyst)
            grazingfile = ncfilepath+str(livestock)+'/base.'+str(prodsyst)+'.'+str(livestock)+\
                        '.grazing.'+str(sim_year)+'.nc'
            grazingds = FUNC.open_ds(grazingfile)
            grazingemiss = grazingds.NH3emiss_urinepatch + \
                            grazingds.NH3emiss_dungpat + \
                            grazingds.NH3emiss_mixed
            if prodsyst == "grassland":
                allyeargrazing_NH3emission[:] = grazingemiss.values
                allyeargrazing_compNH3[:] = grazingemiss.sum(dim="time")
                grasslandds = FUNC.open_ds(infile_path+animal_data_path+animal_file_name)
                allyeargrazinging_comptotalN[:] = grasslandds['Excreted_N'][prodsyst_idx].values*1e3
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
                    description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

        outfilename = str(livestock)+'.NH3emission.'+str(sim_year)+'.timeseries.nc'
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
                    description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

        outfilename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        outds.to_netcdf(outpath+outfilename)
        print("ncfile saved.")


    ## FEEDLOT CATTLE: HOUSING, MMS, LAND application
    for livestock in ["FEEDLOT_CATTLE"]:

        print("=======================================")
        print("Processing: ",livestock)

        infile_path = CONFIG.infile_path
        animal_data_path = 'animal_data/'
        animal_file_name = CONFIG.CONFIG_animal_file_dict[livestock]

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
                    str(sim_year)+'.nc'
        housingds = FUNC.open_ds(housingfile)
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
                            str(mmsphase)+'.'+str(sim_year)+'.nc'
                    mmsds = FUNC.open_ds(mmsfile)
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
                        '.surf.'+str(sim_year)+'.nc'
            manureappds = FUNC.open_ds(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[:] = manureappemiss.values + landspreading_NH3emission[:]
            landspreading_compNH3[mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(sim_year)+'.manureapp.nc'
            manureappNds = FUNC.open_ds(manureappNfile)
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
                    description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

        outfilename = str(livestock)+'.NH3emission.'+str(sim_year)+'.timeseries.nc'
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
                    description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

        outfilename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        outds.to_netcdf(outpath+outfilename)
        print("ncfile saved.")


    ## Pigs: Housing, MMS, Land application
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
                    str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss_slat + housingds.NH3emiss_pit
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission (slat-pit): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (slat-pit): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'_insitu.'+\
                    str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
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
                        str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
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
                            str(mmsphase)+'.'+str(sim_year)+'.nc'
                    mmsds = FUNC.open_ds(mmsfile)
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
                                        str(mmsphase)+'.insitu.'+str(sim_year)+'.nc'
                                mmsds = FUNC.open_ds(mmsfile)
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
                        '.surf.'+str(sim_year)+'.nc'
            manureappds = FUNC.open_ds(manureappfile)
            manureappemiss = manureappds.NH3emiss
            landspreading_NH3emission[prodsyst_idx][:] = manureappemiss.values + landspreading_NH3emission[prodsyst_idx][:]
            landspreading_compNH3[prodsyst_idx][mmsidx][:] = manureappemiss.sum(dim="time").values

            manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                        str(sim_year)+'.manureapp.nc'
            manureappNds = FUNC.open_ds(manureappNfile)
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
                description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

    outfilename = str(livestock)+'.NH3emission.'+str(sim_year)+'.timeseries.nc'
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
                description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

    outfilename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")



    ## Poultry: Housing, MMS, Land application
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
                    str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
            housing_Nexcret = housingds.Nexcret.values 
            housing_emission = housingds.NH3emiss
            housing_compNH3[prodsyst_idx][:] = housing_emission.sum(dim="time").values
            housing_comptotalN[prodsyst_idx][:] = housing_Nexcret
            housing_NH3emission[prodsyst_idx][:] = housing_emission.values
            print("Housing emission (no litter): ",livestock,prodsyst,np.nansum(housing_emission)/1e9," GgN.")
            print("Housing total N (no litter): ",livestock,prodsyst,np.nansum(housing_Nexcret)/1e9," GgN.")
            housingfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(housetype)+'_litter.'+\
                    str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
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
                        str(sim_year)+'.nc'
            housingds = FUNC.open_ds(housingfile)
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
                                str(mmsphase)+'.'+str(sim_year)+'.nc'
                        mmsds = FUNC.open_ds(mmsfile)
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
                                            str(mmsphase)+'.insitu.'+str(sim_year)+'.nc'
                                    mmsds = FUNC.open_ds(mmsfile)
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
            mmsphase = "solid"
            mmscat = "MMS_open"
            mmsidx = 5
            print("MMS files: ",mmsidx,mmsphase,mmscat)
            mmsfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmscat)+'.'+\
                    str(mmsphase)+'.'+str(sim_year)+'.nc'
            mmsds = FUNC.open_ds(mmsfile)
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
                            '.surf.'+str(sim_year)+'.nc'
                manureappds = FUNC.open_ds(manureappfile)
                manureappemiss = manureappds.NH3emiss
                landspreading_NH3emission[prodsyst_idx][:] = manureappemiss.values + landspreading_NH3emission[prodsyst_idx][:]
                landspreading_compNH3[prodsyst_idx][mmsidx][:] = manureappemiss.sum(dim="time").values

                manureappNfile = ncfilepath+str(livestock)+'/'+str(livestock)+'.'+str(prodsyst)+'.'+str(mmstype)+'.'+\
                            str(sim_year)+'.manureapp.nc'
                manureappNds = FUNC.open_ds(manureappNfile)
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
                description="AMCLIM: NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

    outfilename = str(livestock)+'.NH3emission.'+str(sim_year)+'.timeseries.nc'
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
                description="AMCLIM: Annaul NH3 emissions from "+str(livestock)+" farming in " +str(sim_year),
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

    outfilename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")


    #################################
    ## STATISTICS
    #################################

    print("####################################")
    print("## STATISTICS: livestock sector")
    print("####################################")

    livestock_housing_emission = np.zeros((nlat,nlon))
    livestock_housing_N = np.zeros((nlat,nlon))
    livestock_mms_emission = np.zeros((nlat,nlon))
    livestock_mms_N = np.zeros((nlat,nlon))
    livestock_manureapp_emission = np.zeros((nlat,nlon))
    livestock_manureapp_N = np.zeros((nlat,nlon))
    ruminant_grazing_emission = np.zeros((nlat,nlon))
    ruminant_grazing_N = np.zeros((nlat,nlon))

    for livestock in ruminants_list:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        grazing_emiss = livestockds.seasonalgrazing_NH3 + livestockds.allyeargrazing_NH3 
        grazing_N = livestockds.seasonalgrazing_N + livestockds.allyeargrazing_N 
        
        print("housing emission",housing_emiss.sum().values/1e9,"GgN")
        print("housing N",housing_N.sum().values/1e9,"GgN")
        print("mms emission",mms_emiss.sum().values/1e9,"GgN")
        print("mms N",mms_N.sum().values/1e9,"GgN")
        print("manure app emission",app_emiss.sum().values/1e9,"GgN")
        print("manure app N",app_N.sum().values/1e9,"GgN")
        print("grazing emission",grazing_emiss.sum().values/1e9,"GgN")
        print("grazing N",grazing_N.sum().values/1e9,"GgN")
        
        livestock_housing_emission = livestock_housing_emission + housing_emiss
        livestock_housing_N = livestock_housing_N + housing_N.fillna(0)
        livestock_mms_emission = livestock_mms_emission + mms_emiss
        livestock_mms_N = livestock_mms_N + mms_N
        livestock_manureapp_emission = livestock_manureapp_emission + app_emiss
        livestock_manureapp_N = livestock_manureapp_N + app_N
        ruminant_grazing_emission = ruminant_grazing_emission + grazing_emiss
        ruminant_grazing_N = ruminant_grazing_N + grazing_N
        
    for livestock in ['FEEDLOT_CATTLE']:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        
        print("housing emission",housing_emiss.sum().values/1e9,"GgN")
        print("housing N",housing_N.sum().values/1e9,"GgN")
        print("mms emission",mms_emiss.sum().values/1e9,"GgN")
        print("mms N",mms_N.sum().values/1e9,"GgN")
        print("manure app emission",app_emiss.sum().values/1e9,"GgN")
        print("manure app N",app_N.sum().values/1e9,"GgN")

        livestock_housing_emission = livestock_housing_emission + housing_emiss
        livestock_housing_N = livestock_housing_N + housing_N.fillna(0)
        livestock_mms_emission = livestock_mms_emission + mms_emiss
        livestock_mms_N = livestock_mms_N + mms_N
        livestock_manureapp_emission = livestock_manureapp_emission + app_emiss
        livestock_manureapp_N = livestock_manureapp_N + app_N
        
    for livestock in ["PIG","POULTRY"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3.sum(dim="level")
        housing_N = livestockds.housing_N.sum(dim="level")
        mms_emiss = livestockds.mms_NH3.sum(dim=("level","mmslevel"))
        mms_N = livestockds.mms_N.sum(dim=("level","mmslevel"))
        app_emiss = livestockds.manureapp_NH3.sum(dim=("level","applevel"))
        app_N = livestockds.manureapp_N.sum(dim=("level","applevel"))
        
        print("housing emission",housing_emiss.sum().values/1e9,"GgN")
        print("housing N",housing_N.sum().values/1e9,"GgN")
        print("mms emission",mms_emiss.sum().values/1e9,"GgN")
        print("mms N",mms_N.sum().values/1e9,"GgN")
        print("manure app emission",app_emiss.sum().values/1e9,"GgN")
        print("manure app N",app_N.sum().values/1e9,"GgN")
        
        livestock_housing_emission = livestock_housing_emission + housing_emiss
        livestock_housing_N = livestock_housing_N + housing_N.fillna(0)
        livestock_mms_emission = livestock_mms_emission + mms_emiss
        livestock_mms_N = livestock_mms_N + mms_N
        livestock_manureapp_emission = livestock_manureapp_emission + app_emiss
        livestock_manureapp_N = livestock_manureapp_N + app_N
        
    livestock_total_emiss = livestock_housing_emission + livestock_mms_emission + livestock_manureapp_emission + ruminant_grazing_emission
    livestock_total_N = livestock_housing_N + ruminant_grazing_N
    print("livestock total NH3 emission",livestock_total_emiss.sum().values/1e12," TgN")
    # print("livestock total N",livestock_total_N.sum().values/1e12)


    cattle_emission = np.zeros((nlat,nlon))
    cattle_N = np.zeros((nlat,nlon))
    sheep_emission = np.zeros((nlat,nlon))
    sheep_N = np.zeros((nlat,nlon))
    goat_emission = np.zeros((nlat,nlon))
    goat_N = np.zeros((nlat,nlon))
    pig_emission = np.zeros((nlat,nlon))
    pig_N = np.zeros((nlat,nlon))
    poultry_emission = np.zeros((nlat,nlon))
    poultry_N = np.zeros((nlat,nlon))

    for livestock in ["BEEF_CATTLE","DAIRY_CATTLE","OTHER_CATTLE","BUFFALO_BEEF","BUFFALO_DAIRY"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        grazing_emiss = livestockds.seasonalgrazing_NH3 + livestockds.allyeargrazing_NH3 
        grazing_N = livestockds.seasonalgrazing_N + livestockds.allyeargrazing_N 
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing emission",grazing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing N",grazing_N.where(region_mask).sum().values/1e9,"GgN")
        
        cattle_emission = cattle_emission + housing_emiss + mms_emiss + app_emiss + grazing_emiss
        cattle_N = cattle_N + housing_N.fillna(0)+ grazing_N
        
    for livestock in ['FEEDLOT_CATTLE']:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        
        cattle_emission = cattle_emission + housing_emiss + mms_emiss + app_emiss
        cattle_N = cattle_N + housing_N.fillna(0)
        
    for livestock in ["SHEEP"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        grazing_emiss = livestockds.seasonalgrazing_NH3 + livestockds.allyeargrazing_NH3 
        grazing_N = livestockds.seasonalgrazing_N + livestockds.allyeargrazing_N 
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing emission",grazing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing N",grazing_N.where(region_mask).sum().values/1e9,"GgN")
        
        sheep_emission = housing_emiss + mms_emiss + app_emiss + grazing_emiss
        sheep_N = housing_N.fillna(0) + grazing_N

    for livestock in ["GOAT"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3
        housing_N = livestockds.housing_N
        mms_emiss = livestockds.mms_NH3.sum(dim="mmslevel")
        mms_N = livestockds.mms_N.sum(dim="mmslevel")
        app_emiss = livestockds.manureapp_NH3.sum(dim="applevel")
        app_N = livestockds.manureapp_N.sum(dim="applevel")
        grazing_emiss = livestockds.seasonalgrazing_NH3 + livestockds.allyeargrazing_NH3 
        grazing_N = livestockds.seasonalgrazing_N + livestockds.allyeargrazing_N 
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing emission",grazing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("grazing N",grazing_N.where(region_mask).sum().values/1e9,"GgN")
        
        goat_emission = housing_emiss + mms_emiss + app_emiss + grazing_emiss
        goat_N = housing_N.fillna(0) + grazing_N
        
    for livestock in ["PIG"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3.sum(dim="level")
        housing_N = livestockds.housing_N.sum(dim="level")
        mms_emiss = livestockds.mms_NH3.sum(dim=("level","mmslevel"))
        mms_N = livestockds.mms_N.sum(dim=("level","mmslevel"))
        app_emiss = livestockds.manureapp_NH3.sum(dim=("level","applevel"))
        app_N = livestockds.manureapp_N.sum(dim=("level","applevel"))
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        
        pig_emission = housing_emiss + mms_emiss + app_emiss 
        pig_N = housing_N.fillna(0) 
        
    for livestock in ["POULTRY"]:
        print(livestock)
        filename = str(livestock)+'.NH3emission.annual.'+str(sim_year)+'.nc'
        livestockds = xr.open_dataset(outpath+filename)
        housing_emiss = livestockds.housing_NH3.sum(dim="level")
        housing_N = livestockds.housing_N.sum(dim="level")
        mms_emiss = livestockds.mms_NH3.sum(dim=("level","mmslevel"))
        mms_N = livestockds.mms_N.sum(dim=("level","mmslevel"))
        app_emiss = livestockds.manureapp_NH3.sum(dim=("level","applevel"))
        app_N = livestockds.manureapp_N.sum(dim=("level","applevel"))
        
        # print("housing emission",housing_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("housing N",housing_N.where(region_mask).sum().values/1e9,"GgN")
        # print("mms emission",mms_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("mms N",mms_N.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app emission",app_emiss.where(region_mask).sum().values/1e9,"GgN")
        # print("manure app N",app_N.where(region_mask).sum().values/1e9,"GgN")
        
        poultry_emission = housing_emiss + mms_emiss + app_emiss 
        poultry_N = housing_N.fillna(0)
        
    print(cattle_emission.sum().values/1e9,sheep_emission.sum().values/1e9,
        pig_emission.sum().values/1e9,poultry_emission.sum().values/1e9)
    print(cattle_N.sum().values/1e9,sheep_N.sum().values/1e9,
        pig_N.sum().values/1e9,poultry_N.sum().values/1e9)


########################
## Chemical fertilizer
#########################

crop_list = ['barley',
            'cassava',
            'cotton',
            'groundnut',
            'maize',
            'millet',
            'potato',
            'rapeseed',
            'rice',
            'rye',
            'sorghum',
            'soybean',
            'sugarbeet',
            'sugarcane',
            'sunflower',
            'wheat']


ferts = ['ammonium','urea']
techs = ['surf','disk','injection']

years = np.arange(2000,2020).astype(int)

def crop_output():

    NH4appemission = np.zeros((3,len(crop_list),Days,nlat,nlon))
    allcropNH4 = np.zeros((3,len(crop_list),nlat,nlon))
    ureaappemission = np.zeros((3,len(crop_list),Days,nlat,nlon))
    allcropurea = np.zeros((3,len(crop_list),nlat,nlon))
    allcropnitrate = np.zeros((3,len(crop_list),nlat,nlon))

    ## other N pathways
    NH4app_TANwashoff = np.zeros((3,len(crop_list),nlat,nlon))
    NH4app_TANleaching = np.zeros((3,len(crop_list),nlat,nlon))
    NH4app_TANdiffaq = np.zeros((3,len(crop_list),nlat,nlon))
    NH4app_TANdiffgas = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NH4nitrif = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NH4uptake = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NO3uptake = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NO3washoff = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NO3leaching = np.zeros((3,len(crop_list),nlat,nlon)) 
    NH4app_NO3diff = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_TANwashoff = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_TANleaching = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_TANdiffaq = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_TANdiffgas = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_NH4nitrif = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_NH4uptake = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_NO3uptake = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_NO3washoff = np.zeros((3,len(crop_list),nlat,nlon)) 
    ureaapp_NO3leaching = np.zeros((3,len(crop_list),nlat,nlon))
    ureaapp_NO3diff = np.zeros((3,len(crop_list),nlat,nlon))


    for fert in ["ammonium","urea"]:
        for ii in np.arange(3):
            for nn in np.arange(len(crop_list)):
                simdf = FUNC.open_ds(ncfilepath+'chemfert/base.'+crop_list[nn]+'.'+fert+'.'+techs[ii]+'.'+str(yearofstudy)+'.pp.nc')
                cropinfo = FUNC.open_ds(infile_path+'crop_data/AMCLIM_fert_info/chemfertinfo.'+crop_list[nn]+'.21st.nc')
                if fert == 'ammonium':
                    ammN = cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.ammN_area.sel(year=yearofstudy)
                    nitN = cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.nitN_area.sel(year=yearofstudy)
                    totalammN = np.nan_to_num(ammN.values)
                    totalnitN = np.nan_to_num(nitN.values)
                    NH3emiss = np.nan_to_num(simdf.NH3emiss)
                    TANwashoff = np.nan_to_num(simdf.TANwashoff)
                    TANleaching = np.nan_to_num(simdf.TANleaching)
                    TANdiffaq = np.nan_to_num(simdf.TANdiffaq)
                    TANdiffgas = np.nan_to_num(simdf.NH3diffgas)
                    NH4nitrif = np.nan_to_num(simdf.NH4nitrif)
                    NH4uptake = np.nan_to_num(simdf.NH4uptake)
                    NO3uptake = np.nan_to_num(simdf.NO3uptake)
                    NO3washoff = np.nan_to_num(simdf.NO3washoff)
                    NO3leaching = np.nan_to_num(simdf.NO3leaching)
                    NO3diff = np.nan_to_num(simdf.NO3diff)
                    
                    NH4appemission[ii][nn] = NH3emiss
                    NH4app_TANwashoff[ii][nn] = TANwashoff
                    NH4app_TANleaching[ii][nn] = TANleaching
                    NH4app_TANdiffaq[ii][nn] = TANdiffaq
                    NH4app_TANdiffgas[ii][nn] = TANdiffgas
                    NH4app_NH4nitrif[ii][nn] = NH4nitrif
                    NH4app_NH4uptake[ii][nn] = NH4uptake
                    NH4app_NO3uptake[ii][nn] = NO3uptake
                    NH4app_NO3washoff[ii][nn] = NO3washoff
                    NH4app_NO3leaching[ii][nn] = NO3leaching
                    NH4app_NO3diff[ii][nn] = NO3diff
                    allcropNH4[ii][nn] = totalammN
                    allcropnitrate[ii][nn] = totalnitN
                    print(fert,crop_list[nn],techs[ii])
                    print("Crop emission: ",np.nansum(NH3emiss)/1e9,"Crop N app: ",np.nansum(totalammN)/1e9)
                    print("Emission: ",np.nansum(NH4appemission[ii])/1e9,"N app: ",np.nansum(allcropNH4[ii])/1e9)
                elif fert == 'urea':
                    ureaN = cropinfo.Nrate.sel(year=yearofstudy)*cropinfo.ureaN_area.sel(year=yearofstudy)
                    totalureaN = np.nan_to_num(ureaN.values)
                    NH3emiss = np.nan_to_num(simdf.NH3emiss)
                    TANwashoff = np.nan_to_num(simdf.TANwashoff)
                    TANleaching = np.nan_to_num(simdf.TANleaching)
                    TANdiffaq = np.nan_to_num(simdf.TANdiffaq)
                    TANdiffgas = np.nan_to_num(simdf.NH3diffgas)
                    NH4nitrif = np.nan_to_num(simdf.NH4nitrif)
                    NH4uptake = np.nan_to_num(simdf.NH4uptake)
                    NO3uptake = np.nan_to_num(simdf.NO3uptake)
                    NO3washoff = np.nan_to_num(simdf.NO3washoff)
                    NO3leaching = np.nan_to_num(simdf.NO3leaching)
                    NO3diff = np.nan_to_num(simdf.NO3diff)
                    
                    ureaappemission[ii][nn] = NH3emiss
                    ureaapp_TANwashoff[ii][nn] = TANwashoff
                    ureaapp_TANleaching[ii][nn] = TANleaching
                    ureaapp_TANdiffaq[ii][nn] = TANdiffaq
                    ureaapp_TANdiffgas[ii][nn] = TANdiffgas
                    ureaapp_NH4nitrif[ii][nn] = NH4nitrif
                    ureaapp_NH4uptake[ii][nn] = NH4uptake
                    ureaapp_NO3uptake[ii][nn] = NO3uptake
                    ureaapp_NO3washoff[ii][nn] = NO3washoff
                    ureaapp_NO3leaching[ii][nn] = NO3leaching
                    ureaapp_NO3diff[ii][nn] = NO3diff
                    allcropurea[ii][nn] = totalureaN
                    print(fert,crop_list[nn],techs[ii])
                    print("Crop emission: ",np.nansum(NH3emiss)/1e9,"Crop N app: ",np.nansum(totalureaN)/1e9)
                    print("Emission: ",np.nansum(ureaappemission[ii])/1e9,"N app: ",np.nansum(allcropurea[ii])/1e9)


    wbincomelvlds = FUNC.open_ds(infile_path+'crop_data/WB_21st_century_income_classification.nc')
    wb_incomelvl = wbincomelvlds.income_lvl.sel(year=yearofstudy)

    tech_algorithim = PARA.tech_algorithim

    global_tech_used = np.zeros((3,nlat,nlon))
    for ii in np.arange(3):
        for jj in np.arange(1,5):
            global_tech_used[ii][wb_incomelvl.values==jj] = tech_algorithim[jj-1][ii]


    chemfert_total_emiss = np.zeros((2,len(crop_list),Days,nlat,nlon))
    chemfert_Napp = np.zeros((2,len(crop_list),nlat,nlon))

    chemfert_total_emiss[0] = NH4appemission[0]*global_tech_used[0] + \
                                NH4appemission[1]*global_tech_used[1] + \
                                    NH4appemission[2]*global_tech_used[2]
    chemfert_total_emiss[1] = ureaappemission[0]*global_tech_used[0] + \
                                ureaappemission[1]*global_tech_used[1] + \
                                    ureaappemission[2]*global_tech_used[2]

    chemfert_Napp[0] = allcropNH4[0]
    chemfert_Napp[1] = allcropurea[0]

    croplvls = np.arange(16)
    fertlvls = np.arange(2)


    outds = xr.Dataset(
            data_vars=dict(
                chemfert_NH3=(['fertlvl','croplvl','time','lat','lon'],chemfert_total_emiss),
                chemfert_Napp=(['fertlvl','croplvl','lat','lon'],chemfert_Napp/1e3),
                nitapp=(['croplvl','lat','lon'],allcropnitrate[0]/1e3),
            ),
            coords = dict(
                fertlvl=(["fertlvl"], fertlvls),
                croplvl=(["croplvl"], croplvls),
                time=(["time"], times),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: Annaul NH3 emissions from chemical fertilizer application in " +str(sim_year),
                info = "NH3 emission from ammonium and urea fertilizer application.\
                        Fertilizer level: 0-ammonium, 1-urea. Nitrate application is listed separately.\
                        Crop level: 0)barley, 1)cassava, 2)cotton, 3)groundnut, 4)maize, 5)millet, 6)potato, 7)rapeseed,\
                                    8)rice ,9)rye, 10)sorghum, 11)soybean, 12)sugarbeet, 13)sugarcane, 14)sunflower, 15)wheat",
                units= "gN per grid (emission); kgN per grid (application)",
            ),
        )
    outds.chemfert_NH3.attrs["unit"] = 'gN/day'
    outds.chemfert_NH3.attrs["long name"] = 'NH3 emission from chemical fertilizer application'
    outds.chemfert_Napp.attrs["unit"] = 'kgN/year'
    outds.chemfert_Napp.attrs["long name"] = 'Chemical fertilizer application'    
    outds.nitapp.attrs["unit"] = 'kgN/year'
    outds.nitapp.attrs["long name"] = 'Chemical fertilizer (nitrate) application'

    outfilename = 'base.chemfert.NH3emission.'+str(sim_year)+'.timeseries.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

    chemfert_NH3emiss = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_TANwashoff = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_TANleaching = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_TANdiffaq = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_TANdiffgas = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NH4nitrif = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NH4uptake = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NO3uptake = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NO3washoff = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NO3leaching = np.zeros((2,len(crop_list),nlat,nlon))
    chemfert_NO3diff = np.zeros((2,len(crop_list),nlat,nlon))


    chemfert_NH3emiss[0] = np.nansum(chemfert_total_emiss[0],axis=1)
    chemfert_NH3emiss[1] = np.nansum(chemfert_total_emiss[1],axis=1)


    chemfert_TANwashoff[0] = NH4app_TANwashoff[0]*global_tech_used[0] + \
                                NH4app_TANwashoff[1]*global_tech_used[1] + \
                                    NH4app_TANwashoff[2]*global_tech_used[2]
    chemfert_TANwashoff[1] = ureaapp_TANwashoff[0]*global_tech_used[0] + \
                                ureaapp_TANwashoff[1]*global_tech_used[1] + \
                                    ureaapp_TANwashoff[2]*global_tech_used[2]

    chemfert_TANleaching[0] = NH4app_TANleaching[0]*global_tech_used[0] + \
                                NH4app_TANleaching[1]*global_tech_used[1] + \
                                    NH4app_TANleaching[2]*global_tech_used[2]
    chemfert_TANleaching[1] = ureaapp_TANleaching[0]*global_tech_used[0] + \
                                ureaapp_TANleaching[1]*global_tech_used[1] + \
                                    ureaapp_TANleaching[2]*global_tech_used[2]

    chemfert_TANdiffaq[0] = NH4app_TANdiffaq[0]*global_tech_used[0] + \
                                NH4app_TANdiffaq[1]*global_tech_used[1] + \
                                    NH4app_TANdiffaq[2]*global_tech_used[2]
    chemfert_TANdiffaq[1] = ureaapp_TANdiffaq[0]*global_tech_used[0] + \
                                ureaapp_TANdiffaq[1]*global_tech_used[1] + \
                                    ureaapp_TANdiffaq[2]*global_tech_used[2]

    chemfert_TANdiffgas[0] = NH4app_TANdiffgas[0]*global_tech_used[0] + \
                                NH4app_TANdiffgas[1]*global_tech_used[1] + \
                                    NH4app_TANdiffgas[2]*global_tech_used[2]
    chemfert_TANdiffgas[1] = ureaapp_TANdiffgas[0]*global_tech_used[0] + \
                                ureaapp_TANdiffgas[1]*global_tech_used[1] + \
                                    ureaapp_TANdiffgas[2]*global_tech_used[2]

    chemfert_NH4nitrif[0] = NH4app_NH4nitrif[0]*global_tech_used[0] + \
                                NH4app_NH4nitrif[1]*global_tech_used[1] + \
                                    NH4app_NH4nitrif[2]*global_tech_used[2]
    chemfert_NH4nitrif[1] = ureaapp_NH4nitrif[0]*global_tech_used[0] + \
                                ureaapp_NH4nitrif[1]*global_tech_used[1] + \
                                    ureaapp_NH4nitrif[2]*global_tech_used[2]

    chemfert_NH4uptake[0] = NH4app_NH4uptake[0]*global_tech_used[0] + \
                                NH4app_NH4uptake[1]*global_tech_used[1] + \
                                    NH4app_NH4uptake[2]*global_tech_used[2]
    chemfert_NH4uptake[1] = ureaapp_NH4uptake[0]*global_tech_used[0] + \
                                ureaapp_NH4uptake[1]*global_tech_used[1] + \
                                    ureaapp_NH4uptake[2]*global_tech_used[2]

    chemfert_NO3uptake[0] = NH4app_NO3uptake[0]*global_tech_used[0] + \
                                NH4app_NO3uptake[1]*global_tech_used[1] + \
                                    NH4app_NO3uptake[2]*global_tech_used[2]
    chemfert_NO3uptake[1] = ureaapp_NO3uptake[0]*global_tech_used[0] + \
                                ureaapp_NO3uptake[1]*global_tech_used[1] + \
                                    ureaapp_NO3uptake[2]*global_tech_used[2]

    chemfert_NO3washoff[0] = NH4app_NO3washoff[0]*global_tech_used[0] + \
                                NH4app_NO3washoff[1]*global_tech_used[1] + \
                                    NH4app_NO3washoff[2]*global_tech_used[2]
    chemfert_NO3washoff[1] = ureaapp_NO3washoff[0]*global_tech_used[0] + \
                                ureaapp_NO3washoff[1]*global_tech_used[1] + \
                                    ureaapp_NO3washoff[2]*global_tech_used[2]

    chemfert_NO3leaching[0] = NH4app_NO3leaching[0]*global_tech_used[0] + \
                                NH4app_NO3leaching[1]*global_tech_used[1] + \
                                    NH4app_NO3leaching[2]*global_tech_used[2]
    chemfert_NO3leaching[1] = ureaapp_NO3leaching[0]*global_tech_used[0] + \
                                ureaapp_NO3leaching[1]*global_tech_used[1] + \
                                    ureaapp_NO3leaching[2]*global_tech_used[2]

    chemfert_NO3diff[0] = NH4app_NO3diff[0]*global_tech_used[0] + \
                                NH4app_NO3diff[1]*global_tech_used[1] + \
                                    NH4app_NO3diff[2]*global_tech_used[2]
    chemfert_NO3diff[1] = ureaapp_NO3diff[0]*global_tech_used[0] + \
                                ureaapp_NO3diff[1]*global_tech_used[1] + \
                                    ureaapp_NO3diff[2]*global_tech_used[2]

    outds = xr.Dataset(
            data_vars=dict(
                chemfert_Napp=(['fertlvl','lat','lon'],np.nansum(chemfert_Napp,axis=1)/1e3),
                nitapp=(['lat','lon'],np.nansum(allcropnitrate[0],axis=0)/1e3),
                NH3emiss=(['fertlvl','lat','lon'],np.nansum(chemfert_NH3emiss,axis=1)),
                TANwashoff=(['fertlvl','lat','lon'],np.nansum(chemfert_TANwashoff,axis=1)),
                TANleaching=(['fertlvl','lat','lon'],np.nansum(chemfert_TANleaching,axis=1)),
                TANdiffaq=(['fertlvl','lat','lon'],np.nansum(chemfert_TANdiffaq,axis=1)),
                NH3diffgas=(['fertlvl','lat','lon'],np.nansum(chemfert_TANdiffgas,axis=1)),
                NH4nitrif=(['fertlvl','lat','lon'],np.nansum(chemfert_NH4nitrif,axis=1)),
                NH4uptake=(['fertlvl','lat','lon'],np.nansum(chemfert_NH4uptake,axis=1)),
                NO3uptake=(['fertlvl','lat','lon'],np.nansum(chemfert_NO3uptake,axis=1)),
                NO3washoff=(['fertlvl','lat','lon'],np.nansum(chemfert_NO3washoff,axis=1)),
                NO3leaching=(['fertlvl','lat','lon'],np.nansum(chemfert_NO3leaching,axis=1)),
                NO3diff=(['fertlvl','lat','lon'],np.nansum(chemfert_NO3diff,axis=1)),
            ),
            coords = dict(
                fertlvl=(["fertlvl"], fertlvls),
                lon=(["lon"], lons),
                lat=(["lat"], lats),
                        ),
            attrs=dict(
                description="AMCLIM: Annaul NH3 emissions from chemical fertilizer application in " +str(sim_year),
                info = "NH3 emission from ammonium and urea fertilizer application.\
                        Fertilizer level: 0-ammonium, 1-urea. Nitrate application is listed separately.",
                units= "gN per grid (emission); kgN per grid (application)",
            ),
        )

    outds.chemfert_Napp.attrs["unit"] = 'kgN/year'
    outds.chemfert_Napp.attrs["long name"] = 'Chemical fertilizer application'  
    outds.nitapp.attrs["unit"] = 'kgN/year'
    outds.nitapp.attrs["long name"] = 'Chemical fertilizer (nitrate) application'
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
    
    outfilename = 'base.chemfert.NH3emission.Npathways.annual.'+str(sim_year)+'.nc'
    outds.to_netcdf(outpath+outfilename)
    print("ncfile saved.")

    chemfertemissds = xr.open_dataset(ncfilepath+'base.chemfert.NH3emission.Npathways.annual.'+str(sim_year)+'.nc')
    print("N app",chemfertemissds.chemfert_Napp[0].sum().values/1e6,
              chemfertemissds.chemfert_Napp[1].sum().values/1e6,
                 chemfertemissds.nitapp.sum().values/1e6)

    print("NH3 emiss",chemfertemissds.NH3emiss.sum().values/1e9)
    print("TAN washoff",chemfertemissds.TANwashoff.sum().values/1e9)
    print("TAN leaching",chemfertemissds.TANleaching.sum().values/1e9)
    print("TAN diffaq",chemfertemissds.TANdiffaq.sum().values/1e9)
    print("TAN diffgas",chemfertemissds.NH3diffgas.sum().values/1e9)
    print("nitrif",chemfertemissds.NH4nitrif.sum().values/1e9)
    print("NH4 uptake",chemfertemissds.NH4uptake.sum().values/1e9)
    print("NO3 uptake",chemfertemissds.NO3uptake.sum().values/1e9)
    print("NO3 washoff",chemfertemissds.NO3washoff.sum().values/1e9)
    print("NO3 leaching",chemfertemissds.NO3leaching.sum().values/1e9)
    print("NO3 diff",chemfertemissds.NO3diff.sum().values/1e9)

livestock_output()
# crop_output()