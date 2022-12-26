from mpi4py import MPI
import numpy as np
import xarray as xr
import pandas as pd
import time
import sys

import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.LAND_develop as LAND
# import MODULES.LAND_v3 as LAND

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

livestock = str(sys.argv[1])

startday = 0
endday = 730
band = int(360/size)
arr = (730,band,720)

ruminant_prodsysts = ["grassland","mixed"]
sim = "base"

start = time.time()

for prodsyst in ruminant_prodsysts:
# for prodsyst in ["mixed"]:
    print(prodsyst)
    grazing_sim = LAND.LAND_module(prank=rank,psize=band,
                        fert_type='manure',manure_added=np.zeros(arr),urea_added=np.zeros(arr),
                        UA_added=np.zeros(arr),avail_N_added=np.zeros(arr),resist_N_added=np.zeros(arr),
                        unavail_N_added=np.zeros(arr),
                        TAN_added=np.zeros(arr),NO3_added=np.zeros(arr),water_added=np.zeros(arr),
                        pH_value=7.0,grazing=True)

    grazing_sim.grazing_main(livestock_name=livestock,production_system=prodsyst,
                                start_day_idx=0,end_day_idx=CONFIG.Days,stat=True)

    mid = time.time()
    cal_time = mid - start
    print("=======================================================")
    print(livestock,prodsyst)
    print("Rank: ",rank)
    print('sim time: ', np.round(cal_time/60,decimals=2),' mins')

    nlat = int((180.0/CONFIG.CONFIG_dlat)/size)
    nlon = int(360.0/CONFIG.CONFIG_dlon)
    ntime = CONFIG.Days
    lats = 90 - 0.5*np.arange(rank*nlat,(rank+1)*nlat)
    lons = -180 + 0.5*np.arange(nlon)
    yearidx = str(CONFIG.sim_year)+'-01-01'
    times = pd.date_range(yearidx,periods=ntime)

    outds = xr.Dataset(
        data_vars=dict(
            NH3emiss_urinepatch=(['time','lat','lon'],grazing_sim.o_NH3flux),
            TANwashoff_urinepatch=(['time','lat','lon'],grazing_sim.o_washoff),
            NH4nitrif_urinepatch=(['time','lat','lon'],grazing_sim.o_nitrif),
            TANleaching_urinepatch=(['time','lat','lon'],grazing_sim.o_NH4leaching),
            TANdiffaq_urinepatch=(['time','lat','lon'],grazing_sim.o_diffaq),
            NH3diffgas_urinepatch=(['time','lat','lon'],grazing_sim.o_diffgas),
            TANtosoil_urinepatch=(['time','lat','lon'],grazing_sim.o_soil_TAN),
            orgNtosoil_urinepatch=(['time','lat','lon'],grazing_sim.o_soil_orgN),

            NH3emiss_dungpat=(['time','lat','lon'],grazing_sim.o_NH3flux_dung),
            TANwashoff_dungpat=(['time','lat','lon'],grazing_sim.o_washoff_dung),
            NH4nitrif_dungpat=(['time','lat','lon'],grazing_sim.o_nitrif_dung),
            TANleaching_dungpat=(['time','lat','lon'],grazing_sim.o_NH4leaching_dung),
            TANdiffaq_dungpat=(['time','lat','lon'],grazing_sim.o_diffaq_dung),
            TANtosoil_dungpat=(['time','lat','lon'],grazing_sim.o_soil_TAN_dung),
            orgNtosoil_dungpat=(['time','lat','lon'],grazing_sim.o_soil_TAN_dung),

            NH3emiss_mixed=(['time','lat','lon'],grazing_sim.o_NH3flux_FYM),
            TANwashoff_mixed=(['time','lat','lon'],grazing_sim.o_washoff_FYM),
            NH4nitrif_mixed=(['time','lat','lon'],grazing_sim.o_nitrif_FYM),
            TANleaching_mixed=(['time','lat','lon'],grazing_sim.o_NH4leaching_FYM),
            TANdiffaq_mixed=(['time','lat','lon'],grazing_sim.o_diffaq_FYM),
            TANtosoil_mixed=(['time','lat','lon'],grazing_sim.o_soil_TAN_FYM),
            orgNtosoil_mixed=(['time','lat','lon'],grazing_sim.o_soil_orgN_FYM),
                    ),
        coords = dict(
            time=(["time"], times),
            lon=(["lon"], lons),
            lat=(["lat"], lats),
                    ),
        attrs=dict(
            description="AMCLIM-Land_grazing: N pathways in " +str(CONFIG.sim_year),
            units="gN per grid",
        ),
    )

    outds.NH3emiss_urinepatch.attrs["unit"] = 'gN/day'
    outds.NH3emiss_urinepatch.attrs["long name"] = 'NH3 emission from fertilizer application'
    outds.TANwashoff_urinepatch.attrs["unit"] = 'gN/day'
    outds.TANwashoff_urinepatch.attrs["long name"] = 'TAN washoff from fertilizer application'
    outds.NH4nitrif_urinepatch.attrs["unit"] = 'gN/day'
    outds.NH4nitrif_urinepatch.attrs["long name"] = 'Nitrification'
    outds.TANleaching_urinepatch.attrs["unit"] = 'gN/day'
    outds.TANleaching_urinepatch.attrs["long name"] = 'TAN leaching from fertilizer application'
    outds.TANdiffaq_urinepatch.attrs["unit"] = 'gN/day'
    outds.TANdiffaq_urinepatch.attrs["long name"] = 'TAN diffusion to deep soil from fertilizer application'
    outds.NH3diffgas_urinepatch.attrs["unit"] = 'gN/day'
    outds.NH3diffgas_urinepatch.attrs["long name"] = 'NH3 diffusion to deep soil from fertilizer application'
    outds.TANtosoil_urinepatch.attrs["unit"] = 'gN/day'
    outds.TANtosoil_urinepatch.attrs["long name"] = 'TAN incorporated in soil'
    outds.orgNtosoil_urinepatch.attrs["unit"] = 'gN/day'
    outds.orgNtosoil_urinepatch.attrs["long name"] = 'OrgN incorporated in soil'

    outds.NH3emiss_dungpat.attrs["unit"] = 'gN/day'
    outds.NH3emiss_dungpat.attrs["long name"] = 'NH3 emission from fertilizer application'
    outds.TANwashoff_dungpat.attrs["unit"] = 'gN/day'
    outds.TANwashoff_dungpat.attrs["long name"] = 'TAN washoff from fertilizer application'
    outds.NH4nitrif_dungpat.attrs["unit"] = 'gN/day'
    outds.NH4nitrif_dungpat.attrs["long name"] = 'Nitrification'
    outds.TANleaching_dungpat.attrs["unit"] = 'gN/day'
    outds.TANleaching_dungpat.attrs["long name"] = 'TAN leaching from fertilizer application'
    outds.TANdiffaq_dungpat.attrs["unit"] = 'gN/day'
    outds.TANdiffaq_dungpat.attrs["long name"] = 'TAN diffusion to deep soil from fertilizer application'
    outds.TANtosoil_dungpat.attrs["unit"] = 'gN/day'
    outds.TANtosoil_dungpat.attrs["long name"] = 'TAN incorporated in soil'
    outds.orgNtosoil_dungpat.attrs["unit"] = 'gN/day'
    outds.orgNtosoil_dungpat.attrs["long name"] = 'OrgN incorporated in soil'

    outds.NH3emiss_mixed.attrs["unit"] = 'gN/day'
    outds.NH3emiss_mixed.attrs["long name"] = 'NH3 emission from fertilizer application'
    outds.TANwashoff_mixed.attrs["unit"] = 'gN/day'
    outds.TANwashoff_mixed.attrs["long name"] = 'TAN washoff from fertilizer application'
    outds.NH4nitrif_mixed.attrs["unit"] = 'gN/day'
    outds.NH4nitrif_mixed.attrs["long name"] = 'Nitrification'
    outds.TANleaching_mixed.attrs["unit"] = 'gN/day'
    outds.TANleaching_mixed.attrs["long name"] = 'TAN leaching from fertilizer application'
    outds.TANdiffaq_mixed.attrs["unit"] = 'gN/day'
    outds.TANdiffaq_mixed.attrs["long name"] = 'TAN diffusion to deep soil from fertilizer application'
    outds.TANtosoil_mixed.attrs["unit"] = 'gN/day'
    outds.TANtosoil_mixed.attrs["long name"] = 'TAN incorporated in soil'
    outds.orgNtosoil_mixed.attrs["unit"] = 'gN/day'
    outds.orgNtosoil_mixed.attrs["long name"] = 'OrgN incorporated in soil'


    comp = dict(zlib=True, complevel=9)
    encoding = {var: comp for var in outds.data_vars}

    filename = sim+'.'+str(prodsyst)+'.'+\
        str(livestock)+'.grazing.'+str(rank)+'.'+str(CONFIG.sim_year)+'.nc'
    full_path = CONFIG.output_path+'grazing_compnc/'+filename
    print("ncfile",filename)
    print("path",CONFIG.output_path+'grazing_compnc/')
    print("full path",full_path)
    print("xarray saved.")
    outds.to_netcdf(path=str(full_path),mode="w",encoding=encoding)
    print("ncfile saved.")

end = time.time()
runtime = end - start
print(np.round(runtime/60,decimals = 1),' mins')
print("#############################################")