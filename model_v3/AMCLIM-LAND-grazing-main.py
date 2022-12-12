from mpi4py import MPI
import numpy as np
import xarray as xr
import pandas as pd
import time
import sys

import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.LAND as LAND


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

    # grazing_sim.grazing_main(livestock_name=livestock,production_system=prodsyst,
    #                         start_day_idx=0,end_day_idx=1,stat=True)

    mid = time.time()
    cal_time = mid - start
    print("=======================================================")
    print(livestock,prodsyst)
    print("Rank: ",rank)
    print('sim time: ', np.round(cal_time/60,decimals=2),' mins')

    newshape = int(365*720)
    sendNH3flux_urinepatch = np.transpose(grazing_sim.o_NH3flux,(1,0,2)).reshape((band,newshape))
    sendwashoff_urinepatch = np.transpose(grazing_sim.o_washoff,(1,0,2)).reshape((band,newshape))
    sendnitrif_urinepatch = np.transpose(grazing_sim.o_nitrif,(1,0,2)).reshape((band,newshape))
    sendNH4leaching_urinepatch = np.transpose(grazing_sim.o_NH4leaching,(1,0,2)).reshape((band,newshape))
    senddiffaq_urinepatch = np.transpose(grazing_sim.o_diffaq,(1,0,2)).reshape((band,newshape))
    senddiffgas_urinepatch = np.transpose(grazing_sim.o_diffgas,(1,0,2)).reshape((band,newshape))
    sendsoilTAN_urinepatch = np.transpose(grazing_sim.o_soil_TAN,(1,0,2)).reshape((band,newshape))
    sendsoilorgN_urinepatch = np.transpose(grazing_sim.o_soil_orgN,(1,0,2)).reshape((band,newshape))

    sendNH3flux_dungpat = np.transpose(grazing_sim.o_NH3flux_dung,(1,0,2)).reshape((band,newshape))
    sendwashoff_dungpat = np.transpose(grazing_sim.o_washoff_dung,(1,0,2)).reshape((band,newshape))
    sendnitrif_dungpat = np.transpose(grazing_sim.o_nitrif_dung,(1,0,2)).reshape((band,newshape))
    sendNH4leaching_dungpat = np.transpose(grazing_sim.o_NH4leaching_dung,(1,0,2)).reshape((band,newshape))
    senddiffaq_dungpat = np.transpose(grazing_sim.o_diffaq_dung,(1,0,2)).reshape((band,newshape))
    sendsoilTAN_dungpat = np.transpose(grazing_sim.o_soil_TAN_dung,(1,0,2)).reshape((band,newshape))
    sendsoilorgN_dungpat = np.transpose(grazing_sim.o_soil_orgN_dung,(1,0,2)).reshape((band,newshape))

    sendNH3flux_mixed = np.transpose(grazing_sim.o_NH3flux_FYM,(1,0,2)).reshape((band,newshape))
    sendwashoff_mixed = np.transpose(grazing_sim.o_washoff_FYM,(1,0,2)).reshape((band,newshape))
    sendnitrif_mixed = np.transpose(grazing_sim.o_nitrif_FYM,(1,0,2)).reshape((band,newshape))
    sendNH4leaching_mixed = np.transpose(grazing_sim.o_NH4leaching_FYM,(1,0,2)).reshape((band,newshape))
    senddiffaq_mixed = np.transpose(grazing_sim.o_diffaq_FYM,(1,0,2)).reshape((band,newshape))
    senddiffgas_mixed = np.transpose(grazing_sim.o_diffgas_FYM,(1,0,2)).reshape((band,newshape))
    sendsoilTAN_mixed = np.transpose(grazing_sim.o_soil_TAN_FYM,(1,0,2)).reshape((band,newshape))
    sendsoilorgN_mixed = np.transpose(grazing_sim.o_soil_orgN_FYM,(1,0,2)).reshape((band,newshape))

    NH3flux_urinepatch = None
    Nwashoff_urinepatch = None
    TANnitrif_urinepatch = None
    NH4leaching_urinepatch = None
    TANdiff_urinepatch = None
    NH3diff_urinepatch = None
    soilTAN_urinepatch = None
    soilorgN_urinepatch = None

    NH3flux_dungpat = None
    Nwashoff_dungpat = None
    TANnitrif_dungpat = None
    NH4leaching_dungpat = None
    TANdiff_dungpat = None
    NH3diff_dungpat = None
    soilTAN_dungpat = None
    soilorgN_dungpat = None

    NH3flux_mixed = None
    Nwashoff_mixed = None
    TANnitrif_mixed = None
    NH4leaching_mixed = None
    TANdiff_mixed = None
    NH3diff_mixed = None
    soilTAN_mixed = None
    soilorgN_mixed = None

    if rank == 0:
        NH3flux_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        Nwashoff_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANnitrif_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        NH4leaching_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANdiff_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        NH3diff_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilTAN_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilorgN_urinepatch = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))

        NH3flux_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        Nwashoff_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANnitrif_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        NH4leaching_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANdiff_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilTAN_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilorgN_dungpat = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))

        NH3flux_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        Nwashoff_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANnitrif_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        NH4leaching_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        TANdiff_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilTAN_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
        soilorgN_mixed = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))

    comm.Gather(sendNH3flux_urinepatch, NH3flux_urinepatch, root=0)
    comm.Gather(sendwashoff_urinepatch, Nwashoff_urinepatch, root=0)
    comm.Gather(sendnitrif_urinepatch, TANnitrif_urinepatch, root=0)
    comm.Gather(sendNH4leaching_urinepatch, NH4leaching_urinepatch, root=0)
    comm.Gather(senddiffaq_urinepatch, TANdiff_urinepatch, root=0)
    comm.Gather(senddiffgas_urinepatch, NH3diff_urinepatch, root=0)
    comm.Gather(sendsoilTAN_urinepatch, soilTAN_urinepatch, root=0)
    comm.Gather(sendsoilorgN_urinepatch, soilorgN_urinepatch, root=0)

    comm.Gather(sendNH3flux_dungpat, NH3flux_dungpat, root=0)
    comm.Gather(sendwashoff_dungpat, Nwashoff_dungpat, root=0)
    comm.Gather(sendnitrif_dungpat, TANnitrif_dungpat, root=0)
    comm.Gather(sendNH4leaching_dungpat, NH4leaching_dungpat, root=0)
    comm.Gather(senddiffaq_dungpat, TANdiff_dungpat, root=0)
    comm.Gather(sendsoilTAN_dungpat, soilTAN_dungpat, root=0)
    comm.Gather(sendsoilorgN_dungpat, soilorgN_dungpat, root=0)

    comm.Gather(sendNH3flux_mixed, NH3flux_mixed, root=0)
    comm.Gather(sendwashoff_mixed, Nwashoff_mixed, root=0)
    comm.Gather(sendnitrif_mixed, TANnitrif_mixed, root=0)
    comm.Gather(sendNH4leaching_mixed, NH4leaching_mixed, root=0)
    comm.Gather(senddiffaq_mixed, TANdiff_mixed, root=0)
    comm.Gather(sendsoilTAN_mixed, soilTAN_mixed, root=0)
    comm.Gather(sendsoilorgN_mixed, soilorgN_mixed, root=0)

    if rank == 0:
        nlat = int(180.0/CONFIG.CONFIG_dlat)
        nlon = int(360.0/CONFIG.CONFIG_dlon)
        ntime = CONFIG.Days
        lats = 90 - 0.5*np.arange(nlat)
        lons = -180 + 0.5*np.arange(nlon)
        yearidx = str(CONFIG.sim_year)+'-01-01'
        times = pd.date_range(yearidx,periods=ntime)

        outds = xr.Dataset(
            data_vars=dict(
                NH3emiss_urinepatch=(['time','lat','lon'],np.transpose(NH3flux_urinepatch.reshape((360,365,720)),(1,0,2))),
                TANwashoff_urinepatch=(['time','lat','lon'],np.transpose(Nwashoff_urinepatch.reshape((360,365,720)),(1,0,2))),
                NH4nitrif_urinepatch=(['time','lat','lon'],np.transpose(TANnitrif_urinepatch.reshape((360,365,720)),(1,0,2))),
                TANleaching_urinepatch=(['time','lat','lon'],np.transpose(NH4leaching_urinepatch.reshape((360,365,720)),(1,0,2))),
                TANdiffaq_urinepatch=(['time','lat','lon'],np.transpose(TANdiff_urinepatch.reshape((360,365,720)),(1,0,2))),
                NH3diffgas_urinepatch=(['time','lat','lon'],np.transpose(NH3diff_urinepatch.reshape((360,365,720)),(1,0,2))),
                TANtosoil_urinepatch=(['time','lat','lon'],np.transpose(soilTAN_urinepatch.reshape((360,365,720)),(1,0,2))),
                orgNtosoil_urinepatch=(['time','lat','lon'],np.transpose(soilorgN_urinepatch.reshape((360,365,720)),(1,0,2))),

                NH3emiss_dungpat=(['time','lat','lon'],np.transpose(NH3flux_dungpat.reshape((360,365,720)),(1,0,2))),
                TANwashoff_dungpat=(['time','lat','lon'],np.transpose(Nwashoff_dungpat.reshape((360,365,720)),(1,0,2))),
                NH4nitrif_dungpat=(['time','lat','lon'],np.transpose(TANnitrif_dungpat.reshape((360,365,720)),(1,0,2))),
                TANleaching_dungpat=(['time','lat','lon'],np.transpose(NH4leaching_dungpat.reshape((360,365,720)),(1,0,2))),
                TANdiffaq_dungpat=(['time','lat','lon'],np.transpose(TANdiff_dungpat.reshape((360,365,720)),(1,0,2))),
                TANtosoil_dungpat=(['time','lat','lon'],np.transpose(soilTAN_dungpat.reshape((360,365,720)),(1,0,2))),
                orgNtosoil_dungpat=(['time','lat','lon'],np.transpose(soilorgN_dungpat.reshape((360,365,720)),(1,0,2))),

                NH3emiss_mixed=(['time','lat','lon'],np.transpose(NH3flux_mixed.reshape((360,365,720)),(1,0,2))),
                TANwashoff_mixed=(['time','lat','lon'],np.transpose(Nwashoff_mixed.reshape((360,365,720)),(1,0,2))),
                NH4nitrif_mixed=(['time','lat','lon'],np.transpose(TANnitrif_mixed.reshape((360,365,720)),(1,0,2))),
                TANleaching_mixed=(['time','lat','lon'],np.transpose(NH4leaching_mixed.reshape((360,365,720)),(1,0,2))),
                TANdiffaq_mixed=(['time','lat','lon'],np.transpose(TANdiff_mixed.reshape((360,365,720)),(1,0,2))),
                TANtosoil_mixed=(['time','lat','lon'],np.transpose(soilTAN_mixed.reshape((360,365,720)),(1,0,2))),
                orgNtosoil_mixed=(['time','lat','lon'],np.transpose(soilorgN_mixed.reshape((360,365,720)),(1,0,2))),
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
            str(livestock)+'.grazing.'+str(CONFIG.sim_year)+'.nc'
        full_path = CONFIG.output_path+filename
        print("ncfile",filename)
        print("path",CONFIG.output_path)
        print("full path",CONFIG.output_path+filename)
        print("xarray saved.")
        outds.to_netcdf(path=str(full_path),mode="w",encoding=encoding)
        print("ncfile saved.")

end = time.time()
runtime = end - start
print(np.round(runtime/60,decimals = 1),' mins')
print("#############################################")