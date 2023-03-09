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

if livestock == "PIG1":
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict["PIG"][0]
    livestock = "PIG"
elif livestock == "PIG2":
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict["PIG"][1]
    livestock = "PIG"
elif livestock == "PIG3":
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict["PIG"][2]
    livestock = "PIG"

elif livestock == "broiler_CHICKEN":
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict["POULTRY"][0]
    livestock = "POULTRY"
elif livestock == "layer_CHICKEN":
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict["POULTRY"][1]
    livestock = "POULTRY"
else:
    prodsyst = CONFIG.CONFIG_livestock_prodsyst_Nappland_lict[livestock][0]
print(livestock,prodsyst)

startday = 0
endday = CONFIG.Days*2
# methods = ['broadcasting-surf','incorporating-disk','deep injection']
methods = ['broadcasting-surf']

techs = {
    "broadcasting-surf": 'surf',
    "incorporating-disk": 'disk',
    "deep injection": 'injection'
}

manure_types = {
    "liquid": "slurry",
    "solid": "solid manure",
}

mmscats = {
    "liquid":["MMS_indoor","MMS_open","MMS_cover"],
    "solid":["MMS_indoor"],
}

if livestock == "BUFFALO_BEEF":
    mmscats = {
        "liquid":["MMS_indoor","MMS_open"],
        "solid":["MMS_indoor"],
    }
elif livestock == "BUFFALO_DAIRY":
    mmscats = {
        "liquid":["MMS_indoor","MMS_open"],
        "solid":["MMS_indoor"],
    }

sim = 'base'

#sens_vars = ['temp','sw','srf','sub','Napp','ph']
# sens_vars = ['temp']
#sens_tests=['+','-']
# sens_tests=['+']

band = int(360/size)
arr = (730,band,720)


start = time.time()
print("Year of study:",CONFIG.sim_year)

comm.Barrier()

for mmsphase in ["liquid","solid"]:
    for mmscat in mmscats[mmsphase]:
# for mmsphase in ["liquid"]:
#     for mmscat in ["MMS_indoor"]:
        for method_used in methods:
            print("=======================================================")
            manurefertapp_sim = LAND.LAND_module(prank=rank,psize=band,
    fert_type='manure',manure_added=np.zeros(arr),urea_added=np.zeros(arr),
    UA_added=np.zeros(arr),avail_N_added=np.zeros(arr),resist_N_added=np.zeros(arr),
    unavail_N_added=np.zeros(arr),
    TAN_added=np.zeros(arr),NO3_added=np.zeros(arr),water_added=np.zeros(arr),
    pH_value=LAND.pH_info[livestock])
            manurefertapp_sim.main(techs[method_used],None,None,
                        startday,730,manure_types[mmsphase],
                        livestock,prodsyst,
                        mmscat,mmsphase,
                        sim_type='base',
                        output_stat=True,quality_check=True)
            mid = time.time()
            cal_time = mid - start
            print("=======================================================")
            print("Rank: ",rank)
            print('sim time: ', np.round(cal_time/60,decimals=2),' mins')

            comm.Barrier()

            newshape = int(365*720)
            sendNH3flux = np.transpose(manurefertapp_sim.o_NH3flux,(1,0,2)).reshape((band,newshape))
            sendwashoff = np.transpose(manurefertapp_sim.o_washoff,(1,0,2)).reshape((band,newshape))
            sendnitrif = np.transpose(manurefertapp_sim.o_nitrif,(1,0,2)).reshape((band,newshape))
            sendNH4leaching = np.transpose(manurefertapp_sim.o_NH4leaching,(1,0,2)).reshape((band,newshape))
            senddiffaq = np.transpose(manurefertapp_sim.o_diffaq,(1,0,2)).reshape((band,newshape))
            senddiffgas = np.transpose(manurefertapp_sim.o_diffgas,(1,0,2)).reshape((band,newshape))
            sendammNuptake = np.transpose(manurefertapp_sim.o_ammNuptake,(1,0,2)).reshape((band,newshape))
            sendnitNuptake = np.transpose(manurefertapp_sim.o_nitNuptake,(1,0,2)).reshape((band,newshape))
            sendNO3washoff = np.transpose(manurefertapp_sim.o_NO3washoff,(1,0,2)).reshape((band,newshape))
            sendNO3leaching = np.transpose(manurefertapp_sim.o_NO3leaching,(1,0,2)).reshape((band,newshape))
            sendNO3diff = np.transpose(manurefertapp_sim.o_NO3diff,(1,0,2)).reshape((band,newshape))

            NH3flux = None
            Nwashoff = None
            TANnitrif = None
            NH4leaching = None
            TANdiff = None
            NH3diff = None
            Ammuptake = None
            Nituptake = None
            NO3washoff = None
            NO3leaching = None
            NO3diff = None

            if rank == 0:
                NH3flux = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                Nwashoff = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                TANnitrif = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                NH4leaching = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                TANdiff = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                NH3diff = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                Ammuptake = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                Nituptake = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                NO3washoff = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                NO3leaching = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))
                NO3diff = np.transpose(np.zeros((365,360,720), dtype=np.float64),(1,0,2)).reshape((360,newshape))

            comm.Barrier()

            print("Rank", rank, "is sending data")
            comm.Gather(sendNH3flux, NH3flux, root=0)
            print("Rank", rank, "sent data 1")
            comm.Gather(sendwashoff, Nwashoff, root=0)
            print("Rank", rank, "sent data 2")
            comm.Gather(sendnitrif, TANnitrif, root=0)
            print("Rank", rank, "sent data 3")
            comm.Gather(sendNH4leaching, NH4leaching, root=0)
            print("Rank", rank, "sent data 4")
            comm.Gather(senddiffaq, TANdiff, root=0)
            print("Rank", rank, "sent data 5")
            comm.Gather(senddiffgas, NH3diff, root=0)
            print("Rank", rank, "sent data 6")
            comm.Gather(sendammNuptake, Ammuptake, root=0)
            print("Rank", rank, "sent data 7")
            comm.Gather(sendnitNuptake, Nituptake, root=0)
            print("Rank", rank, "sent data 8")
            comm.Gather(sendNO3washoff, NO3washoff, root=0)
            print("Rank", rank, "sent data 9")
            comm.Gather(sendNO3leaching, NO3leaching, root=0)
            print("Rank", rank, "sent data 10")
            comm.Gather(sendNO3diff, NO3diff, root=0)
            print("Rank", rank, "sent data 11")

            comm.Barrier()

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
                        NH3emiss=(['time','lat','lon'],np.transpose(NH3flux.reshape((360,365,720)),(1,0,2))),
                        TANwashoff=(['time','lat','lon'],np.transpose(Nwashoff.reshape((360,365,720)),(1,0,2))),
                        TANleaching=(['time','lat','lon'],np.transpose(NH4leaching.reshape((360,365,720)),(1,0,2))),
                        TANdiffaq=(['time','lat','lon'],np.transpose(TANdiff.reshape((360,365,720)),(1,0,2))),
                        NH3diffgas=(['time','lat','lon'],np.transpose(NH3diff.reshape((360,365,720)),(1,0,2))),
                        NH4nitrif=(['time','lat','lon'],np.transpose(TANnitrif.reshape((360,365,720)),(1,0,2))),
                        NH4uptake=(['time','lat','lon'],np.transpose(Ammuptake.reshape((360,365,720)),(1,0,2))),
                        NO3uptake=(['time','lat','lon'],np.transpose(Nituptake.reshape((360,365,720)),(1,0,2))),
                        NO3washoff=(['time','lat','lon'],np.transpose(NO3washoff.reshape((360,365,720)),(1,0,2))),
                        NO3leaching=(['time','lat','lon'],np.transpose(NO3leaching.reshape((360,365,720)),(1,0,2))),
                        NO3diff=(['time','lat','lon'],np.transpose(NO3diff.reshape((360,365,720)),(1,0,2))),
                                ),
                    coords = dict(
                        time=(["time"], times),
                        lon=(["lon"], lons),
                        lat=(["lat"], lats),
                                ),
                    attrs=dict(
                        description="AMCLIM-Land_manure_fert: \
                            N pathways of manure fertilizer application in " +str(CONFIG.sim_year),
                        units="gN per grid",
                    ),
                )

                outds.NH3emiss.attrs["unit"] = 'gN/day'
                outds.NH3emiss.attrs["long name"] = 'NH3 emission from manure application'
                outds.TANwashoff.attrs["unit"] = 'gN/day'
                outds.TANwashoff.attrs["long name"] = 'TAN washoff from manure application'
                outds.TANleaching.attrs["unit"] = 'gN/day'
                outds.TANleaching.attrs["long name"] = 'TAN leaching from  manure application'
                outds.TANdiffaq.attrs["unit"] = 'gN/day'
                outds.TANdiffaq.attrs["long name"] = 'TAN diffusion to deep soil from manure application'
                outds.NH3diffgas.attrs["unit"] = 'gN/day'
                outds.NH3diffgas.attrs["long name"] = 'NH3 diffusion to deep soil from manure application'
                outds.NH4nitrif.attrs["unit"] = 'gN/day'
                outds.NH4nitrif.attrs["long name"] = 'Nitrification'
                outds.NH4uptake.attrs["unit"] = 'gN/day'
                outds.NH4uptake.attrs["long name"] = 'Uptake of NH4+ by crops'
                outds.NO3uptake.attrs["unit"] = 'gN/day'
                outds.NO3uptake.attrs["long name"] = 'Uptake of NO3- by crops'
                outds.NO3washoff.attrs["unit"] = 'gN/day'
                outds.NO3washoff.attrs["long name"] = 'Nitrate washoff from manure application'
                outds.NO3leaching.attrs["unit"] = 'gN/day'
                outds.NO3leaching.attrs["long name"] = 'Nitrate leaching from manure application'
                outds.NO3diff.attrs["unit"] = 'gN/day'
                outds.NO3diff.attrs["long name"] = 'Nitrate diffusion to deep soil from manure application'

                comp = dict(zlib=True, complevel=9)
                encoding = {var: comp for var in outds.data_vars}

                print("xarray saved.")

                filename = sim+'.'+str(prodsyst)+'.'+\
                    str(livestock)+'.'+str(mmscat)+'.'+str(mmsphase)+'.'+str(techs[method_used])+'.'+str(CONFIG.sim_year)+'.nc'
                
                full_path = CONFIG.output_path+filename
                print("ncfile",filename)
                print("path",CONFIG.output_path)
                print("full path",CONFIG.output_path+filename)
                outds.to_netcdf(path=str(full_path),mode="w",encoding=encoding)
                print("ncfile saved.")

            comm.Barrier()

end = time.time()
runtime = end - start
print(np.round(runtime/60,decimals = 1),' mins')
