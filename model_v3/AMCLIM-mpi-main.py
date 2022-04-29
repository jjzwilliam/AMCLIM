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


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()


crops = [sys.argv[1]]
chem_ferts = ['ammonium','urea']
startday = 0
endday = 730
methods = ['broadcasting-surf','incorporating-disk','deep injection']
# methods = ['broadcasting-surf']

techs = {
    "broadcasting-surf": 'surf',
    "incorporating-disk": 'disk',
    "deep injection": 'injection'
}
sim = 'base'

#sens_vars = ['temp','sw','srf','sub','Napp','ph']
# sens_vars = ['temp']
#sens_tests=['+','-']
# sens_tests=['+']

band = int(360/size)
arr = (731,band,720)


start = time.time()


for chem_fert_use in chem_ferts:
    for crop_item in crops:
        for method_used in methods:
            test_chemfert_sim1 = LAND.LAND_module(prank=rank,psize=band,
fert_type='mineral',manure_added=np.zeros(arr),urea_added=np.zeros(arr),
UA_added=np.zeros(arr),avail_N_added=np.zeros(arr),resist_N_added=np.zeros(arr),
unavail_N_added=np.zeros(arr),
TAN_added=np.zeros(arr),NO3_added=np.zeros(arr),water_added=np.zeros(arr),
pH_value=7.0)
            test_chemfert_sim1.main(techs[method_used],crop_item,chem_fert_use,startday,endday,
sim_type='base',output_stat=True,quality_check=True)
            mid = time.time()
            cal_time = mid - start
            print("=======================================================")
            print("Rank: ",rank)
            print('sim time: ', np.round(cal_time/60,decimals=2),' mins')

            newshape = int(365*720)
            sendNH3flux = np.transpose(test_chemfert_sim1.o_NH3flux,(1,0,2)).reshape((band,newshape))
            sendwashoff = np.transpose(test_chemfert_sim1.o_washoff,(1,0,2)).reshape((band,newshape))
            sendnitrif = np.transpose(test_chemfert_sim1.o_nitrif,(1,0,2)).reshape((band,newshape))
            sendNH4leaching = np.transpose(test_chemfert_sim1.o_NH4leaching,(1,0,2)).reshape((band,newshape))
            senddiffaq = np.transpose(test_chemfert_sim1.o_diffaq,(1,0,2)).reshape((band,newshape))
            senddiffgas = np.transpose(test_chemfert_sim1.o_diffgas,(1,0,2)).reshape((band,newshape))
            sendammNuptake = np.transpose(test_chemfert_sim1.o_ammNuptake,(1,0,2)).reshape((band,newshape))
            sendnitNuptake = np.transpose(test_chemfert_sim1.o_nitNuptake,(1,0,2)).reshape((band,newshape))
            sendNO3washoff = np.transpose(test_chemfert_sim1.o_NO3washoff,(1,0,2)).reshape((band,newshape))
            sendNO3leaching = np.transpose(test_chemfert_sim1.o_NO3leaching,(1,0,2)).reshape((band,newshape))
            sendNO3diff = np.transpose(test_chemfert_sim1.o_NO3diff,(1,0,2)).reshape((band,newshape))

            NH3flux = None
            TANwashoff = None
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
                NH3flux = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                TANwashoff = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                TANnitrif = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                NH4leaching = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                TANdiff = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                NH3diff = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                Ammuptake = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                Nituptake = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                NO3washoff = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                NO3leaching = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))
                NO3diff = np.transpose(np.zeros((365,360,720), dtype=np.float),(1,0,2)).reshape((360,newshape))

            comm.Gather(sendNH3flux, NH3flux, root=0)
            comm.Gather(sendwashoff, TANwashoff, root=0)
            comm.Gather(sendnitrif, TANnitrif, root=0)
            comm.Gather(sendNH4leaching, NH4leaching, root=0)
            comm.Gather(senddiffaq, TANdiff, root=0)
            comm.Gather(senddiffgas, NH3diff, root=0)
            comm.Gather(sendammNuptake, Ammuptake, root=0)
            comm.Gather(sendnitNuptake, Nituptake, root=0)
            comm.Gather(sendNO3washoff, NO3washoff, root=0)
            comm.Gather(sendNO3leaching, NO3leaching, root=0)
            comm.Gather(sendNO3diff, NO3diff, root=0)

            if rank == 0:
                nlat = int(180.0/CONFIG.dlat)
                nlon = int(360.0/CONFIG.dlon)
                ntime = CONFIG.Days
                lats = 90 - 0.5*np.arange(nlat)
                lons = -180 + 0.5*np.arange(nlon)
                yearidx = str(CONFIG.sim_year)+'-01-01'
                times = pd.date_range(yearidx,periods=ntime)

                outds = xr.Dataset(
                    data_vars=dict(
                        NH3emiss=(['time','lat','lon'],np.transpose(NH3flux.reshape((360,365,720)),(1,0,2))),
                        TANwashoff=(['time','lat','lon'],np.transpose(TANwashoff.reshape((360,365,720)),(1,0,2))),
                        TANleaching=(['time','lat','lon'],np.transpose(TANnitrif.reshape((360,365,720)),(1,0,2))),
                        TANdiffaq=(['time','lat','lon'],np.transpose(NH4leaching.reshape((360,365,720)),(1,0,2))),
                        NH3diffgas=(['time','lat','lon'],np.transpose(TANdiff.reshape((360,365,720)),(1,0,2))),
                        NH4nitrif=(['time','lat','lon'],np.transpose(NH3diff.reshape((360,365,720)),(1,0,2))),
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
                        description="AMCLIM-Land_chem_fert: \
                            N pathways of chemical fertilizer application in " +str(CONFIG.sim_year),
                        info = method_used+" of "+chem_fert_use+" for: "+str(crop_item),
                        units="gN per grid",
                    ),
                )

                outds.NH3emiss.attrs["unit"] = 'gN/day'
                outds.NH3emiss.attrs["long name"] = 'NH3 emission from fertilizer application'
                outds.TANwashoff.attrs["unit"] = 'gN/year'
                outds.TANwashoff.attrs["long name"] = 'TAN washoff from fertilizer application'
                outds.TANleaching.attrs["unit"] = 'gN/year'
                outds.TANleaching.attrs["long name"] = 'TAN leaching from fertilizer application'
                outds.TANdiffaq.attrs["unit"] = 'gN/year'
                outds.TANdiffaq.attrs["long name"] = 'TAN diffusion to deep soil from fertilizer application'
                outds.NH3diffgas.attrs["unit"] = 'gN/year'
                outds.NH3diffgas.attrs["long name"] = 'NH3 diffusion to deep soil from fertilizer application'
                outds.NH4nitrif.attrs["unit"] = 'gN/year'
                outds.NH4nitrif.attrs["long name"] = 'Nitrification'
                outds.NH4uptake.attrs["unit"] = 'gN/year'
                outds.NH4uptake.attrs["long name"] = 'Uptake of NH4+ by crops'
                outds.NO3uptake.attrs["unit"] = 'gN/year'
                outds.NO3uptake.attrs["long name"] = 'Uptake of NO3- by crops'
                outds.NO3washoff.attrs["unit"] = 'gN/year'
                outds.NO3washoff.attrs["long name"] = 'Nitrate washoff from fertilizer application'
                outds.NO3leaching.attrs["unit"] = 'gN/year'
                outds.NO3leaching.attrs["long name"] = 'Nitrate leaching from fertilizer application'
                outds.NO3diff.attrs["unit"] = 'gN/year'
                outds.NO3diff.attrs["long name"] = 'Nitrate diffusion to deep soil from fertilizer application'

                comp = dict(zlib=True, complevel=9)
                encoding = {var: comp for var in outds.data_vars}

                outds.to_netcdf(CONFIG.output_path+sim+'.'+str(crop_item)+'.'+\
                    str(chem_fert_use)+'.'+str(techs[method_used])+'.'+str(CONFIG.sim_year)+'.nc',encoding=encoding)
                print("ncfile saved.")

end = time.time()
runtime = end - start
print(np.round(runtime/60,decimals = 1),' mins')
