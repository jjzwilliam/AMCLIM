from mpi4py import MPI
import numpy as np
import xarray as xr
import pandas as pd
import time


import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.LAND_develop as LAND


comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

crops = ['wheat']
chem_ferts = ['ammonium']
startday = 0
endday = 730
# methods = ['broadcasting-surf','incorporating-disk','deep injection']
methods = ['broadcasting-surf']

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

band = int(320/size)
arr = (731,band,720)


start = time.time()


for chem_fert_use in chem_ferts:
    for crop_item in crops:
        for method_used in methods:
            test_chemfert_sim1 = LAND.LAND_module(array_shape=arr,prank=rank,psize=band,
fert_type='mineral',manure_added=np.zeros(arr),urea_added=np.zeros(arr),
UA_added=np.zeros(arr),avail_N_added=np.zeros(arr),resist_N_added=np.zeros(arr),
unavail_N_added=np.zeros(arr),
TAN_added=np.zeros(arr),NO3_added=np.zeros(arr),water_added=np.zeros(arr),
pH_value=7.0)
            test_chemfert_sim1.main(techs[method_used],crop_item,chem_fert_use,startday,endday,sim_type=sim,
sim_stat=True,ncfile_o=False,quality_check=True)
            mid = time.time()
            cal_time = mid - start
            print("=======================================================")
            print("Rank: ",rank)
            print('sim time: ', np.round(cal_time/60,decimals=2),' mins')
'''            if rank == 0:
                NH3flux = np.zeros((365,320,720), dtype=np.float)
                TANwashoff = np.zeros((365,320,720), dtype=np.float)
                TANnitrif = np.zeros((365,320,720), dtype=np.float)
                NH4leaching = np.zeros((365,320,720), dtype=np.float)
                TANdiff = np.zeros((365,320,720), dtype=np.float)
                NH3diff = np.zeros((365,320,720), dtype=np.float)
                Ammuptake = np.zeros((365,320,720), dtype=np.float)
                Nituptake = np.zeros((365,320,720), dtype=np.float)
                NO3washoff = np.zeros((365,320,720), dtype=np.float)
                NO3leaching = np.zeros((365,320,720), dtype=np.float)
                NO3diff = np.zeros((365,320,720), dtype=np.float)

            for dd in np.arange(365):
                dNH3flux = None
                dTANwashoff = None
                dTANnitrif = None
                dNH4leaching = None
                dTANdiff = None
                dNH3diff = None
                dAmmuptake = None
                dNituptake = None
                dNO3washoff = None
                dNO3leaching = None
                dNO3diff = None

                if rank == 0:
                    dNH3flux = np.zeros((320,720), dtype=np.float)
                    dTANwashoff = np.zeros((320,720), dtype=np.float)
                    dTANnitrif = np.zeros((320,720), dtype=np.float)
                    dNH4leaching = np.zeros((320,720), dtype=np.float)
                    dTANdiff = np.zeros((320,720), dtype=np.float)
                    dNH3diff = np.zeros((320,720), dtype=np.float)
                    dAmmuptake = np.zeros((320,720), dtype=np.float)
                    dNituptake = np.zeros((320,720), dtype=np.float)
                    dNO3washoff = np.zeros((320,720), dtype=np.float)
                    dNO3leaching = np.zeros((320,720), dtype=np.float)
                    dNO3diff = np.zeros((320,720), dtype=np.float)

                comm.Gather(test_chemfert_sim1.o_NH3flux[dd], dNH3flux, root=0)
                comm.Gather(test_chemfert_sim1.o_washoff[dd], dTANwashoff, root=0)
                comm.Gather(test_chemfert_sim1.o_nitrif[dd], dTANnitrif, root=0)
                comm.Gather(test_chemfert_sim1.o_NH4leaching[dd], dNH4leaching, root=0)
                comm.Gather(test_chemfert_sim1.o_diffaq[dd], dTANdiff, root=0)
                comm.Gather(test_chemfert_sim1.o_diffgas[dd], dNH3diff, root=0)
                comm.Gather(test_chemfert_sim1.o_ammNuptake[dd], dAmmuptake, root=0)
                comm.Gather(test_chemfert_sim1.o_nitNuptake[dd], dNituptake, root=0)
                comm.Gather(test_chemfert_sim1.o_NO3washoff[dd], dNO3washoff, root=0)
                comm.Gather(test_chemfert_sim1.o_NO3leaching[dd], dNO3leaching, root=0)
                comm.Gather(test_chemfert_sim1.o_NO3diff[dd], dNO3diff, root=0)

                if rank == 0:
                    NH3flux[dd] = dNH3flux
                    TANwashoff[dd] = dTANwashoff
                    TANnitrif[dd] = dTANnitrif
                    NH4leaching[dd] = dNH4leaching
                    TANdiff[dd] = dTANdiff
                    NH3diff[dd] = dNH3diff
                    Ammuptake[dd] = dAmmuptake
                    Nituptake[dd] = dNituptake
                    NO3washoff[dd] = dNO3washoff
                    NO3leaching[dd] = dNO3leaching
                    NO3diff[dd] = dNO3diff'''

end = time.time()
runtime = end - start
print(np.round(runtime/60,decimals = 1),' mins')
