import numpy as np
import xarray as xr
import pandas as pd
import time
import sys

import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.HOUSING_v3 as HOUSING

housing_envs = ["insulated","naturally ventilated"]
housing_types = ["slat/pit house","barn","poultry_house"]

livestock_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE',
                    'PIG',,
                    'SHEEP',
                    'POULTRY']

livestock = sys.argv[1]

# lvl_idx = 1
# type_idx = 1
# barn_cf = 1
# slatpit_cf = 180
cleaning_freq = {"slat/pit house": 180,
                "barn": 1,
                "poultry_house": 60}



##############################
## BEEF CATTLE: Housing & MMS
##############################
if livestock == "BEEF_CATTLE":
    ## BEEF: mixed production system (only barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['BEEF_CATTLE'].index("mixed")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="BEEF_CATTLE",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq["barn"],ncfile_o=True)

    for mmsphase in ["liquid","solid"]:
        if mmsphase == "liquid":
            mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
        elif mmsphase == "solid":
            mmscats = ['MMS_indoor',"MMS_open"]
        for mmscat in mmscats: 
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                resist_N_added=housing_manure.resist_N_pool_to_storage,
                unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                TAN_added=housing_manure.TAN_pool_to_storage,
                water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")

###############################
## DAIRY Cattle: Housing & MMS
###############################
elif livestock == "DAIRY_CATTLE":
    ## DAIRY: mixed production system (slat-pit house, barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['DAIRY_CATTLE'].index("mixed")
    for housetype in ["slat/pit house","barn"]:
        housing_manure = HOUSING.HOUSING_MODULE(livestock_name="DAIRY_CATTLE",
                                                production_system_lvl_idx=prodsyst_idx,
                                                housing_type=housetype)
        if housetype == "slat/pit house":
            housing_manure.slat_pit_sim_main(housing_type=housetype,
                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],ncfile_o=True)
        elif housetype == "barn":
            housing_manure.barn_sim_main(housing_type=housetype,
                                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],
                                            ncfile_o=True)

        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats: 
                mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                    mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                    urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                    resist_N_added=housing_manure.resist_N_pool_to_storage,
                    unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                    TAN_added=housing_manure.TAN_pool_to_storage,
                    water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
                mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                if mmsphase == "solid":
                    if mmscat == "MMS_open":
                        mms_manure.output_MMS_pathway(output_stat=True)
                print("=============================================")
                
###############################
## other CATTLE: Housing & MMS
###############################
elif livestock == "OTHER_CATTLE":
    ## other Cattle: mixed production system (slat-pit house, barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['OTHER_CATTLE'].index("mixed")
    for housetype in ["slat/pit house","barn"]:
        housing_manure = HOUSING.HOUSING_MODULE(livestock_name="OTHER_CATTLE",
                                                production_system_lvl_idx=prodsyst_idx,
                                                housing_type=housetype)
        if housetype == "slat/pit house":
            housing_manure.slat_pit_sim_main(housing_type=housetype,
                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],ncfile_o=True)
        elif housetype == "barn":
            housing_manure.barn_sim_main(housing_type=housetype,
                                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],
                                            ncfile_o=True)

        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats: 
                mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                    mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                    urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                    resist_N_added=housing_manure.resist_N_pool_to_storage,
                    unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                    TAN_added=housing_manure.TAN_pool_to_storage,
                    water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
                mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                if mmsphase == "solid":
                    if mmscat == "MMS_open":
                        mms_manure.output_MMS_pathway(output_stat=True)
                print("=============================================")

########################
## SHEEP: Housing & MMS
########################
elif livestock == "SHEEP":
    ## SHEEP: mixed production system (only barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['SHEEP'].index("mixed")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="SHEEP",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq["barn"],ncfile_o=True)

    for mmsphase in ["liquid","solid"]:
        if mmsphase == "liquid":
            mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
        elif mmsphase == "solid":
            mmscats = ['MMS_indoor',"MMS_open"]
        for mmscat in mmscats: 
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                resist_N_added=housing_manure.resist_N_pool_to_storage,
                unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                TAN_added=housing_manure.TAN_pool_to_storage,
                water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")

#######################
## PIG: HOUSING & MMS
#######################
elif livestock == "PIG":
    ## PIG: industrial, intermediate, backyard (slat-pit house, barn)
    for prodsyst_idx in np.arange(3):
        for housetype in ["slat/pit house","barn"]:
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="OTHER_CATTLE",
                                                    production_system_lvl_idx=prodsyst_idx,
                                                    housing_type=housetype)
            if housetype == "slat/pit house":
                housing_manure.slat_pit_sim_main(housing_type=housetype,
                                start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],ncfile_o=True)
            elif housetype == "barn":
                housing_manure.barn_sim_main(housing_type=housetype,
                                                start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],
                                                ncfile_o=True)
            for mmsphase in ["liquid","solid"]:
                if mmsphase == "liquid":
                    mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
                elif mmsphase == "solid":
                    mmscats = ['MMS_indoor',"MMS_open"]
                for mmscat in mmscats: 
                    mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                        mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                        urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                        resist_N_added=housing_manure.resist_N_pool_to_storage,
                        unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                        TAN_added=housing_manure.TAN_pool_to_storage,
                        water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
                    mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                    if mmsphase == "solid":
                        if mmscat == "MMS_open":
                            mms_manure.output_MMS_pathway(output_stat=True)
                    print("=============================================")


###########################
## POULTRY: HOUSING & MMS
###########################
elif livestock == "POULTRY":
    ## POULTRY (chicken): broiler, layer, backyard (poultry house)
    for prodsyst_idx in np.arange(3):
        housing_manure = HOUSING.HOUSING_MODULE(livestock_name="POULTRY",
                                                    production_system_lvl_idx=prodsyst_idx,
                                                    housing_type=housetype)
        if prodsyst_idx == 2:
            housing_manure.poultry_house_sim_main(housing_type="poultry_house",
                                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],
                                            insitu_storage=False)
        else:
            housing_manure.poultry_house_sim_main(housing_type="poultry_house",
                                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=cleaning_freq[housetype],
                                            insitu_storage=True)

        for mmsphase in ["liquid","solid"]:
            if mmsphase == "liquid":
                mmscats = ['MMS_indoor',"MMS_open","MMS_cover","MMS_lagoon"]
            elif mmsphase == "solid":
                mmscats = ['MMS_indoor',"MMS_open"]
            for mmscat in mmscats: 
                mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                    mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                    urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                    resist_N_added=housing_manure.resist_N_pool_to_storage,
                    unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                    TAN_added=housing_manure.TAN_pool_to_storage,
                    water_added=housing_manure.Total_water_pool_to_storage,pH_value=HOUSING.pH,area_housing=housing_manure.floor_area)
                mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                if mmsphase == "solid":
                    if mmscat == "MMS_open":
                        mms_manure.output_MMS_pathway(output_stat=True)
                print("=============================================")