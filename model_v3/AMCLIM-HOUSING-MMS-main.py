import numpy as np
import xarray as xr
import pandas as pd
import time
import sys

import MODULES.FUNC as FUNC
import INPUT.input as INPUT
import CONFIG.config as CONFIG
import MODULES.HOUSING_v3 as HOUSING
import MODULES.MMS as MMS

housing_envs = ["insulated","naturally ventilated"]
housing_types = ["slat/pit house","barn","poultry_house"]

livestock_list = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE',
                    'PIG',
                    'SHEEP',
                    'POULTRY']

livestock = sys.argv[1]
print(livestock)

##############################
## BEEF CATTLE: Housing & MMS
##############################
if livestock == "BEEF_CATTLE":
    ## BEEF: mixed production system (only barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['BEEF_CATTLE'].index("mixed")
    print("=============================================")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="BEEF_CATTLE",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],ncfile_o=True)
    print("=============================================")
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
                water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")

###############################
## DAIRY Cattle: Housing & MMS
###############################
elif livestock == "DAIRY_CATTLE":
    ## DAIRY: mixed production system (barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['DAIRY_CATTLE'].index("mixed")
    print("=============================================")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="DAIRY_CATTLE",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],
                                    ncfile_o=True)
    print("=============================================")
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
                water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")
                
###############################
## other DAIRY: Housing & MMS
###############################
elif livestock == "OTHER_CATTLE":
    ## other Cattle: mixed production system (barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['OTHER_CATTLE'].index("mixed")
    print("=============================================")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="OTHER_CATTLE",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],
                                    ncfile_o=True)
    print("=============================================")
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
                water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")

#################################
## feedlot CATTLE: Housing & MMS
#################################
elif livestock == "FEEDLOT_CATTLE":
    ## BEEF: mixed production system (only barn)
    prodsyst_idx = CONFIG.CONFIG_production_system_dict['FEEDLOT_CATTLE'].index("feedlot")
    print("=============================================")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="FEEDLOT_CATTLE",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],ncfile_o=True)
    print("=============================================")

    for mmsphase in ["liquid","solid"]:
        if mmsphase == "liquid":
            mmscats = ['MMS_indoor',"MMS_open","MMS_cover"]
        elif mmsphase == "solid":
            mmscats = ['MMS_indoor',"MMS_open"]
        for mmscat in mmscats: 
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                resist_N_added=housing_manure.resist_N_pool_to_storage,
                unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                TAN_added=housing_manure.TAN_pool_to_storage,
                water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
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
    print("=============================================")
    housing_manure = HOUSING.HOUSING_MODULE(livestock_name="SHEEP",
                                            production_system_lvl_idx=prodsyst_idx,
                                            housing_type="barn")
    housing_manure.barn_sim_main(housing_type="barn",
                                    start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],ncfile_o=True)
    print("=============================================")
    for mmsphase in ["liquid","solid"]:
        if mmsphase == "liquid":
            mmscats = ['MMS_indoor',"MMS_open","MMS_cover"]
        elif mmsphase == "solid":
            mmscats = ['MMS_indoor',"MMS_open"]
        for mmscat in mmscats: 
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                resist_N_added=housing_manure.resist_N_pool_to_storage,
                unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                TAN_added=housing_manure.TAN_pool_to_storage,
                water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
            if mmsphase == "solid":
                if mmscat == "MMS_open":
                    mms_manure.output_MMS_pathway(output_stat=True)
            print("=============================================")

#######################
## PIG: HOUSING & MMS
#######################
elif livestock == "PIG":
    ## PIG: industrial, intermediate, backyard (slat-pit_house, barn)
    for prodsyst_idx in np.arange(3):
    # for prodsyst_idx in np.arange(1):
        if prodsyst_idx == 0:
            housetype = "slat-pit_house"
            print("=============================================")
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="PIG",
                                                    production_system_lvl_idx=prodsyst_idx,
                                                    housing_type=housetype)
            housing_manure.slat_pit_sim_main(housing_type=housetype,
                            start_idx=0,end_idx=CONFIG.Days*2,cleaning_frequency=HOUSING.cleaning_freq["slat-pit_house_insitu"],ncfile_o=True)
            print("=============================================")
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                        mms_cat="MMS_open",phase="solid",manure_added=housing_manure.manure_pool_to_storage,
                        urea_added=housing_manure.urea_pool_to_storage,UA_added=np.zeros(CONFIG.CONFIG_mtrx),avail_N_added=housing_manure.avail_N_pool_to_storage,
                        resist_N_added=housing_manure.resist_N_pool_to_storage,
                        unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                        TAN_added=housing_manure.TAN_pool_to_storage,
                        water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat="MMS_open",phase="solid",start_day_idx=0,end_day_idx=CONFIG.Days*2,from_insitu=True,stat=True)
            
            print("=============================================")
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="PIG",
                                                    production_system_lvl_idx=prodsyst_idx,
                                                    housing_type=housetype)
            housing_manure.slat_pit_sim_main(housing_type=housetype,
                            start_idx=0,end_idx=CONFIG.Days*2,cleaning_frequency=HOUSING.cleaning_freq["slat-pit_house"],ncfile_o=True)
            print("=============================================")
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
                        water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
                    mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                    if mmsphase == "solid":
                        if mmscat == "MMS_open":
                            mms_manure.output_MMS_pathway(output_stat=True)
                    print("=============================================")
        else:
            housetype = "barn"
            print("=============================================")
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="PIG",
                                                    production_system_lvl_idx=prodsyst_idx,
                                                    housing_type=housetype)
            housing_manure.barn_sim_main(housing_type=housetype,
                                            start_idx=0,end_idx=CONFIG.Days,cleaning_frequency=HOUSING.cleaning_freq["barn"],
                                            ncfile_o=True)
            print("=============================================")
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
                        water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
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
    # for prodsyst_idx in np.arange(1,2):
        ## broiler 0 and layer 1 (no backyard chicken 2)
        if prodsyst_idx != 2:
            print("=============================================")
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="POULTRY",
                                                        production_system_lvl_idx=prodsyst_idx,
                                                        housing_type="poultry_house")

            housing_manure.poultry_house_sim_main(housing_type="poultry_house",
                                            start_idx=0,end_idx=CONFIG.Days*2,cleaning_frequency=HOUSING.cleaning_freq["poultry_house_with_litter"],
                                            litter=True)
            print("=============================================")
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                            mms_cat="MMS_open",phase="solid",manure_added=housing_manure.manure_pool_to_storage,
                            urea_added=np.zeros(CONFIG.CONFIG_mtrx),UA_added=housing_manure.UA_pool_to_storage,avail_N_added=housing_manure.avail_N_pool_to_storage,
                            resist_N_added=housing_manure.resist_N_pool_to_storage,
                            unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                            TAN_added=housing_manure.TAN_pool_to_storage,
                            water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            mms_manure.MMS_sim_main(mms_cat="MMS_open",phase="solid",start_day_idx=0,end_day_idx=CONFIG.Days*2,from_insitu=True,stat=True)
            print("=============================================")

            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="POULTRY",
                                                        production_system_lvl_idx=prodsyst_idx,
                                                        housing_type="poultry_house")

            housing_manure.poultry_house_sim_main(housing_type="poultry_house",
                                            start_idx=0,end_idx=CONFIG.Days*2,cleaning_frequency=HOUSING.cleaning_freq["poultry_house"],
                                            litter=False)
            print("=============================================")

            
            for mmsphase in ["liquid","solid"]:
                if mmsphase == "liquid":
                    mmscats = ['MMS_indoor',"MMS_open","MMS_cover"]
                elif mmsphase == "solid":
                    mmscats = ['MMS_indoor',"MMS_open"]
                for mmscat in mmscats: 
                    mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                        mms_cat=mmscat,phase=mmsphase,manure_added=housing_manure.manure_pool_to_storage,
                        urea_added=np.zeros(CONFIG.CONFIG_mtrx),UA_added=housing_manure.UA_pool_to_storage,avail_N_added=housing_manure.avail_N_pool_to_storage,
                        resist_N_added=housing_manure.resist_N_pool_to_storage,
                        unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                        TAN_added=housing_manure.TAN_pool_to_storage,
                        water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
                    mms_manure.MMS_sim_main(mms_cat=mmscat,phase=mmsphase,start_day_idx=0,end_day_idx=CONFIG.Days*2,stat=True)
                    if mmsphase == "solid":
                        if mmscat == "MMS_open":
                            mms_manure.output_MMS_pathway(output_stat=True)
                    print("=============================================")

        ## backyard poultry
        elif prodsyst_idx == 2:
            housing_manure = HOUSING.HOUSING_MODULE(livestock_name="POULTRY",
                                                        production_system_lvl_idx=prodsyst_idx,
                                                        housing_type="poultry_house")
            housing_manure.poultry_house_sim_main(housing_type="poultry_house",
                                            start_idx=0,end_idx=CONFIG.Days*2,cleaning_frequency=HOUSING.cleaning_freq["poultry_house"],
                                            litter=False)
            print("=============================================")
            mms_manure = MMS.MMS_module(livestock_name=livestock,production_system_lvl_idx=prodsyst_idx,
                            mms_cat="MMS_open",phase="solid",manure_added=housing_manure.manure_pool_to_storage,
                            urea_added=np.zeros(CONFIG.CONFIG_mtrx),UA_added=housing_manure.UA_pool_to_storage,avail_N_added=housing_manure.avail_N_pool_to_storage,
                            resist_N_added=housing_manure.resist_N_pool_to_storage,
                            unavail_N_added=housing_manure.unavail_N_pool_to_storage,
                            TAN_added=housing_manure.TAN_pool_to_storage,
                            water_added=housing_manure.Total_water_pool_to_storage,area_housing=housing_manure.floor_area)
            ## this is not truly from "is-situ", but from_insitu is set to be True to modify the area of field on which manure is applied
            mms_manure.MMS_sim_main(mms_cat="MMS_open",phase="solid",start_day_idx=0,end_day_idx=CONFIG.Days*2,from_insitu=True,stat=True)
            mms_manure.output_MMS_pathway(output_stat=True)