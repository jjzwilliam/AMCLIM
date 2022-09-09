from ast import Pass
from cgi import test
from cmath import phase, tan
from logging import raiseExceptions
from os import times
from re import sub
from time import time
from wsgiref.simple_server import demo_app

from pandas import array
from INPUT.input import *
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *
import sys

###################################
## Crop info
###################################
## major crops
crop_list = ['barley',
            'cassava',
            'cotton',
            'groundnut',
            'maize',
            'millet',
            'oilpalm',
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
crop_idx = {}

###################################
## LAND parameters
###################################
# ## dry matter (DM) content of solid manure 
# DM_content = solid_m_DM[livestock]
# ## dry matter (DM) content of liquid manure is assumed to be 5%
# f_DM_liquid = 0.1
# ## maximum water content of manure
# f_wcmax = 1 - (DM_content/100)/2
# ## assuming the density of manure; 1t kg/m^3 or 1g/cm^3
# manure_density = rho_m[livestock]
## thickness of the topsoil layer: 7cm; mid point 3.5cm
z_topsoil = 0.07
p_topsoil = 0.035
## thickness of the 2nd soil layer (under topsoil); mid point 17.5cm
z_2ndsoil = 0.21
p_2ndsoil = 0.175
## assuming infiltration of manure water to the soil is 10mm/day (1cm/day; 10 000 g/m^2/day) ref: Vira et al.,2020 GMD (2x d0)
# dailymaxinfil = 0.01
dailymaxinfil = 0.003 ## test: 3mm /day
## assuming the water holding capacity of manure is 3.0 g water/ g manure
absorb_factor = 3.0
## washoff coefficients: 0.1%/mm water for N species, and 0.05%/mm water for non N species (manure)
f_washoff_nonN = 0.0005
f_washoff_N = 0.001
## potential irrigation: 50 mm
irrig_water = 50*1e3
# irrig_water = 0.0
## days for seeds to germinate
day_germinate = 7.0

## thichkness of vertical layer 1, 2, 3, 4: 2cm(1), 5cm(2), 7cm(3), 14cm(4)
zlayers = [0.02, 0.05, 0.07, 0.14]
## midpoints of each layer
pmids = [0.01, 0.045, 0.105, 0.21]

## grazing 
## patch area; 0.25 m2
patch_area = 0.25
## grazing density; 2500 m2/head
grazing_density = {"BEEF_CATTLE":2500,"DAIRY_CATTLE":2500,"SHEEP":400,"GOAT":400} 
## fraction of N_avail, N_resist and N_unavail; 50, 45, 5 per cent
## ref: CLM_FANv1 (Riddick et al., 2016)
f_avail = 0.5
f_resist = 0.45
f_unavail = 0.05
## source layer for NH3 emissions of grazing soils; 4mm (Moring et al., BG 2016)
z_source = 0.004
## fraction of effective source area for NH3 emission; annual urine coverage 17 %, dung coverage 4 %  
f_urine_patch = 0.17
f_dung_pat = 0.04
## fraction of FYM of dung pats
frac_FYM = 0.5

class LAND_module:
    def __init__(self,prank,psize,fert_type,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,NO3_added,water_added,pH_value,grazing=False):
        
        print('LAND Module - current fertilizer application is: '+str(fert_type))

        ## lat idx for parallelization
        self.plat1 = prank*psize
        self.plat2 = (prank+1)*psize

        print(self.plat1,self.plat2)
        
        ## array shape [lats,lons]
        field_shape = (psize,CONFIG_lons)
        ## include the time dimension
        array_shape = (25,) + field_shape
        ## output shape
        outarray_shape = (Days,) + field_shape

        ## shape of 2 levels fields
        dlvl = (2,) + array_shape
        ## shape of 3 levels fields
        tlvl = (3,) + array_shape

        if fert_type == 'manure':
            ## feces input from MMS/HOUSING
            self.manure_added = manure_added
            self.manure = np.zeros(array_shape)
            ## feces pool
            self.manure_pool = np.zeros(array_shape)
            ## washoff flux of manure
            self.manure_washoff = np.zeros(array_shape)
            ## min amount of water in manure (function of T, RH)
            ## equilibrium moisture content, water content cannot be lower than this threshold under "natural" processes,
            ## i.e., without drying processes
            # self.manure_water = np.zeros(array_shape)
            # self.manure_minwc = np.zeros(array_shape)
            ## input of UA from poultry HOUSING/MMS
            self.UA_added = UA_added
            self.UA = np.zeros(array_shape)
            ## input of available N component from MMS/HOUSING
            self.avail_N_added = avail_N_added
            self.avail_N = np.zeros(array_shape)
            ## Nitrogen pool that is available (easily) to form TAN
            self.avail_N_pool = np.zeros(array_shape)
            ## washoff flux of available org N 
            self.avail_N_washoff = np.zeros(array_shape)
            ## input of resistant N component from MMS/HOUSING
            self.resist_N_added = resist_N_added
            self.resist_N = np.zeros(array_shape)
            ## Nitrogen pool that is resistant (slowly) to form TAN
            self.resist_N_pool = np.zeros(array_shape)
            ## washoff flux of resistant org N 
            self.resist_N_washoff = np.zeros(array_shape)
            ## input of unavailable N component from MMS/HOUSING
            self.unavail_N_added = unavail_N_added
            self.unavail_N = np.zeros(array_shape)
            ## Nitrogen pool that is unavailable (cannot) to form TAN
            self.unavail_N_pool = np.zeros(array_shape)
            ## washoff flux of unavailable org N 
            self.unavail_N_washoff = np.zeros(array_shape)
            ## decomposition of org N; 
            self.orgN_decomp = np.zeros(tlvl)
            ## water added from MMS/HOUSING
            self.water_added = water_added
            self.water = np.zeros(array_shape)
            ## total water pool of the system (soil water; manure water+soil water)
            # self.Total_water_pool = np.zeros(array_shape)
            ## TAN pool of the surface slurry
            self.slurry_TAN_pool = np.zeros(array_shape)
            ## TAN concentration of the surface slurry
            self.slurry_TAN_amount = np.zeros(array_shape)
            ## gas concentration of the surface slurry
            self.slurry_NH3_gas = np.zeros(array_shape)
            ## water pool of the surface slurry
            self.slurry_water_pool = np.zeros(array_shape)
            ## TAN infiltration from the surface slurry to the topsoil
            self.slurry_TANinfil = np.zeros(array_shape)
            ## TAN diffusion from the surface slurry to the topsoil
            self.slurry_TANdiffusiondown = np.zeros(array_shape)
            ## slurry pH
            self.slurry_pH = pH_value
            self.slurry_ccH = np.float(10**(-pH_value))
            ## NH3 emissions from grazing FYM
            self.NH3flux_FYM = np.zeros(array_shape)
            ## cropland area
            self.croplandarea = np.zeros(array_shape[1:])
            ## pasture/grassland area
            self.pastarea = np.zeros(array_shape[1:])
            ## area with manure/slurry application
            self.croparea = np.zeros(array_shape[1:])

        else:
            ## cropping area for nitrate N fertilizer
            self.nitN_area = np.zeros(array_shape[1:])
            ## cropping area for ammonium N fertilizer
            self.ammN_area = np.zeros(array_shape[1:])
            ## cropping area for urea N fertilizer
            self.ureaN_area = np.zeros(array_shape[1:])
            ## potential irrigation water
            self.water_added = water_added
            self.water = np.zeros(array_shape)


        ############################
        ## fields with 1 level
        ############################
        #### inputs from fertilizers: urea, ammonium (TAN)
        ## urea input
        self.urea_added = urea_added
        self.urea = np.zeros(array_shape)
        ## TAN added from MMS/HOUSING
        self.TAN_added = TAN_added
        self.TAN = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3 = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3_added = NO3_added
        #### washoff fluxes: urea, TAN, NO3
        ## urea washoff
        self.ureawashoff = np.zeros(array_shape)
        ## TAN washoff 
        self.TANwashoff = np.zeros(array_shape)
        ## NO3 washoff
        self.NO3washoff = np.zeros(array_shape)
        #### volatilization flux: NH3
        ## final NH3 emission
        self.NH3flux = np.zeros(array_shape)

        #### met and soil fields:
        ## atmospheric resistances: 1) aerodynamic resistance, 2) boundary layer resistance
        self.R_atm = np.zeros(array_shape)
        ## evaporation
        self.evap = np.zeros(array_shape)
        ## rain fall (Note the unit)
        self.rainfall = np.zeros(array_shape)
        ## surface runoff rate (m/s)
        self.surfrunoffrate = np.zeros(array_shape)
        ## subsurface runoff rate (m/s)
        self.subrunoffrate = np.zeros(array_shape)
        ## soil pH
        self.soil_pH = np.zeros((Days,) + field_shape)
        self.soil_ccH = np.zeros((Days,) + field_shape)

        ##############################
        ## fields with 2 levels
        ##############################
        #### soil properties: temperature, moisture
        ## soil temperature 0-7cm (surface, topsoil), 7-28cm (midsoil, deepsoil)
        self.soil_temp = np.zeros(dlvl)
        ## soil moisture 0-7cm (surface, topsoil), 7-28cm (midsoil, deepsoil)
        self.soil_moist = np.zeros(dlvl)
        ## soil porosity (moiture at saturation) 0-7cm (surface, topsoil), 7-28cm (midsoil, deepsoil)
        self.soil_satmoist = np.zeros(dlvl)
        #### N uptake: NH4+, NO3-
        ## uptake of NH4+ by plants (topsoil, midsoil)
        self.ammNuptake = np.zeros(dlvl)
        ## uptake of NO3- by plants (topsoil, midsoil)
        self.nitNuptake = np.zeros(dlvl)
        ## water uptake by plants
        self.wateruptake = np.zeros(dlvl)
        #### upwards fluxes: diffusion
        ## urea diffusion (aq) (topsoil->surface; midsoil->topsoil)
        self.ureadiffusionup = np.zeros(tlvl)
        ## TAN diffusion (aq) (topsoil->surface; midsoil->topsoil)
        self.TANdiffusionup = np.zeros(tlvl)
        ## NH3 diffusion (gas) (topsoil->surface; midsoil->topsoil)
        self.NH3diffusionup = np.zeros(tlvl)
        ## NO3 diffusion (aq) (topsoil->surface; midsoil->topsoil)
        self.NO3diffusionup = np.zeros(tlvl)

        ##############################
        ## fields with 3 levels
        ##############################
        #### variables affected by management
        ## soil water content (after irrigation/watering)
        self.theta = np.zeros(tlvl)
        ## soil water drainage (m/s) (Note the difference between infiltration and drainage)
        self.drainagerate = np.zeros(tlvl)
        #### N pools and concentrations: urea, TAN (NH4 + NH3), NO3
        ## urea pools (surface, topsoil, midsoil)
        self.urea_pool = np.zeros(tlvl)
        ## urea concentrations (surface, topsoil, midsoil)
        self.urea_amount = np.zeros(tlvl)
        ## TAN pools (surface, topsoil, midsoil)
        self.TAN_pool = np.zeros(tlvl)
        ## TAN concentrations (surface, topsoil, midsoil)
        self.TAN_amount = np.zeros(tlvl)
        ## NH3 gaseous concentrations (surface, topsoil, midsoil)
        self.NH3_gas = np.zeros(tlvl)
        ## NO3- pools (surface, topsoil, midsoil)
        self.NO3_pool = np.zeros(tlvl)
        ## NO3- concentrations (surface, topsoil, midsoil)
        self.NO3_amount = np.zeros(tlvl)

        #### downwards fluxes: diffusion, infiltration (leaching)
        ## urea diffusion (aq) (surface->topsoil, topsoil->midsoil, midsoil->deep soil)
        self.ureadiffusiondown = np.zeros(tlvl)
        ## TAN diffusion (aq) (surface->topsoil, topsoil->midsoil, midsoil->deep soil)
        self.TANdiffusiondown = np.zeros(tlvl)
        ## NH3 diffusion (gas) (surface->topsoil, topsoil->midsoil, midsoil->deep soil)
        self.NH3diffusiondown = np.zeros(tlvl)
        ## NO3 diffusion (aq) (surface->topsoil, topsoil->midsoil, midsoil->deep soil)
        self.NO3diffusiondown = np.zeros(tlvl)
        ## urea infiltration/leaching
        self.ureainfil = np.zeros(tlvl)
        ## TAN infiltration/leaching
        self.TANinfil = np.zeros(tlvl)
        ## NO3 infiltration/leaching 
        self.NO3infil = np.zeros(tlvl)

        #### chemical transformations: urea hydrolysis, nitrification
        ## urea hydrolysis
        self.ureahydrolysis = np.zeros(tlvl)
        ## nitrification
        self.NH4nitrif = np.zeros(tlvl)

        #### variables: resistances
        ## resistances
        self.Rdiffaq = np.zeros(tlvl)
        self.Rdiffgas = np.zeros(tlvl)
        
        
        ## pH and H+ ions concentration
        self.pH = pH_value
        self.cc_H = np.float(10**(-pH_value))

        ## adsorption 
        self.Kd = np.zeros(array_shape[1:])

        ## irrigation idx
        self.irrigation_idx = np.zeros(array_shape[1:])

        ## output fields
        self.o_NH3flux = np.zeros(outarray_shape)
        self.o_washoff = np.zeros(outarray_shape)
        self.o_nitrif = np.zeros(outarray_shape)
        self.o_NH4leaching = np.zeros(outarray_shape)
        self.o_diffaq = np.zeros(outarray_shape)
        self.o_diffgas = np.zeros(outarray_shape)
        self.o_ammNuptake = np.zeros(outarray_shape)
        self.o_nitNuptake = np.zeros(outarray_shape)
        self.o_NO3washoff = np.zeros(outarray_shape)
        self.o_NO3leaching = np.zeros(outarray_shape)
        self.o_NO3diff = np.zeros(outarray_shape)

        if grazing is True:
            ##
            self.urine_N = np.zeros(array_shape)
            self.manure_N = np.zeros(array_shape)
            self.urine = np.zeros(array_shape)
            self.manure_water = np.zeros(array_shape)
            ## manure pools of FYM and dung
            self.manure_pool_FYM = np.zeros(array_shape)
            self.manure_pool_dung = np.zeros(array_shape)
            ## manure washoff
            self.manure_washoff_FYM = np.zeros(array_shape)
            self.manure_washoff_dung = np.zeros(array_shape)
            ## water pools of FYM and dung
            self.Total_water_pool_FYM = np.zeros(array_shape)
            self.Total_water_pool_dung = np.zeros(array_shape)
            ## org N pools of FYM and dung
            self.avail_N_pool_FYM = np.zeros(array_shape)
            self.resist_N_pool_FYM = np.zeros(array_shape)
            self.unavail_N_pool_FYM = np.zeros(array_shape)
            self.avail_N_pool_dung = np.zeros(array_shape)
            self.resist_N_pool_dung = np.zeros(array_shape)
            self.unavail_N_pool_dung = np.zeros(array_shape)
            ## org N washoff from FYM and dung
            self.avail_N_washoff_FYM = np.zeros(array_shape)
            self.resist_N_washoff_FYM = np.zeros(array_shape)
            self.unavail_N_washoff_FYM = np.zeros(array_shape)
            self.avail_N_washoff_dung = np.zeros(array_shape)
            self.resist_N_washoff_dung = np.zeros(array_shape)
            self.unavail_N_washoff_dung = np.zeros(array_shape)
            ## NH3 emissions from grazing FYM and dung
            self.NH3flux_FYM = np.zeros(array_shape)
            self.ureawashoff_FYM = np.zeros(array_shape)
            self.TANwashoff_FYM = np.zeros(array_shape)
            self.NH3flux_dung = np.zeros(array_shape)
            self.TANwashoff_dung = np.zeros(array_shape)
            self.o_NH3flux_FYM = np.zeros(outarray_shape)
            self.o_washoff_FYM = np.zeros(outarray_shape)
            # self.o_nitrif_FYM = np.zeros(outarray_shape)
            self.o_NH4leaching_FYM = np.zeros(outarray_shape)
            self.o_diffaq_FYM = np.zeros(outarray_shape)
            self.o_diffgas_FYM = np.zeros(outarray_shape)
            self.o_NH3flux_dung = np.zeros(outarray_shape)
            self.o_washoff_dung = np.zeros(outarray_shape)
            # self.o_nitrif_dung = np.zeros(outarray_shape)
            self.o_NH4leaching_dung = np.zeros(outarray_shape)
            self.o_diffaq_dung = np.zeros(outarray_shape)
            self.o_diffgas_dung = np.zeros(outarray_shape)
            ## manure orgN and TAN (incl. urea) enter soils
            self.o_soil_orgN_FYM = np.zeros(outarray_shape)
            self.o_soil_TAN_FYM = np.zeros(outarray_shape)
            self.o_soil_orgN_dung = np.zeros(outarray_shape)
            self.o_soil_TAN_dung = np.zeros(outarray_shape)
            self.o_soil_orgN = np.zeros(outarray_shape)
            self.o_soil_TAN = np.zeros(outarray_shape)
        
    def sim_env(self,dayidx):
        # print('LAND ENV: open env')
        #### environmental conditions
        if CONFIG_machine == "STREAM":
            ## on STREAM, met data are daily average data
            ## soil temperature; K to degC
            self.soil_temp[0,1:] = groundtemp_filelvl1.stl1[dayidx,self.plat1:self.plat2,:] - 273.15
            self.soil_temp[1,1:] = groundtemp_filelvl2.stl2[dayidx,self.plat1:self.plat2,:] - 273.15
            ## soil moisture; m3/m3
            self.soil_moist[0,1:] = soilmoist_filelvl1.swvl1[dayidx,self.plat1:self.plat2,:]
            self.soil_moist[1,1:] = soilmoist_filelvl2.swvl2[dayidx,self.plat1:self.plat2,:]
             ## evaporation from bare soil; m to g/m2/s
            self.evap[1:] = evap_file.evabs[dayidx,self.plat1:self.plat2,:]*(-1e6)/(24*3600)
            ## rainfall; kg/m2 to g/m2/s
            self.rainfall[1:] = rain_file.tcrw[dayidx,self.plat1:self.plat2,:]*1e3/(24*3600)
            ## surface runoff; m to m/s
            self.surfrunoffrate[1:] = runoff_file.sro[dayidx,self.plat1:self.plat2,:]/(24*3600)
            ## sub-surface runoff; m to m/s
            self.subrunoffrate[1:] = subrunoff_file.ssro[dayidx,self.plat1:self.plat2,:]/(24*3600)
            ## atmospheric resistances; s/m
            self.R_atm[1:] = ratm_file.RAM1[dayidx,self.plat1:self.plat2,:]+ratm_file.RB1[dayidx,self.plat1:self.plat2,:]
        
        else:
            ## on Archer and Jasmin, met data are hourly data (higher temperal resolution)
            hhidx = dayidx*24
            ## soil temperature; K to degC
            self.soil_temp[0,1:] = groundtemp_filelvl1.stl1[hhidx:hhidx+24,self.plat1:self.plat2,:] - 273.15
            self.soil_temp[1,1:] = groundtemp_filelvl2.stl2[hhidx:hhidx+24,self.plat1:self.plat2,:] - 273.15
            ## soil moisture; m3/m3
            self.soil_moist[0,1:] = soilmoist_filelvl1.swvl1[hhidx:hhidx+24,self.plat1:self.plat2,:]
            self.soil_moist[1,1:] = soilmoist_filelvl2.swvl2[hhidx:hhidx+24,self.plat1:self.plat2,:]
            self.soil_moist = np.maximum(self.soil_moist,0.0)

            ## evaporation from bare soil; m to g/m2/s
            self.evap[1:] = evap_file.evabs[hhidx:hhidx+24,self.plat1:self.plat2,:]*(-1e6)/(timestep*3600)
            ## rainfall; kg/m2 to g/m2/s
            self.rainfall[1:] = rain_file.tcrw[hhidx:hhidx+24,self.plat1:self.plat2,:]*1e3/(timestep*3600)
            ## surface runoff; m to m/s
            self.surfrunoffrate[1:] = runoff_file.sro[hhidx:hhidx+24,self.plat1:self.plat2,:]/(timestep*3600)
            ## sub-surface runoff; m to m/s
            self.subrunoffrate[1:] = subrunoff_file.ssro[hhidx:hhidx+24,self.plat1:self.plat2,:]/(timestep*3600)
            ## atmospheric resistances; s/m
            self.R_atm[1:] = ratm_file.RAM1[hhidx:hhidx+24,self.plat1:self.plat2,:]+ratm_file.RB1[hhidx:hhidx+24,self.plat1:self.plat2,:]
        
        soilbd_ds = open_ds(infile_path+soil_data_path+soilbdfile)
        ## buld density unit: kg/dm3
        soilbd = soilbd_ds.T_BULK_DEN.values[self.plat1:self.plat2,:]
        soilporosity = 1 - (soilbd/(rho_soil/1000))
        self.soil_satmoist[:] = soilporosity
        self.soil_satmoist = np.minimum(self.soil_satmoist,0.99)
           
    
    def met_input_interp(self,template):
        ##################################
        ## fill land input data
        ##################################
        self.soil_temp[0] = field_var_fill(sd_template=template,
                                            input_field=self.soil_temp[0])  ## degC
        self.soil_temp[1] = field_var_fill(sd_template=template,
                                            input_field=self.soil_temp[1]) ## degC
        # rhum_data = field_var_fill(sd_template=template,input_field=rhum_data)  ## per cent
        # wind_data = field_var_fill(sd_template=template,input_field=wind_data)  ## m/s
        self.evap = field_var_fill(sd_template=template,
                                            input_field=self.evap) ## g/day
        # soilmoist_data = field_var_fill(sd_template=template,input_field=soilmoist_data)  ## m3/m3
        self.soil_moist[0]= field_var_fill(sd_template=template,
                                            input_field=self.soil_moist[0])  ## m3/m3
        self.soil_moist[1] = field_var_fill(sd_template=template,
                                            input_field=self.soil_moist[1])  ## m3/m3
        self.soil_satmoist[0] = field_var_fill(sd_template=template,
                                            input_field=self.soil_satmoist[0]) ## m3/m3
        self.soil_satmoist[1] = field_var_fill(sd_template=template,
                                            input_field=self.soil_satmoist[1]) ## m3/m3
        self.R_atm = field_var_fill(sd_template=template,
                                            input_field=self.R_atm)  ## s/m
        self.surfrunoffrate = field_var_fill(sd_template=template,
                                            input_field=self.surfrunoffrate)  ## m/day
        self.subrunoffrate = field_var_fill(sd_template=template,
                                            input_field=self.subrunoffrate)  ## m/day
        self.soil_moist[self.soil_moist>self.soil_satmoist] = self.soil_satmoist[self.soil_moist>self.soil_satmoist]
        return

    

    def spreading_time(self,fert_type='mineral',
                        crop_type=None,N_input=None,plant_calendar=None,harvest_calendar=None,
                        fert_freq=None,soil_pH=None):
        ## lats and lons
        lats = CONFIG_lats
        lons = CONFIG_lons
        N_app = np.zeros(CONFIG_mtrx2)
        water_app = np.zeros(CONFIG_mtrx2)
        ## mark all application day
        N_app_mark = np.zeros(CONFIG_mtrx2)

        # print(np.where((np.isnan(harvest_calendar))&(~np.isnan(plant_calendar))))

        for dd in np.arange(CONFIG_mtrx2[0]):
            # N_app[dd][dd==plant_calendar] = N_input[dd==plant_calendar]
            N_app[dd][dd==plant_calendar] = N_input[dd==plant_calendar]/2
            water_app[dd][dd==plant_calendar+1] = irrig_water
            N_app[dd][dd==np.floor((plant_calendar+harvest_calendar)/2)] = N_input[dd==np.floor((plant_calendar+harvest_calendar)/2)]/2

        N_app_mark[N_app!=0] = 1
        # water_app[N_app!=0] = irrig_water

        # print(crop_type)
        # for lat in np.arange(lats):
        #     for lon in np.arange(lons):
        #         ## timing of planting and harvesting
        #         plt_time = plant_calendar[lat,lon]
        #         har_time = harvest_calendar[lat,lon]
        #         ## how many times that N is applied on fields in an annual cycle
        #         app_freq = fert_freq[lat,lon]
        #         ## harvesting goes into the next year
        #         if ~np.isnan(app_freq):
        #             if ~np.isnan(plt_time):
        #                 if ~np.isnan(har_time):
                            # if har_time<plt_time:
                            #     har_time = har_time+Days
                            # if app_freq <=1.0:
                            #     N_app[int(plt_time),lat,lon] = N_input[lat,lon]*app_freq
                            #     ## mark up the N application
                            #     N_app_mark[int(plt_time),lat,lon] = 1 
                            #     ## potential irrigation: 50 mm
                            #     water_app[int(plt_time),lat,lon] = irrig_water
                            # elif app_freq>1.0:
                            #     ## get the integer times of N application
                            #     tapp = np.floor(app_freq)
                            #     ## application intervals
                            #     app_int = int(abs(int(har_time)-int(plt_time))/(tapp+1))
                            #     ## index for application
                            #     app_idx = np.arange(int(plt_time),int(har_time)-1,app_int)
                            #     ## in tapp times, the application rate equals the readed applcation rates
                            #     for idx in app_idx[:-1]:
                            #         N_app[int(idx),lat,lon] = N_input[lat,lon]
                            #         N_app_mark[int(idx),lat,lon] = 1
                            #         water_app[int(idx),lat,lon] = irrig_water
                            #     ## the residual of N application rates
                            #     N_app[int(app_idx[-1]),lat,lon] = (app_freq-tapp)*N_input[lat,lon]
                            #     N_app_mark[int(app_idx[-1]),lat,lon] = 1
                            #     water_app[int(app_idx[-1]),lat,lon] = irrig_water
        
        ## determine soil pH after urea application
        if soil_pH is not None:
            # print('check pH over.')
            self.soil_pH = soil_pH_postapp(base_pH=soil_pH,app_timing_map=N_app_mark,fert_pH=8.5)[:,self.plat1:self.plat2,:] 
            self.soil_ccH = 10**(-self.soil_pH)
        return N_app,water_app

    def chem_fert_type(self):
        ## country level chemical fertilizer use: 1) ammonium N, 2) nitrate N, 3) urea N
        fertds = open_ds(infile_path+crop_data_path+fertfilename)
        nitN = fertds.nitrate_fert.sel(year=sim_year).values
        ammN = fertds.ammonium_fert.sel(year=sim_year).values
        ureaN = fertds.urea_fert.sel(year=sim_year).values
        totalN = nitN + ammN + ureaN
        fnitN = nitN/totalN
        fnitNmean = np.nansum(nitN)/np.nansum(totalN)
        fammN = ammN/totalN
        fammNmean = np.nansum(ammN)/np.nansum(totalN)
        fureaN = ureaN/totalN
        fureaNmean = np.nansum(ureaN)/np.nansum(totalN)
        return fnitN, fammN, fureaN, fnitNmean, fammNmean, fureaNmean
    
    def crop_calendar(self,filepath):
        ## crop calendar dataset
        cropcalds = open_ds(filepath)
        ## read planting and harvesting dates
        # plantdate = cropcalds['plant.start'].values
        # harvestdate = cropcalds['harvest.start'].values
        # plantfill = cropcalds['plant'].values
        # harvestfill = cropcalds['harvest'].values
        # plantdate[(np.isnan(plantdate))&(~np.isnan(plantfill))] = plantfill[(np.isnan(plantdate))&(~np.isnan(plantfill))]
        # harvestdate[(np.isnan(harvestdate))&(~np.isnan(harvestfill))] = harvestfill[(np.isnan(harvestdate))&(~np.isnan(harvestfill))]
        # ## harvesting date goes into next year
        # harvestdate[harvestdate<plantdate] = harvestdate[harvestdate<plantdate]+Days

        ## new GGCMI crop calendar with ir,rf classification
        plantdate = cropcalds.planting_day.values
        harvestdate = cropcalds.harvesting_day.values
        harvestdate[harvestdate<plantdate] = harvestdate[harvestdate<plantdate]+Days
        return plantdate,harvestdate

    def chem_fert_input(self,crop):
        ## read N application rates dataset for crops
        fertds = open_ds(infile_path+crop_data_path+crop+cropfileformat)
        ## crop calendar dataset
        # cropcalspath = infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
        ## new GGCMI crop calendar
        cropcalspath = infile_path+crop_data_path+crop_GGCMI_p3_calendar+crop+crop_GGCMI_p3_calendarformat
        # cropcalds = open_ds(infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat)
        ## base soil pH dataset
        soilpHds = open_ds(infile_path+soil_data_path+soilpHfile)

        totalN = fertds.TotalN.sel(year=sim_year).values*1e3
        print("Total N: ",np.nansum(totalN[self.plat1:self.plat2,:]/1e9))
        ## N application rate is interpolated;
        rateN = fertds.Nrate.sel(year=sim_year).values*1e3/1e4
        croparea = fertds.croparea.values*1e4
        croparea[totalN!=0] = totalN[totalN!=0]/rateN[totalN!=0]
        croparea = np.nan_to_num(croparea)
        ## at least, there will be two applications
        ## 1st: at the start of planting, 2nd: between planting and harvesting
        # rateN = rateN/2
        app_freq = totalN/(rateN*croparea)
        # app_freq = np.nan_to_num(app_freq,posinf=0,neginf=0)
        # print("max app freq: ",np.nanmax(app_freq))
        ## the maximum application frequency in a year is 6
        # app_freq = np.minimum(app_freq,2)

        ## read planting and harvesting dates
        plantidx, harvestidx = self.crop_calendar(filepath = cropcalspath)

        soilph = np.zeros(CONFIG_mtrx2)
        soilph[:] = soilpHds.T_PH_H2O.values
        self.pH = soilph[:,self.plat1:self.plat2,:] 
        self.cc_H = 10**(-self.pH)

        chem_N_tocrop,water_to_crop = self.spreading_time(fert_type='mineral',
                        crop_type=crop,N_input=rateN,plant_calendar=plantidx,harvest_calendar=harvestidx,
                        fert_freq=app_freq,soil_pH = soilph)
        fnitN, fammN, fureaN, nitmean, ammmean, ureamean = self.chem_fert_type()
        ## when there is no country-level information for fertilizer type, use global mean value
        fnitN[(np.isnan(fnitN))&(totalN!=0)] = nitmean
        fammN[(np.isnan(fammN))&(totalN!=0)] = ammmean
        fureaN[(np.isnan(fureaN))&(totalN!=0)] = ureamean
        ## N rate of three types of fertilizer;
        ## the fraction of usage is represented by the 
        ## fractional cropping area that uses the corresponding N fertilizer
        self.NO3_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.nitN_area = fnitN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 
        self.TAN_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.ammN_area = fammN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 
        self.urea_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.ureaN_area = fureaN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 
        self.water_added = water_to_crop[:,self.plat1:self.plat2,:]
        print("NO3 N:",np.nansum(np.nansum(self.NO3_added,axis=0)*self.nitN_area)/1e9)
        print("TAN N:",np.nansum(np.nansum(self.TAN_added,axis=0)*self.ammN_area)/1e9)
        print("urea N:",np.nansum(np.nansum(self.urea_added,axis=0)*self.ureaN_area)/1e9)
        ## met data interpolation
        self.met_input_interp(totalN[self.plat1:self.plat2,:])
        return

    def manure_fert_input(self,livestock_name,production_system,mms_cat,phase,crop=None):
        ## read N application rates dataset for crops
        manureds = open_ds(output_path+livestock_name+'.'+production_system+'.'+mms_cat+'.'+phase+\
                                '.'+str(sim_year)+'.manureapp.nc')
        ## crop calendar dataset
        calds = open_ds(infile_path+crop_data_path+manure_appcalendar)
        # cropcalds = open_ds(infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat)

        soilpHds = open_ds(infile_path+soil_data_path+soilpHfile)
        soilph = np.zeros(CONFIG_mtrx2)
        soilph[:] = soilpHds.T_PH_H2O.values
        self.pH = soilph[:,self.plat1:self.plat2,:] 
        self.cc_H = 10**(-self.pH)

        tan = np.zeros(CONFIG_mtrx[1:])
        manure = np.zeros(CONFIG_mtrx[1:])
        water = np.zeros(CONFIG_mtrx[1:])
        availN = np.zeros(CONFIG_mtrx[1:])
        resistN = np.zeros(CONFIG_mtrx[1:])
        unavailN = np.zeros(CONFIG_mtrx[1:])

        spring_TAN = manureds.spring_TAN.values
        spring_manure = manureds.spring_manure.values
        spring_water = manureds.spring_water.values
        spring_avail_N = manureds.spring_availN.values
        spring_resist_N = manureds.spring_resistN.values
        spring_unavail_N = manureds.spring_unavailN.values
        
        winter_TAN = manureds.winter_TAN.values
        winter_manure = manureds.winter_manure.values
        winter_water = manureds.winter_water.values
        winter_avail_N = manureds.winter_availN.values
        winter_resist_N = manureds.winter_resistN.values
        winter_unavail_N = manureds.winter_unavailN.values
        
        spring_plantdate = calds.spring_plant_mean.values
        winter_plantdate = calds.winter_plant_mean.values

        ## test code: area needs to be revised
        ## current test used fixed application rate of manure (5mm slurry) to determine the crop area
        ## the aim is testing/debugging
        croparea = (spring_manure+spring_water)/5e3
        # print(spring_plantdate.shape,tan.shape)
        for dd in np.arange(Days):
            tan[spring_plantdate==dd] = spring_TAN[spring_plantdate==dd]/croparea[spring_plantdate==dd]
            manure[spring_plantdate==dd] = spring_manure[spring_plantdate==dd]/croparea[spring_plantdate==dd]
            water [spring_plantdate==dd] = spring_water[spring_plantdate==dd]/croparea[spring_plantdate==dd]
            availN[spring_plantdate==dd] = spring_avail_N[spring_plantdate==dd]/croparea[spring_plantdate==dd]
            resistN[spring_plantdate==dd] = spring_resist_N[spring_plantdate==dd]/croparea[spring_plantdate==dd]
            unavailN[spring_plantdate==dd] = spring_unavail_N[spring_plantdate==dd]/croparea[spring_plantdate==dd]

            self.TAN_added[dd] = tan[self.plat1:self.plat2,:]
            self.manure_added[dd] = manure[self.plat1:self.plat2,:]
            self.water_added[dd] = water[self.plat1:self.plat2,:]
            self.avail_N_added[dd] = availN[self.plat1:self.plat2,:]
            self.resist_N_added[dd] = resistN[self.plat1:self.plat2,:]
            self.unavail_N_added[dd] = unavailN[self.plat1:self.plat2,:]

            tan[:] = 0.0
            manure[:]  = 0.0
            water[:]  = 0.0
            availN[:] = 0.0
            resistN[:]  = 0.0
            unavailN[:]  = 0.0

            tan[winter_plantdate==dd] = winter_TAN[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            manure[winter_plantdate==dd] = winter_manure[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            water [winter_plantdate==dd] = winter_water[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            availN[winter_plantdate==dd] = winter_avail_N[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            resistN[winter_plantdate==dd] = winter_resist_N[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            unavailN[winter_plantdate==dd] = winter_unavail_N[winter_plantdate==dd]/croparea[winter_plantdate==dd]
            
            self.TAN_added[dd] = tan[self.plat1:self.plat2,:]
            self.manure_added[dd] = manure[self.plat1:self.plat2,:]
            self.water_added[dd] = water[self.plat1:self.plat2,:]
            self.avail_N_added[dd] = availN[self.plat1:self.plat2,:]
            self.resist_N_added[dd] = resistN[self.plat1:self.plat2,:]
            self.unavail_N_added[dd] = unavailN[self.plat1:self.plat2,:]
        return 

    ## irrigation of the land
    def irrigation_event(self,dayidx,plantday,harvestday,irrigation_map,theta_WP,theta_FC):
        ## justify if the soil water content is lower than the wilting point
        soil_water_deficit = (self.theta[1,0] + self.theta[2,0])/2 - theta_WP
        ## only for the periods of growing season
        soil_water_deficit[dayidx<=plantday] = 1.0
        soil_water_deficit[dayidx>=harvestday] = 1.0
        ## only for irrigated croplands (idx=2), excluding rainfed croplands
        ## Note that irrigated cropland index is 2 (rainfed cropland index is 1)
        soil_water_deficit[irrigation_map!=2] = 1.0
        ## update irrigation idx field: soil water is insuficient area plus 1.0
        self.irrigation_idx[soil_water_deficit<0] = self.irrigation_idx[soil_water_deficit<0] + 1.0

        ## water is added when irrigation idx is 1.0 (soil water content is lower than wilting point)
        ## water is added til the soil moisture of the top three layer (0-14cm) is equivalent to the field capacity
        ## m3/m3 change to g/m2
        self.water_added[dayidx][self.irrigation_idx==1.0] = theta_FC[self.irrigation_idx==1.0]*(zlayers[0]+zlayers[1]+zlayers[2])*1e6
        # self.water_added[dayidx][self.irrigation_idx!=1] = 0.0
        self.water_added[dayidx][irrigation_map!=2] = 0.0
        # self.water_added[dayidx][(irrigation_map!=2)&(self.water_added[dayidx]!=0)] = irrig_water/5.0
        self.water[7] = self.water_added[dayidx]
        ## reset irrigation idx
        self.irrigation_idx[:] = 0.0
       
        # self.water_added[dayidx][np.nanmean(self.soil_moist[0,1:],axis=0)<0.1] = self.water_added[dayidx][np.nanmean(self.soil_moist[0,1:],axis=0)<0.1]
        # self.water_added[dayidx][np.nanmean(self.soil_moist[0,1:],axis=0)>=0.1] = self.water_added[dayidx][np.nanmean(self.soil_moist[0,1:],axis=0)>=0.1]/5.0
        # self.water[7] = self.water_added[dayidx]
        # self.water[7][self.water_added[dayidx]!=0] = self.water_added[dayidx][self.water_added[dayidx]!=0]
        ## get rid of the places that do not need irrigation (rain-fed field)
        # self.water[7][np.nanmean(self.soil_moist[:,1:],axis=(0,1))>0.1] = 0.0

 
    ## soil water content of each layer after irrigation
    def post_soil_moist(self,theta1,theta2,theta3,theta_sat,water_added):
        z_added = water_added/1e6
        delta_z1 = (theta_sat-theta1)*zlayers[0]
        delta_z2 = (theta_sat-theta2)*zlayers[1]
        # delta_z3 = (theta_sat-theta3)*zlayers[2]

        post_theta1 = (theta1*zlayers[0]+z_added)/zlayers[0]
        # post_theta1[post_theta1>theta_sat] = theta_sat[post_theta1>theta_sat]
        post_theta1 = np.minimum(post_theta1,theta_sat)
        water_left = z_added - delta_z1
        # water_left[water_left<0] = 0.0
        water_left = np.maximum(water_left,0.0)
        post_theta2 = (theta2*zlayers[1]+water_left)/zlayers[1]
        # post_theta2[post_theta2>theta_sat] = theta_sat[post_theta2>theta_sat]
        post_theta2 = np.minimum(post_theta2,theta_sat)
        water_left = water_left - delta_z2
        # water_left[water_left<0] = 0.0
        water_left = np.maximum(water_left,0.0)
        post_theta3 = (theta3*zlayers[2]+water_left)/zlayers[2]
        # post_theta3[post_theta3>theta_sat] = theta_sat[post_theta3>theta_sat]
        post_theta3 = np.minimum(post_theta3,theta_sat)

        return post_theta1, post_theta2, post_theta3
    
    ## initial infiltration following irrigation events
    def irr_infil(self,theta1,theta2,theta3,theta_sat,water_added,delta_t=24.0):
        z_added = water_added/1e6
        infilrate1 = (z_added-zlayers[0]*(theta_sat-theta1))/(3600*delta_t)
        infilrate1 = np.maximum(infilrate1,0.0)
        infilrate2 = (z_added-zlayers[0]*(theta_sat-theta1)-zlayers[1]*(theta_sat-theta2))/(3600*delta_t)
        infilrate2 = np.maximum(infilrate2,0.0)
        infilrate3 = (z_added-zlayers[0]*(theta_sat-theta1)-zlayers[1]*(theta_sat-theta2)-\
                        zlayers[2]*(theta_sat-theta3))/(3600*delta_t)
        infilrate3 = np.maximum(infilrate3,0.0)
        # infilrate1 = infilrate3
        # infilrate2 = infilrate3
        return infilrate1, infilrate2, infilrate3

    ## crop life cycle stage
    def crop_lifecycle_stage(self,dayidx,plant_calendar,harvest_calendar,stages=6):
        if dayidx >= Days:
            dayidx = int(dayidx - np.floor(dayidx/Days)*Days)   
        growing_length = harvest_calendar - plant_calendar
        stage_idx = np.floor((dayidx - plant_calendar)/(growing_length/stages))+1
        stage_idx[stage_idx<0] = 0
        return stage_idx

    ## layer 0: surface source layer 2cm
    ## layer 1: topsoil layer 5cm
    ## layer 2: intermediate soil layer 7cm
    ## layer 3; deep soil layer 14cm
    ## determine the order for updating pools and fluxes
    def source_layer(self,tech):
        if tech == "surf":
            sourcelayer = 0
        elif tech == "disk":
            sourcelayer = 1
        elif tech == "injection":
            sourcelayer = 2
        return sourcelayer
    
    ## 
    def daily_init(self,manure=False,grazing=False):
        ## N pools
        self.TAN_pool[:,0] = self.TAN_pool[:,-1]
        self.urea_pool[:,0] = self.urea_pool[:,-1]
        self.NO3_pool[:,0] = self.NO3_pool[:,-1]
        ## concentrations
        self.TAN_amount[:,0] = self.TAN_amount[:,-1]
        self.NH3_gas[:,0] = self.NH3_gas[:,-1]
        self.urea_amount[:,0] = self.urea_amount[:,-1]
        self.NO3_amount[:,0] = self.NO3_amount[:,-1]
        ## upward diffusions
        self.TANdiffusionup[:,0] = self.TANdiffusionup[:,-1]
        self.NH3diffusionup[:,0] = self.NH3diffusionup[:,-1]
        self.ureadiffusionup[:,0] = self.ureadiffusionup[:,-1]
        self.NO3diffusionup[:,0] = self.NO3diffusionup[:,-1]
        ## soil water content
        self.theta[:,0] = self.theta[:,-1]
        ## 
        self.soil_moist[:,0] = self.soil_moist[:,-1]
        if manure is True:
            self.avail_N_pool[0] = self.avail_N_pool[-1]
            self.resist_N_pool[0] = self.resist_N_pool[-1]
            self.unavail_N_pool[0] = self.unavail_N_pool[-1]
        if grazing is True:
            self.avail_N_pool[0] = self.avail_N_pool[-1]
            self.resist_N_pool[0] = self.resist_N_pool[-1]
            self.unavail_N_pool[0] = self.unavail_N_pool[-1]
            self.manure_pool_FYM[0] = self.manure_pool_FYM[-1]
            self.avail_N_pool_FYM[0] = self.avail_N_pool_FYM[-1]
            self.resist_N_pool_FYM[0] = self.resist_N_pool_FYM[-1]
            self.unavail_N_pool_FYM[0] = self.unavail_N_pool_FYM[-1]
            self.manure_pool_dung[0] = self.manure_pool_dung[-1]
            self.avail_N_pool_dung[0] = self.avail_N_pool_dung[-1]
            self.resist_N_pool_dung[0] = self.resist_N_pool_dung[-1]
            self.unavail_N_pool_dung[0] = self.unavail_N_pool_dung[-1]
        return

    ##
    def daily_output(self,dayidx,grazing=False):
        if dayidx >= Days:
            dayidx = int(dayidx - np.floor(dayidx/Days)*Days)

        if grazing is True:
            ## FYM; lvl_idx=0
            self.o_NH3flux_FYM[dayidx] = np.nansum(self.NH3flux_FYM[1:],axis=0) + self.o_NH3flux_FYM[dayidx]
            self.o_washoff_FYM[dayidx] = np.nansum(self.TANwashoff_FYM[1:]+self.ureawashoff_FYM[1:]+\
                                            self.avail_N_washoff_FYM[1:]+self.resist_N_washoff_FYM[1:]+self.unavail_N_washoff_FYM[1:],axis=0) + \
                                            self.o_washoff_FYM[dayidx]
            # self.o_nitrif_FYM[dayidx] = np.nansum(self.NH4nitrif[1,1:],axis=0) + self.o_nitrif_FYM[dayidx]
            self.o_NH4leaching_FYM[dayidx] = np.nansum(self.TANinfil[0,1:]+self.ureainfil[0,1:],axis=0) + self.o_NH4leaching_FYM[dayidx]
            self.o_diffaq_FYM[dayidx] = np.nansum(self.TANdiffusiondown[0,1:]+self.ureadiffusiondown[0,1:],axis=0) + self.o_diffaq_FYM[dayidx]
            self.o_diffgas_FYM[dayidx] = np.nansum(self.NH3diffusiondown[0,1:],axis=0) + self.o_diffgas_FYM[dayidx]
            ## dung; lvl_idx=1
            self.o_NH3flux_dung[dayidx] = np.nansum(self.NH3flux_dung[1:],axis=0) + self.o_NH3flux_dung[dayidx]
            self.o_washoff_dung[dayidx] = np.nansum(self.TANwashoff_dung[1:]+\
                                            self.avail_N_washoff_dung[1:]+self.resist_N_washoff_dung[1:]+self.unavail_N_washoff_dung[1:],axis=0) + \
                                            self.o_washoff_dung[dayidx]
            self.o_NH4leaching_dung[dayidx] = np.nansum(self.TANinfil[1,1:],axis=0) + self.o_NH4leaching_dung[dayidx]
            self.o_diffaq_dung[dayidx] = np.nansum(self.TANdiffusiondown[1,1:],axis=0) + self.o_diffaq_dung[dayidx]
            self.o_diffgas_dung[dayidx] = np.nansum(self.NH3diffusiondown[1,1:],axis=0) + self.o_diffgas_dung[dayidx]
            ## urine patch; lvl_idx=2
            self.o_NH3flux[dayidx] = np.nansum(self.NH3flux[1:],axis=0) + self.o_NH3flux[dayidx]
            self.o_washoff[dayidx] = np.nansum(self.TANwashoff[1:]+self.ureawashoff[1:]+\
                                            self.avail_N_washoff[1:]+self.resist_N_washoff[1:]+self.unavail_N_washoff[1:],axis=0) + \
                                            self.o_washoff[dayidx]
            self.o_nitrif[dayidx] = np.nansum(self.NH4nitrif[2,1:],axis=0) + self.o_nitrif[dayidx]
            self.o_NH4leaching[dayidx] = np.nansum(self.TANinfil[2,1:]+self.ureainfil[2,1:],axis=0) + self.o_NH4leaching[dayidx]
            self.o_diffaq[dayidx] = np.nansum(self.TANdiffusiondown[2,1:]+self.ureadiffusiondown[2,1:],axis=0) + self.o_diffaq[dayidx]
            self.o_diffgas[dayidx] = np.nansum(self.NH3diffusiondown[2,1:],axis=0) + self.o_diffgas[dayidx]

        else:
            self.o_NH3flux[dayidx] = np.nansum(self.NH3flux[1:],axis=0)
            self.o_washoff[dayidx] = np.nansum(self.TANwashoff[1:]+self.ureawashoff[1:],axis=0)
            self.o_nitrif[dayidx] = np.nansum(self.NH4nitrif[:,1:],axis=(0,1))
            self.o_NH4leaching[dayidx] = np.nansum(self.TANinfil[-1,1:]+self.ureainfil[-1,1:],axis=0)
            self.o_diffaq[dayidx] = np.nansum(self.TANdiffusiondown[-1,1:]+self.ureadiffusiondown[-1,1:],axis=0)
            self.o_diffgas[dayidx] = np.nansum(self.NH3diffusiondown[-1,1:],axis=0)
            self.o_ammNuptake[dayidx] = np.nansum(self.ammNuptake[:,1:],axis=(0,1))
            self.o_nitNuptake[dayidx] = np.nansum(self.nitNuptake[:,1:],axis=(0,1))
            self.o_NO3washoff[dayidx] = np.nansum(self.NO3washoff[1:],axis=0)
            self.o_NO3leaching[dayidx] = np.nansum(self.NO3infil[-1,1:],axis=0)
            self.o_NO3diff[dayidx] = np.nansum(self.NO3diffusiondown[-1,1:],axis=0)
        # else:
        #     dayidx = dayidx - Days
        #     self.o_NH3flux[dayidx] = self.o_NH3flux[dayidx] + np.nansum(self.NH3flux[1:],axis=0)
        #     self.o_washoff[dayidx] = self.o_washoff[dayidx] + np.nansum(self.TANwashoff[1:]+self.ureawashoff[1:],axis=0)
        #     self.o_nitrif[dayidx] = self.o_nitrif[dayidx] + np.nansum(self.NH4nitrif[:,1:],axis=(0,1))
        #     self.o_NH4leaching[dayidx] = self.o_NH4leaching[dayidx] + np.nansum(self.TANinfil[-1,1:]+self.ureainfil[-1,1:],axis=(0,1))
        #     self.o_diffaq[dayidx] = self.o_diffaq[dayidx] + np.nansum(self.TANdiffusiondown[-1,1:]+self.ureadiffusiondown[-1,1:],axis=(0,1))
        #     self.o_diffgas[dayidx] = self.o_diffgas[dayidx] + np.nansum(self.NH3diffusiondown[-1,1:],axis=(0,1))
        #     self.o_ammNuptake[dayidx] = self.o_ammNuptake[dayidx] + np.nansum(self.ammNuptake[:,1:],axis=(0,1))
        #     self.o_nitNuptake[dayidx] = self.o_nitNuptake[dayidx] + np.nansum(self.nitNuptake[:,1:],axis=(0,1))
        #     self.o_NO3washoff[dayidx] = self.o_NO3washoff[dayidx] + np.nansum(self.NO3washoff[1:],axis=0)
        #     self.o_NO3leaching[dayidx] = self.o_NO3leaching[dayidx] + np.nansum(self.NO3infil[-1,1:],axis=(0,1))
        #     self.o_NO3diff[dayidx] = self.o_NO3diff[dayidx] + np.nansum(self.NO3diffusiondown[-1,1:],axis=(0,1))
        return

    ## sensitivity tests
    def sensitivity_test(self,var,test):
        
        ####################
        ## SOIL TEMPERATURE
        ####################
        if var == 'temp':
            if test == '-':
                ## soil temperature decrease by 2 degC
                self.soil_temp = self.soil_temp - 2.0
                print('Sensitivity tests for soil temp -2 degC')
            elif test == '+':
                ## soil temperature increase by 2 degC
                self.soil_temp = self.soil_temp + 2.0
                print('Sensitivity tests for soil temp +2 degC')
            elif test == '+3':
                ## soil temperature decrease by 3 degC
                self.soil_temp = self.soil_temp + 3.0
                print('Sensitivity tests for soil temp +3 degC')
            elif test == '+5':
                ## soil temperature decrease by 5 degC
                self.soil_temp = self.soil_temp + 5.0
                print('Sensitivity tests for soil temp +5 degC')

        ####################
        ## SOIL MOISTURE
        ####################
        if var == 'sw':
            if test == '-':
                ## soil moisture decrease by 10 per cent
                self.soil_moist = self.soil_moist * 0.9
                print('Sensitivity tests for soil moisture -10 %')
            else:
                ## soil moisture increase by 10 per cent
                self.soil_moist = self.soil_moist * 1.1
                self.soil_moist[self.soil_moist>self.soil_satmoist] = self.soil_satmoist[self.soil_moist>self.soil_satmoist]
                print('Sensitivity tests for soil moisture +10 %')

        #####################
        ## SURFACE RUNOFF
        #####################
        if var == 'srf':
            if test == '-':
                ## surface runoff decrease by 10 per cent
                self.surfrunoffrate = self.surfrunoffrate * 0.9
                print('Sensitivity tests for surface runoff rate -10 %')
            else:
                ## surface runoff increase by 10 per cent
                self.surfrunoffrate = self.surfrunoffrate * 1.1
                print('Sensitivity tests for surface runoff rate +10 %')

        #####################
        ## SUBSURFACE RUNOFF
        #####################
        if var == 'sub':
            if test == '-':
                ## sub-surface runoff decrease by 10 per cent
                self.subrunoffrate = self.subrunoffrate * 0.9
                print('Sensitivity tests for sub-surface runoff rate -10 %')
            else:
                ## sub-surface runoff increase by 10 per cent
                self.subrunoffrate = self.subrunoffrate * 1.1
                print('Sensitivity tests for sub-surface runoff rate +10 %')

        #####################
        ## N APPLICATION RATE
        #####################
        if var == 'Napp':
            if test == '-':
                ## N application rate decrease by 10 per cent
                self.TAN = self.TAN * 0.9
                self.urea = self.urea * 0.9
                print('Sensitivity tests for N application rate -10 %')
            else:
                ## N application rate increase by 10 per cent
                self.TAN = self.TAN * 1.1
                self.urea = self.urea * 1.1
                print('Sensitivity tests for N application rate +10 %')

        ###########
        ## SOIL PH
        ###########
        if var == 'ph':
            if test == '-':
                ## soil pH decrease by 0.5
                self.pH = self.pH - 0.5
                self.soil_pH = self.soil_pH - 0.5
                print('Sensitivity tests for soil pH -0.5')
            else:
                ## soil pH increase by 0.5
                self.pH = self.pH + 0.5
                self.soil_pH = self.soil_pH + 0.5
                print('Sensitivity tests for soil pH +0.5')

        return
    
    ## main sim function of chemical fertilizer application
    ## st: sensitivity test; stvar: sensitivity test value
    def land_sim(self,start_day_idx,end_day_idx,chem_fert_type,tech,crop=None,sim='base',stvar=False,st=False):

        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is: '+str(tech))
        ## soil clay, sand, silt fraction
        soilclayds = open_ds(infile_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        soilsandds = open_ds(infile_path+soil_data_path+soilsandfile)
        soilsand = soilsandds.T_SAND.values
        soilsiltds = open_ds(infile_path+soil_data_path+soilsiltfile)
        soilsilt = soilsiltds.T_SILT.values
        ## soil organic carbon
        soilocds = open_ds(infile_path+soil_data_path+soilorgCfile)
        soiloc = soilocds.T_OC.values
        ## soil bulk density
        soilbd_ds = open_ds(infile_path+soil_data_path+soilbdfile)
        soilbd = soilbd_ds.T_BULK_DEN.values
        ## adsorption constant
        Kd = ammonium_adsorption(clay_content=soilclay)[self.plat1:self.plat2,:]
        ## soil hydraulic conductivity
        Ks_sat = soil_hydraulic_conductivity(fsilt=soilsilt,fclay=soilclay,fsom=soiloc,BD=soilbd)[self.plat1:self.plat2,:]
        ## field capacity
        theta_fc = soil_field_capacity(BD=soilbd)[self.plat1:self.plat2,:]
        ## wilting point
        theta_wp = soil_wilting_point(fsand=soilsand,fclay=soilclay)[self.plat1:self.plat2,:]
        ## irrigated croplands index
        irrigatedcroplandds = open_ds(infile_path+crop_data_path+croplandtypefile)
        ## Note that irrigated cropland index is 2 (rainfed cropland index is 1)
        irrigated_croplandidx = irrigatedcroplandds.croplandsidx.values[self.plat1:self.plat2,:]
        ## crop calendar that detenmines the N uptake by crops
        if crop is not None:
            # cropcalspath = infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
            ## new GGCMI crop calendar
            cropcalspath = infile_path+crop_data_path+crop_GGCMI_p3_calendar+crop+crop_GGCMI_p3_calendarformat
            plantidx, harvestidx = self.crop_calendar(filepath = cropcalspath)
            plantidx = plantidx[self.plat1:self.plat2,:]
            harvestidx = harvestidx[self.plat1:self.plat2,:]
        sourcelayer = self.source_layer(tech)

        for dd in np.arange(start_day_idx,end_day_idx):
            if dd < Days:
                ## read in environmental variables
                self.sim_env(dd)
                ## irrigation events
                self.irrigation_event(dayidx=dd,plantday=plantidx,harvestday=harvestidx,
                                irrigation_map=irrigated_croplandidx,theta_WP=theta_wp,theta_FC=theta_fc)

                infilrate1, infilrate2, infilrate3 = self.irr_infil(theta1=self.theta[0,0],
                                    theta2=self.theta[0,0],theta3=self.theta[1,0],
                                    theta_sat=self.soil_satmoist[0,0],water_added=self.water_added[dd])

                if chem_fert_type == 'nitrate':
                    self.NO3[5] = self.NO3_added[dd]
                    sim_pH = self.pH[dd]
                    sim_ccH = self.cc_H[dd]
                elif chem_fert_type == 'ammonium':
                    self.TAN[5] = self.TAN_added[dd]
                    sim_pH = self.pH[dd]
                    sim_ccH = self.cc_H[dd]
                elif chem_fert_type == 'urea':
                    self.urea[5] = self.urea_added[dd]
                    sim_pH = self.soil_pH[dd]
                    sim_ccH = self.soil_ccH[dd]
            else:
                self.sim_env(dd-Days)
                self.irrigation_event(dayidx=dd-Days,plantday=plantidx,harvestday=harvestidx,
                                irrigation_map=irrigated_croplandidx,theta_WP=theta_wp,theta_FC=theta_fc)

                infilrate1, infilrate2, infilrate3 = self.irr_infil(theta1=self.theta[0,0],
                                    theta2=self.theta[0,0],theta3=self.theta[1,0],
                                    theta_sat=self.soil_satmoist[0,0],
                                    water_added=np.maximum(self.water_added[dd-Days],self.water_added[dd]))

                ## Note that the np.maximum function is to prevent application on the second year being overwritten
                ## by day index dd-Days
                if chem_fert_type == 'nitrate':
                    # print('chemical fertilizer applied: nitrate')
                    self.NO3[5] = np.maximum(self.NO3_added[dd-Days],self.NO3_added[dd])
                    sim_pH = np.maximum(self.pH[dd-Days],self.pH[dd])
                    sim_ccH = np.maximum(self.cc_H[dd-Days],self.cc_H[dd])
                elif chem_fert_type == 'ammonium':
                    self.TAN[5] = np.maximum(self.TAN_added[dd-Days],self.TAN_added[dd])
                    # print('chemical fertilizer applied: ammonium')
                    sim_pH = np.maximum(self.pH[dd-Days],self.pH[dd])
                    sim_ccH = np.maximum(self.cc_H[dd-Days],self.cc_H[dd])
                elif chem_fert_type == 'urea':
                    self.urea[5] = np.maximum(self.urea_added[dd-Days],self.urea_added[dd])
                    # print('chemical fertilizer applied: urea')
                    sim_pH = np.maximum(self.soil_pH[dd-Days],self.soil_pH[dd])
                    sim_ccH = np.maximum(self.soil_ccH[dd-Days],self.soil_ccH[dd])

            if sim != 'base':
                self.sensitivity_test(var=stvar,test=st)

            self.theta[0,1:] = np.copy(self.soil_moist[0,1:])
            self.theta[1,1:] = np.copy(self.soil_moist[0,1:])
            self.theta[2,1:] = np.copy(self.soil_moist[1,1:])

            # print('dd',dd,self.water_added[dd,0,589])

            for hh in np.arange(0,24):
                if CONFIG_machine == "STREAM":
                    surfrunoffrate = self.surfrunoffrate[hh+1]
                    subrunoffrate = self.subrunoffrate[hh+1] 
                    evap = self.evap[hh+1] 
                else:
                    surfrunoffrate = self.surfrunoffrate[hh+1]-self.surfrunoffrate[hh]
                    surfrunoffrate = np.maximum(surfrunoffrate,0.0)
                    subrunoffrate = self.subrunoffrate[hh+1]-self.subrunoffrate[hh]
                    subrunoffrate = np.maximum(subrunoffrate,0.0)
                    evap = self.evap[hh+1]-self.evap[hh]

                delta_sw1 = self.soil_moist[0,hh+1] - self.soil_moist[0,hh]
                delta_sw2 = self.soil_moist[1,hh+1] - self.soil_moist[1,hh]

                self.theta[0,hh+1] = self.theta[0,hh] + delta_sw1
                self.theta[1,hh+1] = self.theta[1,hh] + delta_sw1
                self.theta[2,hh+1] = self.theta[2,hh] + delta_sw2

                self.theta[0,hh+1] = np.maximum(self.theta[0,hh+1],0.01)
                self.theta[1,hh+1] = np.maximum(self.theta[1,hh+1],0.01)
                self.theta[2,hh+1] = np.maximum(self.theta[2,hh+1],0.01)

                # infilrate1, infilrate2, infilrate3 = self.irr_infil(theta1=self.theta[0,hh],
                #                     theta2=self.theta[0,hh],theta3=self.theta[1,hh],
                #                     theta_sat=self.soil_satmoist[0,hh],water_added=self.water[hh+1])
                # infilrate1 = 0.0
                # infilrate2 = 0.0
                # infilrate3 = 0.0
                self.theta[0,hh+1],self.theta[1,hh+1], self.theta[2,hh+1] = self.post_soil_moist(theta1=self.theta[0,hh+1],
                                    theta2=self.theta[1,hh+1],theta3=self.theta[2,hh+1],
                                    theta_sat=self.soil_satmoist[0,hh+1],water_added=self.water[hh+1])

                self.drainagerate[0,hh+1] = np.maximum(infilrate1,water_drainage(theta=self.theta[0,hh+1],
                                            theta_sat=self.soil_satmoist[0,hh+1],Ksat=Ks_sat,fc=theta_fc,layerthickness=zlayers[0]))
                self.drainagerate[1,hh+1] = np.maximum(infilrate2,water_drainage(theta=self.theta[1,hh+1],
                                            theta_sat=self.soil_satmoist[0,hh+1],Ksat=Ks_sat,fc=theta_fc,layerthickness=zlayers[1]))
                self.drainagerate[2,hh+1] = np.maximum(infilrate3,water_drainage(theta=self.theta[2,hh+1],
                                            theta_sat=self.soil_satmoist[1,hh+1],Ksat=Ks_sat,fc=theta_fc,layerthickness=zlayers[2]))
                
                # print('dd',dd,'hh+1',hh+1,"infil 1: ",infilrate1[0,589])
                # print('dd',dd,'hh+1',hh+1,"infil 2: ",infilrate2[0,589])
                # print('dd',dd,'hh+1',hh+1,"infil 3: ",infilrate3[0,589])

                # print('dd',dd,'hh+1',hh+1,"drainage 1: ",self.drainagerate[0,hh+1][0,589])
                # print('dd',dd,'hh+1',hh+1,"drainage 2: ",self.drainagerate[1,hh+1][0,589])
                # print('dd',dd,'hh+1',hh+1,"drainage 3: ",self.drainagerate[2,hh+1][0,589])

                # print('dd',dd,'hh+1',hh+1,"sub : ",subrunoffrate[0,589])

                # print('dd',dd,'hh+1',hh+1,"theta 1: ",self.theta[0,hh+1][0,589],"sw1",self.soil_moist[0,hh+1][0,589])
                # print('dd',dd,'hh+1',hh+1,"theta 2: ",self.theta[1,hh+1][0,589],"sw1",self.soil_moist[0,hh+1][0,589])
                # print('dd',dd,'hh+1',hh+1,"theta 3: ",self.theta[2,hh+1][0,589],"sw2",self.soil_moist[1,hh+1][0,589])
                # print("wp ",theta_wp[0,589])
                # print("irri idx",irrigated_croplandidx[0,589])

                ## resistance for upward diffusion in the surface layer
                Rdiffsrfaq = diff_resistance(distance=pmids[0],phase='aqueous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[0,hh+1])
                Rdiffsrfgas = diff_resistance(distance=pmids[0],phase='gaseous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[0,hh+1])
                
                for ll in np.arange(3):
                    ## lldix: index for soil temp, moisture
                    llidx = int(np.floor(ll/2))
                    
                    # self.theta[ll,hh+1] = np.maximum(self.soil_moist[llidx,hh+1],self.theta[ll,hh+1])
                    # self.theta[ll,hh+1] = np.maximum(theta_wp,self.theta[ll,hh+1])
                    self.drainagerate[ll,hh+1] = np.maximum(self.drainagerate[ll,hh+1],subrunoffrate)
                    # self.drainagerate[ll,hh+1][self.theta[ll,hh+1]==self.soil_moist[llidx,hh+1]] = subrunoffrate[self.theta[ll,hh+1]==self.soil_moist[llidx,hh+1]]
                    ## resistance for diffusion between surface layer (idx0) and topsoil layer (idx1)
                    self.Rdiffaq[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1],temp=self.soil_temp[llidx,hh+1])
                    self.Rdiffgas[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1],temp=self.soil_temp[llidx,hh+1])

                    ## input includes N from fertilizer application
                    if ll == sourcelayer:
                        if sourcelayer == 1:
                            self.urea_pool[0,hh+1] = self.urea_pool[0,hh+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.urea[hh+1]
                            self.urea_pool[1,hh] = self.urea_pool[1,hh] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.urea[hh+1]
                            self.TAN_pool[0,hh+1] = self.TAN_pool[0,hh+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.TAN[hh+1]
                            self.TAN_pool[1,hh] = self.TAN_pool[1,hh] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.TAN[hh+1]
                            self.NO3_pool[0,hh+1] = self.NO3_pool[0,hh+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.NO3[hh+1]
                            self.NO3_pool[1,hh] = self.NO3_pool[1,hh] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.NO3[hh+1]
                        else:
                            self.urea_pool[ll,hh] = self.urea_pool[ll,hh]+self.urea[hh+1]
                            # print("dd",dd,"hh+1",hh+1,"urea, urea pool:",self.urea[hh+1][0,432],self.urea_pool[ll,hh][0,432])
                            self.TAN_pool[ll,hh] = self.TAN_pool[ll,hh]+self.TAN[hh+1]
                            self.NO3_pool[ll,hh] = self.NO3_pool[ll,hh]+self.NO3[hh+1]

                    ## urea scheme
                    if chem_fert_type == 'urea':
                        if ll == 0:
                            ## urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusionup[ll,hh]
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        theta=self.theta[ll,hh+1],delta_t=timestep,k_h=0.03)
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],
                                                                        theta=self.theta[ll,hh+1])
                            ## urea concentration at the compensation point
                            ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,hh+1],rliq=Rdiffsrfaq,
                                                            qrunoff=surfrunoffrate)
                            ## determine the potential of each flux
                            ureawashoffidx = surf_runoff(N_surfcnc=ureasurfamount,
                                                         qrunoff=surfrunoffrate)*timestep*3600
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=self.urea_amount[ll+1,hh],
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
                            correction_ureadiffdown = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                            ureadiffdownidx = ureadiffdownidx*correction_ureadiffdown
                            ureadiffupidx = False
                            srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,hh+1],
                                        flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                            self.ureawashoff[hh+1] = srfrunoffidx
                            self.ureainfil[ll,hh+1] = subsrfleachingidx
                            self.ureadiffusiondown[ll,hh+1] = diffaqdownidx  
                            ## update urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureawashoff[hh+1]-\
                                                    self.ureainfil[ll,hh+1]-self.ureadiffusiondown[ll,hh+1]
                            ## get rid of rounding error
                            self.urea_pool[ll,hh+1] = np.maximum(self.urea_pool[ll,hh+1],0.0)
                        elif ll == 1:
                            ## urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusionup[ll,hh]+\
                                                        self.ureadiffusiondown[ll-1,hh+1]+self.ureainfil[ll-1,hh+1]
                            # print("dd",dd,"hh+1",hh+1,"urea pool",self.urea_pool[ll,hh+1][0,589])
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        theta=self.theta[ll,hh+1],delta_t=timestep,k_h=0.03)
                            
                            # print("dd",dd,"hh+1",hh+1,"soil temp",self.soil_temp[llidx,hh+1][0,589])
                            # print("dd",dd,"hh+1",hh+1,"soil moisture",self.theta[ll,hh+1][0,589])
                            # print("dd",dd,"hh+1",hh+1,"hydro rate",urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                            #                                                     theta=self.theta[ll,hh+1],delta_t=timestep,k_h=0.03)[0,589])
                            # print("dd",dd,"hh+1",hh+1,"urea hydro",self.ureahydrolysis[ll,hh+1][0,589])
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[ll,hh+1])
                            ## determine the potential of each flux
                            ureawashoffidx = False
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=self.urea_amount[ll+1,hh],
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
                            correction_ureadiff = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                            ureadiffdownidx = ureadiffdownidx*correction_ureadiff
                            ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,hh],cnc2=self.urea_amount[ll-1,hh+1],
                                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600
                            ureadiffupidx = ureadiffupidx*correction_ureadiff
                            srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,hh+1],
                                        flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                            self.ureainfil[ll,hh+1] = subsrfleachingidx
                            self.ureadiffusiondown[ll,hh+1] = diffaqdownidx  
                            self.ureadiffusionup[ll-1,hh+1] = diffaqupidx
                            ## update urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureadiffusionup[ll-1,hh+1]-\
                                                        self.ureainfil[ll,hh+1]-self.ureadiffusiondown[ll,hh+1]
                            ## get rid of rounding error
                            self.urea_pool[ll,hh+1] = np.maximum(self.urea_pool[ll,hh+1],0.0)
                        elif ll == 2:
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusiondown[ll-1,hh+1]+self.ureainfil[ll-1,hh+1]
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        theta=self.theta[ll,hh+1],delta_t=timestep,k_h=0.03)
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[ll,hh+1])
                            ## determine the potential of each flux
                            ureawashoffidx = False
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=0.0,
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
                            correction_ureadiff = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                            ureadiffdownidx = ureadiffdownidx*correction_ureadiff
                            ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,hh],cnc2=self.urea_amount[ll-1,hh+1],
                                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600
                            ureadiffupidx = ureadiffupidx*correction_ureadiff
                            srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,hh+1],
                                        flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)       
                            self.ureainfil[ll,hh+1] = subsrfleachingidx
                            self.ureadiffusiondown[ll,hh+1] = diffaqdownidx              
                            self.ureadiffusionup[ll-1,hh+1] = diffaqupidx
                            ## update urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureadiffusionup[ll-1,hh+1]-\
                                                        self.ureainfil[ll,hh+1]-self.ureadiffusiondown[ll,hh+1]
                            ## get rid of rounding error
                            self.urea_pool[ll,hh+1] = np.maximum(self.urea_pool[ll,hh+1],0.0)
                    ## TAN production
                    try:
                        TANprod = self.ureahydrolysis[ll,hh+1]+self.orgN_decomp[hh+1]
                    except:
                        pass
                        TANprod = self.ureahydrolysis[ll,hh+1]

                    if ll == 0:
                        ########################
                        ## TAN sim
                        ########################
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+TANprod
                        # print("dd",dd,"hh+1",hh+1,"TAN diff up",self.TANdiffusionup[ll,hh][0,589])
                        # print("dd",dd,"hh+1",hh+1,"NH3 diffup",self.NH3diffusionup[ll,hh][0,589])
                        # print("dd",dd,"hh+1",hh+1,"TAN prod",TANprod[0,589])
                        # print("dd",dd,"hh+1",hh+1,"TAN pool",self.TAN_pool[ll,hh+1][0,589])

                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        # print("dd",dd,"hh+1",hh+1,"TAN pool",self.TAN_pool[ll,hh+1][0,589])
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1])
                        ## TAN concentration at the compensation surface
                        TANsurfamount = surf_TAN_cnc(tan_cnc=self.TAN_amount[ll,hh+1],rliq=Rdiffsrfaq,rgas=Rdiffsrfgas,
                                            knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=surfrunoffrate)
                        NH3surfamount = TANsurfamount*KNH3
                        ## determining the potential of each flux
                        emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                        TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                    # qrunoff=(self.surfrunoffrate[hh+1]))*timestep*3600 
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                        qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        # TANdiffaqdownidx[self.theta[llidx,hh+1]==0]=0.0
                        correction_diffaq = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                        TANdiffaqdownidx = TANdiffaqdownidx * correction_diffaq
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        # NH3diffgasdownidx[self.theta[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        correction_NH3diffgasdown = np.nan_to_num((self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1])/\
                                                                    (self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1]),posinf=0.0,neginf=0.0)
                        NH3diffgasdownidx = NH3diffgasdownidx*correction_NH3diffgasdown
                        TANuptakeidx = False  ## N uptake
                        nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                                flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                        self.NH3flux[hh+1] = nh3volidx
                        self.TANwashoff[hh+1] = srfrunoffidx
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH3flux[hh+1]-self.TANwashoff[hh+1]-\
                                        self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-self.NH3diffusiondown[ll,hh+1]
                        # print("dd",dd,"hh+1",hh+1,"TAN pool",self.TAN_pool[ll,hh+1][0,589])
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1] = np.maximum(self.TAN_pool[ll,hh+1],0.0)
                        # print("dd",dd,"hh+1",hh+1,"TAN pool",self.TAN_pool[ll,hh+1][0,589])
                        ########################
                        ## NO3- sim
                        ########################
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+self.NO3diffusionup[ll,hh]

                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[ll,hh+1])
                        ## NO3 concentration at the compensation surface
                        NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,hh+1],rliq=Rdiffsrfaq/f_DNO3,
                                                qrunoff=surfrunoffrate)
                        ## determining the potential of each flux
                        NO3washoffidx = surf_runoff(N_surfcnc=NO3surfamount,
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## NO3 washoff
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                       qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=self.NO3_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffaqdownidx = NO3diffaqdownidx*correction_diffaq
                        NO3diffupidx = False
                        srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3washoffidx,flux2=NO3infilidx,flux3=NO3diffaqdownidx,flux4=NO3diffupidx)
                        self.NO3washoff[hh+1] = srfrunoffidx
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3washoff[hh+1]-\
                                                    self.NO3infil[ll,hh+1]-self.NO3diffusiondown[ll,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1] = np.maximum(self.NO3_pool[ll,hh+1],0.0)
                        ########################
                        ## soil water sim
                        ########################
                        ## update soil moisture 
                        update_sw = self.theta[ll,hh+1]*zlayers[0] - self.drainagerate[ll,hh+1]*timestep*3600
                        ## exclude initial inifltration
                        update_sw[self.drainagerate[ll,hh+1]==infilrate1] = self.theta[ll,hh+1][self.drainagerate[ll,hh+1]==infilrate1]*zlayers[0]
                        # self.theta[ll,hh+1] = update_sw/zlayers[0]
                        self.theta[ll,hh+1][self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]] = (update_sw[self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]]-\
                                                                    evap[[self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]]]*timestep*3600/1e6)/zlayers[0]
                        # self.theta[ll,hh+1][self.theta[ll,hh+1]<0] = 0.01
                        # self.theta[ll,hh+1] = np.maximum(self.theta[ll,hh+1],0.01)
                        self.theta[ll,hh+1] = np.maximum(self.soil_moist[llidx,hh+1],self.theta[ll,hh+1])


                    elif ll == 1:
                        ########################
                        ## TAN sim
                        ########################
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+TANprod+\
                                                self.TANinfil[ll-1,hh+1]+self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        # print("dd",dd,"hh+1",hh+1,"TAN diff up",self.TANdiffusionup[ll,hh][0,589])
                        # print("dd",dd,"hh+1",hh+1,"NH3 diffup",self.NH3diffusionup[ll,hh][0,589])
                        # print("dd",dd,"hh+1",hh+1,"TAN prod",TANprod[0,589])
                        # print("dd",dd,"hh+1",hh+1,"TAN pool",self.TAN_pool[ll,hh+1][0,589])
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        # TANdiffaqdownidx[self.theta[llidx,hh+1]==0]=0.0
                        correction_diffaq = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                        TANdiffaqdownidx = TANdiffaqdownidx * correction_diffaq
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        # NH3diffgasdownidx[self.theta[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        correction_diffgas = np.nan_to_num((self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1])/\
                                                                    (self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1]),posinf=0.0,neginf=0.0)
                        NH3diffgasdownidx = NH3diffgasdownidx*correction_diffgas
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        # NH3diffgasupidx[self.theta[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0

                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        # TANdiffaqdownidx[self.theta[llidx,hh+1]==0]=0.0
                        TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='nh4')*timestep*3600  ## N uptake
                        TANuptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                                flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        self.TANdiffusionup[ll-1,hh+1] = diffaqupidx
                        self.NH3diffusionup[ll-1,hh+1] = diffgasupidx
                        self.ammNuptake[ll-1,hh+1] = uptakeidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.TANdiffusionup[ll-1,hh+1]-self.NH3diffusionup[ll-1,hh+1]-\
                            self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-self.NH3diffusiondown[ll,hh+1]-self.ammNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1] = np.maximum(self.TAN_pool[ll,hh+1],0.0)
                        ########################
                        ## NO3- sim
                        ########################
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+self.NO3diffusionup[ll,hh]+\
                                                    self.NO3diffusiondown[ll-1,hh+1]+self.NO3infil[ll-1,hh+1]
                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[ll,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                       qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=self.NO3_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffaqdownidx = NO3diffaqdownidx * correction_diffaq
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = NO3diffupidx * correction_diffaq
                        NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='no3')*timestep*3600  ## N uptake
                        NO3uptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        self.NO3diffusionup[ll-1,hh+1] = diffaqupidx
                        self.nitNuptake[ll-1,hh+1] = uptakeidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3infil[ll,hh+1]-self.NO3diffusionup[ll-1,hh+1]-\
                                                    self.NO3diffusiondown[ll,hh+1]-self.nitNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1] = np.maximum(self.NO3_pool[ll,hh+1],0.0)
                        ########################
                        ## soil water sim
                        ########################
                        ## soil water uptaken by crops; soil volumetric water content change, m3/m3
                        self.wateruptake[ll-1,hh+1] = crop_water_uptake(theta=self.theta[ll,hh+1],theta_WP=theta_wp)*timestep*3600
                        ## water uptake only takes place between plant date plus seedling periods and harvest date
                        self.wateruptake[ll-1,hh+1][dd<=plantidx+day_germinate] = 0.0
                        self.wateruptake[ll-1,hh+1][dd>=harvestidx] = 0.0
                        # print("dd",dd,"hh+1",hh+1,"theta 2: ",self.theta[ll,hh+1][0,589])
                        self.theta[ll,hh+1] = self.theta[ll,hh+1] - self.wateruptake[ll-1,hh+1]
                        ## update soil moisture 
                        update_sw = self.theta[ll,hh+1]*zlayers[1] - (self.drainagerate[ll,hh+1]-self.drainagerate[ll-1,hh+1])*timestep*3600
                        update_sw[self.drainagerate[ll,hh+1]==infilrate2] = self.theta[ll,hh+1][self.drainagerate[ll,hh+1]==infilrate2]*zlayers[1]
                        self.theta[ll,hh+1][self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]] = (update_sw[self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]])/zlayers[1]
                        # print("dd",dd,"hh+1",hh+1,"theta 2: ",self.theta[ll,hh+1][0,589])

                    elif ll == 2:
                        ########################
                        ## TAN sim
                        ########################
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+TANprod+self.TANinfil[ll-1,hh+1]+\
                                                    self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[ll,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[ll,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        # TANdiffaqdownidx[self.theta[llidx,hh+1]==0]=0.0
                        correction_diffaq = np.nan_to_num(self.theta[llidx,hh+1]/self.theta[llidx,hh+1],posinf=0.0,neginf=0.0)
                        TANdiffaqdownidx = TANdiffaqdownidx * correction_diffaq
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        # NH3diffgasdownidx[self.theta[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        correction_diffgas = np.nan_to_num((self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1])/\
                                                                    (self.soil_satmoist[llidx,hh+1]-self.theta[llidx,hh+1]),posinf=0.0,neginf=0.0)
                        NH3diffgasdownidx = NH3diffgasdownidx*correction_diffgas
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        # NH3diffgasupidx[self.theta[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        NH3diffgasupidx = NH3diffgasupidx * correction_diffgas
                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        # TANdiffaqupidx[self.soil_moist[llidx,hh+1]==0]=0.0
                        TANdiffaqupidx = TANdiffaqupidx * correction_diffaq
                        TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='nh4')*timestep*3600  ## N uptake
                        TANuptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                                flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        self.TANdiffusionup[ll-1,hh+1] = diffaqupidx
                        self.NH3diffusionup[ll-1,hh+1] = diffgasupidx
                        self.ammNuptake[ll-1,hh+1] = uptakeidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.TANdiffusionup[ll-1,hh+1]-self.NH3diffusionup[ll-1,hh+1]-\
                            self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-self.NH3diffusiondown[ll,hh+1]-self.ammNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1] = np.maximum(self.TAN_pool[ll,hh+1],0.0)
                        ########################
                        ## NO3 sim
                        ########################
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+\
                                                    self.NO3diffusiondown[ll-1,hh+1]+self.NO3infil[ll-1,hh+1]
                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[ll,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                       qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffaqdownidx = NO3diffaqdownidx * correction_diffaq
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = NO3diffupidx * correction_diffaq
                        NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='no3')*timestep*3600  ## N uptake
                        NO3uptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        self.NO3diffusionup[ll-1,hh+1] = diffaqupidx
                        self.nitNuptake[ll-1,hh+1] = uptakeidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3infil[ll,hh+1]-self.NO3diffusionup[ll-1,hh+1]-\
                                                    self.NO3diffusiondown[ll,hh+1]-self.nitNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1] = np.maximum(self.NO3_pool[ll,hh+1],0.0)
                        ########################
                        ## soil water sim
                        ########################
                        ## soil water uptaken by crops; soil volumetric water content change, m3/m3
                        self.wateruptake[ll-1,hh+1] = crop_water_uptake(theta=self.theta[ll,hh+1],theta_WP=theta_wp)*timestep*3600
                        ## water uptake only takes place between plant date plus seedling periods and harvest date
                        self.wateruptake[ll-1,hh+1][dd<=plantidx+day_germinate] = 0.0
                        self.wateruptake[ll-1,hh+1][dd>=harvestidx] = 0.0
                        self.theta[ll,hh+1] = self.theta[ll,hh+1] - self.wateruptake[ll-1,hh+1]
                        ## update soil moisture 
                        update_sw = self.theta[ll,hh+1]*zlayers[2] - (self.drainagerate[ll,hh+1]-self.drainagerate[ll-1,hh+1])*timestep*3600 
                        update_sw[self.drainagerate[ll,hh+1]==infilrate3] = self.theta[ll,hh+1][self.drainagerate[ll,hh+1]==infilrate3]*zlayers[2]
                        self.theta[ll,hh+1][self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]] = (update_sw[self.theta[ll,hh+1]!=self.soil_moist[llidx,hh+1]])/zlayers[2]

            self.daily_output(dd)
            self.daily_init()            
        return

    ## main sim function of manure fertilizer application
    ## manure_type: "slurry", "solid manure"
    def land_manure_sim(self,start_day_idx,end_day_idx,manure_type,tech,crop=None,sim='base',stvar=False,st=False):
        print('current simulation is for: '+str(manure_type))
        print('technique used is: '+str(tech))
        soilclayds = open_ds(infile_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)[self.plat1:self.plat2,:]
        span_t = 12/timestep
        ## crop calendar that detenmines the N uptake by crops
        if crop is not None:
            cropcalspath = infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
            plantidx, harvestidx = self.crop_calendar(filepath = cropcalspath)
            plantidx = plantidx[self.plat1:self.plat2,:]
            harvestidx = harvestidx[self.plat1:self.plat2,:]
        sourcelayer = self.source_layer(tech)

        for dd in np.arange(start_day_idx,end_day_idx):
            if dd < Days:
                self.sim_env(dd)
                self.TAN[5] = self.TAN_added[dd]
                self.avail_N[5] = self.avail_N_added[dd]
                self.resist_N[5] = self.resist_N_added[dd]
                self.unavail_N[5] = self.unavail_N_added[dd]
                self.water[5] = self.water_added[dd]
                self.manure[5] = self.manure_added[dd]
                sim_pH = self.pH[dd]
                sim_ccH = self.cc_H[dd]

            else:
                self.sim_env(dd-Days)
                self.TAN[5] = self.TAN_added[dd-Days]
                self.avail_N[5] = self.avail_N_added[dd-Days]
                self.resist_N[5] = self.resist_N_added[dd-Days]
                self.unavail_N[5] = self.unavail_N_added[dd-Days]
                self.water[5] = self.water_added[dd-Days]
                self.manure[5] = self.manure_added[dd-Days]
                sim_pH = self.pH[dd-Days]
                sim_ccH = self.cc_H[dd-Days]

            if sim != 'base':
                self.sensitivity_test(var=stvar,test=st)
            
            self.theta[0] = np.copy(self.soil_moist[0])
            self.theta[1] = np.copy(self.soil_moist[0])
            self.theta[2] = np.copy(self.soil_moist[1])
            # self.theta[0,0] = np.copy(self.soil_moist[0,1])
            # self.theta[1,0] = np.copy(self.soil_moist[0,1])
            # self.theta[2,0] = np.copy(self.soil_moist[1,1])

            for hh in np.arange(0,24):
                if CONFIG_machine == "STREAM":
                    surfrunoffrate = self.surfrunoffrate[hh+1]
                    subrunoffrate = self.subrunoffrate[hh+1] 
                    evap = self.evap[hh+1] 
                else:
                    surfrunoffrate = self.surfrunoffrate[hh+1]-self.surfrunoffrate[hh]
                    surfrunoffrate = np.maximum(surfrunoffrate,0.0)
                    subrunoffrate = self.subrunoffrate[hh+1]-self.subrunoffrate[hh]
                    subrunoffrate = np.maximum(subrunoffrate,0.0)
                    evap = self.evap[hh+1]-self.evap[hh]
                for ll in np.arange(3):
                    ## lldix: index for soil temp, moisture
                    llidx = int(np.floor(ll/2))
                    ## input includes N from fertilizer application
                    if ll == sourcelayer:
                        
                        self.avail_N_pool[hh+1] = self.avail_N_pool[hh]+self.avail_N[hh+1]
                        self.resist_N_pool[hh+1] = self.resist_N_pool[hh]+self.resist_N[hh+1]
                        self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh]+self.unavail_N[hh+1]
                        Na_decomp_rate, Nr_decomp_rate = N_pools_decomp_rate(temp=self.soil_temp[llidx,hh+1],delta_t=timestep)
                        ## org N decomposition contributes to TAN pool at a specific soil layer
                        self.orgN_decomp[ll,hh+1] = self.avail_N_pool[hh+1]*Na_decomp_rate + \
                                                    self.resist_N_pool[hh+1]*Nr_decomp_rate
                        if ll == 0:
                            N_washoff_rate = f_washoff_N*surfrunoffrate*3600*timestep*1e3
                            self.avail_N_washoff[hh+1] = N_washoff_rate*self.avail_N_pool[hh+1]
                            self.resist_N_washoff[hh+1] = N_washoff_rate*self.resist_N_pool[hh+1]
                            self.unavail_N_washoff[hh+1] = N_washoff_rate*self.unavail_N_pool[hh+1]
                        self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1]-(self.avail_N_pool[hh+1]*Na_decomp_rate)-\
                                                    self.avail_N_washoff[hh+1]
                        self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1]-(self.resist_N_pool[hh+1]*Nr_decomp_rate)-\
                                                    self.resist_N_washoff[hh+1]
                        self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh+1] - self.unavail_N_washoff[hh+1]
                                        
                        ## water and TAN inital distribution after slurry application
                        if ll == 0:
                            ## washoff rates for orgN and non-N species (manure)
                            nonN_washoff_rate = f_washoff_nonN*surfrunoffrate*3600*timestep*1e3

                            ## slurry TAN pool and water pool
                            self.slurry_TAN_pool[hh+1] = self.slurry_TAN_pool[hh]+self.TAN[hh+1]
                            self.slurry_water_pool[hh+1] = self.slurry_water_pool[hh]+self.water[hh+1]
                            ## manure pool; slurry DM content on the day of application
                            self.manure_pool[hh+1] = self.manure_pool[hh]+self.manure[hh+1]
                            frac_DM = manure_DM(solidmass=self.manure[hh+1],watermass=self.water[hh+1])
                            frac_DM[self.manure[hh+1]==0] = 0.0
                            ## initial "instant" infiltration
                            init_infil = initial_infil(frac_DM)
                            init_infil[frac_DM==0.0] = 0.0
                            init_waterinfil = init_infil * self.water[hh+1]
                            zsat = (init_waterinfil/(self.soil_satmoist[llidx,hh] - self.soil_moist[llidx,hh]))/1e6
                            self.theta[0,hh+1] = (self.theta[0,hh]*zlayers[0] + init_waterinfil/1e6)/zlayers[0]
                            # self.theta[0,hh+1] = (self.theta[0,hh]*zlayers[0] + init_waterinfil/1e6)/zlayers[0]
                            self.theta[0,hh+1][zsat>zlayers[0]] = self.soil_satmoist[0,hh+1][zsat>zlayers[0]]
                            self.theta[1,hh+1] = (self.theta[1,hh]*zlayers[1] + (zsat - zlayers[0]))/zlayers[1]
                            self.theta[1,hh+1][zsat<zlayers[0]] = self.soil_moist[0,hh+1][zsat<zlayers[0]]
                            
                            ## initial inifiltration of TAN from the slurry to the topsoil
                            TAN_infil = self.slurry_TAN_pool[hh+1] * init_infil
                            ## infiltrated TAN to the first soil layer
                            TAN_infil_z0 = TAN_infil*((self.theta[0,hh+1] - self.soil_moist[llidx,hh])*zlayers[0]*1e6/init_waterinfil)
                            TAN_infil_z0[init_infil==0.0] = 0.0
                            self.TAN_pool[0,hh] = self.TAN_pool[0,hh] + TAN_infil_z0
                            self.TAN_pool[1,hh]= self.TAN_pool[1,hh] + (TAN_infil - TAN_infil_z0)
                            self.slurry_TAN_pool[hh+1] = self.slurry_TAN_pool[hh+1] - TAN_infil
                            self.slurry_water_pool[hh+1] = self.slurry_water_pool[hh+1] - init_waterinfil
                            self.manure_washoff[hh+1] = nonN_washoff_rate*self.manure_pool[hh+1]
                            self.manure_pool[hh+1] = self.manure_pool[hh+1] - self.manure_washoff[hh+1]
                            
                        if ll == 1:
                            self.theta[0,hh+1] = (self.theta[0,hh+1]*(zlayers[0]+zlayers[1]) + self.water[hh+1]/1e6)/(zlayers[0]+zlayers[1])
                            self.theta[0,hh+1][self.theta[0,hh+1]>self.soil_satmoist[llidx,hh+1]] = self.soil_satmoist[llidx,hh+1][self.theta[0,hh+1] >self.soil_satmoist[llidx,hh+1]]
                            self.theta[1,hh+1] = (self.theta[1,hh+1]*(zlayers[0]+zlayers[1]) + self.water[hh+1]/1e6)/(zlayers[0]+zlayers[1])
                            self.theta[1,hh+1][self.theta[1,hh+1]>self.soil_satmoist[llidx,hh+1]] = self.soil_satmoist[llidx,hh+1][self.theta[1,hh+1] >self.soil_satmoist[llidx,hh+1]]
                            self.TAN_pool[0,hh+1] = self.TAN_pool[0,hh+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.TAN[hh+1]
                            self.TAN_pool[1,hh] = self.TAN_pool[1,hh] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.TAN[hh+1]
                            self.NO3_pool[0,hh+1] = self.NO3_pool[0,hh+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.NO3[hh+1]
                            self.NO3_pool[1,hh] = self.NO3_pool[1,hh] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.NO3[hh+1]
                            
                        if ll == 2:
                            self.theta[2,hh+1]  = (self.theta[2,hh+1]*zlayers[2] + self.water[hh+1]/1e6)/zlayers[2]
                            self.theta[2,hh+1][self.theta[2,hh+1]>self.soil_satmoist[llidx,hh+1]] = self.soil_satmoist[llidx,hh+1][self.theta[2,hh+1]>self.soil_satmoist[llidx,hh+1]]
                            self.TAN_pool[2,hh] = self.TAN_pool[2,hh]+self.TAN[hh+1]
                            self.NO3_pool[2,hh] = self.NO3_pool[2,hh]+self.NO3[hh+1]
                        
                    # TANadd = self.ureahydrolysis[ll,hh+1] + self.orgN_decomp[ll,hh+1]
                    TANadd = self.ureahydrolysis[ll,hh+1]  
                    
                    if ll == 0:
                        if ll != sourcelayer:
                            # self.theta[0] = self.soil_moist[llidx,hh+1]
                            self.drainagerate[ll,hh+1] = subrunoffrate
                            # subrunoffrate1 = subrunoffrate
                            
                        ## resistance for upward diffusion in the surface layer
                        Rdiffsrfaq = diff_resistance(distance=pmids[0],phase='aqueous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[llidx,hh+1])
                        Rdiffsrfgas = diff_resistance(distance=pmids[0],phase='gaseous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[llidx,hh+1])
                        self.Rdiffaq[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[llidx,hh+1])
                        self.Rdiffgas[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[0,hh+1],temp=self.soil_temp[llidx,hh+1])
                        
                        if sourcelayer == 0:
                            ## slurry pool at the surface
                            frac_DM = manure_DM(solidmass=self.manure_pool[hh+1],watermass=self.slurry_water_pool[hh+1])
                            slurry_infilrate = infil_rate_slurry(frac_DM)
                            slurry_infil = slurry_infilrate*3600*timestep
                            minwater = self.manure_pool[hh+1] * absorb_factor
                            infil_idx = (self.slurry_water_pool[hh+1] - slurry_infil) - minwater
                            slurry_infil[infil_idx<0] = self.slurry_water_pool[hh+1][infil_idx<0] - minwater[infil_idx<0]
                            slurry_infilrate = slurry_infil/(3600*timestep)
                            
                            zslurry = self.slurry_water_pool[hh+1]/1e6
                            
                            KNH3 = NH3_par_coeff(temp=self.soil_temp[0,hh+1],cncH=self.slurry_ccH)
    #                         self.slurry_TAN_pool[hh+1] = self.slurry_TAN_pool[hh]+self.slurry_TANdiffusionup[hh] 
                            self.slurry_TAN_amount[hh+1] = self.slurry_TAN_pool[hh+1]/self.slurry_water_pool[hh+1]*1e6
                            self.slurry_TAN_amount[hh+1][self.slurry_water_pool[hh+1]==0] = 0.0
                            self.slurry_NH3_gas[hh+1] = KNH3*self.slurry_TAN_amount[hh+1]
                            
                            emissidx = NH3_vol(nh3_surfcnc=self.slurry_NH3_gas[hh+1],ratm=self.R_atm[hh+1])*timestep*3600
                            TANwashoffidx = surf_runoff(N_surfcnc=self.slurry_TAN_amount[hh+1],
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                            TANinfilidx = subsurf_leaching(N_cnc=self.slurry_TAN_amount[hh+1],
                                                            qsubrunoff=slurry_infilrate)*timestep*3600  ## TAN infiltration/leaching
                            TANdiffaqdownidx = N_diffusion(cnc1=self.slurry_TAN_amount[hh+1],cnc2=self.TAN_amount[0,hh],
                                                resist=Rdiffsrfaq)*timestep*3600  ## TAN aqueous downwards diffusion
                            
                            nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx = N_pathways(mN=self.slurry_TAN_pool[hh+1],
                                        flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx)
                            
                            self.NH3flux[hh+1] = nh3volidx
                            self.TANwashoff[hh+1] = srfrunoffidx
                            self.slurry_TANinfil[hh+1] = subsrfleachingidx
                            self.slurry_TANdiffusiondown[hh+1] = diffaqdownidx
                        
                            self.slurry_TAN_pool[hh+1] = self.slurry_TAN_pool[hh+1] - self.NH3flux[hh+1] - self.TANwashoff[hh+1] -\
                                                        self.slurry_TANinfil[hh+1] - self.slurry_TANdiffusiondown[hh+1]
                            self.slurry_water_pool[hh+1] = self.slurry_water_pool[hh+1] - evap*timestep*3600 - slurry_infil
                            
                            ## when water pool of the slurry is less than the min water threshold
                            ## slurry TAN is set to be added to the surface soil layer
                            ## slurry water, manure is reset to zero.
                            water_idx = self.slurry_water_pool[hh+1] - minwater
                            self.slurry_TANinfil[hh+1][water_idx<1.0] = self.slurry_TANinfil[hh+1][water_idx<1.0] + \
                                                                            self.slurry_TAN_pool[hh+1][water_idx<1.0]
                            self.slurry_TAN_pool[hh+1][water_idx<1.0] = 0.0
                            self.slurry_water_pool[hh+1][water_idx<1.0] = 0.0
                            self.manure_pool[hh+1][water_idx<1.0] = 0.0
                            
                            ## soil water content is nudged by reanalysis data for soil moisture
                            theta1_recovery = (self.theta[0,hh+1] - self.soil_moist[0,hh+1])/span_t
                            ## subsurface percolation flux includes correction for infiltration  
                            ## and nudged soil moisture to conserve water mass
                            # subrunoffrate1 = subrunoffrate + slurry_infilrate + theta1_recovery/(3600*timestep)
                            self.drainagerate[0,hh+1] = subrunoffrate + slurry_infilrate

                        theta1_recovery = (self.theta[0,hh+1] - self.soil_moist[0,hh+1])/span_t
                        self.theta[0,hh+1] = self.theta[0,hh+1] - theta1_recovery                       
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+\
                                                    self.slurry_TANinfil[hh+1]+self.slurry_TANdiffusiondown[hh+1]+\
                                                        self.orgN_decomp[ll,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[0,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[0,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='manure',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[0,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[0,hh+1])
                        ## TAN concentration at the compensation surface
                        TANsurfamount = surf_TAN_cnc(tan_cnc=self.TAN_amount[ll,hh+1],rliq=Rdiffsrfaq,rgas=Rdiffsrfgas,
                                            knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=surfrunoffrate)
                                            # knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=self.surfrunoffrate[hh+1])
                        NH3surfamount = TANsurfamount*KNH3
                        ## determining the potential of each flux
                        emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                        ## NH3 from soil layer is zero when slurry TAN pool is not empty
                        emissidx[self.slurry_TAN_pool[hh]!=0] = 0.0
                        TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                        ## surface washoff from the soil layer is zero when slurry TAN pool is not empty
                        TANwashoffidx[self.slurry_TAN_pool[hh]!=0] = 0.0
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                        qsubrunoff=self.drainagerate[0,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.theta[0,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.theta[0,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        TANuptakeidx = False  ## N uptake
                        nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                                flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                        ## when slurry TAN pool is depleted to zero, NH3 is released from the top soil layer
                        self.NH3flux[hh+1][self.slurry_TAN_pool[hh+1]==0] = nh3volidx[self.slurry_TAN_pool[hh+1]==0]
                        self.TANwashoff[hh+1][self.slurry_TAN_pool[hh+1]==0] = srfrunoffidx[self.slurry_TAN_pool[hh+1]==0]
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-\
                            self.NH3diffusiondown[ll,hh+1]
                        self.TAN_pool[ll,hh+1][self.slurry_TAN_pool[hh]==0] = self.TAN_pool[ll,hh+1][self.slurry_TAN_pool[hh]==0]-\
                                                self.NH3flux[hh+1][self.slurry_TAN_pool[hh]==0]-self.TANwashoff[hh+1][self.slurry_TAN_pool[hh]==0]
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0
                    elif ll == 1:
                        # if ll != sourcelayer:
                        #     self.theta[1] = self.soil_moist[0,hh+1]
                            
                        theta2_recovery = (self.theta[1,hh+1] - self.soil_moist[llidx,hh+1])/span_t
                        # subrunoffrate2 = subrunoffrate1 + theta2_recovery/(3600*timestep)
                        self.drainagerate[ll,hh+1] = subrunoffrate
                        self.theta[1,hh+1] = self.theta[1,hh+1] - theta2_recovery
                        # print(dd,"midsoil",self.theta[1,hh+1,12,377],"sm",self.soil_moist[llidx,hh+1,12,377])    
                        self.Rdiffaq[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[1,hh+1],temp=self.soil_temp[llidx,hh+1])
                        self.Rdiffgas[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[1,hh+1],temp=self.soil_temp[llidx,hh+1])
                    
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+self.orgN_decomp[ll,hh+1]+\
                                                self.TANinfil[ll-1,hh+1]+self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[1,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[1,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='manure',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[1,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[1,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.theta[1,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.theta[1,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        NH3diffgasupidx[self.theta[1,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        TANdiffaqdownidx[self.theta[1,hh+1]==0]=0.0
                        TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='nh4')*timestep*3600  ## N uptake
                        # TANuptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                                flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        self.TANdiffusionup[ll-1,hh+1] = diffaqupidx
                        self.NH3diffusionup[ll-1,hh+1] = diffgasupidx
                        self.ammNuptake[ll-1,hh+1] = uptakeidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.TANdiffusionup[ll-1,hh+1]-self.NH3diffusionup[ll-1,hh+1]-\
                            self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-self.NH3diffusiondown[ll,hh+1]-self.ammNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0
                        
                        
                    elif ll == 2:
                        theta3_recovery = (self.theta[2,hh+1] - self.soil_moist[llidx,hh+1])/span_t
                        # subrunoffrate3 = subrunoffrate2 + theta3_recovery/(3600*timestep)
                        self.drainagerate[ll,hh+1] = subrunoffrate
                        self.theta[2,hh+1] = self.theta[2,hh+1] - theta3_recovery 
                        self.Rdiffaq[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[2,hh+1],temp=self.soil_temp[llidx,hh+1])
                        self.Rdiffgas[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[2,hh+1],temp=self.soil_temp[llidx,hh+1])
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.orgN_decomp[ll,hh+1]+self.TANinfil[ll-1,hh+1]+\
                                                    self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.theta[2,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.theta[2,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='manure',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[2,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.theta[2,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.theta[2,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.theta[2,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        NH3diffgasupidx[self.theta[2,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        TANdiffaqupidx[self.theta[2,hh+1]==0]=0.0
                        TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='nh4')*timestep*3600  ## N uptake
                        # TANuptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,hh+1],
                                flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                                flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                        self.TANinfil[ll,hh+1] = subsrfleachingidx
                        self.TANdiffusiondown[ll,hh+1] = diffaqdownidx
                        self.NH3diffusiondown[ll,hh+1] = diffgasdownidx
                        self.TANdiffusionup[ll-1,hh+1] = diffaqupidx
                        self.NH3diffusionup[ll-1,hh+1] = diffgasupidx
                        self.ammNuptake[ll-1,hh+1] = uptakeidx
                        ## update TAN pool: subtracting all losses
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.TANdiffusionup[ll-1,hh+1]-self.NH3diffusionup[ll-1,hh+1]-\
                            self.TANinfil[ll,hh+1]-self.TANdiffusiondown[ll,hh+1]-self.NH3diffusiondown[ll,hh+1]-self.ammNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0

                    ## NO3 scheme
                    if ll == 0:
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+self.NO3diffusionup[ll,hh]

                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[0,hh+1])
                        ## NO3 concentration at the compensation surface
                        NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,hh+1],rliq=Rdiffsrfaq/f_DNO3,
                                                qrunoff=surfrunoffrate)
                        # NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,hh+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[hh+1])
                        ## determining the potential of each flux
                        NO3washoffidx = surf_runoff(N_surfcnc=NO3surfamount,
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## NO3 washoff
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                        # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                    qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=self.NO3_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = False
                        srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3washoffidx,flux2=NO3infilidx,flux3=NO3diffaqdownidx,flux4=NO3diffupidx)
                        self.NO3washoff[hh+1] = srfrunoffidx
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3washoff[hh+1]-\
                                                    self.NO3infil[ll,hh+1]-self.NO3diffusiondown[ll,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1][self.NO3_pool[ll,hh+1]<0.0] = 0.0
                    elif ll == 1:
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+self.NO3diffusionup[ll,hh]+\
                                                    self.NO3diffusiondown[ll-1,hh+1]+self.NO3infil[ll-1,hh+1]
                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[1,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                    qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=self.NO3_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='no3')*timestep*3600  ## N uptake
                        # NO3uptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        self.NO3diffusionup[ll-1,hh+1] = diffaqupidx
                        self.nitNuptake[ll-1,hh+1] = uptakeidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3infil[ll,hh+1]-self.NO3diffusionup[ll-1,hh+1]-\
                                                    self.NO3diffusiondown[ll,hh+1]-self.nitNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1][self.NO3_pool[ll,hh+1]<0.0] = 0.0
                    elif ll == 2:
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+\
                                                    self.NO3diffusiondown[ll-1,hh+1]+self.NO3infil[ll-1,hh+1]
                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.theta[2,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                    qsubrunoff=self.drainagerate[ll,hh+1])*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,hh+1]*fNH4,mNO3=self.NO3_pool[ll,hh],
                                        temp=self.soil_temp[llidx,hh+1],uptake='no3')*timestep*3600  ## N uptake
                        # NO3uptakeidx[harvestidx<(hh+1)] = 0.0
                        subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,hh+1],
                                        flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                        self.NO3infil[ll,hh+1] = subsrfleachingidx
                        self.NO3diffusiondown[ll,hh+1] = diffaqdownidx
                        self.NO3diffusionup[ll-1,hh+1] = diffaqupidx
                        self.nitNuptake[ll-1,hh+1] = uptakeidx
                        ## update NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh+1]-self.NO3infil[ll,hh+1]-self.NO3diffusionup[ll-1,hh+1]-\
                                                    self.NO3diffusiondown[ll,hh+1]-self.nitNuptake[ll-1,hh+1]
                        ## get rid of rounding error
                        self.NO3_pool[ll,hh+1][self.NO3_pool[ll,hh+1]<0.0] = 0.0
            self.daily_output(dd)
            self.daily_init(manure=True)            
        return

    ## main sim function for grazing
    def grazing_sim(self,livestock_name,start_day_idx,end_day_idx,span_day=60,sim='base',stvar=False,st=False):
        print('current simulation is for: '+str(livestock_name)+" grazing")
        soilclayds = open_ds(infile_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        soilsandds = open_ds(infile_path+soil_data_path+soilsandfile)
        soilsand = soilsandds.T_SAND.values
        soilsiltds = open_ds(infile_path+soil_data_path+soilsiltfile)
        soilsilt = soilsiltds.T_SILT.values
        soilocds = open_ds(infile_path+soil_data_path+soilorgCfile)
        soiloc = soilocds.T_OC.values
        soilbd_ds = open_ds(infile_path+soil_data_path+soilbdfile)
        soilbd = soilbd_ds.T_BULK_DEN.values
        soilpHds = open_ds(infile_path+soil_data_path+soilpHfile)
        soilph = np.zeros(CONFIG_mtrx2)
        soilph[:] = soilpHds.T_PH_H2O.values
        Kd = ammonium_adsorption(clay_content=soilclay)[self.plat1:self.plat2,:]
        Kd_manure = 1.0
        Ks_sat = soil_hydraulic_conductivity(fsilt=soilsilt,fclay=soilclay,fsom=soiloc,BD=soilbd)[self.plat1:self.plat2,:]
        theta_FC = soil_field_capacity(BD=soilbd)[self.plat1:self.plat2,:]
        theta_wp = soil_wilting_point(fsand=soilsand,fclay=soilclay)[self.plat1:self.plat2,:]

        animal_file_name = CONFIG_animal_file_dict[livestock_name]
        livestockds = open_ds(infile_path+animal_data_path+animal_file_name)
        ## lvl idx 0 is the [Grassland] production system; see config.py
        excretN_info = livestockds['Excreted_N'][0][self.plat1:self.plat2,:]
        print("Total N from grazing "+str(livestock_name)+" is: "+str(np.nansum(excretN_info*1e3/1e9))+" Gg.")
        animal_head = livestockds['Animal_head'][0][self.plat1:self.plat2,:]
        field_area = animal_head * grazing_density[livestock_name]
        # source_area = (1-np.exp(-animal_head*patch_area/field_area))*field_area
        source_area = field_area * (f_urine_patch+f_dung_pat)
        self.pastarea = source_area
        excret_N = excretN_info/source_area
        excret_N = xr_to_np(excret_N)
        durine_N, durea, dmanure_N, durine, dmanure, manure_wc,sim_pH = livestock_waste_info(livestock_type=livestock_name, 
                        waste_N=excret_N)
        test_N = (durine_N+dmanure_N)*source_area
        print("test total N ",np.nansum(test_N)*365*24/1e9)
        self.urea_added = durea*source_area
        self.avail_N_added = ((durine_N - durea) + dmanure_N)*f_avail*source_area
        self.resist_N_added = ((durine_N - durea) + dmanure_N)*f_resist*source_area
        self.unavail_N_added = ((durine_N - durea) + dmanure_N)*f_unavail*source_area
        ## input of manure N pool at each timestep, manure N only falls on dung pats
        durine_N = durine_N/((f_urine_patch+f_dung_pat*frac_FYM)/(f_urine_patch+f_dung_pat))
        durea = durea/((f_urine_patch+f_dung_pat*frac_FYM)/(f_urine_patch+f_dung_pat))
        dmanure_N = dmanure_N/(f_dung_pat/(f_urine_patch+f_dung_pat))
        durine = durine/((f_urine_patch+f_dung_pat*frac_FYM)/(f_urine_patch+f_dung_pat))
        dmanure = dmanure/(f_dung_pat/f_urine_patch+f_dung_pat)
        manure_wc = manure_wc/(f_dung_pat/f_urine_patch+f_dung_pat)
        test_N = durine_N*source_area*(f_urine_patch+frac_FYM*f_dung_pat)/(f_urine_patch+f_dung_pat)+\
                    dmanure_N*source_area*(f_dung_pat)/(f_urine_patch+f_dung_pat)
        print("test total N ",np.nansum(test_N)*365*24/1e9)
        durine = durine * 1000
        manure_wc = manure_wc * 1000
        ## pH and H+ ions concentration
        sim_ccH = np.float(10**(-sim_pH))

        for dd in np.arange(start_day_idx,end_day_idx):
            ## 
            for dd_span in np.arange(dd,dd+span_day):

                if dd_span < Days:
                    self.sim_env(dd_span)
                else:
                    self.sim_env(int(dd_span-np.floor(dd_span/Days)*Days))
                self.theta[0,1:] = np.copy(self.soil_moist[0,1:])
                self.theta[1,1:] = np.copy(self.soil_moist[0,1:]) 

                ## manure and N deposited on the field on the first day of each span loop
                if dd_span == dd:
                    self.manure[:] = dmanure
                    self.urine[:] = durine
                    self.urea[:] = durea
                    self.urine_N[:] = durine_N
                    self.manure_N[:] = dmanure_N
                    self.manure_water[:] = manure_wc

                    ## mark the day when urine deposited
                    N_app_mark = np.zeros(soilph.shape)
                    if dd_span < Days:
                        N_app_mark[dd_span] = 1
                    else:
                        N_app_mark[int(dd_span-np.floor(dd_span/Days)*Days)] = 1
                    self.soil_pH = soil_pH_postapp(base_pH=soilph,app_timing_map=N_app_mark,fert_pH=8.5)[:,self.plat1:self.plat2,:] 
                    self.soil_ccH = 10**(-self.soil_pH) 
                else:
                    self.manure[:] = 0.0
                    self.urine[:] = 0.0
                    self.urea[:] = 0.0
                    self.urine_N[:] = 0.0
                    self.manure_N[:] = 0.0
                    self.manure_water[:] = 0.0

                for hh in np.arange(0,24):
                    if CONFIG_machine == "STREAM":
                        surfrunoffrate = self.surfrunoffrate[hh+1]
                        subrunoffrate = self.subrunoffrate[hh+1] 
                        evap = self.evap[hh+1] 
                    else:
                        surfrunoffrate = self.surfrunoffrate[hh+1]-self.surfrunoffrate[hh]
                        surfrunoffrate = np.maximum(surfrunoffrate,0.0)
                        subrunoffrate = self.subrunoffrate[hh+1]-self.subrunoffrate[hh]
                        subrunoffrate = np.maximum(subrunoffrate,0.0)
                        evap = self.evap[hh+1]-self.evap[hh]
                    ## surface washoff rate of N pools/non N pool
                    nonN_washoff_rate = f_washoff_nonN*surfrunoffrate*3600*timestep*1e3
                    N_washoff_rate = f_washoff_N*surfrunoffrate*3600*timestep*1e3
                    ## urea hydrolysis rate
                    urea_hydro_rate = urea_hydrolysis_rate(temp=self.soil_temp[0,hh+1],theta=1.0,delta_t=timestep)
                    ## org N decomposition rate
                    Na_decomp_rate, Nr_decomp_rate = N_pools_decomp_rate(temp=self.soil_temp[0,hh+1],delta_t=timestep)
                    ## NH3 partition coefficient
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[0,hh+1],cncH=sim_ccH)
                    ################################################################
                    ## dung pat scheme; level index=[0,1] 0-FYM; 1-dung
                    ################################################################
                    ## FYM: urine+feces
                    llidx = 0
                    ## manure fully deposited
                    self.manure_pool_FYM[hh+1] = self.manure_pool_FYM[hh] + self.manure[hh+1]
                    ## a fraction of urine deposited, partially infiltrated to the underlying soil, and partially absorbed by feces
                    init_infil = (self.urine[hh+1]+self.manure_water[hh+1]) - self.manure[hh+1]*absorb_factor
                    init_infil = np.maximum(init_infil,0)
                    init_water = (self.urine[hh+1]+self.manure_water[hh+1]) - init_infil
                    frac_init_infil = init_infil/(durine+manure_wc)
                    self.Total_water_pool_FYM[hh+1] = self.Total_water_pool_FYM[hh] + init_water + (self.rainfall[hh+1]-evap)*timestep*3600
                    mois_coeff = min_manurewc(temp=self.soil_temp[0,hh+1],rhum=50)
                    manure_minwc = self.manure_pool_FYM[hh+1] * mois_coeff
                    self.Total_water_pool_FYM[hh+1] = np.maximum(self.Total_water_pool_FYM[hh+1],manure_minwc)
                    ## thichness of the manure layer
                    z_manure = (self.manure_pool_FYM[hh+1]+self.Total_water_pool_FYM[hh+1])/1e6
                    ## infiltration rate
                    self.drainagerate[llidx,hh+1] = (self.Total_water_pool_FYM[hh+1] - self.manure_pool_FYM[hh+1]*absorb_factor)/(1e6*timestep*3600)
                    self.drainagerate[llidx,hh+1] = np.maximum(self.drainagerate[llidx,hh+1],0)
                    self.Rdiffaq[llidx,hh+1] = ((z_manure+z_source)/2)/diffusivity_NH4(temp=self.soil_temp[0,hh+1],phase="aqueous")
                    ## urea pools of manure
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh] + self.urea[hh+1]*(1-frac_init_infil)
                    ## org N pools
                    avail_N_FYM = (self.manure_N[hh+1] + (self.urine_N[hh+1]-self.urea[hh+1]))*f_avail
                    resist_N_FYM = (self.manure_N[hh+1] + (self.urine_N[hh+1]-self.urea[hh+1]))*f_resist
                    unavail_N_FYM = (self.manure_N[hh+1] + (self.urine_N[hh+1]-self.urea[hh+1]))*f_unavail
                    self.avail_N_pool_FYM[hh+1] = self.avail_N_pool_FYM[hh] + avail_N_FYM
                    self.resist_N_pool_FYM[hh+1] = self.resist_N_pool_FYM[hh] + resist_N_FYM
                    self.unavail_N_pool_FYM[hh+1] = self.unavail_N_pool_FYM[hh] + unavail_N_FYM
                    ## surface washoff of manure pool and org N pools
                    self.manure_washoff_FYM[hh+1] = nonN_washoff_rate*self.manure_pool_FYM[hh+1]
                    self.avail_N_washoff_FYM[hh+1] = N_washoff_rate*self.avail_N_pool_FYM[hh+1]
                    self.resist_N_washoff_FYM[hh+1] = N_washoff_rate*self.resist_N_pool_FYM[hh+1]
                    self.unavail_N_washoff_FYM[hh+1] = N_washoff_rate*self.unavail_N_pool_FYM[hh+1]

                    TANprod_FYM = urea_hydro_rate*self.urea_pool[llidx,hh+1]+ \
                                                Na_decomp_rate*self.avail_N_pool_FYM[hh+1] + \
                                                Nr_decomp_rate*self.resist_N_pool_FYM[hh+1]
                    ## update manure pool and org N pools
                    self.manure_pool_FYM[hh+1] = self.manure_pool_FYM[hh+1] - self.manure_washoff_FYM[hh+1]
                    self.avail_N_pool_FYM[hh+1] = self.avail_N_pool_FYM[hh+1]*(1-Na_decomp_rate) - \
                                                        self.avail_N_washoff_FYM[hh+1]
                    self.resist_N_pool_FYM[hh+1] = self.resist_N_pool_FYM[hh+1]*(1-Nr_decomp_rate) - \
                                                        self.resist_N_washoff_FYM[hh+1]
                    self.unavail_N_pool_FYM[hh+1] = self.unavail_N_pool_FYM[hh+1] - self.unavail_N_washoff_FYM[hh+1]
                    ## update urea and org N pools
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh+1]*(1 - urea_hydro_rate)
                    ## urea concentration
                    self.urea_amount[llidx,hh+1] = self.urea_pool[llidx,hh+1]/(self.Total_water_pool_FYM[hh+1]+self.manure_pool_FYM[hh+1]*Kd_manure/(manure_PD/1e3))
                    self.urea_amount[llidx,hh+1] = self.urea_amount[llidx,hh+1] * 1e6
                    ## urea concentration at the compensation point; set equivalent to the bulk concentration
                    ureasurfamount = self.urea_amount[llidx,hh+1]
                    ## determine the potential of each flux
                    ureawashoffidx = surf_runoff(N_surfcnc=ureasurfamount,
                                                    qrunoff=surfrunoffrate)*timestep*3600
                    ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[llidx,hh+1],
                                                    qsubrunoff=self.drainagerate[llidx,hh+1])*timestep*3600
                    ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[llidx,hh+1],cnc2=0,
                                                    resist=self.Rdiffaq[llidx,hh+1])*timestep*3600
                    ureadiffupidx = False
                    srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[llidx,hh+1],
                                flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                    self.ureawashoff_FYM[hh+1] = srfrunoffidx
                    self.ureainfil[llidx,hh+1] = subsrfleachingidx
                    self.ureadiffusiondown[llidx,hh+1] = diffaqdownidx  
                    ## update urea pool
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh+1]-self.ureawashoff_FYM[hh+1]-\
                                            self.ureainfil[llidx,hh+1]-self.ureadiffusiondown[llidx,hh+1]
                    ## include initial infiltration
                    self.ureainfil[llidx,hh+1] = self.ureainfil[llidx,hh+1] + self.urea[hh+1] * frac_init_infil
                    ## TAN pool
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh] + TANprod_FYM 
                    ## TAN concentration
                    self.TAN_amount[llidx,hh+1] = self.TAN_pool[llidx,hh+1]/(self.Total_water_pool_FYM[hh+1]+self.manure_pool_FYM[hh+1]*Kd_manure/(manure_PD/1e3))
                    self.TAN_amount[llidx,hh+1] = self.TAN_amount[llidx,hh+1]*1e6
                    ## NH3 gas concentration
                    self.NH3_gas[llidx,hh+1] = KNH3*self.TAN_amount[llidx,hh+1]
                    ## TAN concentration at the compensation surface
                    TANsurfamount = self.TAN_amount[llidx,hh+1]
                    NH3surfamount = TANsurfamount*KNH3

                    ## determining the potential of each flux
                    emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                    TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[llidx,hh+1],
                                                    qsubrunoff=self.drainagerate[llidx,hh+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[llidx,hh+1],cnc2=0,
                                        resist=self.Rdiffaq[llidx,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.theta[llidx+1,hh+1]==0]=0.0
                    NH3diffgasdownidx = False
                    TANuptakeidx = False  ## N uptake
                    nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[llidx,hh+1],
                            flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                            flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                    self.NH3flux_FYM[hh+1] = nh3volidx
                    self.TANwashoff_FYM[hh+1] = srfrunoffidx
                    self.TANinfil[llidx,hh+1] = subsrfleachingidx
                    self.TANdiffusiondown[llidx,hh+1] = diffaqdownidx
                    self.NH3diffusiondown[llidx,hh+1] = diffgasdownidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh+1]-self.NH3flux_FYM[hh+1]-self.TANwashoff_FYM[hh+1]-\
                                    self.TANinfil[llidx,hh+1]-self.TANdiffusiondown[llidx,hh+1]-self.NH3diffusiondown[llidx,hh+1]
                    ## get rid of rounding error
                    self.TAN_pool[llidx,hh+1] = np.maximum(self.TAN_pool[llidx,hh+1],0.0)
                    ## updating manure water pool
                    self.Total_water_pool_FYM[hh+1] = self.Total_water_pool_FYM[hh+1] - self.drainagerate[llidx,hh+1]*timestep*3600*1e6
                    ## dung: feces only, no urine
                    llidx = 1
                    ## manure pool
                    self.manure_pool_dung[hh+1] = self.manure_pool[hh] + self.manure[hh+1]
                    init_infil = self.manure_water[hh+1] - self.manure[hh+1]*absorb_factor
                    init_infil = np.maximum(init_infil,0)
                    init_water = self.manure_water[hh+1] - init_infil
                    frac_init_infil = init_infil/manure_wc
                    self.Total_water_pool_dung[hh+1] = self.Total_water_pool_dung[hh] + init_water + (self.rainfall[hh+1]-evap)*timestep*3600
                    manure_minwc = self.manure_pool_dung[hh+1] * mois_coeff
                    self.Total_water_pool_dung[hh+1] = np.maximum(self.Total_water_pool_dung[hh+1],manure_minwc)
                    ## thichness of the manure layer
                    z_manure = (self.manure_pool_dung[hh+1]+self.Total_water_pool_dung[hh+1])/1e6
                    ## infiltration rate
                    self.drainagerate[llidx,hh+1] = (self.Total_water_pool_dung[hh+1] - self.manure_pool_dung[hh+1]*absorb_factor)/(1e6*timestep*3600)
                    self.drainagerate[llidx,hh+1] = np.maximum(self.drainagerate[llidx,hh+1],0)
                    ## resistances
                    self.Rdiffaq[llidx,hh+1] = ((z_manure+z_source)/2)/diffusivity_NH4(temp=self.soil_temp[0,hh+1],phase="aqueous")
                    avail_N_dung = self.manure_N[hh+1] *f_avail
                    resist_N_dung = self.manure_N[hh+1] *f_resist
                    unavail_N_dung = self.manure_N[hh+1] *f_unavail
                    self.avail_N_pool_dung[hh+1] = self.avail_N_pool_dung[hh] + avail_N_dung
                    self.resist_N_pool_dung[hh+1] = self.resist_N_pool_dung[hh] + resist_N_dung
                    self.unavail_N_pool_dung[hh+1] = self.unavail_N_pool_dung[hh] + unavail_N_dung
                    ## surface washoff of manure pool and org N pools
                    self.manure_washoff_dung[hh+1] = nonN_washoff_rate*self.manure_pool_dung[hh+1]
                    self.avail_N_washoff_dung[hh+1] = N_washoff_rate*self.avail_N_pool_dung[hh+1]
                    self.resist_N_washoff_dung[hh+1] = N_washoff_rate*self.resist_N_pool_dung[hh+1]
                    self.unavail_N_washoff_dung[hh+1] = N_washoff_rate*self.unavail_N_pool_dung[hh+1]  
                    TANprod_dung = Na_decomp_rate*self.avail_N_pool_dung[hh+1] + Nr_decomp_rate*self.resist_N_pool_dung[hh+1]
                    ## update manure pool and org N pools
                    self.manure_pool_dung[hh+1] = self.manure_pool_dung[hh+1] - self.manure_washoff_dung[hh+1]
                    self.avail_N_pool_dung[hh+1] = self.avail_N_pool_dung[hh+1]*(1-Na_decomp_rate) - \
                                                        self.avail_N_washoff_dung[hh+1]
                    self.resist_N_pool_dung[hh+1] = self.resist_N_pool_dung[hh+1]*(1-Nr_decomp_rate) - \
                                                        self.resist_N_washoff_dung[hh+1]
                    self.unavail_N_pool_dung[hh+1] = self.unavail_N_pool_dung[hh+1] - self.unavail_N_washoff_dung[hh+1]
                    
                    ## TAN pool of soil
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh]+TANprod_dung
                    ## TAN and NH3 concentration
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[0,hh+1],cncH=sim_ccH)
                    self.TAN_amount[llidx,hh+1] = self.TAN_pool[llidx,hh+1]/(self.Total_water_pool_dung[hh+1]+self.manure_pool_dung[hh+1]*Kd_manure/(manure_PD/1e3))
                    self.TAN_amount[llidx,hh+1] = self.TAN_amount[llidx,hh+1]*1e6
                    self.NH3_gas[llidx,hh+1] = self.TAN_amount[llidx,hh+1]*KNH3
                    TANsurfamount = self.TAN_amount[llidx,hh+1]
                    NH3surfamount = TANsurfamount*KNH3
                    ## determining the potential of each flux
                    emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600
                    TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[llidx,hh+1],
                                                    qsubrunoff=self.drainagerate[llidx,hh+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[llidx,hh+1],cnc2=0,
                                        resist=self.Rdiffaq[llidx,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.theta[llidx+1,hh+1]==0]=0.0
                    NH3diffgasdownidx = False
                    TANuptakeidx = False  ## N uptake
                    nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[llidx,hh+1],
                            flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                            flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                    self.NH3flux_dung[hh+1] = nh3volidx
                    self.TANwashoff_dung[hh+1] = srfrunoffidx
                    self.TANinfil[llidx,hh+1] = subsrfleachingidx
                    self.TANdiffusiondown[llidx,hh+1] = diffaqdownidx
                    self.NH3diffusiondown[llidx,hh+1] = diffgasdownidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh+1]-self.NH3flux_dung[hh+1]-self.TANwashoff_dung[hh+1]-\
                                    self.TANinfil[llidx,hh+1]-self.TANdiffusiondown[llidx,hh+1]-self.NH3diffusiondown[llidx,hh+1]
                    ## get rid of rounding error
                    # self.TAN_pool[llidx,hh+1][self.TAN_pool[llidx,hh+1]<0] = 0.0
                    self.TAN_pool[llidx,hh+1] = np.maximum(self.TAN_pool[llidx,hh+1],0.0)
                    ## updating manure water pool
                    self.Total_water_pool_dung[hh+1] = self.Total_water_pool_dung[hh+1] - self.drainagerate[llidx,hh+1]*timestep*3600*1e6

                    ##########################################
                    ## "urine patch" scheme; level index = 2
                    ##########################################
                    ## urine only, no feces
                    llidx = 2
                    if dd_span < Days:
                        urinepatch_ccH = self.soil_ccH[dd_span]
                    else:
                        urinepatch_ccH = self.soil_ccH[int(dd_span-np.floor(dd_span/Days)*Days)]
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[0,hh+1],cncH=urinepatch_ccH)
                    self.theta[llidx,hh+1] = (self.theta[llidx,hh]*z_source + self.urine[hh+1]/1e6)/z_source
                    self.theta[llidx,hh+1][self.theta[llidx,hh+1]>self.soil_satmoist[0,hh+1]] = self.soil_satmoist[0,hh+1][self.theta[llidx,hh+1]>self.soil_satmoist[0,hh+1]]
                    self.drainagerate[llidx,hh+1] = np.maximum(subrunoffrate,water_drainage(theta=self.theta[llidx,hh+1],theta_sat=self.soil_satmoist[0,hh+1],
                                                        Ksat=Ks_sat,fc=theta_FC,layerthickness=z_source))
                    Rdiffsrfaq = diff_resistance(distance=z_source/2,phase='aqueous',
                                theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1],temp=self.soil_temp[0,hh+1])
                    Rdiffsrfgas = diff_resistance(distance=z_source/2,phase='gaseous',
                                theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1],temp=self.soil_temp[0,hh+1])
                    self.Rdiffaq[llidx,hh+1] = diff_resistance(distance=pmids[0]-z_source/2,phase='aqueous',
                                theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1],temp=self.soil_temp[0,hh+1])
                    self.Rdiffgas[llidx,hh+1] = diff_resistance(distance=pmids[0]-z_source/2,phase='gaseous',
                                theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1],temp=self.soil_temp[0,hh+1])
                    ## urea pool
                    # self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh] + (1-frac_FYM)*self.urea[hh+1]
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh] + self.urea[hh+1]
                    ## org N pools (no fecal N)
                    self.avail_N[hh+1] = (self.urine_N[hh+1]-self.urea[hh+1])*f_avail
                    self.resist_N[hh+1] = (self.urine_N[hh+1]-self.urea[hh+1])*f_resist
                    self.unavail_N[hh+1] = (self.urine_N[hh+1]-self.urea[hh+1])*f_unavail
                    self.avail_N_pool[hh+1] = self.avail_N_pool[hh] + self.avail_N[hh+1]
                    self.resist_N_pool[hh+1] = self.resist_N_pool[hh] + self.resist_N[hh+1]
                    self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh] + self.unavail_N[hh+1]
                    ## surface washoff of org N pools
                    self.avail_N_washoff[hh+1] = N_washoff_rate*self.avail_N_pool[hh+1]
                    self.resist_N_washoff[hh+1] = N_washoff_rate*self.resist_N_pool[hh+1]
                    self.unavail_N_washoff[hh+1] = N_washoff_rate*self.unavail_N_pool[hh+1]
                    ## urea hydrolysis
                    self.ureahydrolysis[llidx,hh+1] = self.urea_pool[llidx,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[0,hh+1],
                                                                                                theta=self.theta[llidx,hh+1],delta_t=timestep,k_h=0.03)
                    ## subtracting chemical losses
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh+1]-self.ureahydrolysis[llidx,hh+1]
                    ## urea concentration
                    self.urea_amount[llidx,hh+1] = N_concentration(mN=self.urea_pool[llidx,hh+1],zlayer=z_source,
                                                                theta=self.theta[llidx,hh+1])
                    ## urea concentration at the compensation point
                    ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[llidx,hh+1],rliq=Rdiffsrfaq,
                                                    qrunoff=surfrunoffrate)
                    ## determine the potential of each flux
                    ureawashoffidx = surf_runoff(N_surfcnc=ureasurfamount,
                                                    qrunoff=surfrunoffrate)*timestep*3600
                    ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[llidx,hh+1],
                                                    qsubrunoff=self.drainagerate[llidx,hh+1])*timestep*3600
                    ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[llidx,hh+1],cnc2=0,
                                                    resist=self.Rdiffaq[llidx,hh+1])*timestep*3600
                    ureadiffupidx = False
                    srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[llidx,hh+1],
                                flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                    self.ureawashoff[hh+1] = srfrunoffidx
                    self.ureainfil[llidx,hh+1] = subsrfleachingidx
                    self.ureadiffusiondown[llidx,hh+1] = diffaqdownidx  
                    ## update urea pool
                    self.urea_pool[llidx,hh+1] = self.urea_pool[llidx,hh+1]-self.ureawashoff[hh+1]-\
                                            self.ureainfil[llidx,hh+1]-self.ureadiffusiondown[llidx,hh+1]
                    ## get rid of rounding error
                    self.urea_pool[llidx,hh+1][self.urea_pool[llidx,hh+1]<0] = 0.0
                    TANprod_soil = self.ureahydrolysis[llidx,hh+1] + \
                                                Na_decomp_rate*self.avail_N_pool[hh+1] + \
                                                Nr_decomp_rate*self.resist_N_pool[hh+1] 
                    ## update org N pools
                    self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1]*(1-Na_decomp_rate)-\
                                                        self.avail_N_washoff[hh+1]
                    self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1]*(1-Nr_decomp_rate)-\
                                                self.resist_N_washoff[hh+1]
                    self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh+1] - self.unavail_N_washoff[hh+1]
                    ## TAN pool
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh] + TANprod_soil
                    fNH4 = frac_NH4(theta=self.theta[llidx,hh+1],theta_sat=self.soil_satmoist[0,hh+1],
                                            temp=self.soil_temp[0,hh+1],cncH=urinepatch_ccH,kd=Kd)
                    self.NH4nitrif[llidx,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[llidx,hh+1],temp=self.soil_temp[0,hh+1],
                                                theta=self.theta[llidx,hh+1],theta_sat=self.soil_satmoist[0,hh+1],
                                                pH=sim_pH,fert_type='manure',frac_nh4=fNH4)*timestep*3600
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh+1] - self.NH4nitrif[llidx,hh+1]
                    self.TAN_amount[llidx,hh+1] = TAN_concentration(mtan=self.TAN_pool[llidx,hh+1],zlayer=z_source,
                                                        theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1],
                                                        knh3=KNH3,kd=Kd)
                    self.NH3_gas[llidx,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[llidx,hh+1],knh3=KNH3,
                                                    theta_sat=self.soil_satmoist[0,hh+1],theta=self.theta[llidx,hh+1])
                    ## TAN concentration at the compensation surface
                    TANsurfamount = surf_TAN_cnc(tan_cnc=self.TAN_amount[llidx,hh+1],rliq=Rdiffsrfaq,rgas=Rdiffsrfgas,
                                        knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=self.drainagerate[llidx,hh+1])
                    NH3surfamount = TANsurfamount*KNH3
                    ## determining the potential of each flux
                    emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                    ## NH3 from soil layer is zero when slurry TAN pool is not empty
                    emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                    TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[llidx,hh+1],
                                                    qsubrunoff=self.drainagerate[llidx,hh+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[llidx,hh+1],cnc2=0,
                                        resist=self.Rdiffaq[llidx,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.theta[llidx,hh+1]==0]=0.0
                    NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[llidx,hh+1],cnc2=0,
                                        resist=self.Rdiffgas[llidx,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                    NH3diffgasdownidx[self.theta[llidx,hh+1]==self.soil_satmoist[0,hh+1]]=0.0
                    TANuptakeidx = False  ## N uptake
                    nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[llidx,hh+1],
                            flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                            flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                    self.NH3flux[hh+1] = nh3volidx
                    self.TANwashoff[hh+1] = srfrunoffidx
                    self.TANinfil[llidx,hh+1] = subsrfleachingidx
                    self.TANdiffusiondown[llidx,hh+1] = diffaqdownidx
                    self.NH3diffusiondown[llidx,hh+1] = diffgasdownidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[llidx,hh+1] = self.TAN_pool[llidx,hh+1]-self.NH3flux[hh+1]-self.TANwashoff[hh+1]-\
                                    self.TANinfil[llidx,hh+1]-self.TANdiffusiondown[llidx,hh+1]-self.NH3diffusiondown[llidx,hh+1]
                    ## get rid of rounding error
                    self.TAN_pool[llidx,hh+1] = np.maximum(self.TAN_pool[llidx,hh+1],0.0)
                    
                    ## update soil moisture 
                    update_sw = (self.theta[llidx,hh+1]*z_source - self.drainagerate[llidx,hh+1]*timestep*3600-evap)/z_source
                    self.theta[llidx,hh+1] = np.maximum(update_sw,self.soil_moist[0,hh+1])
                self.daily_output(dd_span,grazing=True)
                self.daily_init(grazing=True)
                ## manure and N pools re-initialise on the last day of each span loop
                if dd_span == dd+span_day-1:
                    
                    if dd_span < Days:
                        dd_spanidx = dd_span  
                    else:
                        dd_spanidx = int(dd_span-np.floor(dd_span/Days)*Days)
                    ##
                    self.o_soil_orgN_FYM[dd_spanidx] = self.avail_N_pool_FYM[-1]+self.resist_N_pool_FYM[-1]+self.unavail_N_pool_FYM[-1]
                    self.o_soil_TAN_FYM[dd_spanidx] = self.TAN_pool[0,-1]+self.urea_pool[0,-1]
                    self.o_soil_orgN_dung[dd_spanidx] = self.avail_N_pool_dung[-1]+self.resist_N_pool_dung[-1]+self.unavail_N_pool_dung[-1]
                    self.o_soil_TAN_dung[dd_spanidx] = self.TAN_pool[1,-1]
                    self.o_soil_orgN[dd_spanidx] = self.avail_N_pool[-1]+self.resist_N_pool[-1]+self.unavail_N_pool[-1]
                    self.o_soil_TAN[dd_spanidx] = self.TAN_pool[2,-1]+self.urea_pool[2,-1]
                    self.avail_N_pool_FYM[:] = 0.0
                    self.resist_N_pool_FYM[:] = 0.0
                    self.unavail_N_pool_FYM[:] = 0.0
                    self.avail_N_pool_dung[:] = 0.0
                    self.resist_N_pool_dung[:] = 0.0
                    self.unavail_N_pool_dung[:] = 0.0
                    self.avail_N_pool[:] = 0.0
                    self.resist_N_pool[:] = 0.0
                    self.unavail_N_pool[:] = 0.0
                    self.TAN_pool[:] = 0.0
                    self.urea_pool[:] = 0.0
                    ## concentrations
                    self.TAN_amount[:] = 0.0
                    self.NH3_gas[:] = 0.0
                    self.urea_amount[:] = 0.0
                    ## reinitialise soil pH
                    self.soil_pH[:] = 0.0
                    self.soil_ccH[:] = 0.0

        # print(test_annual_N*source_area[34,388]/1e9)   
        return

    def land_sim_reshape(self,sim_result):
        shape = sim_result.shape
        dim1 = int((shape[0]/2))
        output = np.zeros([dim1,shape[1],shape[2]])
        print(dim1,output.shape)
        output = sim_result[:dim1]+sim_result[dim1:]
        return output

    def para_out(self,chem_fert_type,output_stat=False,quality_check=False):

        if chem_fert_type == 'ammonium':
            sim_area = self.ammN_area
            chemfert_Ntotal = self.land_sim_reshape(self.TAN_added)*sim_area
        elif chem_fert_type == 'urea':
            sim_area = self.ureaN_area
            chemfert_Ntotal = self.land_sim_reshape(self.urea_added)*sim_area
        elif chem_fert_type == 'nitrate':
            sim_area = self.nitN_area
            chemfert_Ntotal = self.land_sim_reshape(self.NO3_added)*sim_area   

        self.o_NH3flux = self.o_NH3flux*sim_area
        self.o_washoff = self.o_washoff*sim_area
        self.o_nitrif = self.o_nitrif*sim_area
        self.o_NH4leaching = self.o_NH4leaching*sim_area
        self.o_diffaq = self.o_diffaq*sim_area
        self.o_diffgas = self.o_diffgas*sim_area
        self.o_ammNuptake = self.o_ammNuptake*sim_area
        self.o_nitNuptake = self.o_nitNuptake*sim_area
        self.o_NO3washoff = self.o_NO3washoff*sim_area
        self.o_NO3leaching = self.o_NO3leaching*sim_area
        self.o_NO3diff = self.o_NO3diff*sim_area

        if output_stat is True:
            print('Total N applied: '+str(sum_totalGg(chemfert_Ntotal))+' Gg')
            print('NH3 emission: '+str(sum_totalGg(self.o_NH3flux))+' Gg')
            print('TAN washoff: '+str(sum_totalGg(self.o_washoff))+' Gg')
            print('NH4 nitrification: '+str(sum_totalGg(self.o_nitrif))+' Gg')
            print('NH4 leaching: '+str(sum_totalGg(self.o_NH4leaching))+' Gg')
            print('TAN diffusion (aq) to deeper soil: '+ str(sum_totalGg(self.o_diffaq))+' Gg')
            print('NH3 diffusion (gas) to deeper soil: '+ str(sum_totalGg(self.o_diffgas))+' Gg')
            print('NH4 uptake by plants: '+ str(sum_totalGg(self.o_ammNuptake))+' Gg')
            print('NO3 uptake by plants: '+ str(sum_totalGg(self.o_nitNuptake))+' Gg')
            print('NO3 washoff: '+str(sum_totalGg(self.o_NO3washoff))+' Gg')
            print('NO3 leaching: '+ str(sum_totalGg(self.o_NO3leaching))+' Gg')
            print('NO3 diffusion to deeper soil: '+ str(sum_totalGg(self.o_NO3diff))+' Gg')

        TANleft = (self.TAN_pool[0,-1]+self.TAN_pool[1,-1]+self.TAN_pool[2,-1])*sim_area
        urealeft = (self.urea_pool[0,-1]+self.urea_pool[1,-1]+self.urea_pool[2,-1])*sim_area

        if quality_check is True:
            check = np.nansum(chemfert_Ntotal-self.o_NH3flux-self.o_washoff-self.o_nitrif-self.o_NH4leaching-\
                    self.o_diffaq-self.o_diffgas-self.o_ammNuptake,axis=0)
            result = np.round(sum_totalGg(check)/sum_totalGg(chemfert_Ntotal)*100,decimals=3)
            print("Integrity check result: "+str(result)+' %')
            # result2 = np.nansum(check,axis=0)/np.nansum(chemfert_Ntotal,axis=0)
            # print("max diff:",np.where(result2==np.nanmax(result2)),np.nanmax(result2))      
            # print("min diff:",np.where(result2==np.nanmin(result2)),np.nanmin(result2))     
    
    def main(self,fert_method,crop_item,chem_fert_type,start_day_idx,end_day_idx,
                sim_type='base',senstest_var=False,senstest=False,output_stat=False,quality_check=False):
        self.chem_fert_input(crop=crop_item)
        self.land_sim(start_day_idx,end_day_idx,chem_fert_type,tech=fert_method,crop=crop_item,
                        sim=sim_type,stvar=senstest_var,st=senstest)
        self.para_out(chem_fert_type,output_stat,quality_check)
        return

    def manure_sim_main(self,fert_method,manure_type,livestock_name,production_system,mms_cat,phase,start_day_idx,end_day_idx,
                        crop_item=None,sim_type='base',senstest_var=False,senstest=False,
                        output_stat=False,quality_check=False):
            # manure_fert_input(self,livestock_name,production_system,mms_cat,phase,crop=None):
        self.manure_fert_input(livestock_name,production_system,mms_cat,phase)
        self.land_manure_sim(start_day_idx,end_day_idx,manure_type,tech=fert_method,
                        sim=sim_type,stvar=senstest_var,st=senstest)
        # self.para_out(chem_fert_type,output_stat,quality_check)
        return
