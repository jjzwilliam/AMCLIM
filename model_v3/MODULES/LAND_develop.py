from ast import Pass
from cmath import phase
from logging import raiseExceptions
from os import times
from re import sub

from pandas import array
from INPUT.input import *
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *
from POSTPROC.postprocessing import *
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
## dry matter (DM) content of solid manure 
DM_content = solid_m_DM[livestock]
## dry matter (DM) content of liquid manure is assumed to be 5%
f_DM_liquid = 0.1
## maximum water content of manure
f_wcmax = 1 - (DM_content/100)/2
## assuming the density of manure; 1t kg/m^3 or 1g/cm^3
manure_density = rho_m[livestock]
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
# dailymaxinfil = 0.0 ## test: shut down infiltration
## infiltration flux within manure (m/s)
qpsoil_manure = (dailymaxinfil/1e6)/(24*3600)
## washoff coefficients: 0.1%/mm water for N species, and 0.05%/mm water for non N species (manure)
f_washoff_nonN = 0.0005
f_washoff_N = 0.001

## thichkness of vertical layer 1, 2, 3, 4: 2cm(1), 5cm(2), 7cm(3), 14cm(4)
zlayers = [0.02, 0.05, 0.07, 0.14]
## midpoints of each layer
pmids = [0.01, 0.045, 0.105, 0.21]


class LAND_module:
    def __init__(self,array_shape,fert_type,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,NO3_added,water_added,pH_value):
        
        print('LAND Module - current fertilizer application is: '+str(fert_type))

        ## shape of 2 levels fields
        dlvl = [2] + array_shape
        ## shape of 3 levels fields
        tlvl = [3] + array_shape

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
            self.manure_water = np.zeros(array_shape)
            self.manure_minwc = np.zeros(array_shape)
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
            ## 
            self.orgN_decomp = np.zeros(array_shape)
            ## water added from MMS/HOUSING
            self.water_added = water_added
            self.water = np.zeros(array_shape)
            ## total water pool of the system (soil water; manure water+soil water)
            self.Total_water_pool = np.zeros(array_shape)
            ## total pool of the topsoil layer
            self.water_pool_soil = np.zeros(array_shape)
            ## cropland area
            self.cropland = np.zeros(array_shape[1:])

        else:
            ## cropping area for nitrate N fertilizer
            self.nitN_area = np.zeros(array_shape[1:])
            ## cropping area for ammonium N fertilizer
            self.ammN_area = np.zeros(array_shape[1:])
            ## cropping area for urea N fertilizer
            self.ureaN_area = np.zeros(array_shape[1:])


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

        #### met  and soil fields:
        ## atmospheric resistances: 1) aerodynamic resistance, 2) boundary layer resistance
        self.R_atm = np.zeros(array_shape)
        ## evaporation
        self.evap_sim = np.zeros(array_shape)
        ## rain fall (Note the unit)
        self.rainfall = np.zeros(array_shape)
        ## surface runoff rate (m/s)
        self.surfrunoffrate = np.zeros(array_shape)
        ## subsurface runoff rate (m/s)
        self.subrunoffrate = np.zeros(array_shape)
        ## soil pH
        self.soil_pH = np.zeros(array_shape)
        self.soil_ccH = np.zeros(array_shape)

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
        
    
    def met_input_interp(self,template):
        ##################################
        ## fill land input data
        ##################################
        self.soil_temp[0] = field_var_fill(sd_template=template,input_field=self.soil_temp[0])  ## degC
        self.soil_temp[1] = field_var_fill(sd_template=template,input_field=self.soil_temp[1])  ## degC
        # rhum_data = field_var_fill(sd_template=template,input_field=rhum_data)  ## per cent
        # wind_data = field_var_fill(sd_template=template,input_field=wind_data)  ## m/s
        self.evap_sim = field_var_fill(sd_template=template,input_field=self.evap_sim) ## g/day
        # soilmoist_data = field_var_fill(sd_template=template,input_field=soilmoist_data)  ## m3/m3
        self.soil_moist[0]= field_var_fill(sd_template=template,input_field=self.soil_moist[0])  ## m3/m3
        self.soil_moist[1] = field_var_fill(sd_template=template,input_field=self.soil_moist[1])  ## m3/m3
        self.soil_satmoist[0] = field_var_fill(sd_template=template,input_field=self.soil_satmoist[0])  ## m3/m3
        self.soil_satmoist[1] = field_var_fill(sd_template=template,input_field=self.soil_satmoist[1])  ## m3/m3
        self.R_atm = field_var_fill(sd_template=template,input_field=self.R_atm)  ## s/m
        self.surfrunoffrate = field_var_fill(sd_template=template,input_field=self.surfrunoffrate)  ## m/day
        self.subrunoffrate = field_var_fill(sd_template=template,input_field=self.subrunoffrate)  ## m/day
        self.soil_moist[self.soil_moist>self.soil_satmoist] = self.soil_satmoist[self.soil_moist>self.soil_satmoist]
        return

    def sim_env(self):
        print('LAND ENV: open env')
        #### environmental conditions
        ## soil temperature
        self.soil_temp[0,:366] = groundtemp_datalvl1
        self.soil_temp[0,366:] = groundtemp_datalvl1[1:]
        self.soil_temp[1,:366] = groundtemp_datalvl2
        self.soil_temp[1,366:] = groundtemp_datalvl2[1:]
        ## soil moisture and porosity
        # self.soil_moist[0,:366] = soilmoist_data
        # self.soil_moist[0,366:] = soilmoist_data[1:]
        # self.soil_moist[1,:366] = soilmoist_data
        # self.soil_moist[1,366:] = soilmoist_data[1:]
        self.soil_moist[0,:366] = soilmoist_datalvl1
        self.soil_moist[0,366:] = soilmoist_datalvl1[1:]
        self.soil_moist[1,:366] = soilmoist_datalvl2
        self.soil_moist[1,366:] = soilmoist_datalvl2[1:]
        self.soil_moist[self.soil_moist<0] = 0.0
        soilbd_ds = open_ds(file_path+soil_data_path+soilbdfile)
        ## buld density unit: kg/dm3
        soilbd = soilbd_ds.T_BULK_DEN.values
        soilporosity = 1 - (soilbd/(rho_soil/1000))
        # self.soil_satmoist[0,:366] = soilmoist_data/(persm_data/100)
        # self.soil_satmoist[0,366:] = soilmoist_data[1:]/(persm_data[1:]/100)
        # self.soil_satmoist[1,:366] = soilmoist_data/(persm_data/100)
        # self.soil_satmoist[1,366:] = soilmoist_data[1:]/(persm_data[1:]/100)
        self.soil_satmoist[:] = soilporosity
        self.soil_satmoist[self.soil_satmoist>1.0] = 0.99
        ## evaporation from bare soil
        self.evap_sim[:366] = evap_data
        self.evap_sim[366:] = evap_data[1:]
        ## rainfall
        self.rainfall[:366] = rain_data
        self.rainfall[366:] = rain_data[1:]
        ## convert into m/s
        self.surfrunoffrate[:366] = runoff_data/(timestep*3600)
        self.surfrunoffrate[366:] = runoff_data[1:]/(timestep*3600)
        ## convert into m/s
        self.subrunoffrate[:366] = subrunoff_data/(timestep*3600)
        self.subrunoffrate[366:] = subrunoff_data[1:]/(timestep*3600)
        self.R_atm[:366] = ram1_data+rb1_data
        self.R_atm[366:] = ram1_data[1:]+rb1_data[1:]
        return

    def spreading_time(self,fert_type='manure',
                        crop_type=None,N_input=None,plant_calendar=None,harvest_calendar=None,
                        fert_freq=None,soil_pH=None):
        ## lats and lons
        lats = mtrx2[1]
        lons = mtrx2[2]
        N_app = np.zeros(mtrx2)
        ## mark all application day
        N_app_mark = np.zeros(mtrx2)

        if fert_type == "manure":
            for lat in np.arange(lats):
                for lon in np.arange(lons):
                    ## here
                    app_time = None
        else:
            print(crop_type)
            for lat in np.arange(lats):
                for lon in np.arange(lons):
                    ## timing of planting and harvesting
                    plt_time = plant_calendar[lat,lon]
                    har_time = harvest_calendar[lat,lon]
                    ## how many times that N is applied on fields in an annual cycle
                    app_freq = fert_freq[lat,lon]
                    ## harvesting goes into the next year
                    if ~np.isnan(app_freq):
                        if ~np.isnan(plt_time):
                            if ~np.isnan(har_time):
                                if har_time<plt_time:
                                    har_time = har_time+365
                                if app_freq <=1.0:
                                    N_app[int(plt_time),lat,lon] = N_input[lat,lon]*app_freq
                                    ## mark up the N application
                                    N_app_mark[int(plt_time),lat,lon] = 1 
                                elif app_freq>1.0:
                                    ## get the integer times of N application
                                    tapp = np.floor(app_freq)
                                    ## application intervals
                                    app_int = int(abs(int(har_time)-int(plt_time))/(tapp+1))
                                    ## index for application
                                    app_idx = np.arange(int(plt_time),int(har_time)-1,app_int)
                                    ## in tapp times, the application rate equals the readed applcation rates
                                    for idx in app_idx[:-1]:
                                        N_app[int(idx),lat,lon] = N_input[lat,lon]
                                        N_app_mark[int(idx),lat,lon] = 1
                                    ## the residual of N application rates
                                    N_app[int(app_idx[-1]),lat,lon] = (app_freq-tapp)*N_input[lat,lon]
                                    N_app_mark[int(app_idx[-1]),lat,lon] = 1
        
        ## determine soil pH after urea application
        if soil_pH is not None:
            # print('check pH over.')
            self.soil_pH = soil_pH_postapp(base_pH=soil_pH,app_timing_map=N_app_mark,fert_pH=8.5)
            self.soil_ccH = 10**(-self.soil_pH)
        return N_app

    def chem_fert_type(self):
        ## country level chemical fertilizer use: 1) ammonium N, 2) nitrate N, 3) urea N
        fertds = open_ds(file_path+crop_data_path+fertfilename)
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
        plantdate = cropcalds['plant.start'].values
        harvestdate = cropcalds['harvest.start'].values
        ## harvesting date goes into next year
        harvestdate[harvestdate<plantdate] = harvestdate[harvestdate<plantdate]+365
        return plantdate,harvestdate

    def chem_fert_input(self,crop):
        ## read N application rates dataset for crops
        fertds = open_ds(file_path+crop_data_path+crop+cropfileformat)
        ## crop calendar dataset
        cropcalspath = file_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
        # cropcalds = open_ds(file_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat)
        ## base soil pH dataset
        soilpHds = open_ds(file_path+soil_data_path+soilpHfile)

        totalN = fertds.TotalN.values*1e3
        ## N application rate is interpolated;
        rateN = fertds.Nrate.values*1e3/1e4
        croparea = fertds.croparea.values*1e4
        croparea[totalN!=0] = totalN[totalN!=0]/rateN[totalN!=0]
        app_freq = totalN/(rateN*croparea)
        ## the maximum application frequency in a year is 5
        app_freq[app_freq>5] = 5

        ## read planting and harvesting dates
        plantidx, harvestidx = self.crop_calendar(filepath = cropcalspath)
        # plantidx = cropcalds['plant.start'].values
        # harvestidx = cropcalds['harvest.start'].values

        soilph = np.zeros(mtrx2)
        soilph[:] = soilpHds.T_PH_H2O.values
        self.pH = soilph
        self.cc_H = 10**(-self.pH)

        chem_N_tocrop = self.spreading_time(fert_type='mineral',
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
        self.NO3_added = chem_N_tocrop
        self.nitN_area = fnitN * croparea
        self.TAN_added = chem_N_tocrop
        self.ammN_area = fammN * croparea
        self.urea_added = chem_N_tocrop
        self.ureaN_area = fureaN * croparea

        ## met data interpolation
        self.met_input_interp(totalN)
        return 

    ## determine fertilizer application depth
    ## two basic application techniques:
    ## 1) surface spreading (broadcasting): depth ~ 4mm (0.004m)
    ## 2) deep injection: dpeth ~ 10cm (0.1m)
    ## this will affect the N transport (diffusion to the surface) then affect NH3 emission
    def app_method(self,application_method):
        if application_method is None:
            self.app_depth = 0.004
        else:
            self.app_depth = 0.1

        return

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

    ## main function
    def land_sim(self,start_day_idx,end_day_idx,chem_fert_type,tech,crop=None):

        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is: '+str(tech))
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)

        if chem_fert_type == 'nitrate':
            print('chemical fertilizer applied: nitrate')
            self.NO3 = self.NO3_added
            sim_pH = self.pH
            sim_ccH = self.cc_H
        elif chem_fert_type == 'ammonium':
            self.TAN = self.TAN_added
            print('chemical fertilizer applied: ammonium')
            sim_pH = self.pH
            sim_ccH = self.cc_H
        elif chem_fert_type == 'urea':
            self.urea = self.urea_added
            print('chemical fertilizer applied: urea')
            sim_pH = self.soil_pH
            sim_ccH = self.soil_ccH

        ## crop calendar that detenmines the N uptake by crops
        if crop is not None:
            cropcalspath = file_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
            plantidx, harvestidx = self.crop_calendar(filepath = cropcalspath)

        sourcelayer = self.source_layer(tech)

        for dd in np.arange(start_day_idx,end_day_idx):
            ## resistance for upward diffusion in the surface layer
            Rdiffsrfaq = diff_resistance(distance=pmids[0],phase='aqueous',
                        theta_sat=self.soil_satmoist[0,dd+1],theta=self.soil_moist[0,dd+1],temp=self.soil_temp[0,dd+1])
            Rdiffsrfgas = diff_resistance(distance=pmids[0],phase='gaseous',
                        theta_sat=self.soil_satmoist[0,dd+1],theta=self.soil_moist[0,dd+1],temp=self.soil_temp[0,dd+1])
            for ll in np.arange(3):
                ## lldix: index for soil temp, moisture
                llidx = int(np.floor(ll/2))
                
                ## resistance for diffusion between surface layer (idx0) and topsoil layer (idx1)
                self.Rdiffaq[ll,dd+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                        theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1],temp=self.soil_temp[llidx,dd+1])
                self.Rdiffgas[ll,dd+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                        theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1],temp=self.soil_temp[llidx,dd+1])
                
                ## input includes N from fertilizer application
                if ll == sourcelayer:
                    if sourcelayer == 1:
                        self.urea_pool[0,dd+1] = self.urea_pool[0,dd+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.urea[dd+1]
                        self.urea_pool[1,dd] = self.urea_pool[1,dd] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.urea[dd+1]
                        self.TAN_pool[0,dd+1] = self.TAN_pool[0,dd+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.TAN[dd+1]
                        self.TAN_pool[1,dd] = self.TAN_pool[1,dd] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.TAN[dd+1]
                        self.NO3_pool[0,dd+1] = self.NO3_pool[0,dd+1] + (zlayers[0]/(zlayers[0]+zlayers[1]))*self.NO3[dd+1]
                        self.NO3_pool[1,dd] = self.NO3_pool[1,dd] + (zlayers[1]/(zlayers[0]+zlayers[1]))*self.NO3[dd+1]
                    else:
                        self.urea_pool[ll,dd] = self.urea_pool[ll,dd]+self.urea[dd+1]
                        self.TAN_pool[ll,dd] = self.TAN_pool[ll,dd]+self.TAN[dd+1]
                        self.NO3_pool[ll,dd] = self.NO3_pool[ll,dd]+self.NO3[dd+1]

                ## urea scheme
                if chem_fert_type == 'urea':
                    if ll == 0:
                        ## urea pool
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd]+self.ureadiffusionup[ll,dd]
                        ## urea hydrolysis
                        self.ureahydrolysis[ll,dd+1] = self.urea_pool[ll,dd+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,dd+1],
                                                                                                    delta_t=timestep,k_h=0.03)
                        ## subtracting chemical losses
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureahydrolysis[ll,dd+1]
                        ## urea concentration
                        self.urea_amount[ll,dd+1] = N_concentration(mN=self.urea_pool[ll,dd+1],zlayer=zlayers[ll],
                                                                    theta=self.soil_moist[llidx,dd+1])
                        ## urea concentration at the compensation point
                        ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                        ## determine the potential of each flux
                        ureawashoffidx = surf_runoff(N_surfcnc=ureasurfamount,qrunoff=self.surfrunoffrate[dd+1])*timestep*3600
                        ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600
                        ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,dd+1],cnc2=self.urea_amount[ll+1,dd],
                                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600
                        ureadiffupidx = False
                        srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,dd+1],
                                    flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                        self.ureawashoff[dd+1] = srfrunoffidx
                        self.ureainfil[ll,dd+1] = subsrfleachingidx
                        self.ureadiffusiondown[ll,dd+1] = diffaqdownidx  
                        ## update urea pool
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureawashoff[dd+1]-\
                                                self.ureainfil[ll,dd+1]-self.ureadiffusiondown[ll,dd+1]
                        ## get rid of rounding error
                        self.urea_pool[ll,dd+1][self.urea_pool[ll,dd+1]<0] = 0.0
                    elif ll == 1:
                        ## urea pool
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd]+self.ureadiffusionup[ll,dd]+\
                                                    self.ureadiffusiondown[ll-1,dd+1]+self.ureainfil[ll-1,dd+1]
                        ## urea hydrolysis
                        self.ureahydrolysis[ll,dd+1] = self.urea_pool[ll,dd+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,dd+1],
                                                                                                    delta_t=timestep,k_h=0.03)
                        ## subtracting chemical losses
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureahydrolysis[ll,dd+1]
                        ## urea concentration
                        self.urea_amount[ll,dd+1] = N_concentration(mN=self.urea_pool[ll,dd+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,dd+1])
                        ## urea concentration at the compensation point
                        ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                        ## determine the potential of each flux
                        ureawashoffidx = False
                        ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600
                        ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,dd+1],cnc2=self.urea_amount[ll+1,dd],
                                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600
                        ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,dd],cnc2=self.urea_amount[ll-1,dd+1],
                                                        resist=self.Rdiffaq[ll-1,dd+1])*timestep*3600
                        srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,dd+1],
                                    flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                        self.ureainfil[ll,dd+1] = subsrfleachingidx
                        self.ureadiffusiondown[ll,dd+1] = diffaqdownidx  
                        self.ureadiffusionup[ll-1,dd+1] = diffaqupidx
                        ## update urea pool
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureadiffusionup[ll-1,dd+1]-\
                                                    self.ureainfil[ll,dd+1]-self.ureadiffusiondown[ll,dd+1]
                        ## get rid of rounding error
                        self.urea_pool[ll,dd+1][self.urea_pool[ll,dd+1]<0] = 0.0
                    elif ll == 2:
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd]+self.ureadiffusiondown[ll-1,dd+1]+self.ureainfil[ll-1,dd+1]
                        ## urea hydrolysis
                        self.ureahydrolysis[ll,dd+1] = self.urea_pool[ll,dd+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,dd+1],
                                                                                                    delta_t=timestep,k_h=0.03)
                        ## subtracting chemical losses
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureahydrolysis[ll,dd+1]
                        ## urea concentration
                        self.urea_amount[ll,dd+1] = N_concentration(mN=self.urea_pool[ll,dd+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,dd+1])
                        ## urea concentration at the compensation point
                        ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                        ## determine the potential of each flux
                        ureawashoffidx = False
                        ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600
                        ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,dd+1],cnc2=0.0,
                                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600
                        ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,dd],cnc2=self.urea_amount[ll-1,dd+1],
                                                        resist=self.Rdiffaq[ll-1,dd+1])*timestep*3600
                        srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,dd+1],
                                    flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)       
                        self.ureainfil[ll,dd+1] = subsrfleachingidx
                        self.ureadiffusiondown[ll,dd+1] = diffaqdownidx              
                        self.ureadiffusionup[ll-1,dd+1] = diffaqupidx
                        ## update urea pool
                        self.urea_pool[ll,dd+1] = self.urea_pool[ll,dd+1]-self.ureadiffusionup[ll-1,dd+1]-\
                                                    self.ureainfil[ll,dd+1]-self.ureadiffusiondown[ll,dd+1]
                        ## get rid of rounding error
                        self.urea_pool[ll,dd+1][self.urea_pool[ll,dd+1]<0] = 0.0
                ## TAN scheme
                try:
                    TANprod = self.ureahydrolysis[ll,dd+1]+self.orgN_decomp[dd+1]
                except:
                    Pass
                    TANprod = self.ureahydrolysis[ll,dd+1]
                
                if ll == 0:
                    ## TAN pool
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd]+self.TANdiffusionup[ll,dd]+self.NH3diffusionup[ll,dd]+TANprod
                    ## fraction of aqueous NH4
                    fNH4 = frac_NH4(theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                    temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1],kd=Kd)
                    self.NH4nitrif[ll,dd+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,dd+1],temp=self.soil_temp[llidx,dd+1],
                                                theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                                pH=sim_pH[dd+1],fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                    ## subtracting chemical transformation
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.NH4nitrif[ll,dd+1]
                    ## TAN and NH3 concentration
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1])
                    self.TAN_amount[ll,dd+1] = TAN_concentration(mtan=self.TAN_pool[ll,dd+1],zlayer=zlayers[ll],
                                                theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1],
                                                knh3=KNH3,kd=Kd)
                    self.NH3_gas[ll,dd+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,dd+1],knh3=KNH3,
                                            theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1])
                    ## TAN concentration at the compensation surface
                    TANsurfamount = surf_TAN_cnc(tan_cnc=self.TAN_amount[ll,dd+1],rliq=Rdiffsrfaq,rgas=Rdiffsrfgas,
                                        knh3=KNH3,ratm=self.R_atm[dd+1],qrunoff=self.surfrunoffrate[dd+1])
                    NH3surfamount = TANsurfamount*KNH3
                    ## determining the potential of each flux
                    emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[dd+1])*timestep*3600  ## NH3 volatlization
                    TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,qrunoff=self.surfrunoffrate[dd+1])*timestep*3600  ## TAN washoff
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,dd+1],
                                                        qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,dd+1],cnc2=self.TAN_amount[ll+1,dd],
                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.soil_moist[llidx,dd+1]==0]=0.0
                    NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,dd+1],cnc2=self.NH3_gas[ll+1,dd],
                                        resist=self.Rdiffgas[ll,dd+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                    NH3diffgasdownidx[self.soil_moist[llidx,dd+1]==self.soil_satmoist[llidx,dd+1]]=0.0
                    TANuptakeidx = False  ## N uptake
                    nh3volidx,srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffgasdownidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,dd+1],
                            flux1=emissidx,flux2=TANwashoffidx,flux3=TANinfilidx,flux4=TANdiffaqdownidx,
                            flux5=NH3diffgasdownidx,flux6=TANuptakeidx)
                    self.NH3flux[dd+1] = nh3volidx
                    self.TANwashoff[dd+1] = srfrunoffidx
                    self.TANinfil[ll,dd+1] = subsrfleachingidx
                    self.TANdiffusiondown[ll,dd+1] = diffaqdownidx
                    self.NH3diffusiondown[ll,dd+1] = diffgasdownidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.NH3flux[dd+1]-self.TANwashoff[dd+1]-\
                                    self.TANinfil[ll,dd+1]-self.TANdiffusiondown[ll,dd+1]-self.NH3diffusiondown[ll,dd+1]
                    ## get rid of rounding error
                    self.TAN_pool[ll,dd+1][self.TAN_pool[ll,dd+1]<0] = 0.0

                elif ll == 1:
                    ## TAN pool
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd]+self.TANdiffusionup[ll,dd]+self.NH3diffusionup[ll,dd]+TANprod+\
                                            self.TANinfil[ll-1,dd+1]+self.TANdiffusiondown[ll-1,dd+1]+self.NH3diffusiondown[ll-1,dd+1]
                    ## fraction of aqueous NH4
                    fNH4 = frac_NH4(theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                    temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1],kd=Kd)
                    self.NH4nitrif[ll,dd+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,dd+1],temp=self.soil_temp[llidx,dd+1],
                                                theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                                pH=sim_pH[dd+1],fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                    ## subtracting chemical transformation
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.NH4nitrif[ll,dd+1]
                    ## TAN and NH3 concentration
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1])
                    self.TAN_amount[ll,dd+1] = TAN_concentration(mtan=self.TAN_pool[ll,dd+1],zlayer=zlayers[ll],
                                                theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1],
                                                knh3=KNH3,kd=Kd)
                    self.NH3_gas[ll,dd+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,dd+1],knh3=KNH3,
                                            theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1])
                    ## determining the potential of each flux
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,dd+1],
                                                        qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,dd+1],cnc2=self.TAN_amount[ll+1,dd],
                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.soil_moist[llidx,dd+1]==0]=0.0
                    NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,dd+1],cnc2=self.NH3_gas[ll+1,dd],
                                        resist=self.Rdiffgas[ll,dd+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                    NH3diffgasdownidx[self.soil_moist[llidx,dd+1]==self.soil_satmoist[llidx,dd+1]]=0.0
                    NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,dd],cnc2=self.NH3_gas[ll-1,dd+1],
                                        resist=self.Rdiffgas[ll-1,dd+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                    NH3diffgasupidx[self.soil_moist[llidx,dd+1]==self.soil_satmoist[llidx,dd+1]]=0.0
                    TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,dd],cnc2=self.TAN_amount[ll-1,dd+1],
                                        resist=self.Rdiffaq[ll-1,dd+1])*timestep*3600  ## TAN aqueous diffusion
                    TANdiffaqdownidx[self.soil_moist[llidx,dd+1]==0]=0.0
                    TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,dd+1],mNO3=self.NO3_pool[ll,dd],
                                    temp=self.soil_temp[llidx,dd+1],uptake='nh4')*timestep*3600  ## N uptake
                    TANuptakeidx[harvestidx<(dd+1)] = 0.0
                    subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,dd+1],
                            flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                            flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                    self.TANinfil[ll,dd+1] = subsrfleachingidx
                    self.TANdiffusiondown[ll,dd+1] = diffaqdownidx
                    self.NH3diffusiondown[ll,dd+1] = diffgasdownidx
                    self.TANdiffusionup[ll-1,dd+1] = diffaqupidx
                    self.NH3diffusionup[ll-1,dd+1] = diffgasupidx
                    self.ammNuptake[ll-1,dd+1] = uptakeidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.TANdiffusionup[ll-1,dd+1]-self.NH3diffusionup[ll-1,dd+1]-\
                        self.TANinfil[ll,dd+1]-self.TANdiffusiondown[ll,dd+1]-self.NH3diffusiondown[ll,dd+1]-self.ammNuptake[ll-1,dd+1]
                    ## get rid of rounding error
                    self.TAN_pool[ll,dd+1][self.TAN_pool[ll,dd+1]<0] = 0.0
                elif ll == 2:
                    ## TAN pool
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd]+TANprod+self.TANinfil[ll-1,dd+1]+\
                                                self.TANdiffusiondown[ll-1,dd+1]+self.NH3diffusiondown[ll-1,dd+1]
                    ## fraction of aqueous NH4
                    fNH4 = frac_NH4(theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                    temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1],kd=Kd)
                    self.NH4nitrif[ll,dd+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,dd+1],temp=self.soil_temp[llidx,dd+1],
                                                theta=self.soil_moist[llidx,dd+1],theta_sat=self.soil_satmoist[llidx,dd+1],
                                                pH=sim_pH[dd+1],fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                    ## subtracting chemical transformation
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.NH4nitrif[ll,dd+1]
                    ## TAN and NH3 concentration
                    KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,dd+1],cncH=sim_ccH[dd+1])
                    self.TAN_amount[ll,dd+1] = TAN_concentration(mtan=self.TAN_pool[ll,dd+1],zlayer=zlayers[ll],
                                                theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1],
                                                knh3=KNH3,kd=Kd)
                    self.NH3_gas[ll,dd+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,dd+1],knh3=KNH3,
                                            theta_sat=self.soil_satmoist[llidx,dd+1],theta=self.soil_moist[llidx,dd+1])
                    ## determining the potential of each flux
                    TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,dd+1],
                                                        qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## TAN infiltration/leaching
                    TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,dd+1],cnc2=0,
                                        resist=self.Rdiffaq[ll,dd+1])*timestep*3600  ## TAN aqueous downwards diffusion
                    TANdiffaqdownidx[self.soil_moist[llidx,dd+1]==0]=0.0
                    NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,dd+1],cnc2=0,
                                        resist=self.Rdiffgas[ll,dd+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                    NH3diffgasdownidx[self.soil_moist[llidx,dd+1]==self.soil_satmoist[llidx,dd+1]]=0.0
                    NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,dd],cnc2=self.NH3_gas[ll-1,dd+1],
                                        resist=self.Rdiffgas[ll-1,dd+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                    NH3diffgasupidx[self.soil_moist[llidx,dd+1]==self.soil_satmoist[llidx,dd+1]]=0.0
                    TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,dd],cnc2=self.TAN_amount[ll-1,dd+1],
                                        resist=self.Rdiffaq[ll-1,dd+1])*timestep*3600  ## TAN aqueous diffusion
                    TANdiffaqupidx[self.soil_moist[llidx,dd+1]==0]=0.0
                    TANuptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,dd+1],mNO3=self.NO3_pool[ll,dd],
                                    temp=self.soil_temp[llidx,dd+1],uptake='nh4')*timestep*3600  ## N uptake
                    TANuptakeidx[harvestidx<(dd+1)] = 0.0
                    subsrfleachingidx,diffaqdownidx,diffgasdownidx,diffaqupidx,diffgasupidx,uptakeidx = TAN_pathways(mN=self.TAN_pool[ll,dd+1],
                            flux1=TANinfilidx,flux2=TANdiffaqdownidx,flux3=NH3diffgasdownidx,flux4=TANdiffaqupidx,
                            flux5=NH3diffgasupidx,flux6=TANuptakeidx)
                    self.TANinfil[ll,dd+1] = subsrfleachingidx
                    self.TANdiffusiondown[ll,dd+1] = diffaqdownidx
                    self.NH3diffusiondown[ll,dd+1] = diffgasdownidx
                    self.TANdiffusionup[ll-1,dd+1] = diffaqupidx
                    self.NH3diffusionup[ll-1,dd+1] = diffgasupidx
                    self.ammNuptake[ll-1,dd+1] = uptakeidx
                    ## update TAN pool: subtracting all losses
                    self.TAN_pool[ll,dd+1] = self.TAN_pool[ll,dd+1]-self.TANdiffusionup[ll-1,dd+1]-self.NH3diffusionup[ll-1,dd+1]-\
                        self.TANinfil[ll,dd+1]-self.TANdiffusiondown[ll,dd+1]-self.NH3diffusiondown[ll,dd+1]-self.ammNuptake[ll-1,dd+1]
                    ## get rid of rounding error
                    self.TAN_pool[ll,dd+1][self.TAN_pool[ll,dd+1]<0] = 0.0

                ## NO3 scheme
                if ll == 0:
                    ## NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd]+self.NH4nitrif[ll,dd+1]+self.NO3diffusionup[ll,dd]
                                                
                    ## NO3 concentration
                    self.NO3_amount[ll,dd+1] = N_concentration(mN=self.NO3_pool[ll,dd+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,dd+1])
                    ## NO3 concentration at the compensation surface
                    NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                    ## determining the potential of each flux
                    NO3washoffidx = surf_runoff(N_surfcnc=NO3surfamount, qrunoff=self.surfrunoffrate[dd+1])*timestep*3600  ## NO3 washoff
                    NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## NO3 leaching
                    NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,dd+1],cnc2=self.NO3_amount[ll+1,dd],
                                        resist=self.Rdiffaq[ll,dd+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                    NO3diffupidx = False
                    srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.NO3_pool[ll,dd+1],
                                    flux1=NO3washoffidx,flux2=NO3infilidx,flux3=NO3diffaqdownidx,flux4=NO3diffupidx)
                    self.NO3washoff[dd+1] = srfrunoffidx
                    self.NO3infil[ll,dd+1] = subsrfleachingidx
                    self.NO3diffusiondown[ll,dd+1] = diffaqdownidx
                    ## update NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd+1]-self.NO3washoff[dd+1]-\
                                                self.NO3infil[ll,dd+1]-self.NO3diffusiondown[ll,dd+1]
                    ## get rid of rounding error
                    self.NO3_pool[ll,dd+1][self.NO3_pool[ll,dd+1]<0.0] = 0.0
                elif ll == 1:
                    ## NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd]+self.NH4nitrif[ll,dd+1]+self.NO3diffusionup[ll,dd]+\
                                                self.NO3diffusiondown[ll-1,dd+1]+self.NO3infil[ll-1,dd+1]
                    ## NO3 concentration
                    self.NO3_amount[ll,dd+1] = N_concentration(mN=self.NO3_pool[ll,dd+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,dd+1])
                    ## NO3 concentration at the compensation surface
                    NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                    ## determining the potential of each flux
                    NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## NO3 leaching
                    NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,dd+1],cnc2=self.NO3_amount[ll+1,dd],
                                        resist=self.Rdiffaq[ll,dd+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                    NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,dd],cnc2=self.NO3_amount[ll-1,dd+1],
                                        resist=self.Rdiffaq[ll-1,dd+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                    NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,dd+1],mNO3=self.NO3_pool[ll,dd],
                                    temp=self.soil_temp[llidx,dd+1],uptake='no3')*timestep*3600  ## N uptake
                    NO3uptakeidx[harvestidx<(dd+1)] = 0.0
                    subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,dd+1],
                                    flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                    self.NO3infil[ll,dd+1] = subsrfleachingidx
                    self.NO3diffusiondown[ll,dd+1] = diffaqdownidx
                    self.NO3diffusionup[ll-1,dd+1] = diffaqupidx
                    self.nitNuptake[ll-1,dd+1] = uptakeidx
                    ## update NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd+1]-self.NO3infil[ll,dd+1]-self.NO3diffusionup[ll-1,dd+1]-\
                                                self.NO3diffusiondown[ll,dd+1]-self.nitNuptake[ll-1,dd+1]
                    ## get rid of rounding error
                    self.NO3_pool[ll,dd+1][self.NO3_pool[ll,dd+1]<0.0] = 0.0
                elif ll == 2:
                    ## NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd]+self.NH4nitrif[ll,dd+1]+\
                                                self.NO3diffusiondown[ll-1,dd+1]+self.NO3infil[ll-1,dd+1]
                    ## NO3 concentration
                    self.NO3_amount[ll,dd+1] = N_concentration(mN=self.NO3_pool[ll,dd+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,dd+1])
                    ## NO3 concentration at the compensation surface
                    NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,dd+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[dd+1])
                    ## determining the potential of each flux
                    NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,dd+1],qsubrunoff=self.subrunoffrate[dd+1])*timestep*3600  ## NO3 leaching
                    NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,dd+1],cnc2=0,
                                        resist=self.Rdiffaq[ll,dd+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                    NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,dd],cnc2=self.NO3_amount[ll-1,dd+1],
                                        resist=self.Rdiffaq[ll-1,dd+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                    NO3uptakeidx = plant_N_uptake(mNH4=self.TAN_pool[ll,dd+1],mNO3=self.NO3_pool[ll,dd],
                                    temp=self.soil_temp[llidx,dd+1],uptake='no3')*timestep*3600  ## N uptake
                    NO3uptakeidx[harvestidx<(dd+1)] = 0.0
                    subsrfleachingidx,diffaqdownidx,diffaqupidx,uptakeidx = N_pathways(mN=self.NO3_pool[ll,dd+1],
                                    flux1=NO3infilidx,flux2=NO3diffaqdownidx,flux3=NO3diffupidx,flux4=NO3uptakeidx)
                    self.NO3infil[ll,dd+1] = subsrfleachingidx
                    self.NO3diffusiondown[ll,dd+1] = diffaqdownidx
                    self.NO3diffusionup[ll-1,dd+1] = diffaqupidx
                    self.nitNuptake[ll-1,dd+1] = uptakeidx
                    ## update NO3 pool
                    self.NO3_pool[ll,dd+1] = self.NO3_pool[ll,dd+1]-self.NO3infil[ll,dd]-self.NO3diffusionup[ll-1,dd+1]-\
                                                self.NO3diffusiondown[ll,dd+1]-self.nitNuptake[ll-1,dd+1]
                    ## get rid of rounding error
                    self.NO3_pool[ll,dd+1][self.NO3_pool[ll,dd+1]<0.0] = 0.0
        return

    def land_sim_reshape(self,sim_result):
        shape = sim_result.shape
        dim1 = int((shape[0]-1)/2+1)
        output = np.zeros([dim1,shape[1],shape[2]])
        output = sim_result[1:dim1]+sim_result[dim1:]
        return output
    
    def N_stat(self,crop_item,chem_fert_type,fert_method,ncfile_o=False,quality_check=False):
        
        if chem_fert_type == 'ammonium':
            sim_area = self.ammN_area
            chemfert_Ntotal = self.land_sim_reshape(self.TAN_added)*sim_area

        elif chem_fert_type == 'urea':
            sim_area = self.ureaN_area
            chemfert_Ntotal = self.land_sim_reshape(self.urea_added)*sim_area
            
        elif chem_fert_type == 'nitrate':
            sim_area = self.nitN_area
            chemfert_Ntotal = self.land_sim_reshape(self.NO3_added)*sim_area     
        print('Total N applied: '+str(sum_totalGg(chemfert_Ntotal))+' Gg')
        
        chemfert_NH3emiss = self.land_sim_reshape(self.NH3flux)*sim_area
        print('NH3 emission: '+str(sum_totalGg(chemfert_NH3emiss))+' Gg')
        chemfert_TANwashoff = self.land_sim_reshape(self.TANwashoff+self.ureawashoff)*sim_area
        print('TAN washoff: '+str(sum_totalGg(chemfert_TANwashoff))+' Gg')
        chemfert_nitrif = self.land_sim_reshape(np.nansum(self.NH4nitrif,axis=0))*sim_area
        print('NH4 nitrification: '+str(sum_totalGg(chemfert_nitrif))+' Gg')
        chemfert_leaching = self.land_sim_reshape(self.TANinfil[-1]+self.ureainfil[-1])*sim_area
        print('NH4 leaching: '+str(sum_totalGg(chemfert_leaching))+' Gg')
        chemfert_diffaq = self.land_sim_reshape(self.TANdiffusiondown[-1]+self.ureadiffusiondown[-1])*sim_area
        print('TAN diffusion (aq) to deeper soil: '+ str(sum_totalGg(chemfert_diffaq))+' Gg')
        chemfert_diffgas = self.land_sim_reshape(self.NH3diffusiondown[-1])*sim_area
        print('NH3 diffusion (gas) to deeper soil: '+ str(sum_totalGg(chemfert_diffgas))+' Gg')
        chemfert_ammN_uptake = self.land_sim_reshape(np.nansum(self.ammNuptake,axis=0))*sim_area
        print('NH4 uptake by plants: '+ str(sum_totalGg(chemfert_ammN_uptake))+' Gg')
        chemfert_nitN_uptake = self.land_sim_reshape(np.nansum(self.nitNuptake,axis=0))*sim_area
        print('NO3 uptake by plants: '+ str(sum_totalGg(chemfert_nitN_uptake))+' Gg')
        chemfert_NO3washoff = self.land_sim_reshape(self.NO3washoff)*sim_area
        print('NO3 washoff: '+str(sum_totalGg(chemfert_NO3washoff))+' Gg')
        chemfert_NO3leaching = self.land_sim_reshape(self.NO3infil[-1])*sim_area
        print('NO3 leaching: '+ str(sum_totalGg(chemfert_NO3leaching))+' Gg')
        chemfert_NO3diffusion = self.land_sim_reshape(self.NO3diffusiondown[-1])*sim_area
        print('NO3 diffusion to deeper soil: '+ str(sum_totalGg(chemfert_NO3diffusion))+' Gg')

        if quality_check is True:
            check = chemfert_Ntotal-chemfert_NH3emiss-chemfert_TANwashoff-chemfert_nitrif-chemfert_leaching-\
                    chemfert_diffaq-chemfert_diffgas-chemfert_ammN_uptake
            result = np.round(sum_totalGg(check)/sum_totalGg(chemfert_Ntotal)*100,decimals=3)
            print("Integrity check result: "+str(result)+' %')

        ## generate an ncfile that contains the N pathways
        if ncfile_o is True:
            ## define output dims
            nlat = int(180.0/dlat)
            nlon = int(360.0/dlon)
            ntime = Days
            lats = 90 - 0.5*np.arange(nlat)
            lons = -180 + 0.5*np.arange(nlon)
            yearidx = str(sim_year)+'-01-01'
            times = pd.date_range(yearidx,periods=ntime)

            outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss=(['time','lat','lon'],chemfert_NH3emiss),
                    TANwashoff=(['lat','lon'],annual_total(chemfert_TANwashoff)),
                    TANleaching=(['lat','lon'],annual_total(chemfert_leaching)),
                    TANdiffaq=(['lat','lon'],annual_total(chemfert_diffaq)),
                    NH3diffgas=(['lat','lon'],annual_total(chemfert_diffgas)),
                    NH4nitrif=(['lat','lon'],annual_total(chemfert_nitrif)),
                    NH4uptake=(['lat','lon'],annual_total(chemfert_ammN_uptake)),
                    NO3uptake=(['lat','lon'],annual_total(chemfert_nitN_uptake)),
                    NO3washoff=(['lat','lon'],annual_total(chemfert_NO3washoff)),
                    NO3leaching=(['lat','lon'],annual_total(chemfert_NO3leaching)),
                    NO3diff=(['lat','lon'],annual_total(chemfert_NO3diffusion)),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-Land_chem_fert: \
                        N pathways of chemical fertilizer application in " +str(sim_year),
                    info = fert_method+" of "+chem_fert_type+" for: "+str(crop_item),
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

            outds.to_netcdf(output_path+str(crop_item)+'.'+\
                str(chem_fert_type)+'.'+str(fert_method)+'.'+str(sim_year)+'.nc',encoding=encoding)
            print("ncfile saved.")
        return
    
    def quality_check(self):
        ptest_stat = {'urea pool': 0,
             'urea conc': 0,
             'urea washoff': 0,
             'urea diffusion': 0,
             'urea hydrolysis': 0,
             'TAN pool': 0,
             'TAN conc': 0,
             'NH3 flux': 0,
             'TAN washoff': 0,
             'TAN leaching': 0,
             'TAN diffusion': 0,
             'NH3 diffusion': 0,
             'NH4 nitrif': 0,
             'NH4 uptake': 0,
             'NO3 pool': 0,
             'NO3 conc': 0,
             'NO3 washoff': 0,
             'NO3 leaching': 0,
             'NO3 diffusion': 0,
             'NO3 uptake': 0}

        if np.where(self.urea_pool<0)[0].size!=0:
            ptest_stat['urea pool'] = 1
        if np.where(self.urea_amount<0)[0].size!=0:
            ptest_stat['urea conc'] = 1
        if np.where(self.ureawashoff<0)[0].size!=0:
            ptest_stat['urea washoff'] = 1
        if np.where(self.ureadiffusiondown<0)[0].size!=0:
            ptest_stat['urea diffusion'] = 1
        if np.where(self.ureadiffusionup<0)[0].size!=0:
            ptest_stat['urea diffusion'] = 1
        if np.where(self.ureahydrolysis<0)[0].size!=0:
            ptest_stat['urea hydrolysis'] = 1
            
        if np.where(self.TAN_pool<0)[0].size!=0:
            ptest_stat['TAN pool'] = 1
        if np.where(self.TAN_amount<0)[0].size!=0:
            ptest_stat['TAN conc'] = 1
        if np.where(self.NH3flux<0)[0].size!=0:
            ptest_stat['NH3 flux'] = 1
        if np.where(self.TANwashoff<0)[0].size!=0:
            ptest_stat['TAN washoff'] = 1
        if np.where(self.TANinfil<0)[0].size!=0:
            ptest_stat['TAN leaching'] = 1
        if np.where(self.TANdiffusiondown<0)[0].size!=0:
            ptest_stat['TAN diffusion'] = 1
        if np.where(self.TANdiffusionup<0)[0].size!=0:
            ptest_stat['TAN diffusion'] = 1
        if np.where(self.NH3diffusiondown<0)[0].size!=0:
            ptest_stat['NH3 diffusion'] = 1
        if np.where(self.NH3diffusionup<0)[0].size!=0:
            ptest_stat['NH3 diffusion'] = 1
        if np.where(self.NH4nitrif<0)[0].size!=0:
            ptest_stat['NH4 nitrif'] = 1
        if np.where(self.ammNuptake<0)[0].size!=0:
            ptest_stat['NH4 uptake'] = 1

        if np.where(self.NO3_pool<0)[0].size!=0:
            ptest_stat['NO3 pool'] = 1
        if np.where(self.NO3_amount<0)[0].size!=0:
            ptest_stat['NO3 conc'] = 1
        if np.where(self.NO3washoff<0)[0].size!=0:
            ptest_stat['NO3 washoff'] = 1
        if np.where(self.NO3infil<0)[0].size!=0:
            ptest_stat['NO3 leaching'] = 1
        if np.where(self.NO3diffusiondown<0)[0].size!=0:
            ptest_stat['NO3 diffusion'] = 1
        if np.where(self.NO3diffusionup<0)[0].size!=0:
            ptest_stat['NO3 diffusion'] = 1
        if np.where(self.nitNuptake<0)[0].size!=0:
            ptest_stat['NO3 uptake'] = 1

        if all(value == 0 for value in ptest_stat.values()) is True:
            print("Basic test passsed.")
        else:
            print("Basic test failed! Please check!")
            print(ptest_stat)
    
    def main(self,fert_method,crop_item,chem_fert_type,start_day_idx,end_day_idx,sim_stat=False,ncfile_o=False,quality_check=False):
        self.sim_env()
        self.chem_fert_input(crop=crop_item)
        self.land_sim(start_day_idx,end_day_idx,chem_fert_type,tech=fert_method,crop=crop_item)
        if quality_check is True:
            self.quality_check()
        if sim_stat is True:
            ## output ncfile
            if ncfile_o is True:
                self.N_stat(crop_item,chem_fert_type,fert_method,ncfile_o=True,quality_check=quality_check)
            else:
                self.N_stat(crop_item,chem_fert_type,fert_method,ncfile_o=False,quality_check=quality_check)

        return
