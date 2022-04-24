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
    def __init__(self,prank,psize,fert_type,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,NO3_added,water_added,pH_value):
        
        print('LAND Module - current fertilizer application is: '+str(fert_type))

        ## lat idx for parallelization
        self.plat1 = prank*psize
        self.plat2 = (prank+1)*psize

        print(self.plat1,self.plat2)
        
        ## array shape [lats,lons]
        field_shape = (psize,lons)
        ## include the time dimension
        array_shape = (25,) + field_shape
        ## intermediate field shape
#         intarray_shape = mtrx2
        ## output shape
        outarray_shape = (365,) + field_shape

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
        self.soil_pH = np.zeros((365,) + field_shape)
        self.soil_ccH = np.zeros((365,) + field_shape)

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

        ## adsorption 
        self.Kd = np.zeros(array_shape[1:])

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
        
    def sim_env(self,dayidx):
        # print('LAND ENV: open env')
        #### environmental conditions
        ## soil temperature
        hhidx = dayidx*24
        self.soil_temp[0,1:] = groundtemp_filelvl1.stl1[hhidx:hhidx+24,self.plat1:self.plat2,:] - 273.15
        self.soil_temp[1,1:] = groundtemp_filelvl2.stl2[hhidx:hhidx+24,self.plat1:self.plat2,:] - 273.15
        # self.soil_temp[0,1:] = groundtemp_filelvl1.stl1[dayidx,self.plat1:self.plat2,:] - 273.15
        # self.soil_temp[1,1:] = groundtemp_filelvl2.stl2[dayidx,self.plat1:self.plat2,:] - 273.15

        self.soil_moist[0,1:] = soilmoist_filelvl1.swvl1[hhidx:hhidx+24,self.plat1:self.plat2,:]
        self.soil_moist[1,1:] = soilmoist_filelvl2.swvl2[hhidx:hhidx+24,self.plat1:self.plat2,:]
        # self.soil_moist[0,1:] = soilmoist_filelvl1.swvl1[dayidx,self.plat1:self.plat2,:]
        # self.soil_moist[1,1:] = soilmoist_filelvl2.swvl2[dayidx,self.plat1:self.plat2,:]
        self.soil_moist[self.soil_moist<0] = 0.0

        soilbd_ds = open_ds(infile_path+soil_data_path+soilbdfile)
        ## buld density unit: kg/dm3
        soilbd = soilbd_ds.T_BULK_DEN.values[self.plat1:self.plat2,:]
        soilporosity = 1 - (soilbd/(rho_soil/1000))
        self.soil_satmoist[:] = soilporosity
        self.soil_satmoist[self.soil_satmoist>1.0] = 0.99

        ## evaporation from bare soil
        # self.evap_sim[1:] = evap_file.evabs[hhidx:hhidx+24,self.plat1:self.plat2,:]*(-1e6)
        ## rainfall
        # self.rainfall[1:] = rain_file.tcrw[hhidx:hhidx+24,self.plat1:self.plat2,:]*1e3
        ## convert into m/s
        self.surfrunoffrate[1:] = runoff_file.sro[hhidx:hhidx+24,self.plat1:self.plat2,:]/(timestep*3600)
        # self.surfrunoffrate[1:] = runoff_file.sro[dayidx,self.plat1:self.plat2,:]/(24*3600)
        ## convert into m/s
        self.subrunoffrate[1:] = subrunoff_file.ssro[hhidx:hhidx+24,self.plat1:self.plat2,:]/(timestep*3600)
        # self.subrunoffrate[1:] = subrunoff_file.ssro[dayidx,self.plat1:self.plat2,:]/(24*3600)
        ## atmospheric resistances
        self.R_atm[1:] = ratm_file.RAM1[hhidx:hhidx+24,self.plat1:self.plat2,:]+ratm_file.RB1[hhidx:hhidx+24,self.plat1:self.plat2,:]
        # self.R_atm[1:] = ratm_file.RAM1[dayidx,self.plat1:self.plat2,:]+ratm_file.RB1[dayidx,self.plat1:self.plat2,:]
           
    
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
        self.evap_sim = field_var_fill(sd_template=template,
                                            input_field=self.evap_sim) ## g/day
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
            self.soil_pH = soil_pH_postapp(base_pH=soil_pH,app_timing_map=N_app_mark,fert_pH=8.5)[:,self.plat1:self.plat2,:] 
            self.soil_ccH = 10**(-self.soil_pH)
        return N_app

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
        plantdate = cropcalds['plant.start'].values
        harvestdate = cropcalds['harvest.start'].values
        plantfill = cropcalds['plant'].values
        harvestfill = cropcalds['harvest'].values
        plantdate[(np.isnan(plantdate))&(~np.isnan(plantfill))] = plantfill[(np.isnan(plantdate))&(~np.isnan(plantfill))]
        harvestdate[(np.isnan(harvestdate))&(~np.isnan(harvestfill))] = harvestfill[(np.isnan(harvestdate))&(~np.isnan(harvestfill))]
        ## harvesting date goes into next year
        harvestdate[harvestdate<plantdate] = harvestdate[harvestdate<plantdate]+365
        return plantdate,harvestdate

    def chem_fert_input(self,crop):
        ## read N application rates dataset for crops
        fertds = open_ds(infile_path+crop_data_path+crop+cropfileformat)
        ## crop calendar dataset
        cropcalspath = infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat
        # cropcalds = open_ds(infile_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat)
        ## base soil pH dataset
        soilpHds = open_ds(infile_path+soil_data_path+soilpHfile)

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
        self.pH = soilph[:,self.plat1:self.plat2,:] 
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
        self.NO3_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.nitN_area = fnitN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 
        self.TAN_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.ammN_area = fammN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 
        self.urea_added = chem_N_tocrop[:,self.plat1:self.plat2,:] 
        self.ureaN_area = fureaN[self.plat1:self.plat2,:]  * croparea[self.plat1:self.plat2,:] 

        ## met data interpolation
        self.met_input_interp(totalN[self.plat1:self.plat2,:])
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
    
    ## 
    def daily_init(self):
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
        return

    ##
    def daily_output(self,dayidx):
        if dayidx >= Days:
            dayidx = dayidx - Days
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
    
    ## main function
    def land_sim(self,start_day_idx,end_day_idx,chem_fert_type,tech,crop=None,sim='base',stvar=False,st=False):

        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is: '+str(tech))
        soilclayds = open_ds(infile_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)[self.plat1:self.plat2,:]

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
                if chem_fert_type == 'nitrate':
                    # print('chemical fertilizer applied: nitrate')
                    self.NO3[5] = self.NO3_added[dd]
                    sim_pH = self.pH[dd]
                    sim_ccH = self.cc_H[dd]
                elif chem_fert_type == 'ammonium':
                    self.TAN[5] = self.TAN_added[dd]
                    # print('chemical fertilizer applied: ammonium')
                    sim_pH = self.pH[dd]
                    sim_ccH = self.cc_H[dd]
                elif chem_fert_type == 'urea':
                    self.urea[5] = self.urea_added[dd]
                    # print('chemical fertilizer applied: urea')
                    sim_pH = self.soil_pH[dd]
                    sim_ccH = self.soil_ccH[dd]
            else:
                self.sim_env(dd-365)
                if chem_fert_type == 'nitrate':
                    # print('chemical fertilizer applied: nitrate')
                    self.NO3[5] = self.NO3_added[dd-365]
                    sim_pH = self.pH[dd-365]
                    sim_ccH = self.cc_H[dd-365]
                elif chem_fert_type == 'ammonium':
                    self.TAN[5] = self.TAN_added[dd-365]
                    # print('chemical fertilizer applied: ammonium')
                    sim_pH = self.pH[dd-365]
                    sim_ccH = self.cc_H[dd]
                elif chem_fert_type == 'urea':
                    self.urea[5] = self.urea_added[dd-365]
                    # print('chemical fertilizer applied: urea')
                    sim_pH = self.soil_pH[dd-365]
                    sim_ccH = self.soil_ccH[dd-365]

            if sim != 'base':
                self.sensitivity_test(var=stvar,test=st)

            for hh in np.arange(0,24):
                ## resistance for upward diffusion in the surface layer
                Rdiffsrfaq = diff_resistance(distance=pmids[0],phase='aqueous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.soil_moist[0,hh+1],temp=self.soil_temp[0,hh+1])
                Rdiffsrfgas = diff_resistance(distance=pmids[0],phase='gaseous',
                            theta_sat=self.soil_satmoist[0,hh+1],theta=self.soil_moist[0,hh+1],temp=self.soil_temp[0,hh+1])
                surfrunoffrate = self.surfrunoffrate[hh+1]-self.surfrunoffrate[hh]
                surfrunoffrate[surfrunoffrate<0] = 0.0
                subrunoffrate = self.subrunoffrate[hh+1]-self.subrunoffrate[hh]
                subrunoffrate[subrunoffrate<0] = 0.0
                for ll in np.arange(3):
                    ## lldix: index for soil temp, moisture
                    llidx = int(np.floor(ll/2))

                    ## resistance for diffusion between surface layer (idx0) and topsoil layer (idx1)
                    self.Rdiffaq[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='aqueous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1],temp=self.soil_temp[llidx,hh+1])
                    self.Rdiffgas[ll,hh+1] = diff_resistance(distance=pmids[ll+1]-pmids[ll],phase='gaseous',
                            theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1],temp=self.soil_temp[llidx,hh+1])

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
                            self.TAN_pool[ll,hh] = self.TAN_pool[ll,hh]+self.TAN[hh+1]
                            self.NO3_pool[ll,hh] = self.NO3_pool[ll,hh]+self.NO3[hh+1]

                    ## urea scheme
                    if chem_fert_type == 'urea':
                        if ll == 0:
                            ## urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusionup[ll,hh]
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        delta_t=timestep,k_h=0.03)
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],
                                                                        theta=self.soil_moist[llidx,hh+1])
                            ## urea concentration at the compensation point
                            ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,hh+1],rliq=Rdiffsrfaq,
                                                            qrunoff=surfrunoffrate)
                            # ureasurfamount = surf_Ncnc(N_cnc=self.urea_amount[ll,hh+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[hh+1])
                            ## determine the potential of each flux
                            ureawashoffidx = surf_runoff(N_surfcnc=ureasurfamount,
                                                        # qrunoff=(self.surfrunoffrate[hh+1]))*timestep*3600
                                                         qrunoff=surfrunoffrate)*timestep*3600
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                            qsubrunoff=subrunoffrate)*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=self.urea_amount[ll+1,hh],
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
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
                            self.urea_pool[ll,hh+1][self.urea_pool[ll,hh+1]<0] = 0.0
                        elif ll == 1:
                            ## urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusionup[ll,hh]+\
                                                        self.ureadiffusiondown[ll-1,hh+1]+self.ureainfil[ll-1,hh+1]
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        delta_t=timestep,k_h=0.03)
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,hh+1])
                            ## determine the potential of each flux
                            ureawashoffidx = False
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                            qsubrunoff=subrunoffrate)*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=self.urea_amount[ll+1,hh],
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
                            ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,hh],cnc2=self.urea_amount[ll-1,hh+1],
                                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600
                            srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,hh+1],
                                        flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)
                            self.ureainfil[ll,hh+1] = subsrfleachingidx
                            self.ureadiffusiondown[ll,hh+1] = diffaqdownidx  
                            self.ureadiffusionup[ll-1,hh+1] = diffaqupidx
                            ## update urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureadiffusionup[ll-1,hh+1]-\
                                                        self.ureainfil[ll,hh+1]-self.ureadiffusiondown[ll,hh+1]
                            ## get rid of rounding error
                            self.urea_pool[ll,hh+1][self.urea_pool[ll,hh+1]<0] = 0.0
                        elif ll == 2:
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh]+self.ureadiffusiondown[ll-1,hh+1]+self.ureainfil[ll-1,hh+1]
                            ## urea hydrolysis
                            self.ureahydrolysis[ll,hh+1] = self.urea_pool[ll,hh+1]*urea_hydrolysis_rate(temp=self.soil_temp[llidx,hh+1],
                                                                                                        delta_t=timestep,k_h=0.03)
                            ## subtracting chemical losses
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureahydrolysis[ll,hh+1]
                            ## urea concentration
                            self.urea_amount[ll,hh+1] = N_concentration(mN=self.urea_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,hh+1])
                            ## determine the potential of each flux
                            ureawashoffidx = False
                            ureainfilidx = subsurf_leaching(N_cnc=self.urea_amount[ll,hh+1],
                                                            # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                            qsubrunoff=subrunoffrate)*timestep*3600
                            ureadiffdownidx = N_diffusion(cnc1=self.urea_amount[ll,hh+1],cnc2=0.0,
                                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600
                            ureadiffupidx = N_diffusion(cnc1=self.urea_amount[ll,hh],cnc2=self.urea_amount[ll-1,hh+1],
                                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600
                            srfrunoffidx,subsrfleachingidx,diffaqdownidx,diffaqupidx = N_pathways(mN=self.urea_pool[ll,hh+1],
                                        flux1=ureawashoffidx,flux2=ureainfilidx,flux3=ureadiffdownidx,flux4=ureadiffupidx)       
                            self.ureainfil[ll,hh+1] = subsrfleachingidx
                            self.ureadiffusiondown[ll,hh+1] = diffaqdownidx              
                            self.ureadiffusionup[ll-1,hh+1] = diffaqupidx
                            ## update urea pool
                            self.urea_pool[ll,hh+1] = self.urea_pool[ll,hh+1]-self.ureadiffusionup[ll-1,hh+1]-\
                                                        self.ureainfil[ll,hh+1]-self.ureadiffusiondown[ll,hh+1]
                            ## get rid of rounding error
                            self.urea_pool[ll,hh+1][self.urea_pool[ll,hh+1]<0] = 0.0
                    ## TAN scheme
                    try:
                        TANprod = self.ureahydrolysis[ll,hh+1]+self.orgN_decomp[hh+1]
                    except:
                        Pass
                        TANprod = self.ureahydrolysis[ll,hh+1]

                    if ll == 0:
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+TANprod
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1])
                        ## TAN concentration at the compensation surface
                        TANsurfamount = surf_TAN_cnc(tan_cnc=self.TAN_amount[ll,hh+1],rliq=Rdiffsrfaq,rgas=Rdiffsrfgas,
                                            knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=surfrunoffrate)
                                            # knh3=KNH3,ratm=self.R_atm[hh+1],qrunoff=self.surfrunoffrate[hh+1])
                        NH3surfamount = TANsurfamount*KNH3
                        ## determining the potential of each flux
                        emissidx = NH3_vol(nh3_surfcnc=NH3surfamount,ratm=self.R_atm[hh+1])*timestep*3600  ## NH3 volatlization
                        TANwashoffidx = surf_runoff(N_surfcnc=TANsurfamount,
                                                    # qrunoff=(self.surfrunoffrate[hh+1]))*timestep*3600 
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## TAN washoff
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                        # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                        qsubrunoff=subrunoffrate)*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.soil_moist[llidx,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.soil_moist[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
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
                        ## get rid of rounding error
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0

                    elif ll == 1:
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+self.TANdiffusionup[ll,hh]+self.NH3diffusionup[ll,hh]+TANprod+\
                                                self.TANinfil[ll-1,hh+1]+self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                            qsubrunoff=subrunoffrate)*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=self.TAN_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.soil_moist[llidx,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=self.NH3_gas[ll+1,hh],
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.soil_moist[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        NH3diffgasupidx[self.soil_moist[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        TANdiffaqdownidx[self.soil_moist[llidx,hh+1]==0]=0.0
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
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0
                    elif ll == 2:
                        ## TAN pool
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh]+TANprod+self.TANinfil[ll-1,hh+1]+\
                                                    self.TANdiffusiondown[ll-1,hh+1]+self.NH3diffusiondown[ll-1,hh+1]
                        ## fraction of aqueous NH4
                        fNH4 = frac_NH4(theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                        temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH,kd=Kd)
                        self.NH4nitrif[ll,hh+1] = TAN_nitrif(tan_pool=self.TAN_pool[ll,hh+1],temp=self.soil_temp[llidx,hh+1],
                                                    theta=self.soil_moist[llidx,hh+1],theta_sat=self.soil_satmoist[llidx,hh+1],
                                                    pH=sim_pH,fert_type='mineral',frac_nh4=fNH4)*timestep*3600
                        ## subtracting chemical transformation
                        self.TAN_pool[ll,hh+1] = self.TAN_pool[ll,hh+1]-self.NH4nitrif[ll,hh+1]
                        ## TAN and NH3 concentration
                        KNH3 = NH3_par_coeff(temp=self.soil_temp[llidx,hh+1],cncH=sim_ccH)
                        self.TAN_amount[ll,hh+1] = TAN_concentration(mtan=self.TAN_pool[ll,hh+1],zlayer=zlayers[ll],
                                                    theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1],
                                                    knh3=KNH3,kd=Kd)
                        self.NH3_gas[ll,hh+1] = NH3_concentration(tan_cnc=self.TAN_amount[ll,hh+1],knh3=KNH3,
                                                theta_sat=self.soil_satmoist[llidx,hh+1],theta=self.soil_moist[llidx,hh+1])
                        ## determining the potential of each flux
                        TANinfilidx = subsurf_leaching(N_cnc=self.TAN_amount[ll,hh+1],
                                                            # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                            qsubrunoff=subrunoffrate)*timestep*3600  ## TAN infiltration/leaching
                        TANdiffaqdownidx = N_diffusion(cnc1=self.TAN_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1])*timestep*3600  ## TAN aqueous downwards diffusion
                        TANdiffaqdownidx[self.soil_moist[llidx,hh+1]==0]=0.0
                        NH3diffgasdownidx = N_diffusion(cnc1=self.NH3_gas[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffgas[ll,hh+1])*timestep*3600  ## NH3 gaseous downwards diffusion
                        NH3diffgasdownidx[self.soil_moist[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        NH3diffgasupidx = N_diffusion(cnc1=self.NH3_gas[ll,hh],cnc2=self.NH3_gas[ll-1,hh+1],
                                            resist=self.Rdiffgas[ll-1,hh+1])*timestep*3600  ## NH3 gaseous upwards diffusion
                        NH3diffgasupidx[self.soil_moist[llidx,hh+1]==self.soil_satmoist[llidx,hh+1]]=0.0
                        TANdiffaqupidx = N_diffusion(cnc1=self.TAN_amount[ll,hh],cnc2=self.TAN_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1])*timestep*3600  ## TAN aqueous diffusion
                        TANdiffaqupidx[self.soil_moist[llidx,hh+1]==0]=0.0
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
                        self.TAN_pool[ll,hh+1][self.TAN_pool[ll,hh+1]<0] = 0.0

                    ## NO3 scheme
                    if ll == 0:
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+self.NO3diffusionup[ll,hh]

                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,hh+1])
                        ## NO3 concentration at the compensation surface
                        NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,hh+1],rliq=Rdiffsrfaq,
                                                qrunoff=surfrunoffrate)
                        # NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[ll,hh+1],rliq=Rdiffsrfaq,qrunoff=self.surfrunoffrate[hh+1])
                        ## determining the potential of each flux
                        NO3washoffidx = surf_runoff(N_surfcnc=NO3surfamount,
                                                    # qrunoff=(self.surfrunoffrate[hh+1]))*timestep*3600
                                                    qrunoff=surfrunoffrate)*timestep*3600  ## NO3 washoff
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                        # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                       qsubrunoff=subrunoffrate)*timestep*3600  ## NO3 leaching
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
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                        # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600 
                                                       qsubrunoff=subrunoffrate)*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=self.NO3_amount[ll+1,hh],
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
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
                        self.NO3_pool[ll,hh+1][self.NO3_pool[ll,hh+1]<0.0] = 0.0
                    elif ll == 2:
                        ## NO3 pool
                        self.NO3_pool[ll,hh+1] = self.NO3_pool[ll,hh]+self.NH4nitrif[ll,hh+1]+\
                                                    self.NO3diffusiondown[ll-1,hh+1]+self.NO3infil[ll-1,hh+1]
                        ## NO3 concentration
                        self.NO3_amount[ll,hh+1] = N_concentration(mN=self.NO3_pool[ll,hh+1],zlayer=zlayers[ll],theta=self.soil_moist[llidx,hh+1])
                        ## determining the potential of each flux
                        NO3infilidx = subsurf_leaching(N_cnc=self.NO3_amount[ll,hh+1],
                                                        # qsubrunoff=(self.subrunoffrate[hh+1]))*timestep*3600
                                                       qsubrunoff=subrunoffrate)*timestep*3600  ## NO3 leaching
                        NO3diffaqdownidx = N_diffusion(cnc1=self.NO3_amount[ll,hh+1],cnc2=0,
                                            resist=self.Rdiffaq[ll,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
                        NO3diffupidx = N_diffusion(cnc1=self.NO3_amount[ll,hh],cnc2=self.NO3_amount[ll-1,hh+1],
                                            resist=self.Rdiffaq[ll-1,hh+1]/f_DNO3)*timestep*3600  ## NO3 aqueous diffusion
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
                        self.NO3_pool[ll,hh+1][self.NO3_pool[ll,hh+1]<0.0] = 0.0
            self.daily_output(dd)
            self.daily_init()            
        return

    def land_sim_reshape(self,sim_result):
        shape = sim_result.shape
        dim1 = int((shape[0]-1)/2+1)
        output = np.zeros([dim1,shape[1],shape[2]])
        output = sim_result[1:dim1]+sim_result[dim1:]
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

        if quality_check is True:
            check = chemfert_Ntotal-self.o_NH3flux-self.o_washoff-self.o_nitrif-self.o_NH4leaching-\
                    self.o_diffaq-self.o_diffgas-self.o_ammNuptake
            result = np.round(sum_totalGg(check)/sum_totalGg(chemfert_Ntotal)*100,decimals=3)
            print("Integrity check result: "+str(result)+' %')        
    
    def main(self,fert_method,crop_item,chem_fert_type,start_day_idx,end_day_idx,
                sim_type='base',senstest_var=False,senstest=False,output_stat=False,quality_check=False):
        self.chem_fert_input(crop=crop_item)
        self.land_sim(start_day_idx,end_day_idx,chem_fert_type,tech=fert_method,crop=crop_item,
                        sim=sim_type,stvar=senstest_var,st=senstest)
        self.para_out(chem_fert_type,output_stat,quality_check)
        return
