from logging import raiseExceptions
from os import times

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

class LAND_module:
    def __init__(self,array_shape,fert_type,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,NO3_added,water_added,pH_value,cropping_area,application_method_index):
        
        print('LAND Module - current fertilizer application is: '+str(fert_type))
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
            ## water added from MMS/HOUSING
            self.water_added = water_added
            self.water = np.zeros(array_shape)
            ## total water pool of the system (soil water; manure water+soil water)
            self.Total_water_pool = np.zeros(array_shape)
            ## total pool of the topsoil layer
            self.water_pool_soil = np.zeros(array_shape)

        else:
            ## cropping area for nitrate N fertilizer
            self.nitN_area = np.zeros(array_shape[1:])
            ## cropping area for ammonium N fertilizer
            self.ammN_area = np.zeros(array_shape[1:])
            ## cropping area for urea N fertilizer
            self.ureaN_area = np.zeros(array_shape[1:])

        ## urea input
        self.urea_added = urea_added
        self.urea = np.zeros(array_shape)
        ## urea pool
        self.urea_pool = np.zeros(array_shape)
        ## urea pool of the topsoil layer
        self.urea_pool_soil = np.zeros(array_shape)
        ## urea concentration (layer where major N input is; source layer)
        self.urea_amount = np.zeros(array_shape)
        ## urea concentration of the topsoil layer
        self.urea_soil_amount = np.zeros(array_shape)
        ## washoff flux of urea
        self.urea_washoff = np.zeros(array_shape)
        ## urea diffusion to deeper soil (aq)
        self.urea_diff = np.zeros(array_shape)
        ## urea infiltration from the 1st layer to the underlying layer
        self.urea_infil = np.zeros(array_shape)
        ## urea leaching to deeper soil
        self.urealeaching = np.zeros(array_shape)
        ## uric acid input from MMS/HOUSING
        self.UA_added = UA_added
        self.UA = np.zeros(array_shape)
        ## uric acid pool
        self.UA_pool = np.zeros(array_shape)
        ## washoff flux of UA
        self.UA_washoff = np.zeros(array_shape)
        ## TAN added from MMS/HOUSING
        self.TAN_added = TAN_added
        self.TAN = np.zeros(array_shape)
        ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool (source layer)
        self.TAN_prod = np.zeros(array_shape)
        ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool (topsoil layer)
        self.TAN_prod_soil = np.zeros(array_shape)
        ## TAN pool (layer where major N input is; source layer)
        self.TAN_pool = np.zeros(array_shape)
        ## TAN pool of the topsoil layer
        self.TAN_pool_soil = np.zeros(array_shape)
        ## surface washoff flux of TAN 
        self.TAN_washoff = np.zeros(array_shape)
        ## TAN pool conc (aqueous phase)
        self.TAN_amount = np.zeros(array_shape)
        ## TAN pool in molar concentration
        self.TAN_amount_M = np.zeros(array_shape)
        ## TAN conc at the surface; 
        self.TAN_surf_amount_M = np.zeros(array_shape)
        ## TAN conc at the topsoil layer
        self.TAN_soil_amount = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3 = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3_added = NO3_added
        ## NO3- from nitrification in the source layer
        self.nitrif_NO3_sourcelayer = np.zeros(array_shape)
        ## NO3- washoff at the surface
        self.NO3_washoff = np.zeros(array_shape)
        ## NO3- pool in the source layer
        self.NO3_pool = np.zeros(array_shape)
        ## NO3- conc in the source layer
        self.NO3_amount = np.zeros(array_shape)
        ## NO3- from nitrification in the topsoil layer
        self.nitrif_NO3_soil = np.zeros(array_shape)
        ## NO3- pool in the topsoil layer
        self.NO3_pool_soil = np.zeros(array_shape)
        ## NO3- conc in the topsoil layer
        self.NO3_soil_amount = np.zeros(array_shape)
        ## infiltration of NO3 from the source layer to the topsoil layer
        self.NO3_infilsourcelayer = np.zeros(array_shape)
        ## diffusive aqueous NO3- from the source layer to the topsoil layer
        self.NO3_diffusivesourcelayer = np.zeros(array_shape)
        ## NO3- leaching to deeper soil
        self.NO3_leaching = np.zeros(array_shape)
        ## diffusive aqueous NO3- to deeper soil
        self.NO3_diffusivedeep = np.zeros(array_shape)
        ## ratio of [NH4+]/[H+] in the system
        self.Gamma_manure = np.zeros(array_shape)
        ## surface NH3 concentrtion at equilirium (in molar mass)
        self.NH3_gas_M = np.zeros(array_shape)
        ## NH3 concentration in the source layer
        self.NH3_gas_bulk = np.zeros(array_shape)
        ## topsoil layer NH3 concentration
        self.NH3_gas_soil = np.zeros(array_shape)

        ## emission potential
        self.modelled_emiss = np.zeros(array_shape)
        ## final emission
        ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3_flux]
        ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to
        ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
        ## during land spreading stages, such as canopy recapture, deap injection...
        self.NH3_flux = np.zeros(array_shape)
        ## leaching of aqueous TAN to deeper soil (not diffusive)
        self.leachingflux = np.zeros(array_shape)
        ## diffusion of aqueous TAN from the topsoil layer to the deeper soil
        self.diffusivefluxsoil_aq = np.zeros(array_shape)
        ## diffusion of gaseous NH3 from the topsoil layer to the deeper soil
        self.diffusivefluxsoil_gas = np.zeros(array_shape)
        ## infiltration of aqueous TAN from the 1st layer to the underlying layer; very similar to the leaching flux
        self.infilflux = np.zeros(array_shape)
        ## Note that when there is a thin source layer above (within) the topsoil layer, 
        ## there is no N uptake taken place in the source layer,
        ## because N uptake occurrs in the soil layer where the seedling zone is
        ## uptake of ammonium nitrogen by plants in the topsoil layer
        self.ammN_uptake = np.zeros(array_shape)
        ## uptake of nitrate nitrogen by plants in the topsoil layer
        self.nitN_uptake = np.zeros(array_shape)

        ## [broadcasting-disl] technique:
        ## the following fluxes are valid when source layer refers to the whole topsoil layer (broadcasting-disk)   
        if application_method_index == 'surf':
            ## diffusion of aqueous TAN from the source layer to the topsoil layer
            self.diffusivefluxsourcelayer_aq = np.zeros(array_shape)
            ## diffusion of gaseous NH3 from the source layer to the topsoil layer
            self.diffusivefluxsourcelayer_gas = np.zeros(array_shape)
            ## diffusion of urea from the source layer to the topsoil layer
            self.ureadiffusivefluxsourcelayer = np.zeros(array_shape)
            ## upward diffusion of aqueous TAN: from the topsoil layer to the source layer 
            self.diffusivefluxup_aq = np.zeros(array_shape)
            ## upward diffusion of gaseous NH3: from the topsoil layer to the  source layer
            self.diffusivefluxup_gas = np.zeros(array_shape)
            ## upward diffusion of urea: from the topsoil layer to the source layer
            # self.ureadiffusivefluxup = np.zeros(array_shape)

        ## [deep injection] technique:
        ## vertical soil profile in the model (from top to the bottom):
        ## surface (0 cm) - topsoil layer (0 - 7 cm) - source layer (7 - 13 cm) - deeper soil
        elif application_method_index == 'deep injection':
            ## diffusion of aqueous TAN from the source layer to the topsoil layer
            self.diffusivefluxsourcelayer_aq = np.zeros(array_shape)
            ## diffusion of gaseous NH3 from the source layer to the topsoil layer
            self.diffusivefluxsourcelayer_gas = np.zeros(array_shape)
            ## diffusion of urea from the source layer to the topsoil layer
            self.ureadiffusivefluxsourcelayer = np.zeros(array_shape)
            ## downward diffusion of aqueous TAN: from the topsoil layer to the source layer 
            self.diffusivefluxdown_aq = np.zeros(array_shape)
            ## downward diffusion of gaseous NH3: from the topsoil layer to the  source layer
            self.diffusivefluxdown_gas = np.zeros(array_shape)
            ## downward diffusion of urea: from topsoil layer to the source layer
            # self.ureadiffusivefluxdown = np.zeros(array_shape)
            ## infiltration of NO3 from the topsoillayer to the source layer; downwards
            self.NO3_infilsoil = np.zeros(array_shape)
            ## diffusive aqueous NO3- from the source layer to the topsoil layer; upwards
            self.NO3_diffusiveup = np.zeros(array_shape)
            ## ammonium N uptake by plants in the source layer
            self.ammN_uptake_sourcelayer = np.zeros(array_shape)
            ## nitrate N uptake by plants in the source layer
            self.nitN_uptake_sourcelayer = np.zeros(array_shape)

        ## temperature
        self.T_sim = np.zeros(array_shape)
        ## wind/ventilation
        self.u_sim = np.zeros(array_shape)
        ## RH
        self.RH_sim = np.zeros(array_shape)
        ## evaporation
        self.evap_sim = np.zeros(array_shape)
        ## volumetric soil moisture
        self.soilmoist = np.zeros(array_shape)
        ## percentage saturation soil moisture
        self.persm = np.zeros(array_shape)
        ## rain fall (Note the unit)
        self.rainfall = np.zeros(array_shape)
        ## rain available for "washoff"
        self.rain_avail_washoff = np.zeros(array_shape)
        ## atmospheric resistances: 1) aerodynamic resistance, 2) boundary layer resistance
        self.R_atm = np.zeros(array_shape)
        ## soil resistance of the source layer; distances may vary
        self.R_sourcelayer_aq = np.zeros(array_shape)
        ## soil resistance of the source layer; distances may vary
        self.R_sourcelayer_gas = np.zeros(array_shape)
        ## soil resistance for aqueous diffusion; distances may vary
        self.R_soilaq = np.zeros(array_shape)
        ## soil resistance for gaseous diffusion; distances may vary
        self.R_soilg = np.zeros(array_shape)
        ## soil resistance for downward aqueous diffusion; distances may vary
        self.R_soilaq_down = np.zeros(array_shape)
        ## soil resistance for downward gaseous diffusion; distances may vary
        self.R_soilg_down = np.zeros(array_shape)
        ## mositure equilirium, mositure content of manure
        self.mois_coeff = np.zeros(array_shape)
        ## daily UA hydrolysis rate
        self.daily_ua_conv_factor = np.zeros(array_shape)
        ## daily urea hydrolysis rate
        self.daily_urea_hydro_rate = np.zeros(array_shape)
        ## daily decomposition rate of available and resistant N components
        self.daily_Na_decomp_rate = np.zeros(array_shape)
        self.daily_Nr_decomp_rate = np.zeros(array_shape)
        ## daily nitrification rate
        self.daily_KNO3 = np.zeros(array_shape)
        ## Henry's law constant and dissociation equilibria; Eq.6
        self.Henry_constant = np.zeros(array_shape)
        ## dissociation constant of NH4+; 298.15 K is 25 degC (room temperature)
        self.k_NH4 = np.zeros(array_shape)
        ## molecular diffusivity of NH4+ in the water
        self.D_aq_NH4 = np.zeros(array_shape)
        ## molecular diffusivity of NH3 in the air
        self.D_air_NH3 = np.zeros(array_shape)
        ## percolation flux; m/s
        self.qpsoil = np.zeros(array_shape)
        ## pH and H+ ions concentration
        self.pH = pH_value
        self.cc_H = np.float(10**(-pH_value))
        ## this refers to soil pH after urea application;
        ## urea hydrolysis causes soil pH increases (consumes H+ ions) 
        self.soil_pH = np.zeros(array_shape)
        self.soil_ccH = np.zeros(array_shape)
        ## housing area that is used to determine MMS area
        self.cropland = cropping_area
        ## fertilzier application depth
        self.app_depth = np.zeros(array_shape[1:])
        ## MMS fractions used for determining fertilizer application method
        self.fMMS = application_method_index
        

    def sim_env(self):
        ## environmental conditions
        self.T_sim[:366] = groundtemp_data
        self.T_sim[366:] = groundtemp_data[1:]
        self.u_sim[:366] = wind_data
        self.u_sim[366:] = wind_data[1:]
        self.RH_sim[:366] = rhum_data
        self.RH_sim[366:] = rhum_data[1:]
        self.evap_sim[:366] = evap_data
        self.evap_sim[366:] = evap_data[1:]
        self.soilmoist[:366] = soilmoist_data
        self.soilmoist[366:] = soilmoist_data[1:]
        self.persm[:366] = soilmoist_data/(persm_data/100)
        self.persm[366:] = soilmoist_data[1:]/(persm_data[1:]/100)
        self.rainfall[:366] = rain_data
        self.rainfall[366:] = rain_data[1:]
        ## convert into m/s
        self.rain_avail_washoff[:366] = runoff_data/(timestep*3600)
        self.rain_avail_washoff[366:] = runoff_data[1:]/(timestep*3600)
        ## convert into m/s
        self.qpsoil[:366] = subrunoff_data/(timestep*3600)
        self.qpsoil[366:] = subrunoff_data[1:]/(timestep*3600)
        self.R_atm[:366] = ram1_data+rb1_data
        self.R_atm[366:] = ram1_data[1:]+rb1_data[1:]
        ## numpy arrays (not xarrays)
        self.T_sim = xr_to_np(self.T_sim)
        self.RH_sim = xr_to_np(self.RH_sim)
        self.u_sim = xr_to_np(self.u_sim)
        self.evap_sim = xr_to_np(self.evap_sim)
        self.soilmoist = xr_to_np(self.soilmoist)
        self.persm = xr_to_np(self.persm)
        self.persm[self.persm>1.0] = 0.9999
        self.rainfall = xr_to_np(self.rainfall)
        self.R_atm = xr_to_np(self.R_atm)
        print('LAND ENV: open env')

        ## mositure equilirium, mositure content of manure
        self.mois_coeff = (-np.log(1.00001-(self.RH_sim/100))/(0.0000534*(self.T_sim+273.15)))**(1/1.41)
        if livestock.lower()=="poultry":
            self.daily_ua_conv_factor = ua_hydrolysis_rate(temp=self.T_sim,rhum=self.RH_sim,ph=self.pH)
        else:
            ## daily urea hydrolysis rate
            self.daily_urea_hydro_rate = urea_hydrolysis_rate(temp=self.T_sim,delta_t=timestep)
            ## daily decomposition rate of available and resistant N components
            self.daily_Na_decomp_rate, self.daily_Nr_decomp_rate = N_pools_decomp_rate(temp=self.T_sim, delta_t=timestep)
        ## Henry's law constant and dissociation equilibria; Eq.6
        self.Henry_constant = (161500/(self.T_sim + 273.15)) * np.exp(-10380/(self.T_sim + 273.15))
        ## dissociation constant of NH4+; 298.15 K is 25 degC (room temperature)
        self.k_NH4 = 5.67e-10*np.exp(-6286*(1/(self.T_sim + 273.15)-1/298.15))
        ## molecular diffusivity of NH4+ in water
        self.D_aq_NH4 = diffusivity_NH4(temp=self.T_sim,phase='aqueous')
        ## molecular diffusivity of NH3 in the aire
        self.D_air_NH3 = diffusivity_NH4(temp=self.T_sim,phase='gaseous')
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
        nitN = fertds.nitrate_fert.values
        ammN = fertds.ammonium_fert.values
        ureaN = fertds.urea_fert.values
        totalN = nitN + ammN
        fnitN = nitN/totalN
        ## urea N is a subcategory of ammonium N in the readed dataset
        fammN = (ammN-ureaN)/totalN
        fureaN = ureaN/totalN
        return fnitN, fammN, fureaN
    
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
        print(soilph.shape)
        self.pH = soilph
        self.cc_H = 10**(-self.pH)

        chem_N_tocrop = self.spreading_time(fert_type='mineral',
                        crop_type=crop,N_input=rateN,plant_calendar=plantidx,harvest_calendar=harvestidx,
                        fert_freq=app_freq,soil_pH = soilph)
        fnitN, fammN, fureaN = self.chem_fert_type()
        ## when there is no country-level information for fertilizer type, use global mean value
        fnitN[(np.isnan(fnitN))&(totalN!=0)] = 0.27
        fammN[(np.isnan(fammN))&(totalN!=0)] = 0.24
        fureaN[(np.isnan(fureaN))&(totalN!=0)] = 0.49
        ## N rate of three types of fertilizer;
        ## the fraction of usage is represented by the 
        ## fractional cropping area that uses the corresponding N fertilizer
        self.NO3_added = chem_N_tocrop
        self.nitN_area = fnitN * croparea
        self.TAN_added = chem_N_tocrop
        self.ammN_area = fammN * croparea
        self.urea_added = chem_N_tocrop
        self.ureaN_area = fureaN * croparea
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

    ## Simulation: Cat C manure stored in open environment (as solid) //// under development 07/Sep
    ## manure/animal waste is left on open land to reduce water content (evaporation+infiltration); 
    ## water pool is determined by:
    ## 1) mositure equilibrium (minimum amount of manure moisture content) 2) water transferred from housing
    ## NEW) incorporate a [surface compensation point] scheme in the model as Cat B [MMS barn solid]
    ## NEW) incorporate a [soil pool] in the model; with additional TAN loss through aqueous/gaseous diffusion to soil
    ##    and infiltration (subsurface leaching) to soil
    ##       gaseouos concentrations (chi): chi_atm; chi_surface; chi_bulk; chi_soil
    ##       aqueous concrntrations: [TAN]_bulk; [TAN]_surf; [TAN]_soil (TAN conc at the interface between manure and land/soil)
    ##       equilibrium between aqueous and gaseous conc: chi = KNH3 * [TAN]
    ##       fluxes (F_upwards): F_atm, surface to atmosphere; F_tosurf: manure to surface;;; 
    ##                   F_atm=(chi_surf-chi_atm)/R_ab;  F_runoff = qr*[TAN]_surf
    ##                   F_tosurfaq=([TAN]_bulk-[TAN]_surf)/R_manureaq; F_tosurfgas = (chi_bulk - chi_surf)/R_manuregas
    ##                   F_tosurf(aq+gas) = F_atm + F_runoff
    ##       fluxes (F_downwards): F_difftosoilaq: TAN from manure to interface layer between manure and soil (aq diffusion);
    ##                             F_difftosoilgas: TAN from manure to interface layer between manure and soil (gas diffusion);
    ##                             F_infiltosoil: TAN from manure to interface layer between manure and soil (infiltration)  
    ##                             F_soildiffusionaq: aqueous diffusion to deeper soil; 
    ##                             F_soildiffusiongas: gaseous diffusion to deeper soil;
    ##                             F_soilinfiltration: infiltration/subsurface leaching of TAN to deeper soil
    ##                   F_difftosoilaq=([TAN]_bulk-[TAN]_soil)/R_manureaq; F_difftosoilgas=(chi_bulk-chi_soil)/R_manuregas
    ##                   F_infiltosoil=[TAN]_bulk*kinfil;
    ##                   F_soildiffusion=[TAN]_soil/R_soil; F_soilinfil=[TAN]_soil*qpsoil; 
    ## Note) soil pool is prognostic pools, each pool is calucalted dependent on its value fro mthe previous time step, 
    ##       rather than solved (assume to be equilibrium) like the [surface compensation point]
    ##                   kinfil and qpsoil are different for both infiltration processes (manrue to interface; interface to deeper soil)
    ##    ([TAN]_bulk is the prognostic variable and is determined by source and loss based on the mass balance approach)
    ##    (solve [TAN]_surf: [TAN]_surf = ((1/R_manureaq+KNH3/R_manuregas)*[TAN]_bulk + chi_atm/R_ab)/
    ##                                     (qr+KNH3/R_ab+1/R_manureaq+KNH3/R_manuregas)
    ## Note) NO3- pools in bulk manure and soil have similar processes as TAN, which do not include 1) adsorption, 2) gaseous diffustion
    ##       NO3- aqueous diffusion is moderated by a scaling factor regarding the different diffusivity of NO3- from NH4+    
    '''def slurry_app_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: ')
        MMS_area["mms_open_solid_area"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_open_solid*MMS_area_factor["mms_open_solid"]
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):

                ## note the equations for the pools are similar to "barn liquid sim"; see comments above
                ## daily manure and urine on per unit area
                self.manure[dd+1] = self.manure_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
                ## manure pool: with manure input (washoff has not been taken into account here)
                self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]

                ## water amount in "solid" manure
                ## self.manure refers to the DM mass; therefore, total manure mass = DM mass/DM%
                ## water in the "solid" manure = water% x total manure mass
                self.manure_water[dd+1] = self.water_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd]*self.mois_coeff[dd+1]/100

                ## water pool
                ## water pool of [MMS open solid] is determined by:
                ##   source: manure water, rain (precipitation)
                ##   loss: evaporation, infiltration to soil/interface
                ## minimum water amount is equivalent to the mositure equilibrium content of the manure  
                ## Note: infiltration of manure water to soil is assumed to be 10mm/day (Vira et al., 2020GMD)
                ##       and this should be differentiate with TAN infiltration to soil
                ## justify the amount of water that is available for infiltration; daily maximum infiltration is 10mm;
                ## convert g/m^2 to m
                # self.qinfil[dd+1] = (self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd])/1e6
                # self.qinfil[dd+1] = (self.Total_water_pool[dd]-absorb_factor*self.manure_pool[dd])/1e6
                # self.qinfil[dd+1][self.qinfil[dd+1]>dailymaxinfil]=dailymaxinfil
                # self.qinfil[dd+1][self.qinfil[dd+1]<0] = 0.0
                ## covert m/day to m/s
                self.qinfil[dd+1] = dailymaxinfil/(timestep*3600)
                ## soil infiltration flux is given in m/s
                self.qpsoil[dd+1] = infiltration_rate_method(dailyinfil=self.qinfil[dd+1]*timestep*3600,
                                                                theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1])
                self.qpsoil[dd+1][self.qpsoil[dd+1]<0] = 0.0
                ## justify the water amount; infil flux is the infiltration within the manure (between 0-10mm/day)
                # water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd]-self.qinfil[dd+1]*timestep*3600*1e6
                # self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] + self.rainfall[dd+1][water_idx>0] + \
                #                                                 self.manure_water[dd+1][water_idx>0] - \
                #                                                     self.evap_sim[dd][water_idx>0] - \
                #                                                         self.qinfil[dd+1][water_idx>0]*timestep*3600*1e6
                # self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd][water_idx<=0] + self.rainfall[dd+1][water_idx<=0] + \
                #                                                 self.manure_water[dd+1][water_idx<=0] 
                ## maximum water capcity of the manure pool; assuming manure water holding capcity + maximum infiltration
                # max_wc = absorb_factor*self.manure_pool[dd] + dailymaxinfil*1e6
                ## justipy whether current water pool exceeds the maximum water holidng capacity
                # water_idx2 =  max_wc - self.Total_water_pool[dd+1]
                # water_idx2 = self.rainfall[dd+1]/1e6 +self.water_added[dd+1]/1e6 - (self.persm[dd+1]-self.soilmoist_data[dd+1])*z_thickness
                water_idx2 = self.rainfall[dd+1]/1e6 - self.qinfil[dd+1]
                ## if exceeds: the surplus amount of water acts of "washoff" water, and the water pool equals the maximum wc
                ## rain available for washoff has the unit of mm (as an accumulation of a day's rainfall in mm)
                self.rain_avail_washoff[dd+1][water_idx2>0] = water_idx2
                # self.Total_water_pool[dd+1][water_idx2<0] = max_wc[water_idx2<0]
                ## if not exceeds: "washoff" water is 0 as water is absorbed by the manure
                self.rain_avail_washoff[dd+1][water_idx2<=0] = 0.0
                self.rain_avail_washoff[dd+1] = self.rain_avail_washoff[dd+1]/(timestep*3600)

                ## water pool of soil interface
                # water_soil_idx = self.qinfil[dd+1]*timestep*3600+z_soil*self.soilmoist[dd+1] - z_soil*self.persm[dd+1]
                # self.water_pool_soil[dd+1][water_soil_idx>0] = (z_soil*self.persm[dd+1][water_soil_idx>0])*1e6
                # self.water_pool_soil[dd+1][water_soil_idx<=0] = (self.qinfil[dd+1][water_soil_idx<=0]*timestep*3600 + \
                #                                                 z_soil*self.soilmoist[dd+1][water_soil_idx<=0])*1e6
                self.soilmoist[dd+1] = self.rainfall[dd+1]/1e6 + self.water_added[dd+1]/1e6 + self.soilmoist_data[dd+1]*z_thickness
                self.soilmoist[dd+1][self.soilmoist[dd+1]>self.persm[dd+1]] = self.persm[dd+1][self.soilmoist[dd+1]>self.persm[dd+1]]

                ## washoff: 1) manure, 2) urea, 3) available org N, 4) reistant org N, 5) TAN
                ## washoff coefficient (m) = washoff water (mm) * washoff (%/mm)
                nonN_washoff_rate = self.rain_avail_washoff[dd+1]*f_washoff_nonN
                N_washoff_rate = self.rain_avail_washoff[dd+1]*f_washoff_N
                self.manure_washoff[dd+1] = nonN_washoff_rate*self.manure_pool[dd]*timestep*3600
                self.urea_washoff[dd+1] =  N_washoff_rate*self.urea_pool[dd]*timestep*3600
                self.avail_N_washoff[dd+1] = N_washoff_rate*self.avail_N_pool[dd]*timestep*3600
                self.resist_N_pool[dd+1] = N_washoff_rate*self.resist_N_pool[dd]*timestep*3600
                self.unavail_N_washoff[dd+1] = N_washoff_rate*self.unavail_N_pool[dd]*timestep*3600
                self.NO3_washoff[dd+1] = self.rain_avail_washoff[dd+1]*self.NO3_pool[dd]*timestep*3600

                ## manure pool: subtracting washoff
                manure_idx = self.manure_pool[dd+1] - self.manure_washoff[dd+1]
                self.manure_pool[dd+1][manure_idx>0] = manure_idx[manure_idx>0]
                self.manure_pool[dd+1][manure_idx<=0] = 0.0  ## manure has been washed off 
                ## water content of the manure
                # vtotal,manurewc,manure_WFPS = manure_properties(solidmass=self.manure_pool[dd+1],
                #                                         watermass=self.Total_water_pool[dd+1])
                # manurewc[manurewc>manure_porosity] = manure_porosity
                # manure_WFPS[manure_WFPS>1.0] = 1.0
                # manure_torl = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='aqueous')
                # manure_torg = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='gaseous')

                ## manure resistance
                ## manure resistance is determined by: R = z/(2*tor_manure*D); 
                ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not considered here
                ## layer thickness of manure (DM+water) in meter: z = Vmanure
                # z_total = vtotal
                # self.R_manurel[dd+1] = z_total/(2*self.D_aq_NH4[dd+1]*manure_torl)
                # self.R_manureg[dd+1] = z_total/(2*self.D_air_NH3[dd+1]*manure_torg)
                # self.R_manureg[dd+1][manure_torg!=0] = z_total[manure_torg!=0]/\
                #                                     (2*self.D_air_NH3[dd+1][manure_torg!=0]*manure_torg[manure_torg!=0])  
                # ## when water content is zero, gaseous diffusion is ceased by infinite resistance
                # self.R_manureg[dd+1][manure_torg==0] = np.inf
                
                ## soil resistance
                ## soil resistance = thickness of the source layer (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                # tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                # tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                self.R_soilaq[dd+1] = self.app_depth/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg[dd+1] = self.app_depth/(tor_soil_gas*self.D_air_NH3[dd+1])
                self.R_soilaq_down[dd+1] = d_deepsoil/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg_down[dd+1] = d_deepsoil/(tor_soil_gas*self.D_air_NH3[dd+1])
                # self.R_soilg[dd+1][tor_soil_gas!=0] = d_deepsoil/\
                #                                 (tor_soil_gas[tor_soil_gas!=0]*self.D_air_NH3[dd+1][tor_soil_gas!=0])
                # ## when water content is zero, gaseous diffusion is ceased by infinite resistance
                # self.R_soilg[dd+1][tor_soil_gas==0] = np.inf

                ## N input in multiple forms
                self.urea[dd+1] = self.urea_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
                self.avail_N[dd+1] = self.avail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
                self.resist_N[dd+1] = self.resist_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
                self.unavail_N[dd+1] = self.unavail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                        self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                        self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
                ## TAN from housing to storage
                self.TAN[dd+1] = self.TAN_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1]) - self.urea_washoff[dd+1]
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]
                ## Org N pools in various forms
                self.avail_N_pool[dd+1] = (self.avail_N_pool[dd] - self.avail_N_washoff[dd+1])*(1 - self.daily_Na_decomp_rate[dd+1]) + \
                    self.avail_N[dd+1] 
                self.resist_N_pool[dd+1] = (self.resist_N_pool[dd] - self.resist_N_washoff[dd+1])* (1 - self.daily_Nr_decomp_rate[dd+1]) + \
                    self.resist_N[dd+1] 
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] - self.unavail_N_washoff[dd+1]
                
                ## TAN pool (different from [MMS barn (liquid,solid)])
                ## Note: source of TAN pool: TAN production from 1) urea hydrolysis, 2) decomposition of org N and 3) input of TAN from housing
                ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil (aq+gas), 
                ##                         3) infiltration (subsurface leaching) to soil, and 4) nitrification
                ##       only aqueous phase diffusion is considered, and gaseous diffusion is not considered in this study
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusiveflux_aq[dd] - self.diffusiveflux_gas[dd]-\
                    self.infilflux[dd] - self.nitrif_NO3[dd] - self.TAN_washoff[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]
                ## TAN pool in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6

                ## TAN conc; g/m3
                ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to manure) NH4+
                ## gaseous NH3 in air-filled pore space; epsilon(manure porosity) - theta
                ## aqueous TAN in water-filled pore space; theta
                ## solid phase TAN adsorbed on solid particles; 1 - epsilon
                ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd = 1.0 m^3/m^3 
                ## Kd = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                ## manure density varies, ~ 0.3-1.9 g/cm^3, we assmume the porosity of mamnure is 40%
                ## the volume (thickness) is determined by solid mass and bulk density (derived from porosity)
                ## NH3(g) = KNH3*[TAN(aq)]
                ## MTAN = Vmanure*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
                ## so, [TAN(aq)] = MTAN/Vmanure * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
                KNH3 = self.Henry_constant[dd+1]/(self.cc_H + self.k_NH4[dd+1])
                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_thickness[self.soilmoist[dd+1]!=0]*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd))
                ## TAN molar conc; mol/L
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
                ## NH3 conc in bulk manure
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when manure water content equals to porosity
                self.NH3_gas_bulk[dd+1][manurewc==manure_porosity] = 0.0

                ## TAN conc at the surface (solved); atmospheric NH3 is ignored
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_manurel[dd+1]+KNH3/self.R_manureg[dd+1])/\
                        (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_manureg[dd+1])+1/self.R_manurel[dd+1]))/\
                            (14*1000)
                ## TAN surface conc in g/m3
                TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
                ## in ug
                NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

                ## determining the maximum emission; 
                emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
                # if dd+1 == 143:
                #     print("emiss_idx: ",emiss_idx[130,363]) 
                ## final emission flux
                # self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
                ## determining the maximum TAN runoff;
                runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
                # if dd+1 == 143:
                #     print("runoff_idx: ",runoff_idx[130,363]) 
                ## determining the maximum TAN aqueous diffusion to soil interface;
                ## diffusion is considered to be unidirectional - fro bulk manure to soil interface
                diffaq_idx = self.TAN_amount[dd+1]/self.R_soilaq_down[dd+1]*timestep*3600
                diffaq_idx[diffaq_idx<=0] = 0.0
                ## when soil mositure is 0, aqueous diffusion stops
                diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
                # if dd+1 == 143:
                #     print("diffaq_idx: ",diffaq_idx[130,363]) 
                ## determining the maximum TAN gaseous diffusion to soil interface
                diffgas_idx = self.NH3_gas_bulk[dd+1]/self.R_soilg_down[dd+1]*timestep*3600
                diffgas_idx[diffgas_idx<0] = 0.0
                ## when the soil moisture reaches the saturation, gaseous diffusion stops
                diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
                # if dd+1 == 143:
                #     print("diffgas_idx: ",diffgas_idx[130,363]) 
                ## determining the maximum infiltration of TAN
                leaching_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600
                # if dd+1 == 143:
                #     print("infil_idx: ",infil_idx[130,363]) 

                ## nirification rate of TAN in bulk manure; daily maximum nitrification rate is 0.1 per day
                KNO3_soil= nitrification_rate_soil(ground_temp=self.T_sim[dd+1],theta=self.soilmoist[dd+1],
                                            theta_sat=self.persm[dd+1],fer_type='manure')*timestep*3600
                KNO3_soil[KNO3_soil>0.1] = 0.1
                ## correction for WPFS response
                KNO3_soil[KNO3_manure<0.0] = 0.0
                KNO3_soil[np.isnan(KNO3_soil)] = 0.0
                ## determining the maximum nitrification of TAN
                # nitrif_idx = KNO3_manure*self.TAN_pool[dd]
                nitrif_idx = KNO3_soil*self.TAN_pool[dd]

                ## fluxes from bulk soil to deeper soil and chemical loss
                ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                all_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + leaching_idx + nitrif_idx
                # if dd+1 == 143:
                #     print("manureall_loss: ",manureall_loss[130,363]) 
                loss_idx = self.TAN_pool[dd+1] - all_loss
                # if dd+1 == 143:
                #     print("manureloss_idx: ",manureloss_idx[130,363]) 
                self.NH3_flux[dd+1][loss_idx>=0] = emiss_idx[loss_idx>=0]
                self.TAN_washoff[dd+1][loss_idx>=0] = runoff_idx[loss_idx>=0]
                self.diffusiveflux_aq[dd+1][loss_idx>=0] = diffaq_idx[loss_idx>=0]
                self.diffusiveflux_gas[dd+1][loss_idx>=0] = diffgas_idx[loss_idx>=0]
                self.leachingflux[dd+1][loss_idx>=0] = leaching_idx[loss_idx>=0]
                self.nitrif_NO3[dd+1][loss_idx>=0] = nitrif_idx[loss_idx>=0]
                ## otherwise, we use a weighted distribution
                self.NH3_flux[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (emiss_idx[loss_idx<0]/all_loss[loss_idx<0])
                self.TAN_washoff[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (runoff_idx[loss_idx<0]/all_loss[loss_idx<0])
                self.diffusiveflux_aq[dd+1][loss_idx<0]= self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (diffaq_idx[loss_idx<0]/all_loss[loss_idx<0])
                self.diffusiveflux_gas[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (diffgas_idx[loss_idx<0]/all_loss[loss_idx<0])
                self.leachingflux[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (leaching_idx[loss_idx<0]/all_loss[loss_idx<0])
                self.nitrif_NO3[dd+1][loss_idx<0]= self.TAN_pool[dd+1][loss_idx<0]*\
                                                            (nitrif_idx[loss_idx<0]/all_loss[loss_idx<0])

                ## NO3- pool of soil layer
                NO3_idx = self.NO3_pool[dd] - self.NO3_leaching[dd] - self.NO3_diffusive[dd] - self.NO3_washoff[dd+1]
                self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3[dd+1][NO3_idx>0] + \
                                                    self.NO3_added[dd+1][NO3_idx>0]
                self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3[dd+1][NO3_idx<=0] + self.NO3_added[dd+1][NO3_idx<=0]

                ## NO3- conc of soil layer; g/m3
                self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/(z_thickness*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
                self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                
                ## diffusive aquous NO3 and infiltration of NO3 from bulk manure to soil interface
                NO3_diffidx = self.NO3_amount[dd+1]*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
                NO3_diffidx[NO3_diffidx<0] = 0.0
                NO3_diffidx[self.soilmoist[dd+1]==0] = 0.0
                NO3_leachingidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
                NO3_lossall = NO3_diffidx + NO3_leachingidx 
                loss_idx = self.NO3_pool[dd+1] - NO3_lossall
                self.NO3_diffusive[dd+1][loss_idx>=0] = NO3_diffidx[loss_idx>=0]
                self.NO3_leaching[dd+1][loss_idx>=0] = NO3_leachingidx[loss_idx>=0]
                self.NO3_diffusive[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                        NO3_diffidx[loss_idx<0]/NO3_lossall[loss_idx<0]                                           
                self.NO3_leaching[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                        NO3_leachingidx[loss_idx<0]/NO3_lossall[loss_idx<0]                                         
                
        return'''

    ## Simulation: BROADCASTING - Incorporated disk scheme
    ## vertical profile: 1 - source layer
    ## fluxes: 1. surface runoff  2. NH3 volatilization
    ##         3. percolation flux (subsurface runoff/infiltration)  4. diffusion (aq,gas)
    ## pathways/processes: 1. nitrification  2. plant N uptake inc. ammonium and nitrate   
    def chem_fert_bcdisk_sim(self,start_day_idx,end_day_idx,chem_fert_type,disk_depth,crop=None):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [broadcasting - incorporated disk], fertilizer placement depth is: '+str(disk_depth*100)+' cm')
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        ## NH4+ adsorption coefficient is an emperical function of the clay content of soils
        Kd = ammonium_adsorption(clay_content=soilclay)
        ## thickness of the source layer is the depth of the incorporated disk
        z_sourcelayer = disk_depth
        ## different types of fertilizers: 1) nitrate, 2) ammonium, 3) urea N
        ## urea hydrolysis affects soil pH (soil pH increases as urea hydrolysis consumed H+)
        ## ammonium/nitrate has little effects on changing soil pH
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

        
        for dd in np.arange(start_day_idx,end_day_idx):
            
            ## soil resistance
            ## soil resistance = half of the thickness of the corresponding soil layer (z) / (tortuosity for diffusion x diffusivity of the species)
            ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
            tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
            tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
            self.R_soilaq[dd+1] = (disk_depth/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg[dd+1] = (disk_depth/2)/(tor_soil_gas*self.D_air_NH3[dd+1])
            ## the distance of diffusion to the deeper soil is set to be from the middle point of the disk layer 
            ## to the middle point of the second soil layer
            self.R_soilaq_down[dd+1] = (p_2ndsoil-disk_depth/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg_down[dd+1] = (p_2ndsoil-disk_depth/2)/(tor_soil_gas*self.D_air_NH3[dd+1])

            if chem_fert_type == 'urea':
                ## Urea pool
                urea_idx = self.urea_pool[dd] - self.TAN_prod[dd] - self.urea_washoff[dd] - \
                            self.urea_diff[dd] - self.urealeaching[dd]
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]
                ## urea concentration of the source layer and at the surface
                self.urea_amount[dd+1] = self.urea_pool[dd+1]/(z_sourcelayer*self.soilmoist[dd+1])
                urea_surf_amount = self.urea_amount[dd+1]/(self.R_soilaq[dd+1]*self.rain_avail_washoff[dd+1]+1)
                ## pathways
                urea_washoffidx = self.rain_avail_washoff[dd+1] * urea_surf_amount*timestep*3600
                urea_diffidx = self.urea_amount[dd+1]/self.R_soilaq_down[dd+1] * timestep*3600
                urea_leachingidx = self.qpsoil[dd+1] * self.urea_amount[dd+1]*timestep*3600
                urea_toTANidx = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd+1]
                ## determining the fluxes: washoff, diffusion, leaching, hydrolysis
                all_loss = urea_washoffidx + urea_diffidx + urea_leachingidx + urea_toTANidx
                loss_idx = self.urea_pool[dd+1] - all_loss
                self.urea_washoff[dd+1][loss_idx>=0] = urea_washoffidx[loss_idx>=0]
                self.urea_diff[dd+1][loss_idx>=0] = urea_diffidx[loss_idx>=0]
                self.urealeaching[dd+1][loss_idx>=0] = urea_leachingidx[loss_idx>=0]
                self.TAN_prod[dd+1][loss_idx>=0] = urea_toTANidx[loss_idx>=0]
                self.urea_washoff[dd+1][loss_idx<0] = urea_washoffidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.urea_diff[dd+1][loss_idx<0] = urea_diffidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.urealeaching[dd+1][loss_idx<0]= urea_leachingidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.TAN_prod[dd+1][loss_idx<0] = urea_toTANidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
            
            ## TAN pool 
            ## Note: source of TAN pool: input of TAN from fertilizer (ammonium or urea hydrolysis)
            ##                              nitrate does not contribute to the TAN pool
            ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil (aq+gas), 
            ##                         3) infiltration (subsurface leaching) to soil, 4) surface runoff
            ##                         5) nitrification, 6) plant uptake
            TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsoil_aq[dd] - self.diffusivefluxsoil_gas[dd]-\
                self.leachingflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd] - self.ammN_uptake[dd]
            self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
            self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

            ## TAN conc; g/m3
            ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to manure) NH4+
            ## gaseous NH3 in air-filled pore space; epsilon(soil porosity) - theta
            ## aqueous TAN in water-filled pore space; theta
            ## solid phase TAN adsorbed on solid particles; 1 - epsilon
            ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd is derived from an emperical relationship with soil clay content (ref: DNDC model)
            ## NH3(g) = KNH3*[TAN(aq)]
            ## MTAN = Vsoil*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
            ## so, [TAN(aq)] = MTAN/Vsoil * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
            KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
            self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## TAN molar conc; mol/L
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
            ## NH3 conc in soil pool
            self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
            ## NH3 conc is 0 when soil water content equals to porosity (saturated soil water content)
            self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## TAN conc at the surface (solved); atmospheric NH3 is ignored
            self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_soilaq[dd+1]+KNH3/self.R_soilg[dd+1])/\
                    (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_soilg[dd+1])+1/self.R_soilaq[dd+1]))/\
                        (14*1000)
            ## TAN surface conc in g/m3
            TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
            ## Gaseous NH3 at the surface
            self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
            ## in ug
            NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

            ## determining the maximum emission; 
            emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
            ## determining the maximum TAN runoff at the surface;
            runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
            ## determining the maximum TAN aqueous diffusion to deeper soil;
            diffaq_idx = self.TAN_amount[dd+1]/self.R_soilaq_down[dd+1]*timestep*3600
            diffaq_idx[diffaq_idx<=0] = 0.0
            ## when soil mositure is 0, aqueous diffusion stops
            diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
            ## determining the maximum TAN gaseous diffusion to deeper soil
            diffgas_idx = self.NH3_gas_bulk[dd+1]/self.R_soilg_down[dd+1]*timestep*3600
            diffgas_idx[diffgas_idx<0] = 0.0
            ## when the soil moisture reaches the saturation, gaseous diffusion stops
            diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
            ## determining the maximum infiltration/subsurface leaching of TAN
            leaching_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600

            ## nirification rate of TAN in the soil layer; daily maximum nitrification rate is 0.1 per day
            KNO3_soil= nitrification_rate_soil(ground_temp=self.T_sim[dd+1],theta=self.soilmoist[dd+1],
                                            theta_sat=self.persm[dd+1],pH=sim_pH[dd+1],fer_type='mineral')*timestep*3600
            KNO3_soil[KNO3_soil>0.1] = 0.1
            ## correction for WPFS response
            KNO3_soil[KNO3_soil<0.0] = 0.0
            KNO3_soil[np.isnan(KNO3_soil)] = 0.0
            ## determine the aqueous NH4+ fraction of TAN
            f_NH4 = self.soilmoist[dd+1]/(self.soilmoist[dd+1]+\
                KNH3*(self.persm[dd+1]-self.soilmoist[dd+1])+\
                (1-self.persm[dd+1])*Kd)*(sim_ccH[dd+1]/(sim_ccH[dd+1]+self.k_NH4[dd+1]))
            ## determining the maximum nitrification of TAN
            nitrif_idx = KNO3_soil*self.TAN_pool[dd+1]*f_NH4
            
            ## plant N uptake inc. ammonium, nitrate
            soilammNuptake_idx, soilnitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool[dd]*f_NH4,\
                                    Nnit=self.NO3_pool[dd],temp=self.T_sim[dd+1])
            
            ## no uptake after harvesting
            soilammNuptake_idx[harvestidx<(dd+1)] = 0.0
            soilnitNuptake_idx[harvestidx<(dd+1)] = 0.0

            ## all loss pathways
            all_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + leaching_idx + nitrif_idx + soilammNuptake_idx
            loss_idx = self.TAN_pool[dd+1] - all_loss
            ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
            self.NH3_flux[dd+1][loss_idx>=0] = emiss_idx[loss_idx>=0]
            self.TAN_washoff[dd+1][loss_idx>=0] = runoff_idx[loss_idx>=0]
            self.diffusivefluxsoil_aq[dd+1][loss_idx>=0] = diffaq_idx[loss_idx>=0]
            self.diffusivefluxsoil_gas[dd+1][loss_idx>=0] = diffgas_idx[loss_idx>=0]
            self.leachingflux[dd+1][loss_idx>=0] = leaching_idx[loss_idx>=0]
            self.nitrif_NO3_sourcelayer[dd+1][loss_idx>=0] = nitrif_idx[loss_idx>=0]
            self.ammN_uptake[dd+1][loss_idx>=0] = soilammNuptake_idx[loss_idx>=0]
            ## otherwise, we use a weighted distribution
            self.NH3_flux[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (emiss_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.TAN_washoff[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (runoff_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.diffusivefluxsoil_aq[dd+1][loss_idx<0]= self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (diffaq_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.diffusivefluxsoil_gas[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (diffgas_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.leachingflux[dd+1][loss_idx<0] = self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (leaching_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.nitrif_NO3_sourcelayer[dd+1][loss_idx<0]= self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (nitrif_idx[loss_idx<0]/all_loss[loss_idx<0])
            self.ammN_uptake[dd+1][loss_idx<0]= self.TAN_pool[dd+1][loss_idx<0]*\
                                                        (soilammNuptake_idx[loss_idx<0]/all_loss[loss_idx<0])

            ## NO3- pool of soil layer
            ## sources: 1) nitrification of TAN, and 2) input of nitrate fertilizer
            ## losses: 1) leaching to deeper soil, 2) diffusion to deeper soil
            ##         3) surface runoff, 4) plant uptake
            NO3_idx = self.NO3_pool[dd] - self.NO3_leaching[dd] - self.NO3_diffusivedeep[dd] - self.NO3_washoff[dd+1] - \
                        self.nitN_uptake[dd]
            self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + \
                                                self.NO3[dd+1][NO3_idx>0]
            self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0]

            ## NO3- conc of soil layer; g/m3
            self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                        (z_sourcelayer*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
            self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
            
            ## diffusive aquous NO3 and infiltration of NO3 to deeper soil
            NO3_diffidx = self.NO3_amount[dd+1]*(f_DNO3/self.R_soilaq_down[dd+1])*timestep*3600
            NO3_diffidx[NO3_diffidx<0] = 0.0
            NO3_diffidx[self.soilmoist[dd+1]==0] = 0.0
            NO3_leachingidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
            NO3_washoffidx = self.rain_avail_washoff[dd+1]*self.NO3_amount[dd+1]*timestep*3600
            NO3_lossall = NO3_diffidx + NO3_leachingidx + NO3_washoffidx + soilnitNuptake_idx 
            loss_idx = self.NO3_pool[dd+1] - NO3_lossall 
            self.NO3_diffusivedeep[dd+1][loss_idx>=0] = NO3_diffidx[loss_idx>=0]
            self.NO3_leaching[dd+1][loss_idx>=0] = NO3_leachingidx[loss_idx>=0]
            self.NO3_washoff[dd+1][loss_idx>=0] = NO3_washoffidx[loss_idx>=0]
            self.nitN_uptake[dd+1][loss_idx>=0] = soilnitNuptake_idx[loss_idx>=0]
            self.NO3_diffusivedeep[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                    NO3_diffidx[loss_idx<0]/NO3_lossall[loss_idx<0]                                           
            self.NO3_leaching[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                    NO3_leachingidx[loss_idx<0]/NO3_lossall[loss_idx<0]
            self.NO3_washoff[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                    NO3_washoffidx[loss_idx<0]/NO3_lossall[loss_idx<0]   
            self.nitN_uptake[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                    soilnitNuptake_idx[loss_idx<0]/NO3_lossall[loss_idx<0]                                      
                
        return

    ## Simulation: BROADCASTING - topdressed (surface spreading with very shallow incroporation of ~ 2cm)
    ## vertical profile: 1 - source layer 2 - topsoil layer 3 - deep soil
    ## fluxes: 1. surface runoff  2. NH3 volatilization
    ##         3. percolation flux (subsurface runoff/infiltration)  4. diffusion (aq,gas)
    ## pathways/processes: 1. nitrification  2. plant N uptake inc. ammonium and nitrate
    def chem_fert_bcsurf_sim(self,start_day_idx,end_day_idx,chem_fert_type,crop=None):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [broadcasting - topdressed], fertilizer applied on the field surface')
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)
        ## thickness of the default source layer: 2cm; mid point 1cm
        z_sourcelayer = 0.02
        p_sourcelayer = 0.01
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

        
        for dd in np.arange(start_day_idx,end_day_idx):

            ## soil resistance
            tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
            tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
            ## resistance for diffusion from the surface source layer to the underlying topsoil layer
            self.R_soilaq[dd+1] = (p_topsoil-p_sourcelayer)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg[dd+1] = (p_topsoil-p_sourcelayer)/(tor_soil_gas*self.D_air_NH3[dd+1])
            ## resistance for diffusion from the topsoil layer to the deeper soil
            self.R_soilaq_down[dd+1] = (p_2ndsoil-p_topsoil)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg_down[dd+1] = (p_2ndsoil-p_topsoil)/(tor_soil_gas*self.D_air_NH3[dd+1])

            ## sourcelayer resistance
            ## sourcelayer resistance is determined by: R = (z/2)/(tor_sourcelayer*D); 
            ##        z is the layer thickness of sourcelayer; D is molecular diffusivity of NH4+ in water
            ##        tortuosity dependence for aqueous diffusion is not considered here
            ## layer thickness of sourcelayer (DM+water) in meter: z = Vsourcelayer
            self.R_sourcelayer_aq[dd+1] = (z_sourcelayer/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_sourcelayer_gas[dd+1] = (z_sourcelayer/2)/(tor_soil_gas*self.D_air_NH3[dd+1])


            '''## N input in multiple forms
            self.urea[dd+1] = self.urea_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.avail_N[dd+1] = self.avail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.resist_N[dd+1] = self.resist_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.unavail_N[dd+1] = self.unavail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])

            ## TAN production from urea hydrolysis and the N decomposition rate from dung
            self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                    self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                    self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
            ## TAN from housing to storage
            self.TAN[dd+1] = self.TAN_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])'''

            if chem_fert_type == 'urea':

                ## Urea pool of the surface source layer
                self.urea_pool[dd+1] = self.urea_pool[dd] + self.urea[dd+1]
                ## urea concentration of the source layer and at the surface
                self.urea_amount[dd+1] = self.urea_pool[dd+1]/(z_sourcelayer*self.soilmoist[dd+1])
                urea_surf_amount = self.urea_amount[dd+1]/(self.R_sourcelayer_aq[dd+1]*self.rain_avail_washoff[dd+1]+1)
                ## pathways
                urea_washoffidx = self.rain_avail_washoff[dd+1] * urea_surf_amount*timestep*3600
                urea_diffidx = (self.urea_amount[dd+1]-self.urea_soil_amount[dd])/self.R_soilaq[dd+1] * timestep*3600
                urea_infilidx = self.qpsoil[dd+1] * self.urea_amount[dd+1]*timestep*3600
                urea_toTANidx = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd+1]
                ## determining the fluxes: washoff, diffusion, infiltration, hydrolysis
                all_loss = urea_washoffidx + urea_diffidx + urea_infilidx + urea_toTANidx
                loss_idx = self.urea_pool[dd+1] - all_loss
                self.urea_washoff[dd+1][loss_idx>=0] = urea_washoffidx[loss_idx>=0]
                self.ureadiffusivefluxsourcelayer[dd+1][loss_idx>=0] = urea_diffidx[loss_idx>=0]
                self.urea_infil[dd+1][loss_idx>=0] = urea_infilidx[loss_idx>=0]
                self.TAN_prod[dd+1][loss_idx>=0] = urea_toTANidx[loss_idx>=0]
                self.urea_washoff[dd+1][loss_idx<0] = urea_washoffidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.ureadiffusivefluxsourcelayer[dd+1][loss_idx<0] = urea_diffidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.urealeaching[dd+1][loss_idx<0]= urea_infilidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.TAN_prod[dd+1][loss_idx<0] = urea_toTANidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                ## determining urea pool after subtracting all loss pathways
                self.urea_pool[dd+1] = self.urea_pool[dd+1] - self.TAN_prod[dd+1] - self.urea_washoff[dd+1] - \
                            self.ureadiffusivefluxsourcelayer[dd+1] - self.urea_infil[dd+1]
                ## getting rid of rounding error
                self.urea_pool[dd+1][self.urea_pool[dd+1]<0.0] = 0.0
                self.urea_amount[dd+1] = self.urea_pool[dd+1]/(z_sourcelayer*self.soilmoist[dd+1])

                ## urea pool of the topsoil layer
                self.urea_pool_soil[dd+1] = self.urea_pool_soil[dd] + self.urea_infil[dd+1] + self.ureadiffusivefluxsourcelayer[dd+1]
                ## urea concentration of the topsoil layer
                self.urea_soil_amount[dd+1] = self.urea_pool_soil[dd+1]/(z_topsoil*self.soilmoist[dd+1])
                ## pathways
                soilurea_diffidx = self.urea_soil_amount[dd+1]/self.R_soilaq_down[dd+1]*timestep*3600
                soilurea_leachingidx = self.qpsoil[dd+1]*self.urea_soil_amount[dd+1]*timestep*3600
                soilurea_toTANidx = self.daily_urea_hydro_rate[dd+1]*self.urea_pool_soil[dd+1]
                ## determining the fluxes: washoff, diffusion, leaching, hydrolysis
                soilall_loss = soilurea_diffidx + soilurea_leachingidx + soilurea_toTANidx
                soilloss_idx = self.urea_pool_soil[dd+1] - soilall_loss
                self.urea_diff[dd+1][soilloss_idx>=0] = soilurea_diffidx[soilloss_idx>=0]
                self.urealeaching[dd+1][soilloss_idx>=0] = soilurea_leachingidx[soilloss_idx>=0]
                self.TAN_prod_soil[dd+1][soilloss_idx>=0] = soilurea_toTANidx[soilloss_idx>=0]
                self.urea_diff[dd+1][soilloss_idx<0] = soilurea_diffidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                self.urealeaching[dd+1][soilloss_idx<0]= soilurea_leachingidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                self.TAN_prod_soil[dd+1][soilloss_idx<0] = soilurea_toTANidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                ## deterining urea pool of the topsoil layer after subtracting all losses
                self.urea_pool_soil[dd+1] = self.urea_pool_soil[dd+1] - self.urea_diff[dd+1] - self.urealeaching[dd+1] - \
                                            self.TAN_prod_soil[dd+1]
                self.urea_pool_soil[dd+1][self.urea_pool_soil[dd+1]<0.0] = 0.0
                self.urea_soil_amount[dd+1] = self.urea_pool_soil[dd+1]/(z_topsoil*self.soilmoist[dd+1])

            '''## Org N pools in various forms
            self.avail_N_pool[dd+1] = (self.avail_N_pool[dd] - self.avail_N_washoff[dd+1])*(1 - self.daily_Na_decomp_rate[dd+1]) + \
                self.avail_N[dd+1] 
            self.resist_N_pool[dd+1] = (self.resist_N_pool[dd] - self.resist_N_washoff[dd+1])* (1 - self.daily_Nr_decomp_rate[dd+1]) + \
                self.resist_N[dd+1] 
            self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] - self.unavail_N_washoff[dd+1]'''

            ## TAN pool (source layer)
            ## Note: source of TAN pool: 1) input of TAN from fertilizer (ammonium or urea hydrolysis)
            ##                              nitrate does not contribute to the TAN pool
            ##                           2) diffusion of TAN from the underlying topsoil layer to the source layer
            ##                                 this is a bi-directional diffusion scheme as the N is removed from the 
            ##                                  the source layer through surface runoff and the NH3 volatilization, which
            ##                                   will lead to decreasing of the TAN concentration of the source layer. 
            ##                                  When the TAN concentration of the underlying topsoil layer is higher than
            ##                                   the source layer TAN concentration, diffusion takes place from the topsoil layer.
            ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) surface runoff
            ##                         3) diffusion to the topsoil (aq+gas) layer, 
            ##                         4) infiltration to the topsoil layer, and 5) nitrification
            ##       Note that there is no plant N uptake at the surface source layer in this scheme
            '''TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
            self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
            self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]
            ## TAN pool in ug
            TAN_pool_ug = self.TAN_pool[dd+1] * 1e6'''

            ## TAN pool with inputs
            self.TAN_pool[dd+1] = self.TAN_pool[dd]+self.TAN_prod[dd+1]+self.TAN[dd+1] + \
                    self.diffusivefluxup_aq[dd] + self.diffusivefluxup_gas[dd]

            ## TAN conc; g/m3
            ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to manure) NH4+
            ## gaseous NH3 in air-filled pore space; epsilon(soil porosity) - theta
            ## aqueous TAN in water-filled pore space; theta
            ## solid phase TAN adsorbed on solid particles; 1 - epsilon
            ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd is derived from an emperical relationship with soil clay content (ref: DNDC model)
            ## NH3(g) = KNH3*[TAN(aq)]
            ## MTAN = Vsoil*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
            ## so, [TAN(aq)] = MTAN/Vsoil * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
            KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
            self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## TAN molar conc; mol/L
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
            ## NH3 conc in the source layer
            self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
            ## NH3 conc is 0 when sourcelayer water content equals to porosity
            self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## TAN conc at the surface (solved); atmospheric NH3 is ignored
            self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_sourcelayer_aq[dd+1]+KNH3/self.R_sourcelayer_gas[dd+1])/\
                                    (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_sourcelayer_gas[dd+1])+1/self.R_sourcelayer_aq[dd+1]))/\
                                        (14*1000)

            ## TAN surface conc in g/m3
            TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
            ## Gaseous NH3 at the surface
            self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
            ## in g
            NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

            ## determining the maximum emission; 
            emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
            ## determining the maximum TAN runoff;
            ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
            ## else need to multiply by timestep*3600
            runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
            ## determining the maximum TAN aqueous diffusion to the topsoil layer;
            ## diffusion is considered to be bi-directional - either from source layer to topsoil or the other way round
            diffaq_idx = (self.TAN_amount[dd+1]-self.TAN_soil_amount[dd])/self.R_soilaq[dd+1]*timestep*3600
            ## if this term <0, diffusion takes place from the underlying topsoil layer
            diffaq_idx[diffaq_idx<=0] = 0.0
            ## when soil mositure is 0, aqueous diffusion stops
            diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
            ## determining the maximum TAN gaseous diffusion to the topsoil layer
            diffgas_idx = (self.NH3_gas_bulk[dd+1]-self.NH3_gas_soil[dd])/self.R_soilg[dd+1]*timestep*3600
            diffgas_idx[diffgas_idx<0] = 0.0
            ## when the soil moisture reaches the saturation, gaseous diffusion stops
            diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
            ## determining the maximum infiltration of TAN from the source layer to the topsoil layer
            infil_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600
            # infil_idx = 0*self.TAN_amount[dd+1]*timestep*3600
            ## nirification rate of TAN in the source layer; daily maximum nitrification rate is 0.1 per day
            KNO3_soil= nitrification_rate_soil(ground_temp=self.T_sim[dd+1],theta=self.soilmoist[dd+1],
                                            theta_sat=self.persm[dd+1],pH=sim_pH[dd+1],fer_type='mineral')*timestep*3600
            KNO3_soil[KNO3_soil>0.1] = 0.1
            ## correction for WPFS response
            KNO3_soil[KNO3_soil<0.0] = 0.0
            KNO3_soil[np.isnan(KNO3_soil)] = 0.0
            ## determining the maximum nitrification of TAN
            # nitrif_idx = KNO3_sourcelayer*self.TAN_pool[dd]
            ## fraction of [NH4+(aq)]
            f_NH4 = self.soilmoist[dd+1]/(self.soilmoist[dd+1]+\
                KNH3*(self.persm[dd+1]-self.soilmoist[dd+1])+\
                (1-self.persm[dd+1])*Kd)*(sim_ccH[dd+1]/(sim_ccH[dd+1]+self.k_NH4[dd+1]))
            nitrif_idx = KNO3_soil*self.TAN_pool[dd+1]*f_NH4

            ## fluxes from sourcelayer to the topsoil layer
            ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
            sourcelayerall_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + infil_idx + nitrif_idx
            sourcelayerloss_idx = self.TAN_pool[dd+1] - sourcelayerall_loss
            
            self.NH3_flux[dd+1][sourcelayerloss_idx>=0] = emiss_idx[sourcelayerloss_idx>=0]
            self.TAN_washoff[dd+1][sourcelayerloss_idx>=0] = runoff_idx[sourcelayerloss_idx>=0]
            self.diffusivefluxsourcelayer_aq[dd+1][sourcelayerloss_idx>=0] = diffaq_idx[sourcelayerloss_idx>=0]
            self.diffusivefluxsourcelayer_gas[dd+1][sourcelayerloss_idx>=0] = diffgas_idx[sourcelayerloss_idx>=0]
            self.infilflux[dd+1][sourcelayerloss_idx>=0] = infil_idx[sourcelayerloss_idx>=0]
            self.nitrif_NO3_sourcelayer[dd+1][sourcelayerloss_idx>=0] = nitrif_idx[sourcelayerloss_idx>=0]
            ## otherwise, we use a weighted distribution
            self.NH3_flux[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (emiss_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])
            self.TAN_washoff[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (runoff_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])
            self.diffusivefluxsourcelayer_aq[dd+1][sourcelayerloss_idx<0]= self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (diffaq_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])
            self.diffusivefluxsourcelayer_gas[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (diffgas_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])
            self.infilflux[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (infil_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])
            self.nitrif_NO3_sourcelayer[dd+1][sourcelayerloss_idx<0]= self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                                        (nitrif_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0])

            ## determine TAN pool after subtracting all loss pathways
            self.TAN_pool[dd+1] = self.TAN_pool[dd+1] - self.NH3_flux[dd+1] - self.diffusivefluxsourcelayer_aq[dd+1] - self.diffusivefluxsourcelayer_gas[dd+1]-\
                self.infilflux[dd+1] - self.nitrif_NO3_sourcelayer[dd+1] - self.TAN_washoff[dd+1]
            ## get rid of rounding error (TAN pool occasionaly become <-1e-13)
            self.TAN_pool[dd+1][self.TAN_pool[dd+1]<0.0] = 0.0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 conc in the source layer
            self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
            ## NH3 conc is 0 when source layer water content reaches saturation
            self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## NO3- pool of the source layer
            ## sources: 1) nitrification of TAN, and 2) input of nitrate fertilizer
            ## losses: 1) infiltration to the underlying topsoil layer, 
            ##         2) diffusion (aq only) to the underlying topsoil layer, 3) surface washoff
            NO3_idx = self.NO3_pool[dd] - self.NO3_infilsourcelayer[dd] - self.NO3_diffusivesourcelayer[dd] - self.NO3_washoff[dd+1]
            self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + self.NO3[dd+1][NO3_idx>0]
            self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0]

            ## NO3- conc of the source layer; g/m3
            self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                        (z_sourcelayer*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
            self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0

            ## diffusive aquous NO3 and infiltration of NO3 from the source layer to the topsoil layer
            NO3_diffidx = (self.NO3_amount[dd+1] - self.NO3_soil_amount[dd])*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
            NO3_diffidx[NO3_diffidx<0] = 0.0
            NO3_diffidx[self.soilmoist[dd+1]==0] = 0.0
            NO3_infilidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
            NO3_washoffidx = self.rain_avail_washoff[dd+1]*self.NO3_amount[dd+1]*timestep*3600
            NO3_lossall = NO3_diffidx + NO3_infilidx + NO3_washoffidx
            sourcelayerloss_idx = self.NO3_pool[dd+1] - NO3_lossall
            self.NO3_diffusivesourcelayer[dd+1][sourcelayerloss_idx>=0] = NO3_diffidx[sourcelayerloss_idx>=0]
            self.NO3_infilsourcelayer[dd+1][sourcelayerloss_idx>=0] = NO3_infilidx[sourcelayerloss_idx>=0]
            self.NO3_washoff[dd+1][sourcelayerloss_idx>=0] = NO3_washoffidx[sourcelayerloss_idx>=0]
            self.NO3_diffusivesourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                    NO3_diffidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                           
            self.NO3_infilsourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                    NO3_infilidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]    
            self.NO3_washoff[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                    NO3_washoffidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                     

            ## TAN pool of the topsoil layer
            ## sources: 1) infiltration from the above source layer, 2) diffusion (aq+gas) source layer to the topsoi layer
            ## losses: 1) diffusion (aq+gas) to the deeper soil, 2) diffusion (aq+gas) to the source layer if possible (bi-directional)
            ##         3) subsurface leaching, 4) nitrification, 5) plant uptake 
            ## TAN pool with inputs
            self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd] + self.infilflux[dd+1] +\
                self.diffusivefluxsourcelayer_aq[dd+1] + self.diffusivefluxsourcelayer_gas[dd+1] + self.TAN_prod_soil[dd+1]

            ## TAN conc of the topsoil layer; g/m3
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                ((z_topsoil-z_sourcelayer)*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 concentration in the soil pore space
            self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
            ## NH3 conc is zero when soil moisture content reaches the saturation
            self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0 

            ## fluxes to deeper soil
            ## TAN loss through aqueous diffusion and leaching to deeper soil
            soildiffaq_idx = self.TAN_soil_amount[dd+1]/(self.R_soilaq_down[dd+1]*timestep*3600)
            soildiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
            soildiffgas_idx = self.NH3_gas_soil[dd+1]/(self.R_soilg_down[dd+1]*timestep*3600)
            soildiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
            
            ## diffusive fluxes to the above source layer when the topsoil TAN concentration is higher
            soilupdiffaq_idx = (self.TAN_soil_amount[dd+1] - self.TAN_amount[dd+1])/self.R_soilaq[dd+1]*timestep*3600
            soilupdiffaq_idx[soilupdiffaq_idx<0] = 0.0 
            soilupdiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
            soilupdiffgas_idx = (self.NH3_gas_soil[dd+1] - self.NH3_gas_bulk[dd+1])/self.R_soilg[dd+1]*timestep*3600
            soilupdiffgas_idx[soilupdiffgas_idx<0] = 0.0 
            soilupdiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            soilleaching_idx = self.qpsoil[dd+1]*self.TAN_soil_amount[dd+1]*timestep*3600
            soilnitrif_idx =  KNO3_soil*self.TAN_pool_soil[dd+1]*f_NH4

            soilammNuptake_idx, soilnitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool_soil[dd+1]*f_NH4,\
                                    Nnit=self.NO3_pool_soil[dd],temp=self.T_sim[dd+1])
            
            ## no uptake after harvesting
            soilammNuptake_idx[harvestidx<(dd+1)] = 0.0
            soilnitNuptake_idx[harvestidx<(dd+1)] = 0.0

            soilall_loss = soildiffaq_idx + soildiffgas_idx + soilleaching_idx + soilnitrif_idx +\
                            soilammNuptake_idx + soilupdiffaq_idx + soilupdiffgas_idx
            soilloss_idx = self.TAN_pool_soil[dd+1] - soilall_loss
            self.diffusivefluxsoil_aq[dd+1][soilloss_idx>=0] = soildiffaq_idx[soilloss_idx>=0]
            self.diffusivefluxsoil_gas[dd+1][soilloss_idx>=0] = soildiffgas_idx[soilloss_idx>=0]
            self.diffusivefluxup_aq[dd+1][soilloss_idx>=0] = soilupdiffaq_idx[soilloss_idx>=0]
            self.diffusivefluxup_gas[dd+1][soilloss_idx>=0] = soilupdiffgas_idx[soilloss_idx>=0]
            self.leachingflux[dd+1][soilloss_idx>=0] = soilleaching_idx[soilloss_idx>=0]
            self.nitrif_NO3_soil[dd+1][soilloss_idx>=0] = soilnitrif_idx[soilloss_idx>=0]
            self.ammN_uptake[dd+1][soilloss_idx>=0] = soilammNuptake_idx[soilloss_idx>=0]
            self.diffusivefluxsoil_aq[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soildiffaq_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.diffusivefluxsoil_gas[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soildiffgas_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.diffusivefluxup_aq[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilupdiffaq_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.diffusivefluxup_gas[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilupdiffgas_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.leachingflux[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilleaching_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.nitrif_NO3_soil[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilnitrif_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
            self.ammN_uptake[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilammNuptake_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]

            ## update TAN pool of the topsoil layer
            self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd+1] - self.diffusivefluxsoil_aq[dd+1] - self.diffusivefluxsoil_gas[dd+1] - \
                - self.diffusivefluxup_aq[dd+1] - self.diffusivefluxup_gas[dd+1] - \
                self.leachingflux[dd+1] - self.nitrif_NO3_soil[dd+1] - self.ammN_uptake[dd+1]
            ## get rid of rounding error (TAN pool occasionaly become <-1e-13)
            self.TAN_pool_soil[dd+1][self.TAN_pool_soil[dd+1]<0.0] = 0.0

            ## TAN conc of the topsoil layer; g/m3
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                ((z_topsoil-z_sourcelayer)*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 concentration in the soil pore space
            self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
            ## NH3 conc is zero when soil moisture content reaches the saturation
            self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0 


            ## NO3- pool of the topsoil layer
            ## sources: 1) nitrification of TAN, 2) infiltration from the above source layer, 
            ##          3) diffusion from the above source layer
            ## losses: 1) leaching to the deeper soil, 2) diffusion (aq only) to the depper soil, 3) plant uptake
            NO3_soil_idx = self.NO3_pool_soil[dd] - self.NO3_leaching[dd] - self.NO3_diffusivedeep[dd] -\
                            self.nitN_uptake[dd] 
            self.NO3_pool_soil[dd+1][NO3_soil_idx>0] = NO3_soil_idx[NO3_soil_idx>0] + self.nitrif_NO3_soil[dd+1][NO3_soil_idx>0] + \
                    self.NO3_infilsourcelayer[dd+1][NO3_soil_idx>0] + self.NO3_diffusivesourcelayer[dd+1][NO3_soil_idx>0]
            self.NO3_pool_soil[dd+1][NO3_soil_idx<=0] = self.nitrif_NO3_soil[dd+1][NO3_soil_idx<=0] + \
                    self.NO3_infilsourcelayer[dd+1][NO3_soil_idx<=0] + self.NO3_diffusivesourcelayer[dd+1][NO3_soil_idx<=0]

            ## NO3- conc of the topsoil layer in g/m3
            self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                                                    ((z_topsoil-z_sourcelayer)*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])   
            self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0

            ## NO3 loss through aqueous diffusion and leaching to deeper soil
            NO3_soildiffidx = self.NO3_soil_amount[dd+1]*(f_DNO3/self.R_soilaq_down[dd+1])*timestep*3600
            NO3_soilleachingidx = self.qpsoil[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
            
            NO3_soilall_loss = NO3_soildiffidx + NO3_soilleachingidx + soilnitNuptake_idx
            soilloss_idx = self.NO3_pool_soil[dd+1] - NO3_soilall_loss
            self.NO3_diffusivedeep[dd+1][soilloss_idx>=0] = NO3_soildiffidx[soilloss_idx>=0]
            self.NO3_leaching[dd+1][soilloss_idx>=0] = NO3_soilleachingidx[soilloss_idx>=0]
            self.nitN_uptake[dd+1][soilloss_idx>=0] = soilnitNuptake_idx[soilloss_idx>=0]
            self.NO3_diffusivedeep[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                            NO3_soildiffidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
            self.NO3_leaching[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                        NO3_soilleachingidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
            self.nitN_uptake[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilnitNuptake_idx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
            
            # if dd+1 == 347:
            #     print('Day: '+str(dd+1))
            #     print('nit uptake idx',soilnitNuptake_idx[183,201])
            #     print('nit uptake',self.nitN_uptake[dd+1,183,201])
            #     print('NO3 soil pool all loss idx',NO3_soilall_loss[183,201])
            #     print('NO3 soil pool idx',soilloss_idx[183,201])
            # elif dd+1 == 348:
            #     print('Day: '+str(dd+1))
            #     print('nit uptake idx',soilnitNuptake_idx[183,201])
            #     print('nit uptake',self.nitN_uptake[dd+1,183,201])
            #     print('NO3 soil pool all loss idx',NO3_soilall_loss[183,201])
            #     print('NO3 soil pool idx',soilloss_idx[183,201])
        return

    ## Simulation: DEEP INJECTION 
    ## vertical profile: 1 - topsoil layer 2 - source layer (where the fertilizer is placed) 3 - deep soil
    ## fluxes: 1. surface runoff  2. NH3 volatilization
    ##         3. percolation flux (subsurface runoff/infiltration)  4. diffusion (aq,gas)
    ## pathways/processes: 1. nitrification  2. plant N uptake inc. ammonium and nitrate
    def chem_fert_deepinjec_sim(self,start_day_idx,end_day_idx,chem_fert_type,injection_depth=0.1,crop=None):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [deep injection], injection depth is: '+str(injection_depth*100)+' cm')
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)
        ## thickness of the source layer; default 10 cm
        # z_sourcelayer = 0.1
        z_sourcelayer = (injection_depth - z_topsoil)*2
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

        
        for dd in np.arange(start_day_idx,end_day_idx):

            ## soil resistance
            tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
            tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
            ## resistance for diffusion from the topsoil layer to the emitting surface
            self.R_soilaq[dd+1] = (z_topsoil/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg[dd+1] = (z_topsoil/2)/(tor_soil_gas*self.D_air_NH3[dd+1])
            ## resistance to the deeper soil from the source layer (deeper soil)
            self.R_soilaq_down[dd+1] = (p_2ndsoil-injection_depth)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_soilg_down[dd+1] = (p_2ndsoil-injection_depth)/(tor_soil_gas*self.D_air_NH3[dd+1])
            # self.R_soilaq_down[dd+1] = (0.14)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            # self.R_soilg_down[dd+1] = (0.14)/(tor_soil_gas*self.D_air_NH3[dd+1])

            ## sourcelayer resistance
            ## resistance from the source layer to the above topsoil layer
            self.R_sourcelayer_aq[dd+1] = (injection_depth-p_topsoil)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            self.R_sourcelayer_gas[dd+1] = (injection_depth-p_topsoil)/(tor_soil_gas*self.D_air_NH3[dd+1])
            # self.R_sourcelayer_aq[dd+1] = ((injection_depth+z_sourcelayer/2)/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
            # self.R_sourcelayer_gas[dd+1] = ((injection_depth+z_sourcelayer/2)/2)/(tor_soil_gas*self.D_air_NH3[dd+1])

            '''## N input in multiple forms
            self.urea[dd+1] = self.urea_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.avail_N[dd+1] = self.avail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.resist_N[dd+1] = self.resist_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])
            self.unavail_N[dd+1] = self.unavail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])

            ## TAN production from urea hydrolysis and the N decomposition rate from dung
            self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                    self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                    self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
            ## TAN from housing to storage
            self.TAN[dd+1] = self.TAN_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])'''

            if chem_fert_type == 'urea':

                ## Urea pool of the source layer (deep injection layer)
                self.urea_pool[dd+1] = self.urea_pool[dd] + self.urea[dd+1] + self.urea_infil[dd]
                ## urea concentration of the source layer and at the surface
                self.urea_amount[dd+1] = self.urea_pool[dd+1]/(z_sourcelayer*self.soilmoist[dd+1])
                ## pathways
                urea_diffidx = self.urea_amount[dd+1]/self.R_soilaq_down[dd+1] * timestep*3600
                urea_diffupidx = (self.urea_amount[dd+1]-self.urea_soil_amount[dd])/self.R_sourcelayer_aq[dd+1]*timestep*3600
                urea_leachingidx = self.qpsoil[dd+1] * self.urea_amount[dd+1]*timestep*3600
                urea_toTANidx = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd+1]
                ## determining the fluxes: washoff, diffusion, leaching, hydrolysis
                all_loss =  urea_diffidx + urea_diffupidx + urea_leachingidx + urea_toTANidx
                loss_idx = self.urea_pool[dd+1] - all_loss
                self.urea_diff[dd+1][loss_idx>=0] = urea_diffidx[loss_idx>=0]
                self.ureadiffusivefluxsourcelayer[dd+1][loss_idx>0] = urea_diffupidx[loss_idx>=0]
                self.urealeaching[dd+1][loss_idx>=0] = urea_leachingidx[loss_idx>=0]
                self.TAN_prod[dd+1][loss_idx>=0] = urea_toTANidx[loss_idx>=0]
                self.urea_diff[dd+1][loss_idx<0] = urea_diffidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.ureadiffusivefluxsourcelayer[dd+1][loss_idx<0] = urea_diffupidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.urealeaching[dd+1][loss_idx<0]= urea_leachingidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                self.TAN_prod[dd+1][loss_idx<0] = urea_toTANidx[loss_idx<0]/\
                                                        all_loss[loss_idx<0]*self.urea_pool[dd+1][loss_idx<0]
                ## determining urea pool after subtracting all loss pathways
                self.urea_pool[dd+1] = self.urea_pool[dd+1] - self.TAN_prod[dd+1] - self.urea_diff[dd+1] - \
                            self.ureadiffusivefluxsourcelayer[dd+1] - self.urealeaching[dd+1]
                ## getting rid of rounding error
                self.urea_pool[dd+1][self.urea_pool[dd+1]<0.0] = 0.0
                self.urea_amount[dd+1] = self.urea_pool[dd+1]/(z_sourcelayer*self.soilmoist[dd+1])

                ## urea pool of the topsoil layer
                self.urea_pool_soil[dd+1] = self.urea_pool_soil[dd] + self.ureadiffusivefluxsourcelayer[dd+1]
                ## urea concentration of the topsoil layer
                self.urea_soil_amount[dd+1] = self.urea_pool_soil[dd+1]/(z_topsoil*self.soilmoist[dd+1])
                urea_surf_amount = self.urea_soil_amount[dd+1]/(self.R_soilaq[dd+1]*self.rain_avail_washoff[dd+1]+1)
                ## pathways
                urea_washoffidx = self.rain_avail_washoff[dd+1] * urea_surf_amount*timestep*3600
                soilurea_infilidx = self.qpsoil[dd+1]*self.urea_soil_amount[dd+1]*timestep*3600
                soilurea_toTANidx = self.daily_urea_hydro_rate[dd+1]*self.urea_pool_soil[dd+1]
                ## determining the fluxes: washoff, diffusion, leaching, hydrolysis
                soilall_loss = urea_washoffidx + soilurea_infilidx + soilurea_toTANidx
                soilloss_idx = self.urea_pool_soil[dd+1] - soilall_loss
                self.urea_washoff[dd+1][soilloss_idx>=0] = urea_washoffidx[soilloss_idx>=0]
                self.urea_infil[dd+1][soilloss_idx>=0] = soilurea_infilidx[soilloss_idx>=0]
                self.TAN_prod_soil[dd+1][soilloss_idx>=0] = soilurea_toTANidx[soilloss_idx>=0]
                self.urea_washoff[dd+1][soilloss_idx<0] = urea_washoffidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                self.urea_infil[dd+1][soilloss_idx<0]= soilurea_infilidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                self.TAN_prod_soil[dd+1][soilloss_idx<0] = soilurea_toTANidx[soilloss_idx<0]/\
                                            soilall_loss[soilloss_idx<0]*self.urea_pool_soil[dd+1][soilloss_idx<0]
                ## deterining urea pool of the topsoil layer after subtracting all losses
                self.urea_pool_soil[dd+1] = self.urea_pool_soil[dd+1] - self.urea_washoff[dd+1] - self.urea_infil[dd+1] - \
                                            self.TAN_prod_soil[dd+1]
                self.urea_pool_soil[dd+1][self.urea_pool_soil[dd+1]<0.0] = 0.0
                self.urea_soil_amount[dd+1] = self.urea_pool_soil[dd+1]/(z_topsoil*self.soilmoist[dd+1])

            '''## Org N pools in various forms
            self.avail_N_pool[dd+1] = (self.avail_N_pool[dd] - self.avail_N_washoff[dd+1])*(1 - self.daily_Na_decomp_rate[dd+1]) + \
                self.avail_N[dd+1] 
            self.resist_N_pool[dd+1] = (self.resist_N_pool[dd] - self.resist_N_washoff[dd+1])* (1 - self.daily_Nr_decomp_rate[dd+1]) + \
                self.resist_N[dd+1] 
            self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] - self.unavail_N_washoff[dd+1]'''

            ## TAN pool (source layer) - beneath the topsoil layer
            ## Note: source of TAN pool: 1) input of TAN from fertilizer (ammonium or urea hydrolysis)
            ##                              nitrate does not contribute to the TAN pool
            ##                           2) diffusion (aq+gas) from the above topsoil layer to the source layer (bidirectional)                                
            ##       loss of TAN pool: 1) diffusion (aq+gas) to the above topsoil layer, 
            ##                         2) diffusion (aq+gas) to the deep soil, 
            ##                         3) subsurface leaching to deep soil, 
            ##                         4) nitrification, 5) plant uptake
            ##       Note that plant N uptake exists in this scheme, compared to the BC-topdressed scheme
            '''TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
            self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
            self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]
            ## TAN pool in ug
            TAN_pool_ug = self.TAN_pool[dd+1] * 1e6'''

            ## TAN pool with input
            self.TAN_pool[dd+1] = self.TAN_pool[dd]+self.TAN_prod[dd+1]+self.TAN[dd+1] + \
                self.diffusivefluxdown_aq[dd] + self.diffusivefluxdown_gas[dd]

            ## TAN conc; g/m3
            KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
            self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## TAN molar conc; mol/L
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
            ## NH3 conc in the source layer
            self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
            ## NH3 conc is 0 when sourcelayer water content equals to porosity
            self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## fluxes to deeper soil
            ## TAN loss through aqueous diffusion and leaching to deeper soil
            sourcelayerdiffaq_idx = self.TAN_amount[dd+1]/(self.R_soilaq_down[dd+1]*timestep*3600)
            sourcelayerdiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
            sourcelayerdiffgas_idx = self.NH3_gas_bulk[dd+1]/(self.R_soilg_down[dd+1]*timestep*3600)
            sourcelayerdiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
            
            ## diffusive fluxes to the above topsoil layer
            sourcelayerupdiffaq_idx = (self.TAN_amount[dd+1] - self.TAN_soil_amount[dd+1])/self.R_sourcelayer_aq[dd+1]*timestep*3600
            sourcelayerupdiffaq_idx[sourcelayerupdiffaq_idx<0] = 0.0 
            sourcelayerupdiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
            sourcelayerupdiffgas_idx = (self.NH3_gas_bulk[dd+1] - self.NH3_gas_soil[dd+1])/self.R_sourcelayer_gas[dd+1]*timestep*3600
            sourcelayerupdiffgas_idx[sourcelayerupdiffgas_idx<0] = 0.0 
            sourcelayerupdiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## leaching
            sourcelayerleaching_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600
            ## nirification rate of TAN in the source layer; daily maximum nitrification rate is 0.1 per day
            KNO3_soil= nitrification_rate_soil(ground_temp=self.T_sim[dd+1],theta=self.soilmoist[dd+1],
                                            theta_sat=self.persm[dd+1],pH=sim_pH[dd+1],fer_type='mineral')*timestep*3600
            KNO3_soil[KNO3_soil>0.1] = 0.1
            ## correction for WPFS response
            KNO3_soil[KNO3_soil<0.0] = 0.0
            KNO3_soil[np.isnan(KNO3_soil)] = 0.0
            f_NH4 = self.soilmoist[dd+1]/(self.soilmoist[dd+1]+\
                KNH3*(self.persm[dd+1]-self.soilmoist[dd+1])+\
                (1-self.persm[dd+1])*Kd)*(sim_ccH[dd+1]/(sim_ccH[dd+1]+self.k_NH4[dd+1]))
            sourcelayernitrif_idx =  KNO3_soil*self.TAN_pool[dd+1]*f_NH4
            sourcelayerammNuptake_idx, sourcelayernitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool[dd+1]*f_NH4,\
                                    Nnit=self.NO3_pool[dd],temp=self.T_sim[dd+1])
            ## no uptake after harvesting
            sourcelayerammNuptake_idx[harvestidx<(dd+1)] = 0.0
            sourcelayernitNuptake_idx[harvestidx<(dd+1)] = 0.0    

            ## fluxes from sourcelayer to the topsoil layer
            ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
            sourcelayerall_loss = sourcelayerdiffaq_idx + sourcelayerdiffgas_idx + sourcelayerupdiffaq_idx + \
                                sourcelayerupdiffgas_idx + sourcelayerleaching_idx + sourcelayernitrif_idx + \
                                sourcelayerammNuptake_idx 
            sourcelayerloss_idx = self.TAN_pool[dd+1] - sourcelayerall_loss

            self.diffusivefluxsoil_aq[dd+1][sourcelayerloss_idx>=0] = sourcelayerdiffaq_idx[sourcelayerloss_idx>=0]
            self.diffusivefluxsoil_gas[dd+1][sourcelayerloss_idx>=0] = sourcelayerdiffgas_idx[sourcelayerloss_idx>=0]
            self.diffusivefluxsourcelayer_aq[dd+1][sourcelayerloss_idx>=0] = sourcelayerupdiffaq_idx[sourcelayerloss_idx>=0]
            self.diffusivefluxsourcelayer_gas[dd+1][sourcelayerloss_idx>=0] = sourcelayerupdiffgas_idx[sourcelayerloss_idx>=0]
            self.leachingflux[dd+1][sourcelayerloss_idx>=0] = sourcelayerleaching_idx[sourcelayerloss_idx>=0]
            self.nitrif_NO3_sourcelayer[dd+1][sourcelayerloss_idx>=0] = sourcelayernitrif_idx[sourcelayerloss_idx>=0]
            self.ammN_uptake_sourcelayer[dd+1][sourcelayerloss_idx>=0] = sourcelayerammNuptake_idx[sourcelayerloss_idx>=0]

            self.diffusivefluxsoil_aq[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerdiffaq_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.diffusivefluxsoil_gas[dd+1][sourcelayerloss_idx<0] = self.TAN_pool_soil[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerdiffgas_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.diffusivefluxsourcelayer_aq[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerupdiffaq_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.diffusivefluxsourcelayer_gas[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerupdiffgas_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.leachingflux[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerleaching_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.nitrif_NO3_sourcelayer[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayernitrif_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]
            self.ammN_uptake_sourcelayer[dd+1][sourcelayerloss_idx<0] = self.TAN_pool[dd+1][sourcelayerloss_idx<0]*\
                                        sourcelayerammNuptake_idx[sourcelayerloss_idx<0]/sourcelayerall_loss[sourcelayerloss_idx<0]

            self.TAN_pool[dd+1] = self.TAN_pool[dd+1] - self.diffusivefluxsoil_aq[dd+1] - self.diffusivefluxsoil_gas[dd+1] -\
                                self.diffusivefluxsourcelayer_aq[dd+1] - self.diffusivefluxsourcelayer_gas[dd+1] - self.leachingflux[dd+1] -\
                                self.nitrif_NO3_sourcelayer[dd+1] - self.ammN_uptake_sourcelayer[dd+1]
            ## get rid of rounding error (TAN pool occasionaly become <-1e-13)
            self.TAN_pool[dd+1][self.TAN_pool[dd+1]<0.0] = 0.0

            self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 conc of the source layer
            self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
            ## NH3 conc is 0 when source layer water content equals to porosity
            self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## NO3- pool of source layer
            ## sources: 1) input of nitrate fertilizer
            ##          2) nitrification, 3) infiltration from the above topsoil layer
            ## losses: 1) leaching to the deep soil, 2) diffusion (aq only) to the deep soil
            ##         3) diffusion (aq only) to the above topsoil layer, 4) plant uptake
            NO3_idx = self.NO3_pool[dd] - self.NO3_leaching[dd] - self.NO3_diffusiveup[dd] - \
                        self.NO3_diffusivedeep[dd] - self.nitN_uptake_sourcelayer[dd]
            self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + \
                        self.NO3[dd+1][NO3_idx>0] + self.NO3_infilsoil[dd][NO3_idx>0]
            self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0] +\
                        self.NO3_infilsoil[dd][NO3_idx<=0]

            ## NO3- conc of the source layer; g/m3
            self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                        (z_sourcelayer*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
            self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0

            ## diffusive aquous NO3 and infiltration of NO3 from the source layer to the deep soil
            NO3_diffdownidx = (self.NO3_amount[dd+1])*(f_DNO3/self.R_soilaq_down[dd+1])*timestep*3600
            NO3_diffdownidx[NO3_diffdownidx<0] = 0.0
            NO3_diffdownidx[self.soilmoist[dd+1]==0] = 0.0
            NO3_diffupidx = (self.NO3_amount[dd+1] - self.NO3_soil_amount[dd])*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
            NO3_diffupidx[NO3_diffupidx<0] = 0.0
            NO3_diffupidx[self.soilmoist[dd+1]==0] = 0.0
            NO3_leachingidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600

            NO3_lossall = NO3_diffdownidx + NO3_diffupidx + NO3_leachingidx + sourcelayernitNuptake_idx
            sourcelayerloss_idx = self.NO3_pool[dd+1] - NO3_lossall
            self.NO3_diffusivedeep[dd+1][sourcelayerloss_idx>=0] = NO3_diffdownidx[sourcelayerloss_idx>=0]
            self.NO3_diffusiveup[dd+1][sourcelayerloss_idx>=0] = NO3_diffupidx[sourcelayerloss_idx>=0]
            self.NO3_leaching[dd+1][sourcelayerloss_idx>=0] = NO3_leachingidx[sourcelayerloss_idx>=0]
            self.nitN_uptake_sourcelayer[dd+1][sourcelayerloss_idx>=0] = sourcelayernitNuptake_idx[sourcelayerloss_idx>=0]
            self.NO3_diffusivedeep[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                            NO3_diffdownidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]   
            self.NO3_diffusiveup[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                            NO3_diffupidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                        
            self.NO3_leaching[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                            NO3_leachingidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]
            self.nitN_uptake_sourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                            sourcelayernitNuptake_idx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]

            ## TAN pool of the topsoil layer
            ## sources: diffusion (aq+gas) from the underlying source layer
            ## losses: 1) NH3 volatilization, 2) surface runoff
            ##         3) infiltration to the underlying source layer
            ##         4) diffusion (aq+gas) to the underlying source layer (bi-directional)
            ##         5) nitrification, 6) plant N uptake
            ## TAN pool with inputs
            self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd] + self.diffusivefluxsourcelayer_aq[dd+1] +\
                                    self.diffusivefluxsourcelayer_gas[dd+1] + self.TAN_prod_soil[dd+1]

            ## TAN conc of the topsoil layer; g/m3
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                ((injection_depth-z_sourcelayer/2)*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 concentration in the soil pore space
            self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
            ## NH3 conc is zero when soil moisture content reaches the saturation
            self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## TAN conc at the surface (solved); atmospheric NH3 is ignored
            self.TAN_surf_amount_M[dd+1] = (self.TAN_soil_amount[dd+1]*(1/self.R_soilaq[dd+1]+KNH3/self.R_soilg[dd+1])/\
                                    (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_soilg[dd+1])+1/self.R_soilaq[dd+1]))/\
                                        (14*1000)
            ## TAN surface conc in g/m3
            TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
            ## Gaseous NH3 at the surface
            self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
            ## in g
            NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

            ## determining the maximum emission; 
            emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
            ## determining the maximum TAN runoff;
            ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
            ## else need to multiply by timestep*3600
            runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
            ## determining the maximum TAN aqueous diffusion to the underlying source layer;
            ## diffusion is considered to be bidirectional
            diffaq_idx = (self.TAN_soil_amount[dd+1]-self.TAN_amount[dd+1])/self.R_sourcelayer_aq[dd+1]*timestep*3600
            diffaq_idx[diffaq_idx<=0] = 0.0
            ## when soil mositure is 0, aqueous diffusion stops
            diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
            ## determining the maximum TAN gaseous diffusion to the underlying source layer;
            diffgas_idx = (self.NH3_gas_soil[dd+1]-self.NH3_gas_bulk[dd+1])/self.R_sourcelayer_gas[dd+1]*timestep*3600
            diffgas_idx[diffgas_idx<0] = 0.0
            ## when the soil moisture reaches the saturation, gaseous diffusion stops
            diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
            ## determining the maximum infiltration of TAN
            infil_idx = self.qpsoil[dd+1]*self.TAN_soil_amount[dd+1]*timestep*3600
            ## nirification rate of TAN in the source layer; daily maximum nitrification rate is 0.1 per day
            nitrif_idx = KNO3_soil*self.TAN_pool_soil[dd+1]*f_NH4
            ## plant N uptake
            soilammNuptake_idx, soilnitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool_soil[dd+1]*f_NH4,\
                                    Nnit=self.NO3_pool_soil[dd],temp=self.T_sim[dd+1])
            ## no uptake after harvesting
            soilammNuptake_idx[harvestidx<(dd+1)] = 0.0
            soilnitNuptake_idx[harvestidx<(dd+1)] = 0.0

            ## fluxes from sourcelayer to the topsoil layer
            ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
            soilall_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + infil_idx + nitrif_idx + soilammNuptake_idx
            soilloss_idx = self.TAN_pool_soil[dd+1] - soilall_loss
            
            self.NH3_flux[dd+1][soilloss_idx>=0] = emiss_idx[soilloss_idx>=0]
            self.TAN_washoff[dd+1][soilloss_idx>=0] = runoff_idx[soilloss_idx>=0]
            self.diffusivefluxdown_aq[dd+1][soilloss_idx>=0] = diffaq_idx[soilloss_idx>=0]
            self.diffusivefluxdown_gas[dd+1][soilloss_idx>=0] = diffgas_idx[soilloss_idx>=0]
            self.infilflux[dd+1][soilloss_idx>=0] = infil_idx[soilloss_idx>=0]
            self.nitrif_NO3_soil[dd+1][soilloss_idx>=0] = nitrif_idx[soilloss_idx>=0]
            self.ammN_uptake[dd+1][soilloss_idx>=0] = soilammNuptake_idx[soilloss_idx>=0]
            ## otherwise, we use a weighted distribution
            self.NH3_flux[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (emiss_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.TAN_washoff[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (runoff_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.diffusivefluxdown_aq[dd+1][soilloss_idx<0]= self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (diffaq_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.diffusivefluxdown_gas[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (diffgas_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.infilflux[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (infil_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.nitrif_NO3_soil[dd+1][soilloss_idx<0]= self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (nitrif_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])
            self.ammN_uptake[dd+1][soilloss_idx<0]= self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                        (soilammNuptake_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0])                             

            ## determine TAN pool after subtracting all loss pathways
            self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd+1] - self.NH3_flux[dd+1] - self.diffusivefluxdown_aq[dd+1] -\
                    self.diffusivefluxdown_gas[dd+1]- self.infilflux[dd+1] - \
                        self.nitrif_NO3_soil[dd+1] - self.TAN_washoff[dd+1] - self.ammN_uptake[dd+1]
            ## get rid of rounding error (TAN pool occasionaly become <-1e-13)
            self.TAN_pool_soil[dd+1][self.TAN_pool_soil[dd+1]<0.0] = 0.0
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0
            self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                ((injection_depth-z_sourcelayer/2)*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
            ## NH3 concentration in the soil pore space
            self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
            ## NH3 conc is zero when soil moisture content reaches the saturation
            self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

            ## NO3- pool of topsoil 
            ## sources: 1) nitrification, 2) diffusion (aq only) from the underlying source layer
            ## losses: 1) infiltration to the underlying source layer, 2) plant uptake, 3) surface washoff
            NO3_soil_idx = self.NO3_pool_soil[dd] - self.NO3_infilsoil[dd] - self.nitN_uptake[dd] - self.NO3_washoff[dd]
            self.NO3_pool_soil[dd+1][NO3_soil_idx>=0] = NO3_soil_idx[NO3_soil_idx>=0] + self.nitrif_NO3_soil[dd+1][NO3_soil_idx>=0] +\
                self.NO3_diffusiveup[dd+1][NO3_soil_idx>=0]
            self.NO3_pool_soil[dd+1][NO3_soil_idx<0] = self.nitrif_NO3_soil[dd+1][NO3_soil_idx<0] + \
                self.NO3_diffusiveup[dd+1][NO3_soil_idx<0]

            ## NO3- conc of the topsoil layer; g/m3
            self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                        ((injection_depth-z_sourcelayer/2)*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])   
            self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
            
            ## NO3 loss through aqueous diffusion and plant uptake
            NO3_soilinfilidx = self.qpsoil[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
            ## NO3 surface runoff
            NO3_washoffidx = self.rain_avail_washoff[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
            
            NO3_soilall_loss = NO3_soilinfilidx + NO3_washoffidx +soilnitNuptake_idx
            soilloss_idx = self.NO3_pool_soil[dd+1] - NO3_soilall_loss
            self.NO3_infilsoil[dd+1][soilloss_idx>=0] = NO3_soilinfilidx[soilloss_idx>=0]
            self.NO3_washoff[dd+1][soilloss_idx>=0] = NO3_washoffidx[soilloss_idx>=0]
            self.nitN_uptake[dd+1][soilloss_idx>=0] = soilnitNuptake_idx[soilloss_idx>=0]
            self.NO3_infilsoil[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                        NO3_soilinfilidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
            self.NO3_washoff[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                        NO3_washoffidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
            self.nitN_uptake[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                                soilnitNuptake_idx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]

        return

    def land_sim_reshape(self,sim_result):
        shape = sim_result.shape
        dim1 = int((shape[0]-1)/2+1)
        output = np.zeros([dim1,shape[1],shape[2]])
        output = sim_result[1:dim1]+sim_result[dim1:]
        return output

    def N_stat(self,crop_item,fert_method,chem_fert_type,ncfile_o=False):
        
        if chem_fert_type == 'ammonium':
            sim_area = self.ammN_area
            chemfert_Ntotal = self.land_sim_reshape(self.TAN_added)*sim_area
            chemfert_NH3emiss = self.land_sim_reshape(self.NH3_flux)*sim_area
            # chemfert_Ntotal = self.TAN_added*sim_area
            # chemfert_NH3emiss = self.NH3_flux*sim_area
            #chemfert_ammN = np.nansum(chemfert_Ntotal,axis=0)
            #chemfert_ammNH3 = np.nansum(chemfert_NH3emiss,axis=0)

        elif chem_fert_type == 'urea':
            sim_area = self.ureaN_area
            chemfert_Ntotal = self.land_sim_reshape(self.urea_added)*sim_area
            chemfert_NH3emiss = self.land_sim_reshape(self.NH3_flux)*sim_area
            # chemfert_Ntotal = self.urea_added*sim_area
            # chemfert_NH3emiss = self.NH3_flux*sim_area
            #chemfert_ureaN = np.nansum(chemfert_Ntotal,axis=0)
            #chemfert_ureaNH3 = np.nansum(chemfert_NH3emiss,axis=0)
            
        elif chem_fert_type == 'nitrate':
            sim_area = self.nitN_area
            chemfert_Ntotal = self.land_sim_reshape(self.NO3_added)*sim_area
            chemfert_NH3emiss = self.land_sim_reshape(self.NH3_flux)*sim_area
            # chemfert_Ntotal = self.NO3_added*sim_area
            # chemfert_NH3emiss = self.NH3_flux*sim_area
            #chemfert_nitN = np.nansum(chemfert_Ntotal,axis=0)
            #chemfert_nitNH3 = np.nansum(chemfert_NH3emiss,axis=0)
        
        print('Total N applied: '+str(sum_totalGg(chemfert_Ntotal)))
        print('NH3 emission: '+str(sum_totalGg(chemfert_NH3emiss)))

        chemfert_TANwashoff = self.land_sim_reshape(self.TAN_washoff+self.urea_washoff)*sim_area
        print('TAN washoff: '+ str(sum_totalGg(chemfert_TANwashoff))+' Gg')
        chemfert_leaching = self.land_sim_reshape(self.leachingflux+self.urealeaching)*sim_area
        print('NH4 leaching: '+ str(sum_totalGg(chemfert_leaching))+' Gg')

        if fert_method == 'broadcasting-surf':
            chemfert_diffaq = self.land_sim_reshape(self.diffusivefluxsourcelayer_aq+self.urea_diff)*sim_area
            print('TAN diff aq to topsoil: '+ str(sum_totalGg(chemfert_diffaq))+' Gg')
            chemfert_diffgas = self.land_sim_reshape(self.diffusivefluxsourcelayer_gas)*sim_area
            print('TAN diff gas to topsoil: '+ str(sum_totalGg(chemfert_diffgas))+' Gg')
            chemfert_updiffaq = self.land_sim_reshape(self.diffusivefluxup_aq)*sim_area
            print('TAN diff aq upwards to source layer: '+ str(sum_totalGg(chemfert_updiffaq))+' Gg')
            chemfert_updiffgas = self.land_sim_reshape(self.diffusivefluxup_gas)*sim_area
            print('TAN diff gas upwards to source layer: '+ str(sum_totalGg(chemfert_updiffgas))+' Gg')
            chemfert_infil = self.land_sim_reshape(self.infilflux)*sim_area
            print('NH4 infiltration to topsoil: '+ str(sum_totalGg(chemfert_infil))+' Gg')
            chemfert_nitrif = self.land_sim_reshape(self.nitrif_NO3_sourcelayer)*sim_area
            print('NH4 nitrification (surf source layer): '+ str(sum_totalGg(chemfert_nitrif))+' Gg')
            chemfert_diffaqdown = self.land_sim_reshape(self.diffusivefluxsoil_aq)*sim_area
            print('TAN diff aq to deeper soil: '+ str(sum_totalGg(chemfert_diffaqdown))+' Gg')
            chemfert_diffgasdown = self.land_sim_reshape(self.diffusivefluxsoil_gas)*sim_area
            print('TAN diff gas to deeper soil: '+ str(sum_totalGg(chemfert_diffgasdown))+' Gg')
            chemfert_nitrif_soil = self.land_sim_reshape(self.nitrif_NO3_soil)*sim_area
            print('NH4 nitrif (topsoil layer): '+ str(sum_totalGg(chemfert_nitrif_soil))+' Gg')
            chemfert_ammN_uptake = self.land_sim_reshape(self.ammN_uptake)*sim_area
            print('NH4 uptake by plants: '+ str(sum_totalGg(chemfert_ammN_uptake))+' Gg')
            chemfert_nitN_uptake = self.land_sim_reshape(self.nitN_uptake)*sim_area
            print('NO3 uptake by plants: '+ str(sum_totalGg(chemfert_nitN_uptake))+' Gg')
            chemfert_NO3diffusive = self.land_sim_reshape(self.NO3_diffusivesourcelayer)*sim_area
            print('NO3 diffusion from source layer to soil pool: '+ str(sum_totalGg(chemfert_NO3diffusive))+' Gg')
            chemfert_NO3leaching = self.land_sim_reshape(self.NO3_leaching)*sim_area
            print('NO3 leaching: '+ str(sum_totalGg(chemfert_NO3leaching))+' Gg')
            chemfert_NO3diffusionsoil = self.land_sim_reshape(self.NO3_diffusivedeep)*sim_area
            print('NO3 diffusion to deeper soil: '+ str(sum_totalGg(chemfert_NO3diffusionsoil))+' Gg')
                
        elif fert_method == 'broadcasting-disk':

            chemfert_nitrif = self.land_sim_reshape(self.nitrif_NO3_sourcelayer)*sim_area
            print('NH4 nitrification: '+ str(sum_totalGg(chemfert_nitrif))+' Gg')
            chemfert_diffaqdown = self.land_sim_reshape(self.diffusivefluxsoil_aq+self.urea_diff)*sim_area
            print('TAN diff aq to deeper soil: '+ str(sum_totalGg(chemfert_diffaqdown))+' Gg')
            chemfert_diffgasdown = self.land_sim_reshape(self.diffusivefluxsoil_gas)*sim_area
            print('TAN diff gas to deeper soil: '+ str(sum_totalGg(chemfert_diffgasdown))+' Gg')
            chemfert_ammN_uptake = self.land_sim_reshape(self.ammN_uptake)*sim_area
            print('NH4 uptake by plants: '+ str(sum_totalGg(chemfert_ammN_uptake))+' Gg')
            chemfert_nitN_uptake = self.land_sim_reshape(self.nitN_uptake)*sim_area
            print('NO3 uptake by plants: '+ str(sum_totalGg(chemfert_nitN_uptake))+' Gg')
            chemfert_NO3leaching = self.land_sim_reshape(self.NO3_leaching)*sim_area
            print('NO3 leaching: '+ str(sum_totalGg(chemfert_NO3leaching))+' Gg')
            chemfert_NO3diffusionsoil = self.land_sim_reshape(self.NO3_diffusivedeep)*sim_area
            print('NO3 diffusion to deeper soil: '+ str(sum_totalGg(chemfert_NO3diffusionsoil))+' Gg')

        elif fert_method == 'deep injection':

            chemfert_nitrif = self.land_sim_reshape(self.nitrif_NO3_sourcelayer)*sim_area
            print('NH4 nitrification (source layer): '+ str(sum_totalGg(chemfert_nitrif))+' Gg')
            chemfert_diffaqdown = self.land_sim_reshape(self.diffusivefluxsoil_aq+self.urea_diff)*sim_area
            print('TAN diff aq to deeper soil: '+ str(sum_totalGg(chemfert_diffaqdown))+' Gg')
            chemfert_diffgasdown = self.land_sim_reshape(self.diffusivefluxsoil_gas)*sim_area
            print('TAN diff gas to deeper soil: '+ str(sum_totalGg(chemfert_diffgasdown))+' Gg')
            chemfert_nitrif_soil = self.land_sim_reshape(self.nitrif_NO3_soil)*sim_area
            print('NH4 nitrif (topsoil layer): '+ str(sum_totalGg(chemfert_nitrif_soil))+' Gg')
            chemfert_ammN_uptakesl = self.land_sim_reshape(self.ammN_uptake_sourcelayer)*sim_area
            print('NH4 uptake by plants (source layer): '+ str(sum_totalGg(chemfert_ammN_uptakesl))+' Gg')
            chemfert_ammN_uptake = self.land_sim_reshape(self.ammN_uptake)*sim_area
            print('NH4 uptake by plants (topsoil layer): '+ str(sum_totalGg(chemfert_ammN_uptake))+' Gg')
            chemfert_nitN_uptakesl = self.land_sim_reshape(self.nitN_uptake_sourcelayer)*sim_area
            print('NO3 uptake by plants (source layer): '+ str(sum_totalGg(chemfert_nitN_uptakesl))+' Gg')
            chemfert_nitN_uptake = self.land_sim_reshape(self.nitN_uptake)*sim_area
            print('NO3 uptake by plants (topsoil layer): '+ str(sum_totalGg(chemfert_nitN_uptake))+' Gg')
            chemfert_NO3leaching = self.land_sim_reshape(self.NO3_leaching)*sim_area
            print('NO3 leaching: '+ str(sum_totalGg(chemfert_NO3leaching))+' Gg')
            chemfert_NO3diffusionsoil = self.land_sim_reshape(self.NO3_diffusivedeep)*sim_area
            print('NO3 diffusion to deeper soil: '+ str(sum_totalGg(chemfert_NO3diffusionsoil))+' Gg')

            chemfert_ammN_uptake = chemfert_ammN_uptake + chemfert_ammN_uptakesl
            chemfert_nitN_uptake = chemfert_nitN_uptake + chemfert_nitN_uptakesl
        
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
                    TANdiffaq=(['lat','lon'],annual_total(chemfert_diffaqdown)),
                    NH3diffgas=(['lat','lon'],annual_total(chemfert_diffgasdown)),
                    NH4nitrif=(['lat','lon'],annual_total(chemfert_nitrif)),
                    NH4uptake=(['lat','lon'],annual_total(chemfert_ammN_uptake)),
                    NO3uptake=(['lat','lon'],annual_total(chemfert_nitN_uptake)),
                    NO3leaching=(['lat','lon'],annual_total(chemfert_NO3leaching)),
                    NO3diff=(['lat','lon'],annual_total(chemfert_NO3diffusionsoil)),
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
            outds.NO3leaching.attrs["unit"] = 'gN/year'
            outds.NO3leaching.attrs["long name"] = 'Nitrate leaching from fertilizer application'
            outds.NO3diff.attrs["unit"] = 'gN/year'
            outds.NO3diff.attrs["long name"] = 'Nitrate diffusion to deep soil from fertilizer application'

            comp = dict(zlib=True, complevel=9)
            encoding = {var: comp for var in outds.data_vars}

            outds.to_netcdf(output_path+str(crop_item)+'.'+\
                str(chem_fert_type)+'.'+str(fert_method)+'.'+str(sim_year)+'.nc',encoding=encoding)
        return

    def chem_fert_main(self,fert_method,crop_item,chem_fert_type,start_day_idx,end_day_idx,fert_depth,sim_stat=False,ncfile_o=False):
        
        self.sim_env()
        
        if fert_method == 'broadcasting-surf':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_bcsurf_sim(start_day_idx,end_day_idx,chem_fert_type,crop=crop_item)              
        elif fert_method == 'broadcasting-disk':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_bcdisk_sim(start_day_idx,end_day_idx,chem_fert_type,fert_depth,crop=crop_item)
        elif fert_method == 'deep injection':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_deepinjec_sim(start_day_idx,end_day_idx,chem_fert_type,fert_depth,crop=crop_item)
        else:
            print('Error: fertilization method')

        if sim_stat is True:
            ## output ncfile
            if ncfile_o is True:
                self.N_stat(crop_item,fert_method,chem_fert_type,ncfile_o=True)
            else:
                self.N_stat(crop_item,fert_method,chem_fert_type,ncfile_o=False)

        return

