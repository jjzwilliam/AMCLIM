from logging import raiseExceptions
from INPUT.input import *
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *
import sys

###################################
## Crop info
###################################
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
## assuming the roughness height of xxx is ~ xxxm (<ref height of 2m)
# zo_xxx = xxx
## assuming the roughness height of manure pile (open land) is ~ 1.0m (<ref height of 2m)
zo_manureland = 1.0
## manure surface thickness is set to be 2 cm
z_manuresurf = 0.02
## dry matter (DM) content of solid manure 
DM_content = solid_m_DM[livestock]
## dry matter (DM) content of liquid manure is assumed to be 5%
f_DM_liquid = 0.1
## maximum water content of manure
f_wcmax = 1 - (DM_content/100)/2
## assuming the density of manure; 1t kg/m^3 or 1g/cm^3
manure_density = rho_m[livestock]
## thickness of the source layer: 2cm; mid point 1cm
z_sourcelayer = 0.02
p_sourcelayer = 0.01
## thickness of the topsoil layer: 7cm; mide point 3.5cm
z_topsoil = 0.07
p_topsoil = 0.035
## thickness of the 2nd soil layer (under topsoil); mid point 17.5cm
z_2ndsoil = 0.21
p_2ndsoil = 0.175
## assuming soil interface (source layer) of 4mm
# z_soil = 0.004
z_soil = 0.1 ## test
## distance to deeper soil; 2cm 
# d_deepsoil = 0.02
d_deepsoil = 0.1 ## test
## assuming infiltration of manure water to the soil is 10mm/day (1cm/day; 10 000 g/m^2/day) ref: Vira et al.,2020 GMD (2x d0)
# dailymaxinfil = 0.01
dailymaxinfil = 0.003 ## test: 3mm /day
## assuming the water holding capacity of manure is 3.0 g water/ g manure
absorb_factor = 3.0
# dailymaxinfil = 0.0 ## test: shut down infiltration
## infiltration flux within manure (m/s)
qpsoil_manure = (dailymaxinfil/1e6)/(24*3600)
## assuming soil characteristics: 1) sand (%), 2) clay (%), 3) bulk density (g/cm^3), 4) particle density (g/cm^3)
soil_sand = 80
soil_clay = 8
soil_bd = 1.5
soil_pd = 2.66
## washoff coefficients: 0.1%/mm water for N species, and 0.05%/mm water for non N species (manure)
f_washoff_nonN = 0.0005
f_washoff_N = 0.001
## anoxic fraction of manure; average of swine (75%) and cow (50%); ref: Wang et al., Sci.Report.2015
f_manure_anoxic = 0.6
## lagoon TAN conc; g/mL
lagoon_TAN_conc = 630.0e-6

class LAND_module:
    def __init__(self,array_shape,fert_type,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,NO3_added,water_added,pH_value,cropping_area,application_method_index):
        
        print('LAND Module - current fertilizer application is: '+str(fert_type))
        if fert_type == 'manure':
            ## feces input from MMS
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
            ## water amount in fresh feces
            # self.manure_initwc = np.zeros(array_shape)
            ## urine input from housing
            # self.urine_added = np.zeros(array_shape)
            # self.urine = np.zeros(array_shape)
            ## urea input from housing
        else:
            ## cropping area for nitrate N fertilizer
            self.nitN_area = np.zeros(array_shape[1:])
            ## cropping area for ammonium N fertilizer
            self.ammN_area = np.zeros(array_shape[1:])
            ## cropping area for urea N fertilizer
            self.ureaN_area = np.zeros(array_shape[1:])

        if application_method_index == 'deep injection':
            ## upward diffusion of aqueous TAN: from the top soil layer to the source layer 
            self.diffusivefluxdown_aq = np.zeros(array_shape)
            ## upward diffusion of gaseous NH3: from the top soil layer to the  source layer
            self.diffusivefluxdown_gas = np.zeros(array_shape)
            ## infiltration of NO3 from the top soil layer to the source layer; downwards
            self.NO3_infilsoil = np.zeros(array_shape)
            ## diffusive aqueous NO3- from the source layer to the top soil layer; upwards
            self.NO3_diffusiveup = np.zeros(array_shape)
            ## diffusive aqueous NO3- from the source layer to the deeper soil; downwards
            self.NO3_diffusivedown = np.zeros(array_shape)
            ##
            self.ammN_uptake_sourcelayer = np.zeros(array_shape)
            ##
            self.nitN_uptake_sourcelayer = np.zeros(array_shape)

        self.urea_added = urea_added
        self.urea = np.zeros(array_shape)
        ## urea pool
        self.urea_pool = np.zeros(array_shape)
        ## washoff flux of urea
        self.urea_washoff = np.zeros(array_shape)
        ## uric acid input from housing
        self.UA_added = UA_added
        self.UA = np.zeros(array_shape)
        ## uric acid pool
        self.UA_pool = np.zeros(array_shape)
        ## washoff flux of UA
        self.UA_washoff = np.zeros(array_shape)
        ## input of available N component from housing
        self.avail_N_added = avail_N_added
        self.avail_N = np.zeros(array_shape)
        ## Nitrogen pool that is available (easily) to form TAN
        self.avail_N_pool = np.zeros(array_shape)
        ## washoff flux of available org N 
        self.avail_N_washoff = np.zeros(array_shape)
        ## input of resistant N component from housing
        self.resist_N_added = resist_N_added
        self.resist_N = np.zeros(array_shape)
        ## Nitrogen pool that is resistant (slowly) to form TAN
        self.resist_N_pool = np.zeros(array_shape)
        ## washoff flux of resistant org N 
        self.resist_N_washoff = np.zeros(array_shape)
        ## input of unavailable N component from housing
        self.unavail_N_added = unavail_N_added
        self.unavail_N = np.zeros(array_shape)
        ## Nitrogen pool that is unavailable (cannot) to form TAN
        self.unavail_N_pool = np.zeros(array_shape)
        ## washoff flux of unavailable org N 
        self.unavail_N_washoff = np.zeros(array_shape)
        ## TAN added from housing
        self.TAN_added = TAN_added
        self.TAN = np.zeros(array_shape)
        ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool
        self.TAN_prod = np.zeros(array_shape)
        ## TAN pool (main soil layer)
        self.TAN_pool = np.zeros(array_shape)
        ## TAN pool of soil interface
        self.TAN_pool_soil = np.zeros(array_shape)
        ## washoff flux of TAN 
        self.TAN_washoff = np.zeros(array_shape)
        ## TAN pool in ug/m2
        # self.TAN_pool_ug = np.zeros(array_shape)
        ## TAN pool conc (aqueous phase)
        self.TAN_amount = np.zeros(array_shape)
        ## TAN pool in molar concentration
        self.TAN_amount_M = np.zeros(array_shape)
        ## TAN conc at the surface; for [MMS (barn,open) solid]
        self.TAN_surf_amount_M = np.zeros(array_shape)
        ## TAN conc at the soil surface/interface between manure and land (soil); for [MMS open solid]
        self.TAN_soil_amount = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3 = np.zeros(array_shape)
        ## NO3 from manure/mineral fertilizer
        self.NO3_added = NO3_added
        ## NO3- from nitrification in bulk manrue 
        self.nitrif_NO3_sourcelayer = np.zeros(array_shape)
        ## NO3- washoff 
        self.NO3_washoff = np.zeros(array_shape)
        ## NO3- pool in soil layer
        self.NO3_pool = np.zeros(array_shape)
        ## NO3- conc in soil layer
        self.NO3_amount = np.zeros(array_shape)
        ## NO3- from nitrification in soil interface
        self.nitrif_NO3_soil = np.zeros(array_shape)
        ## NO3- pool in soil interface
        self.NO3_pool_soil = np.zeros(array_shape)
        ## NO3- conc in soil interface
        self.NO3_soil_amount = np.zeros(array_shape)
        ## infiltration of NO3 from the source layer to the top soil layer
        self.NO3_infilsourcelayer = np.zeros(array_shape)
        ## diffusive aqueous NO3- from the source layer to the top soil layer
        self.NO3_diffusivesourcelayer = np.zeros(array_shape)
        ## NO3- leaching to deeper soil
        self.NO3_leaching = np.zeros(array_shape)
        ## diffusive aqueous NO3- to deeper soil
        self.NO3_diffusivesoil = np.zeros(array_shape)
        ## water added from housing
        self.water_added = water_added
        self.water = np.zeros(array_shape)
        ## total water pool of the system (manure water; urine+manure water+[washing water])
        self.Total_water_pool = np.zeros(array_shape)
        ## total pool of the soil interface
        self.water_pool_soil = np.zeros(array_shape)
        ## ratio of [NH4+]/[H+] in the system
        self.Gamma_manure = np.zeros(array_shape)
        ## surface NH3 concentrtion at equilirium (in molar mass)
        self.NH3_gas_M = np.zeros(array_shape)
        ## NH3 concentration in soil layer
        self.NH3_gas_bulk = np.zeros(array_shape)
        ## soil NH3 concentration
        self.NH3_gas_soil = np.zeros(array_shape)
        ## surface NH3 concentrtion at equilirium (in ug)
        # self.NH3_gas_ug = np.zeros(array_shape)
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
        ## diffusion of aqueous TAN from the source layer to the top soil layer
        self.diffusivefluxsourcelayer_aq = np.zeros(array_shape)
        ## diffusion of gaseous NH3 from the source layer to the top soil layer
        self.diffusivefluxsourcelayer_gas = np.zeros(array_shape)
        ## diffusion of aqueous TAN from the top soil layer to the deeper soil
        self.diffusivefluxsoil_aq = np.zeros(array_shape)
        ## diffusion of gaseous NH3 from the top soil layer to the deeper soil
        self.diffusivefluxsoil_gas = np.zeros(array_shape)
        ## upward diffusion of aqueous TAN: from the top soil layer to the source layer 
        self.diffusivefluxup_aq = np.zeros(array_shape)
        ## upward diffusion of gaseous NH3: from the top soil layer to the  source layer
        self.diffusivefluxup_gas = np.zeros(array_shape)

        ## infiltration of aqueous TAN to soil interface
        self.infilflux = np.zeros(array_shape)
        ## uptake of ammonium nitrogen by plants
        self.ammN_uptake = np.zeros(array_shape)
        ## uptake of nitrate nitrogen by plants
        self.nitN_uptake = np.zeros(array_shape)

        ## temp
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
        ## soil resistance of the source layer
        self.R_sourcelayer_aq = np.zeros(array_shape)
        ## soil resistance of the source layer
        self.R_sourcelayer_gas = np.zeros(array_shape)
        ## soil resistance for aqueous diffusion
        self.R_soilaq = np.zeros(array_shape)
        ## soil resistance for gaseous diffusion
        self.R_soilg = np.zeros(array_shape)
        ## soil resistance for downward aqueous diffusion
        self.R_soilaq_down = np.zeros(array_shape)
        ## soil resistance for downward gaseous diffusion
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
        ## tortuosity for diffusion
        # self.tor_soil = np.zeros(array_shape)
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
        self.T_sim = xr_to_np(self.T_sim)
        self.RH_sim = xr_to_np(self.RH_sim)
        self.u_sim = xr_to_np(self.u_sim)
        self.evap_sim = xr_to_np(self.evap_sim)
        self.soilmoist = xr_to_np(self.soilmoist)
        self.persm = xr_to_np(self.persm)
        self.persm[self.persm>1.0] = 0.9999
        self.rainfall = xr_to_np(self.rainfall)
        self.R_atm = xr_to_np(self.R_atm)
        # self.daily_KNO3 = nitrification_rate_soil(ground_temp=self.T_sim,theta=self.soilmoist,theta_sat=self.persm,
        #                                             fer_type="manure")*timestep*3600
        # self.daily_KNO3[self.daily_KNO3<0] = 0.0
        # self.daily_KNO3[np.isnan(self.daily_KNO3)] = 0.0
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
        # times = mtrx[0]
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
                    plt_time = plant_calendar[lat,lon]
                    har_time = harvest_calendar[lat,lon]
                    app_freq = fert_freq[lat,lon]
            #         print(lat,lon)
                    ## harvesting goes into the next year
                    if ~np.isnan(app_freq):
                        if ~np.isnan(plt_time):
                            if ~np.isnan(har_time):
                                if har_time<plt_time:
                                    har_time = har_time+365
                                if app_freq <=1.0:
                                    N_app[int(plt_time),lat,lon] = N_input[lat,lon]*app_freq
                                    N_app_mark[int(plt_time),lat,lon] = 1 
                                elif app_freq>1.0:
                                    tapp = np.floor(app_freq)
                                    app_int = int(abs(int(har_time)-int(plt_time))/(tapp+1))
                                    app_idx = np.arange(int(plt_time),int(har_time)-1,app_int)
                                    for idx in app_idx[:-1]:
                                        N_app[int(idx),lat,lon] = N_input[lat,lon]
                                        N_app_mark[int(idx),lat,lon] = 1
                                    N_app[int(app_idx[-1]),lat,lon] = (app_freq-tapp)*N_input[lat,lon]
                                    N_app_mark[int(app_idx[-1]),lat,lon] = 1
        
        ## determine soil pH after urea application
        if soil_pH is not None:
            # print('check pH over.')
            self.soil_pH = soil_pH_postapp(base_pH=soil_pH,app_timing_map=N_app_mark,fert_pH=8.5)
            self.soil_ccH = 10**(-self.soil_pH)
        return N_app

    def chem_fert_type(self):
        fertds = open_ds(file_path+crop_data_path+fertfilename)
        nitN = fertds.nitrate_fert.values
        ammN = fertds.ammonium_fert.values
        ureaN = fertds.urea_fert.values
        totalN = nitN + ammN
        fnitN = nitN/totalN
        fammN = (ammN-ureaN)/totalN
        fureaN = ureaN/totalN
        return fnitN, fammN, fureaN

    def chem_fert_input(self,crop):
        fertds = open_ds(file_path+crop_data_path+crop+cropfileformat)
        cropcalds = open_ds(file_path+crop_data_path+crop_filledcalendar+crop+crop_filledcalendarformat)
        soilpHds = open_ds(file_path+soil_data_path+soilpHfile)

        totalN = fertds.TotalN.values*1e3
        ## N application rate is interpolated;
        ##  
        rateN = fertds.Nrate.values*1e3/1e4
        croparea = fertds.croparea.values*1e4
        croparea[totalN!=0] = totalN[totalN!=0]/rateN[totalN!=0]
        app_freq = totalN/(rateN*croparea)
        ## the maximum application frequency in a year is 5
        app_freq[app_freq>5] = 5

        plantidx = cropcalds['plant.start'].values
        harvestidx = cropcalds['harvest.start'].values

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

    
    def chem_fert_bcdisk_sim(self,start_day_idx,end_day_idx,chem_fert_type,disk_depth):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [broadcasting - incorporated disk], fertilizer placement depth is: '+str(disk_depth*100)+' cm')
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
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):
                
                ## soil resistance
                ## soil resistance = thickness of the source layer (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                # tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                # tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                self.R_soilaq[dd+1] = (disk_depth/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg[dd+1] = (disk_depth/2)/(tor_soil_gas*self.D_air_NH3[dd+1])
                self.R_soilaq_down[dd+1] = (p_2ndsoil-disk_depth/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg_down[dd+1] = (p_2ndsoil-disk_depth/2)/(tor_soil_gas*self.D_air_NH3[dd+1])

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1])
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]
                
                ## TAN pool (different from [MMS barn (liquid,solid)])
                ## Note: source of TAN pool: TAN production from 1) urea hydrolysis, 2) decomposition of org N and 3) input of TAN from housing
                ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil (aq+gas), 
                ##                         3) infiltration (subsurface leaching) to soil, and 4) nitrification
                ##       only aqueous phase diffusion is considered, and gaseous diffusion is not considered in this study
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsoil_aq[dd] - self.diffusivefluxsoil_gas[dd]-\
                    self.leachingflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd] - self.ammN_uptake[dd]
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
                KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_soil*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])*Kd))
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
                                                theta_sat=self.persm[dd+1],pH=sim_pH[dd+1],fer_type='mineral')*timestep*3600
                KNO3_soil[KNO3_soil>0.1] = 0.1
                ## correction for WPFS response
                KNO3_soil[KNO3_soil<0.0] = 0.0
                KNO3_soil[np.isnan(KNO3_soil)] = 0.0
                ## determining the maximum nitrification of TAN
                # nitrif_idx = KNO3_manure*self.TAN_pool[dd]
                f_NH4 = self.soilmoist[dd+1]/(self.soilmoist[dd+1]+\
                    KNH3*(self.persm[dd+1]-self.soilmoist[dd+1])+\
                    (1-self.persm[dd+1])*Kd)*(sim_ccH[dd+1]/(sim_ccH[dd+1]+self.k_NH4[dd+1]))
                nitrif_idx = KNO3_soil*self.TAN_pool[dd]*f_NH4
                soilammNuptake_idx, soilnitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool[dd]*f_NH4,\
                                        Nnit=self.NO3_pool[dd],temp=self.T_sim[dd+1])

                ## fluxes from bulk soil to deeper soil and chemical loss
                ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                all_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + leaching_idx + nitrif_idx + soilammNuptake_idx
                # if dd+1 == 143:
                #     print("manureall_loss: ",manureall_loss[130,363]) 
                loss_idx = self.TAN_pool[dd+1] - all_loss
                # if dd+1 == 143:
                #     print("manureloss_idx: ",manureloss_idx[130,363]) 
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
                NO3_idx = self.NO3_pool[dd] - self.NO3_leaching[dd] - self.NO3_diffusivesoil[dd] - self.NO3_washoff[dd+1] - \
                            self.nitN_uptake[dd]
                self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + \
                                                    self.NO3[dd+1][NO3_idx>0]
                self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0]

                ## NO3- conc of soil layer; g/m3
                self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                            (z_soil*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
                self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                
                ## diffusive aquous NO3 and infiltration of NO3 from bulk manure to soil interface
                NO3_diffidx = self.NO3_amount[dd+1]*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
                NO3_diffidx[NO3_diffidx<0] = 0.0
                NO3_diffidx[self.soilmoist[dd+1]==0] = 0.0
                NO3_leachingidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
                NO3_lossall = NO3_diffidx + NO3_leachingidx + soilnitNuptake_idx
                loss_idx = self.NO3_pool[dd+1] - NO3_lossall 
                self.NO3_diffusivesoil[dd+1][loss_idx>=0] = NO3_diffidx[loss_idx>=0]
                self.NO3_leaching[dd+1][loss_idx>=0] = NO3_leachingidx[loss_idx>=0]
                self.nitN_uptake[dd+1][loss_idx>=0] = soilnitNuptake_idx[loss_idx>=0]
                self.NO3_diffusivesoil[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                        NO3_diffidx[loss_idx<0]/NO3_lossall[loss_idx<0]                                           
                self.NO3_leaching[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                        NO3_leachingidx[loss_idx<0]/NO3_lossall[loss_idx<0]   
                self.nitN_uptake[dd+1][loss_idx<0] = self.NO3_pool[dd+1][loss_idx<0]*\
                                                                       soilnitNuptake_idx[loss_idx<0]/NO3_lossall[loss_idx<0]                                      
                
        return



    def chem_fert_bcsurf_sim(self,start_day_idx,end_day_idx,chem_fert_type):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [broadcasting - topdressed], fertilizer applied on the field surface')
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
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):

                ## soil resistance
                ## soil resistance = distance (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                ## resistance for diffusion from the surface source layer to the underlying topsoil layer
                self.R_soilaq[dd+1] = (p_topsoil-p_sourcelayer)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg[dd+1] = (p_topsoil-p_sourcelayer)/(tor_soil_gas*self.D_air_NH3[dd+1])
                ## resistance to the deeper soil from the topsoil layer to the underlying 2nd soil layer
                ## N is transported the to the deeper soil
                self.R_soilaq_down[dd+1] = (p_2ndsoil-p_topsoil)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg_down[dd+1] = (p_2ndsoil-p_topsoil)/(tor_soil_gas*self.D_air_NH3[dd+1])

                ## sourcelayer resistance
                ## sourcelayer resistance is determined by: R = z/(2*tor_sourcelayer*D); 
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

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1]) - self.urea_washoff[dd+1]
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

                '''## Org N pools in various forms
                self.avail_N_pool[dd+1] = (self.avail_N_pool[dd] - self.avail_N_washoff[dd+1])*(1 - self.daily_Na_decomp_rate[dd+1]) + \
                    self.avail_N[dd+1] 
                self.resist_N_pool[dd+1] = (self.resist_N_pool[dd] - self.resist_N_washoff[dd+1])* (1 - self.daily_Nr_decomp_rate[dd+1]) + \
                    self.resist_N[dd+1] 
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] - self.unavail_N_washoff[dd+1]'''

                ## TAN pool (not the top soil layer)
                ## Note: source of TAN pool: TAN production from 1) urea hydrolysis, 2) decomposition of org N and 3) input of TAN from housing
                ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil (aq+gas), 
                ##                         3) infiltration (subsurface leaching) to soil, and 4) nitrification
                ##       only aqueous phase diffusion is considered, and gaseous diffusion is not considered in this study
                '''TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                    self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]
                ## TAN pool in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6'''

                # TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                #     self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
                # self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                # self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

                ## TAN pool with input
                self.TAN_pool[dd+1] = self.TAN_pool[dd]+self.TAN_prod[dd+1]+self.TAN[dd+1] + \
                        self.diffusivefluxup_aq[dd] + self.diffusivefluxup_gas[dd]
                ## TAN pool of source layer in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6


                ## TAN conc; g/m3
                ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to sourcelayer) NH4+
                ## gaseous NH3 in air-filled pore space; epsilon(sourcelayer porosity) - theta
                ## aqueous TAN in water-filled pore space; theta
                ## solid phase TAN adsorbed on solid particles; 1 - epsilon
                ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd = 1.0 m^3/m^3 
                ## Kd = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                ## sourcelayer density varies, ~ 0.3-1.9 g/cm^3, we assmume the porosity of mamnure is 40%
                ## the volume (thickness) is determined by solid mass and bulk density (derived from porosity)
                ## NH3(g) = KNH3*[TAN(aq)]
                ## MTAN = Vsourcelayer*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
                ## so, [TAN(aq)] = MTAN/Vsourcelayer * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
                KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## TAN molar conc; mol/L
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
                ## NH3 conc in bulk sourcelayer
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when sourcelayer water content equals to porosity
                self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                ## TAN conc at the surface (solved); atmospheric NH3 is ignored
                # self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_soilaq[dd+1]+KNH3/self.R_soilg[dd+1])/\
                #                         (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_soilg[dd+1])+1/self.R_soilaq[dd+1]))/\
                #                             (14*1000)
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_sourcelayer_aq[dd+1]+KNH3/self.R_sourcelayer_gas[dd+1])/\
                                        (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_atm[dd+1]+1/self.R_sourcelayer_gas[dd+1])+1/self.R_sourcelayer_aq[dd+1]))/\
                                            (14*1000)

                # self.TAN_surf_amount_M[dd+1] = self.TAN_amount_M[dd+1]
                ## TAN surface conc in g/m3
                TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
                ## in ug
                NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

                ## determining the maximum emission; 
                emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
                ## determining the maximum TAN runoff;
                ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
                ## else need to multiply by timestep*3600
                runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
                ## determining the maximum TAN aqueous diffusion to topsoil layer;
                ## diffusion is considered to be unidirectional - from source layer to top soil
                diffaq_idx = (self.TAN_amount[dd+1]-self.TAN_soil_amount[dd])/self.R_soilaq[dd+1]*timestep*3600
                diffaq_idx[diffaq_idx<=0] = 0.0
                ## when soil mositure is 0, aqueous diffusion stops
                diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
                ## determining the maximum TAN gaseous diffusion to soil interface
                diffgas_idx = (self.NH3_gas_bulk[dd+1]-self.NH3_gas_soil[dd])/self.R_soilg[dd+1]*timestep*3600
                diffgas_idx[diffgas_idx<0] = 0.0
                ## when the soil moisture reaches the saturation, gaseous diffusion stops
                diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
                ## determining the maximum infiltration of TAN
                infil_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600
                # infil_idx = 0*self.TAN_amount[dd+1]*timestep*3600
                ## nirification rate of TAN in bulk sourcelayer; daily maximum nitrification rate is 0.1 per day
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

                ## fluxes from sourcelayer to the top soil layer
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
                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## NH3 conc in bulk sourcelayer
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when sourcelayer water content equals to porosity
                self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                ## NO3- pool of sourcelayer
                NO3_idx = self.NO3_pool[dd] - self.NO3_infilsourcelayer[dd] - self.NO3_diffusivesourcelayer[dd] - self.NO3_washoff[dd+1]
                self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + self.NO3[dd+1][NO3_idx>0]
                self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0]

                ## NO3- conc of bulk sourcelayer; g/mL
                self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                            (z_sourcelayer*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
                self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                ## NO3- conc of bulk sourcelayer in g/m3
                self.NO3_amount[dd+1] = self.NO3_amount[dd+1]*1e6

                ## diffusive aquous NO3 and infiltration of NO3 from bulk sourcelayer to soil interface
                NO3_diffidx = (self.NO3_amount[dd+1] - self.NO3_soil_amount[dd])*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
                NO3_diffidx[NO3_diffidx<0] = 0.0
                NO3_diffidx[self.soilmoist[dd+1]==0] = 0.0
                NO3_infilidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
                NO3_lossall = NO3_diffidx + NO3_infilidx 
                sourcelayerloss_idx = self.NO3_pool[dd+1] - NO3_lossall
                self.NO3_diffusivesourcelayer[dd+1][sourcelayerloss_idx>=0] = NO3_diffidx[sourcelayerloss_idx>=0]
                self.NO3_infilsourcelayer[dd+1][sourcelayerloss_idx>=0] = NO3_infilidx[sourcelayerloss_idx>=0]
                self.NO3_diffusivesourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                        NO3_diffidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                           
                self.NO3_infilsourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                        NO3_infilidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                         

                ## TAN pool in top soil layer
                # TAN_soil_idx = self.TAN_pool_soil[dd] - self.diffusivefluxsoil_aq[dd] - self.diffusivefluxsoil_gas[dd] - \
                #     self.leachingflux[dd] - self.nitrif_NO3_soil[dd] - self.ammN_uptake[dd]
                # self.TAN_pool_soil[dd+1][TAN_soil_idx>0] = TAN_soil_idx[TAN_soil_idx>0] + self.infilflux[dd+1][TAN_soil_idx>0] + \
                #                                             self.diffusivefluxsourcelayer_aq[dd+1][TAN_soil_idx>0] + \
                #                                             self.diffusivefluxsourcelayer_gas[dd+1][TAN_soil_idx>0]
                # self.TAN_pool_soil[dd+1][TAN_soil_idx<=0] =  self.infilflux[dd+1][TAN_soil_idx<=0]+ \
                #                                                 self.diffusivefluxsourcelayer_aq[dd+1][TAN_soil_idx<=0] +\
                #                                                 self.diffusivefluxsourcelayer_gas[dd+1][TAN_soil_idx<=0]

                ## TAN pool of the top soil layer
                self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd] + self.infilflux[dd+1] +\
                    self.diffusivefluxsourcelayer_aq[dd+1] + self.diffusivefluxsourcelayer_gas[dd+1]


                ## TAN conc at the soil surface/interface between sourcelayer and land (soil); g/m3
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_topsoil*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
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

                ## update TAN pool of the top soil layer
                self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd+1] - self.diffusivefluxsoil_aq[dd+1] - self.diffusivefluxsoil_gas[dd+1] - \
                    - self.diffusivefluxup_aq[dd+1] - self.diffusivefluxup_gas[dd+1] - \
                    self.leachingflux[dd+1] - self.nitrif_NO3_soil[dd+1] - self.ammN_uptake[dd+1]

                ## TAN conc at the soil surface/interface between sourcelayer and land (soil); g/m3
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_topsoil*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                    self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## NH3 concentration in the soil pore space
                self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
                ## NH3 conc is zero when soil moisture content reaches the saturation
                self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0 


                ## NO3- pool of soil interface
                NO3_soil_idx = self.NO3_pool_soil[dd] - self.NO3_leaching[dd] - self.NO3_diffusivesoil[dd] -\
                                self.nitN_uptake[dd] 
                self.NO3_pool_soil[dd+1][NO3_soil_idx>0] = NO3_soil_idx[NO3_soil_idx>0] + self.nitrif_NO3_soil[dd+1][NO3_soil_idx>0] + \
                        self.NO3_infilsourcelayer[dd+1][NO3_soil_idx>0] + self.NO3_diffusivesourcelayer[dd+1][NO3_soil_idx>0]
                self.NO3_pool_soil[dd+1][NO3_soil_idx<=0] = self.nitrif_NO3_soil[dd+1][NO3_soil_idx<=0] + \
                        self.NO3_infilsourcelayer[dd+1][NO3_soil_idx<=0] + self.NO3_diffusivesourcelayer[dd+1][NO3_soil_idx<=0]

                ## NO3- conc of soil interface; g/mL
                self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                            (z_topsoil*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])   
                self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                ## NO3- conc of soil interface in g/m3
                self.NO3_soil_amount[dd+1] = self.NO3_soil_amount[dd+1]*1e6

                ## NO3 loss through aqueous diffusion and leaching to deeper soil
                NO3_soildiffidx = self.NO3_soil_amount[dd+1]*(f_DNO3/self.R_soilaq_down[dd+1])*timestep*3600
                NO3_soilleachingidx = self.qpsoil[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
                
                NO3_soilall_loss = NO3_soildiffidx + NO3_soilleachingidx + soilnitNuptake_idx
                soilloss_idx = self.NO3_pool_soil[dd+1] - NO3_soilall_loss
                self.NO3_diffusivesoil[dd+1][soilloss_idx>=0] = NO3_soildiffidx[soilloss_idx>=0]
                self.NO3_leaching[dd+1][soilloss_idx>=0] = NO3_soilleachingidx[soilloss_idx>=0]
                self.nitN_uptake[dd+1][soilloss_idx>=0] = soilnitNuptake_idx[soilloss_idx>=0]
                self.NO3_diffusivesoil[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
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


    def chem_fert_deepinjec_sim(self,start_day_idx,end_day_idx,chem_fert_type,injection_depth):
        print('current simulation is for: '+str(chem_fert_type))
        print('technique used is [deep injection], injection depth is: '+str(injection_depth*100)+' cm')
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kd = ammonium_adsorption(clay_content=soilclay)
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
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):

                ## soil resistance
                ## soil resistance = distance (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="aqueous")
                tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1],phase="gaseous")
                ## resistance for diffusion from the surface source layer to the underlying topsoil layer
                self.R_soilaq[dd+1] = (z_topsoil/2)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg[dd+1] = (z_topsoil/2)/(tor_soil_gas*self.D_air_NH3[dd+1])
                ## resistance to the deeper soil from the topsoil layer to the underlying 2nd soil layer
                ## N is transported the to the deeper soil
                self.R_soilaq_down[dd+1] = (p_2ndsoil-injection_depth)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg_down[dd+1] = (p_2ndsoil-injection_depth)/(tor_soil_gas*self.D_air_NH3[dd+1])

                ## sourcelayer resistance
                ## sourcelayer resistance is determined by: R = z/(2*tor_sourcelayer*D); 
                ##        z is the layer thickness of sourcelayer; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not considered here
                ## layer thickness of sourcelayer (DM+water) in meter: z = Vsourcelayer
                self.R_sourcelayer_aq[dd+1] = (injection_depth-p_topsoil)/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_sourcelayer_gas[dd+1] = (injection_depth-p_topsoil)/(tor_soil_gas*self.D_air_NH3[dd+1])

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

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1]) - self.urea_washoff[dd+1]
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

                '''## Org N pools in various forms
                self.avail_N_pool[dd+1] = (self.avail_N_pool[dd] - self.avail_N_washoff[dd+1])*(1 - self.daily_Na_decomp_rate[dd+1]) + \
                    self.avail_N[dd+1] 
                self.resist_N_pool[dd+1] = (self.resist_N_pool[dd] - self.resist_N_washoff[dd+1])* (1 - self.daily_Nr_decomp_rate[dd+1]) + \
                    self.resist_N[dd+1] 
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] - self.unavail_N_washoff[dd+1]'''

                ## TAN pool (not the top soil layer)
                ## Note: source of TAN pool: TAN production from 1) urea hydrolysis, 2) decomposition of org N and 3) input of TAN from housing
                ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil (aq+gas), 
                ##                         3) infiltration (subsurface leaching) to soil, and 4) nitrification
                ##       only aqueous phase diffusion is considered, and gaseous diffusion is not considered in this study
                '''TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                    self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]
                ## TAN pool in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6'''

                # TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxsourcelayer_aq[dd] - self.diffusivefluxsourcelayer_gas[dd]-\
                #     self.infilflux[dd] - self.nitrif_NO3_sourcelayer[dd] - self.TAN_washoff[dd]
                # self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                # self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

                ## TAN pool with input
                self.TAN_pool[dd+1] = self.TAN_pool[dd]+self.TAN_prod[dd+1]+self.TAN[dd+1] + \
                    self.diffusivefluxdown_aq[dd] + self.diffusivefluxdown_gas[dd]
                ## TAN pool of source layer in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6


                ## TAN conc; g/m3
                ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to sourcelayer) NH4+
                ## gaseous NH3 in air-filled pore space; epsilon(sourcelayer porosity) - theta
                ## aqueous TAN in water-filled pore space; theta
                ## solid phase TAN adsorbed on solid particles; 1 - epsilon
                ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd = 1.0 m^3/m^3 
                ## Kd = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                ## sourcelayer density varies, ~ 0.3-1.9 g/cm^3, we assmume the porosity of mamnure is 40%
                ## the volume (thickness) is determined by solid mass and bulk density (derived from porosity)
                ## NH3(g) = KNH3*[TAN(aq)]
                ## MTAN = Vsourcelayer*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
                ## so, [TAN(aq)] = MTAN/Vsourcelayer * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
                KNH3 = self.Henry_constant[dd+1]/(sim_ccH[dd+1] + self.k_NH4[dd+1])
                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## TAN molar conc; mol/L
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
                ## NH3 conc in bulk sourcelayer
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when sourcelayer water content equals to porosity
                self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                ## fluxes to deeper soil
                ## TAN loss through aqueous diffusion and leaching to deeper soil
                sourcelayerdiffaq_idx = self.TAN_amount[dd+1]/(self.R_soilaq_down[dd+1]*timestep*3600)
                sourcelayerdiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
                sourcelayerdiffgas_idx = self.NH3_gas_bulk[dd+1]/(self.R_soilg_down[dd+1]*timestep*3600)
                sourcelayerdiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                sourcelayerupdiffaq_idx = (self.TAN_amount[dd+1] - self.TAN_soil_amount[dd+1])/self.R_sourcelayer_aq[dd+1]*timestep*3600
                sourcelayerupdiffaq_idx[sourcelayerupdiffaq_idx<0] = 0.0 
                sourcelayerupdiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
                sourcelayerupdiffgas_idx = (self.NH3_gas_bulk[dd+1] - self.NH3_gas_soil[dd+1])/self.R_sourcelayer_gas[dd+1]*timestep*3600
                sourcelayerupdiffgas_idx[sourcelayerupdiffgas_idx<0] = 0.0 
                sourcelayerupdiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                sourcelayerleaching_idx = self.qpsoil[dd+1]*self.TAN_amount[dd+1]*timestep*3600
                ## nirification rate of TAN in bulk sourcelayer; daily maximum nitrification rate is 0.1 per day
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

                ## fluxes from sourcelayer to the top soil layer
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

                self.TAN_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_sourcelayer*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## NH3 conc in bulk sourcelayer
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when sourcelayer water content equals to porosity
                self.NH3_gas_bulk[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                ## NO3- pool of sourcelayer
                NO3_idx = self.NO3_pool[dd] - self.NO3_leaching[dd] - self.NO3_diffusiveup[dd] - \
                            self.NO3_diffusivedown[dd] - self.nitN_uptake_sourcelayer[dd]
                self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_sourcelayer[dd+1][NO3_idx>0] + \
                            self.NO3[dd+1][NO3_idx>0] + self.NO3_infilsoil[dd][NO3_idx>0]
                self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_sourcelayer[dd+1][NO3_idx<=0] + self.NO3[dd+1][NO3_idx<=0] +\
                            self.NO3_infilsoil[dd][NO3_idx<=0]

                ## NO3- conc of bulk sourcelayer; g/mL
                self.NO3_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                            (z_sourcelayer*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])
                self.NO3_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                ## NO3- conc of bulk sourcelayer in g/m3
                self.NO3_amount[dd+1] = self.NO3_amount[dd+1]*1e6

                ## diffusive aquous NO3 and infiltration of NO3 from bulk sourcelayer to soil interface
                NO3_diffdownidx = (self.NO3_amount[dd+1])*(f_DNO3/self.R_sourcelayer_aq[dd+1])*timestep*3600
                NO3_diffdownidx[NO3_diffdownidx<0] = 0.0
                NO3_diffdownidx[self.soilmoist[dd+1]==0] = 0.0
                NO3_diffupidx = (self.NO3_amount[dd+1] - self.NO3_soil_amount[dd])*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
                NO3_diffupidx[NO3_diffupidx<0] = 0.0
                NO3_diffupidx[self.soilmoist[dd+1]==0] = 0.0
                NO3_leachingidx = self.qpsoil[dd+1]*self.NO3_amount[dd+1]*timestep*3600

                NO3_lossall = NO3_diffdownidx + NO3_diffupidx + NO3_leachingidx + sourcelayernitNuptake_idx
                sourcelayerloss_idx = self.NO3_pool[dd+1] - NO3_lossall
                self.NO3_diffusivedown[dd+1][sourcelayerloss_idx>=0] = NO3_diffdownidx[sourcelayerloss_idx>=0]
                self.NO3_diffusiveup[dd+1][sourcelayerloss_idx>=0] = NO3_diffupidx[sourcelayerloss_idx>=0]
                self.NO3_leaching[dd+1][sourcelayerloss_idx>=0] = NO3_leachingidx[sourcelayerloss_idx>=0]
                self.nitN_uptake_sourcelayer[dd+1][sourcelayerloss_idx>=0] = sourcelayernitNuptake_idx[sourcelayerloss_idx>=0]
                self.NO3_diffusivedown[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                NO3_diffdownidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]   
                self.NO3_diffusiveup[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                NO3_diffupidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]                                        
                self.NO3_leaching[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                NO3_leachingidx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]
                self.nitN_uptake_sourcelayer[dd+1][sourcelayerloss_idx<0] = self.NO3_pool[dd+1][sourcelayerloss_idx<0]*\
                                                                sourcelayernitNuptake_idx[sourcelayerloss_idx<0]/NO3_lossall[sourcelayerloss_idx<0]

                ## TAN pool of the top soil layer
                self.TAN_pool_soil[dd+1] = self.TAN_pool_soil[dd] + self.diffusivefluxsourcelayer_aq[dd+1] +\
                                                 self.diffusivefluxsourcelayer_gas[dd+1]

                ## TAN conc at the soil surface/interface between sourcelayer and land (soil); g/m3
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_topsoil*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
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
                # self.TAN_surf_amount_M[dd+1] = self.TAN_amount_M[dd+1]
                ## TAN surface conc in g/m3
                TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
                ## in ug
                NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

                ## determining the maximum emission; 
                emiss_idx = (NH3_gas_g*3600*timestep/self.R_atm[dd+1])
                ## determining the maximum TAN runoff;
                ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
                ## else need to multiply by timestep*3600
                runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
                ## determining the maximum TAN aqueous diffusion to soil interface;
                ## diffusion is considered to be unidirectional - fro bulk sourcelayer to soil interface
                diffaq_idx = (self.TAN_soil_amount[dd+1]-self.TAN_amount[dd+1])/self.R_sourcelayer_aq[dd+1]*timestep*3600
                diffaq_idx[diffaq_idx<=0] = 0.0
                ## when soil mositure is 0, aqueous diffusion stops
                diffaq_idx[self.soilmoist[dd+1]==0] = 0.0
                ## determining the maximum TAN gaseous diffusion to soil interface
                diffgas_idx = (self.NH3_gas_soil[dd+1]-self.NH3_gas_bulk[dd+1])/self.R_sourcelayer_gas[dd+1]*timestep*3600
                diffgas_idx[diffgas_idx<0] = 0.0
                ## when the soil moisture reaches the saturation, gaseous diffusion stops
                diffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
                ## determining the maximum infiltration of TAN
                infil_idx = self.qpsoil[dd+1]*self.TAN_soil_amount[dd+1]*timestep*3600
                ## nirification rate of TAN in bulk sourcelayer; daily maximum nitrification rate is 0.1 per day
                nitrif_idx = KNO3_soil*self.TAN_pool_soil[dd+1]*f_NH4
                ## plant N uptake
                soilammNuptake_idx, soilnitNuptake_idx = plant_N_uptake(Namm=self.TAN_pool_soil[dd+1]*f_NH4,\
                                        Nnit=self.NO3_pool_soil[dd],temp=self.T_sim[dd+1])

                ## fluxes from sourcelayer to the top soil layer
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
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0
                self.TAN_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                    (z_topsoil*(self.soilmoist[dd+1][self.soilmoist[dd+1]!=0]+\
                    KNH3[self.soilmoist[dd+1]!=0]*(self.persm[dd+1][self.soilmoist[dd+1]!=0]-\
                        self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.soilmoist[dd+1]!=0])*Kd[self.soilmoist[dd+1]!=0]))
                ## NH3 concentration in the soil pore space
                self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
                ## NH3 conc is zero when soil moisture content reaches the saturation
                self.NH3_gas_soil[dd+1][self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0

                ## NO3- pool of top soil 
                NO3_soil_idx = self.NO3_pool_soil[dd] - self.NO3_infilsoil[dd] - self.nitN_uptake[dd]
                self.NO3_pool_soil[dd+1][NO3_soil_idx>=0] = NO3_soil_idx[NO3_soil_idx>=0] + self.nitrif_NO3_soil[dd+1][NO3_soil_idx>=0] +\
                    self.NO3_diffusiveup[dd+1][NO3_soil_idx>=0]
                self.NO3_pool_soil[dd+1][NO3_soil_idx<0] = self.nitrif_NO3_soil[dd+1][NO3_soil_idx<0] + \
                    self.NO3_diffusiveup[dd+1][NO3_soil_idx<0]

                ## NO3- conc of soil interface; g/mL
                self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]!=0] = self.NO3_pool_soil[dd+1][self.soilmoist[dd+1]!=0]/\
                                                                            (z_topsoil*self.soilmoist[dd+1][self.soilmoist[dd+1]!=0])   
                self.NO3_soil_amount[dd+1][self.soilmoist[dd+1]==0] = 0.0
                ## NO3- conc of soil interface in g/m3
                self.NO3_soil_amount[dd+1] = self.NO3_soil_amount[dd+1]*1e6

                ## NO3 loss through aqueous diffusion and plant uptake
                NO3_soilinfilidx = self.qpsoil[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
                
                NO3_soilall_loss = NO3_soilinfilidx + soilnitNuptake_idx
                soilloss_idx = self.NO3_pool_soil[dd+1] - NO3_soilall_loss
                self.NO3_infilsoil[dd+1][soilloss_idx>=0] = NO3_soilinfilidx[soilloss_idx>=0]
                self.nitN_uptake[dd+1][soilloss_idx>=0] = soilnitNuptake_idx[soilloss_idx>=0]
                self.NO3_infilsoil[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                            NO3_soilinfilidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
                self.nitN_uptake[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                                    soilnitNuptake_idx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]

        return

    def N_stat(self,fert_method,chem_fert_type):

        if chem_fert_type == 'ammonium':
            sim_area = self.ammN_area
            chemfert_Ntotal = self.TAN_added*sim_area
            chemfert_NH3emiss = self.NH3_flux*sim_area
            chemfert_ammN = np.nansum(chemfert_Ntotal,axis=0)
            chemfert_ammNH3 = np.nansum(chemfert_NH3emiss,axis=0)

        elif chem_fert_type == 'urea':
            sim_area = self.ureaN_area
            chemfert_Ntotal = self.urea_added*sim_area
            chemfert_NH3emiss = self.NH3_flux*sim_area
            chemfert_ureaN = np.nansum(chemfert_Ntotal,axis=0)
            chemfert_ureaNH3 = np.nansum(chemfert_NH3emiss,axis=0)
            
        elif chem_fert_type == 'nitrate':
            sim_area = self.nitN_area
            chemfert_Ntotal = self.NO3_added*sim_area
            chemfert_NH3emiss = self.NH3_flux*sim_area
            chemfert_nitN = np.nansum(chemfert_Ntotal,axis=0)
            chemfert_nitNH3 = np.nansum(chemfert_NH3emiss,axis=0)
        
        print('Total N applied: '+str(np.nansum(chemfert_Ntotal/1e9)))
        print('NH3 emission: '+str(np.nansum(chemfert_NH3emiss/1e9)))

        chemfert_TANwashoff = np.nansum(self.TAN_washoff*sim_area)/1e9
        print('TAN washoff: '+ str(chemfert_TANwashoff))

        chemfert_leaching = np.nansum(self.leachingflux*sim_area)/1e9
        print('NH4 leaching: '+ str(chemfert_leaching))

        if fert_method == 'broadcasting-surf':

            chemfert_diffaq = np.nansum(self.diffusivefluxsourcelayer_aq*sim_area)/1e9
            print('TAN diff aq to top soil: '+ str(chemfert_diffaq))

            chemfert_diffgas = np.nansum(self.diffusivefluxsourcelayer_gas*sim_area)/1e9
            print('TAN diff gas to top soil: '+ str(chemfert_diffgas))

            chemfert_updiffaq = np.nansum(self.diffusivefluxup_aq*sim_area)/1e9
            print('TAN diff aq upwards to source layer: '+ str(chemfert_updiffaq))

            chemfert_updiffgas = np.nansum(self.diffusivefluxup_gas*sim_area)/1e9
            print('TAN diff gas upwards to source layer: '+ str(chemfert_updiffgas))

            chemfert_infil = np.nansum(self.infilflux*sim_area)/1e9
            print('NH4 infiltration to top soil: '+ str(chemfert_infil))

            chemfert_nitrif = np.nansum(self.nitrif_NO3_sourcelayer*sim_area)/1e9
            print('TAN nitrification: '+ str(chemfert_nitrif))

            chemfert_diffaqdown = np.nansum(self.diffusivefluxsoil_aq*sim_area)/1e9
            print('TAN diff aq to deeper soil: '+ str(chemfert_diffaqdown))

            chemfert_diffgasdown = np.nansum(self.diffusivefluxsoil_gas*sim_area)/1e9
            print('TAN diff gas to deeper soil: '+ str(chemfert_diffgasdown))

            chemfert_nitrif_soil = np.nansum(self.nitrif_NO3_soil*sim_area)/1e9
            print('NH4 nitrif: '+ str(chemfert_nitrif_soil))

            chemfert_ammN_uptake = np.nansum(self.ammN_uptake*sim_area)/1e9
            print('NH4 uptake by plants: '+ str(chemfert_ammN_uptake))

            chemfert_nitN_uptake = np.nansum(self.nitN_uptake*sim_area)/1e9
            print('NO3 uptake by plants: '+ str(chemfert_nitN_uptake))

            chemfert_NO3diffusive = np.nansum(self.NO3_diffusivesourcelayer*sim_area)/1e9
            print('NO3 diffusion from source layer to soil pool: '+ str(chemfert_NO3diffusive))

            chemfert_NO3leaching = np.nansum(self.NO3_leaching*sim_area)/1e9
            print('NO3 leaching: '+str(chemfert_NO3leaching))

            chemfert_NO3diffusionsoil = np.nansum(self.NO3_diffusivesoil*sim_area)/1e9
            print('NO3 diffusion to deeper soil: '+str(chemfert_NO3diffusionsoil))
                
        elif fert_method == 'broadcasting-disk':

            chemfert_nitrif = np.nansum(self.nitrif_NO3_sourcelayer*sim_area)/1e9
            print('TAN nitrification: '+ str(chemfert_nitrif))

            chemfert_diffaqdown = np.nansum(self.diffusivefluxsoil_aq*sim_area)/1e9
            print('TAN diff aq to deeper soil: '+ str(chemfert_diffaqdown))

            chemfert_diffgasdown = np.nansum(self.diffusivefluxsoil_gas*sim_area)/1e9
            print('TAN diff gas to deeper soil: '+ str(chemfert_diffgasdown))

            chemfert_ammN_uptake = np.nansum(self.ammN_uptake*sim_area)/1e9
            print('NH4 uptake by plants: '+ str(chemfert_ammN_uptake))

            chemfert_nitN_uptake = np.nansum(self.nitN_uptake*sim_area)/1e9
            print('NO3 uptake by plants: '+ str(chemfert_nitN_uptake))

            chemfert_NO3diffusive = np.nansum(self.NO3_diffusivesourcelayer*sim_area)/1e9
            print('NO3 diffusion from source layer to soil pool: '+ str(chemfert_NO3diffusive))

            chemfert_NO3leaching = np.nansum(self.NO3_leaching*sim_area)/1e9
            print('NO3 leaching: '+str(chemfert_NO3leaching))

            chemfert_NO3diffusionsoil = np.nansum(self.NO3_diffusivesoil*sim_area)/1e9
            print('NO3 diffusion to deeper soil: '+str(chemfert_NO3diffusionsoil))    

        elif fert_method == 'deep injection':

            chemfert_nitrif = np.nansum(self.nitrif_NO3_sourcelayer*sim_area)/1e9
            print('NH4 nitrification (source layer): '+ str(chemfert_nitrif))

            chemfert_diffaqdown = np.nansum(self.diffusivefluxsoil_aq*sim_area)/1e9
            print('TAN diff aq to deeper soil: '+ str(chemfert_diffaqdown))

            chemfert_diffgasdown = np.nansum(self.diffusivefluxsoil_gas*sim_area)/1e9
            print('TAN diff gas to deeper soil: '+ str(chemfert_diffgasdown))

            chemfert_nitrif_soil = np.nansum(self.nitrif_NO3_soil*sim_area)/1e9
            print('NH4 nitrif (topsoil layer): '+ str(chemfert_nitrif_soil))

            chemfert_ammN_uptakesl = np.nansum(self.ammN_uptake_sourcelayer*sim_area)/1e9
            print('NH4 uptake by plants (source layer): '+ str(chemfert_ammN_uptakesl))

            chemfert_ammN_uptake = np.nansum(self.ammN_uptake*sim_area)/1e9
            print('NH4 uptake by plants (topsoil layer): '+ str(chemfert_ammN_uptake))

            chemfert_nitN_uptakesl = np.nansum(self.nitN_uptake_sourcelayer*sim_area)/1e9
            print('NO3 uptake by plants (source layer): '+ str(chemfert_nitN_uptakesl))

            chemfert_nitN_uptake = np.nansum(self.nitN_uptake*sim_area)/1e9
            print('NO3 uptake by plants (topsoil layer): '+ str(chemfert_nitN_uptake))

            chemfert_NO3leaching = np.nansum(self.NO3_leaching*sim_area)/1e9
            print('NO3 leaching: '+str(chemfert_NO3leaching))

            chemfert_NO3diffusionsoil = np.nansum(self.NO3_diffusivedown*sim_area)/1e9
            print('NO3 diffusion to deeper soil: '+str(chemfert_NO3diffusionsoil))

        return

    def chem_fert_main(self,fert_method,crop_item,chem_fert_type,start_day_idx,end_day_idx,fert_depth,sim_stat=False):
        self.sim_env()
        if fert_method == 'broadcasting-surf':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_bcsurf_sim(start_day_idx,end_day_idx,chem_fert_type)
            
                
        elif fert_method == 'broadcasting-disk':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_bcdisk_sim(start_day_idx,end_day_idx,chem_fert_type,fert_depth)

        elif fert_method == 'deep injection':
            self.chem_fert_input(crop=crop_item)
            self.app_method(application_method=None)
            self.chem_fert_deepinjec_sim(start_day_idx,end_day_idx,chem_fert_type,fert_depth)

        if sim_stat is True:
            self.N_stat(fert_method,chem_fert_type)

        return

