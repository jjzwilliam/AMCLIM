
####################################
## import essential AMCLIM modules
####################################
'''file_path = os.getcwd()
module_path = str(Path(file_path).parent)
print(module_path)

if module_path not in sys.path:
    sys.path.append(module_path)'''
from asyncio import format_helpers
from INPUT.input import *
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *

## this refers to the storage surface area in farms relative to the housing area
## e.g., mms_indoor_solid = 0.2, 
##       which means the area of a barn that stores solid manure has 1/5 of housing area in the local farm;
##       housing area = 10 km^2 in a grid; indoor_solid area = 1 km^2 (1/5) * f_MMS_indoor_solid
## these values need to be reviewed (?)
MMS_area_factor = {
    "mms_indoor_solid":0.25,
    "mms_indoor_liquid":0.5,
    "mms_open_solid":0.25,
    "mms_open_liquid":0.5,
    "mms_open_lagoon":2.5}

## areas of each MMS, initial values are None
MMS_area = {
    "mms_indoor_solid_area":None,
    "mms_indoor_liquid_area":None,
    "mms_open_solid_area":None,
    "mms_open_liquid_area":None}

###################################
## MMS parameters
###################################
## adsorption constant for manure; m3/m3
Kd_manure = 1.0
## surface roughness height of water surface
zo_water = 0.002
## assuming the roughness height of manure storage barn is ~ 0.5m (<ref height of 2m)
zo_barn = 0.5  
## assuming the roughness height of manure pile (open land) is ~ 1.0m (<ref height of 2m)
zo_manureland = 1.0
## manure surface thickness is set to be 2 cm
z_manuresurf = 0.02
## source layer of manure for NH3 emission; 1cm thickness
z_topmanure = 0.01
## dry matter (DM) content of liquid manure is assumed to be 5%
f_DM_liquid = 0.05
## assuming the soil interface (source layer) of 4mm
# z_soil = 0.004
z_soil = 0.02  ## test 2cm
## distance to deeper soil; 2cm 
# d_deepsoil = 0.02
d_deepsoil = 0.025  ## test 2.5cm
## assuming infiltration of manure water to the soil is 10mm/day (1cm/day; 10 000 g/m^2/day) ref: Vira et al.,2020 GMD (2x d0)
# dailymaxinfil = 0.01
dailymaxinfil = 0.003 ## test: 3mm /day
## assuming the water holding capacity of manure is 3.0 g water/ g manure
absorb_factor = 3.0
## washoff coefficients: 0.1%/mm water for N species, and 0.05%/mm water for non N species (manure)
f_washoff_nonN = 0.0005
f_washoff_N = 0.001
## anoxic fraction of manure; average of swine (75%) and cow (50%); ref: Wang et al., Sci.Report.2015
f_manure_anoxic = 0.6
## lagoon TAN conc; g/mL
lagoon_TAN_conc = 630.0e-6


##################################
## MMS module
##################################
class MMS_module:
    def __init__(self,livestock_name,production_system_lvl_idx,mms_cat,phase,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,water_added,pH_value,area_housing):
        ## livestock and production system info
        self.livestock = livestock_name
        self.lvl_idx = production_system_lvl_idx
        self.production_system = CONFIG_production_system_dict[self.livestock][self.lvl_idx]
        ## show current config settings, e.g., MMS type
        print('MMS Module - current MMS is for production system: '+str(self.livestock)+', '+str(self.production_system))
        #####################################################
        ## livestock info and MMS info
        #####################################################
        ## read livestock and the corresponding MMS datasets
        self.animal_file_name = CONFIG_animal_file_dict[self.livestock]
        self.MMS_file_name = CONFIG_MMS_file_dict[self.livestock] 
        self.animal_file = xr.open_dataset(infile_path+animal_data_path+self.animal_file_name)
        self.MMS_file = xr.open_dataset(infile_path+animal_data_path+self.MMS_file_name)

        self.MMS_type_list = []
        for MMS in self.MMS_file.data_vars:
            self.MMS_type_list.append(str(MMS))

        self.f_loss = np.zeros(CONFIG_mtrx[1:])
        self.f_sold = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_fuel = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_preserve_solid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_preserve_liquid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_indoor_solid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_indoor_liquid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_open_solid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_open_liquid = np.zeros(CONFIG_mtrx[1:])
        self.f_MMS_open_lagoon = np.zeros(CONFIG_mtrx[1:]) 

        for mms in loss_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_loss = self.f_loss + f_mms 
            except:pass
        for mms in sold_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_sold = self.f_sold + f_mms
            except:pass
        for mms in MMS_fuel_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_fuel = self.f_MMS_fuel + f_mms
            except:pass
        for mms in MMS_preserve_solid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                f_MMS_preserve_solid = f_MMS_preserve_solid + f_mms
            except:pass
        for mms in MMS_preserve_liquid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_preserve_liquid = self.f_MMS_preserve_liquid + f_mms
            except:pass
        for mms in MMS_indoor_solid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_indoor_solid = self.f_MMS_indoor_solid + f_mms
            except:pass
        for mms in MMS_indoor_liquid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_indoor_liquid = self.f_MMS_indoor_liquid + f_mms
            except:pass
        for mms in MMS_outdoor_solid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_open_solid = self.f_MMS_open_solid + f_mms
            except:pass
        for mms in MMS_outdoor_liquid_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_open_liquid = self.f_MMS_open_liquid + f_mms
            except:pass
        for mms in MMS_outdoor_lagoon_list:
            try:
                f_mms = self.MMS_file[mms][self.lvl_idx].values
                f_mms = np.nan_to_num(f_mms)
                self.f_MMS_open_lagoon = self.f_MMS_open_lagoon + f_mms
            except:pass

        ## MMS info excluding f_HOUSING and f_spreading
        self.f_MMS = self.f_MMS_fuel+self.f_MMS_preserve_solid+self.f_MMS_preserve_liquid+self.f_MMS_indoor_solid+\
                        self.f_MMS_indoor_liquid+self.f_MMS_open_solid+self.f_MMS_open_liquid+self.f_MMS_open_lagoon

        field_shape = (lats,lons)
        array_shape = (25,lats,lons)
        ## output shape
        outarray_shape = (Days,lats,lons)

        ## feces input from housing
        self.manure_added = manure_added
        self.manure = np.zeros(array_shape)
        ## feces pool
        self.manure_pool = np.zeros(array_shape)
        ## urea input from housing
        self.urea_added = urea_added
        self.urea = np.zeros(array_shape)
        ## urea pool
        self.urea_pool = np.zeros(array_shape)
        ## uric acid input from housing
        self.UA_added = UA_added
        self.UA = np.zeros(array_shape)
        ## uric acid pool
        self.UA_pool = np.zeros(array_shape)
        ## input of available N component from housing
        self.avail_N_added = avail_N_added
        self.avail_N = np.zeros(array_shape)
        ## Nitrogen pool that is available (easily) to form TAN
        self.avail_N_pool = np.zeros(array_shape)
        ## input of resistant N component from housing
        self.resist_N_added = resist_N_added
        self.resist_N = np.zeros(array_shape)
        ## Nitrogen pool that is resistant (slowly) to form TAN
        self.resist_N_pool = np.zeros(array_shape)
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
        ## TAN pool (the bulk manure)
        self.TAN_pool = np.zeros(array_shape)
        ## TAN pool conc (aqueous phase)
        self.TAN_amount = np.zeros(array_shape)
        ## water added from housing
        self.water_added = water_added
        self.water = np.zeros(array_shape)
        ## total water pool of the system (manure water; urine+manure water+[washing water])
        self.Total_water_pool = np.zeros(array_shape) 
        ## NH3 concentration in the bulk manure
        self.NH3_gas_bulk = np.zeros(array_shape)
        ## surface NH3 concentrtion at equilirium (in ug)
        # self.NH3_gas_ug = np.zeros(array_shape)
        ## emission potential
        self.modelled_emiss = np.zeros(array_shape)
        ## final emission
        ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3flux]
        ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to
        ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
        ## during land spreading stages, such as canopy recapture, deap injection...
        self.NH3flux = np.zeros(array_shape)

        ## temp
        self.T_sim = np.zeros(array_shape)
        self.T_gnd = np.zeros(array_shape)
        ## wind/ventilation
        self.u_sim = np.zeros(array_shape)
        ## RH
        self.RH_sim = np.zeros(array_shape)
        ## evaporation
        self.evap_sim = np.zeros(array_shape)
        ## barn/buidling resistance
        self.R_star = np.zeros(array_shape)
        ## manure resistance for aqueous diffusion
        self.R_manurel = np.zeros(array_shape)
        ## manure resistance for gaseous diffusion
        self.R_manureg = np.zeros(array_shape)
        ## daily UA hydrolysis rate
        self.ua_conv_factor = np.zeros(array_shape)
        ## daily urea hydrolysis rate
        self.urea_hydro_rate = np.zeros(array_shape)
        ## daily decomposition rate of available and resistant N components
        self.Na_decomp_rate = np.zeros(array_shape)
        self.Nr_decomp_rate = np.zeros(array_shape)
        ## pH and H+ ions concentration
        self.pH = pH_value
        self.cc_H = np.float(10**(-pH_value))
        ## housing area that is used to determine MMS area
        self.housingarea = area_housing
        self.mmsarea = area_housing.shape
        ## output of NH3 flux
        self.o_NH3flux = np.zeros(outarray_shape)

        ## spring/winter manure application calendar
        self.spring_app_cal = np.zeros(field_shape)
        self.winter_app_cal = np.zeros(field_shape)

        ## N app: including TAN and urea
        self.spring_Napp = np.zeros(field_shape)
        ## UAN app: UA-N (from poultry litter)
        self.spring_UANapp = np.zeros(field_shape)
        self.spring_manureapp = np.zeros(field_shape)
        self.spring_waterapp = np.zeros(field_shape)
        self.spring_availNapp = np.zeros(field_shape)
        self.spring_resistNapp = np.zeros(field_shape)
        self.spring_unavailNapp = np.zeros(field_shape)

        self.winter_Napp = np.zeros(field_shape)
        self.winter_UANapp = np.zeros(field_shape)
        self.winter_manureapp = np.zeros(field_shape)
        self.winter_waterapp = np.zeros(field_shape)
        self.winter_availNapp = np.zeros(field_shape)
        self.winter_resistNapp = np.zeros(field_shape)
        self.winter_unavailNapp = np.zeros(field_shape)

        ## total N of this MMS
        self.total_N_MMS = np.nansum(self.TAN_added+self.urea_added+self.UA_added+\
                            self.avail_N_added+self.resist_N_added+self.unavail_N_added,axis=0)

        ## manure N to land spreading from MMS
        # self.manure_to_land = np.zeros(outarray_shape)
        # self.avail_N_to_land = np.zeros(outarray_shape)
        # self.resist_N_to_land = np.zeros(outarray_shape)
        # self.unavail_N_to_land = np.zeros(outarray_shape)
        # self.TAN_to_land = np.zeros(outarray_shape)
        # self.water_to_land = np.zeros(outarray_shape) 

        if phase == "solid":
            ## TAN concentration at surface (compensation point)
            self.TAN_surf_amount = np.zeros(array_shape)
            ## NH3 concentration at the surface (compensation point)
            self.NH3_gas_surf = np.zeros(array_shape)
            ## NO3- from nitrification in bulk manrue 
            self.TANnitrif_manure = np.zeros(array_shape)
            ## NO3- pool in the bulk manure
            self.NO3_pool = np.zeros(array_shape)
            ## output of nitrification
            self.o_NH4nitrif = np.zeros(outarray_shape)

            if mms_cat == "MMS_open":
                ## rain fall (Note the unit)
                self.rainfall = np.zeros(array_shape)
                ## atmospheric resistances
                self.R_atm = np.zeros(array_shape)
                ## rain available for "washoff"
                self.rain_avail_washoff = np.zeros(array_shape)
                ## washoff flux of manure
                self.manure_washoff = np.zeros(array_shape)
                ## washoff flux of urea
                self.ureawashoff = np.zeros(array_shape)
                ## washoff flux of UA
                self.UA_washoff = np.zeros(array_shape)
                ## washoff flux of available org N 
                self.avail_N_washoff = np.zeros(array_shape)
                ## washoff flux of resistant org N 
                self.resist_N_washoff = np.zeros(array_shape)
                ## washoff flux of TAN 
                self.TANwashoff = np.zeros(array_shape)
                ## NO3- washoff 
                self.NO3washoff = np.zeros(array_shape)
                ## NO3- conc in the bulk manure
                self.NO3_amount = np.zeros(array_shape)
                ## output of N fluxes (pathways)
                self.o_Nwashoff = np.zeros(outarray_shape)

        elif  phase == "liquid":
            if mms_cat == "MMS_open":
                ## rain fall (Note the unit)
                self.rainfall = np.zeros(array_shape)
                ## atmospheric resistances
                self.R_atm = np.zeros(array_shape)
        
        ## dry matter (DM) content of solid manure 
        self.DM_content = solid_m_DM[self.livestock]
        ## the density of manure; 1t kg/m^3 or 1g/cm^3
        self.manure_density = rho_m[self.livestock]
        ## maximum water content of manure
        self.f_wcmax = 1 - (self.DM_content/100)/2

    def met_input_interp(self,template):
        ##################################
        ## fill land input data
        ##################################
        try:
            self.T_gnd = field_var_fill(sd_template=template,input_field=self.T_gnd)  ## degC
        except:pass
        self.evap_sim = field_var_fill(sd_template=template,input_field=self.evap_sim) ## g/day
        self.u_sim = field_var_fill(sd_template=template,input_field=self.u_sim) ## m/s
        self.R_star = field_var_fill(sd_template=template,input_field=self.R_star)  ## s/m
        return

    def sim_env(self,mms_type,mms_phase,dayidx):
        
        if CONFIG_machine == "STREAM":
            temp_data = temp_file.t2m[dayidx] - 273.15
            rhum_data = rhum_file.Relative_Humidity_2m_06h[dayidx]
            wind_data = wind_file.Wind_Speed_10m_Mean[dayidx]
        else:
            hhidx = dayidx*24
            temp_data = temp_file.t2m[hhidx:hhidx+24] - 273.15
            rhum_data = rhum_file.rh2m[hhidx:hhidx+24]
            wind_data = wind_file.wind10m[hhidx:hhidx+24]
        if mms_type == 'MMS_indoor':
            self.T_sim[1:],self.T_gnd[1:],self.u_sim[1:] = barn_env(temp_data,
                wind_profile(uref=wind_data,height_ref=wind_data_height,height_out=ref_height,zo=zo_barn))
            self.T_sim[1:] = temp_data
            self.RH_sim[1:] = rhum_data
            self.T_sim = xr_to_np(self.T_sim)
            self.T_gnd = xr_to_np(self.T_gnd)
            self.RH_sim = xr_to_np(self.RH_sim)
            self.u_sim = xr_to_np(self.u_sim)
            u2m= wind_profile(uref=self.u_sim,height_ref=10,height_out=2,zo=zo_barn)

            # self.evap_sim = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=u2m,zo=zo_barn)
            # self.R_star = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap_sim)
            # ## convert evap from m/s to g/m2/day
            # self.evap_sim = self.evap_sim*1e6*timestep*3600

            if mms_phase == "liquid":
                ## daily evaporation; aerodynamic method; (m/s)
                ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
                self.evap_sim = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=u2m,zo=zo_barn)
                ## convert evap from m/s to g/m2/day
                self.evap_sim = self.evap_sim*1e6*timestep*3600
                # print('MMS ENV: mms barn (liquid)')
            elif mms_phase == "solid":
                # self.evap_sim = self.evap_sim*1e6*timestep*3600
                self.evap_sim = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=u2m,zo=zo_barn)
                ## convert evap from m/s to g/m2/day
                f_evap_correct = 1.0
                self.evap_sim = self.evap_sim*1e6*timestep*3600*f_evap_correct
                # print('MMS ENV: mms barn (solid)')

        elif mms_type == "MMS_open":
            evap_data = evap_file.evabs[dayidx]*(-1e6)
            rain_data = rain_file.tcrw[dayidx]*1e3
            ram1_data = ratm_file.RAM1[dayidx]
            rb1_data = ratm_file.RB1[dayidx]
            if mms_phase == "liquid":
                self.T_sim[1:] = temp_data
                self.u_sim[1:] = wind_profile(uref=wind_data,height_ref=wind_data_height,height_out=ref_height,zo=zo_barn)
                self.evap_sim[1:] = evap_data/24
                self.rainfall[1:] = rain_data/24
            elif mms_phase == "solid":
                groundtemp_datalvl1 = groundtemp_filelvl1.stl1[dayidx] - 273.15
                self.T_sim[1:] = temp_data
                self.T_gnd[1:] = groundtemp_datalvl1
                self.u_sim[1:] = wind_data
                self.RH_sim[1:] = rhum_data
                self.evap_sim[1:] = evap_data/24
                self.rainfall[1:] = rain_data/24
                self.R_atm[1:] = ram1_data+rb1_data
            self.met_input_interp(template=self.animal_file['Excreted_N'][self.lvl_idx])
            # print('MMS ENV: mms open env')
        return
    
    ## Simulation: Cat A.2 and B.2 manure stored as liquid
    ## water pool is transfered from housing to slurry tank/barn
    def MMS_liquid_sim(self,mms_cat,start_day_idx,end_day_idx):
        print('current simulation is for: ['+str(mms_cat)+', liquid], '+\
                                str(self.livestock)+', '+str(self.production_system))
        if mms_cat == "MMS_indoor":
            f_MMSarea = self.f_MMS_indoor_liquid/self.f_MMS
            ## area of MMS_indoor_liquid
            MMS_area["mms_indoor_liquid_area"] = self.housingarea*(1.0-self.f_loss-self.f_sold)*f_MMSarea*MMS_area_factor['mms_indoor_liquid']
            self.mmsarea = MMS_area["mms_indoor_liquid_area"].values
        elif mms_cat == "MMS_open":
            f_MMSarea = self.f_MMS_open_liquid/self.f_MMS
            ## area of MMS_open_liquid
            MMS_area["mms_open_liquid_area"] = self.housingarea*(1.0-self.f_loss-self.f_sold)*f_MMSarea*MMS_area_factor['mms_open_liquid']
            self.mmsarea = MMS_area["mms_open_liquid_area"].values
        
        if self.livestock.lower()=="poultry":
            print("Simulation for poultry.")

        for dd in np.arange(start_day_idx,end_day_idx):
            ## environment
            if dd<Days:
                self.sim_env(mms_type=mms_cat,mms_phase="liquid",dayidx=dd)
            if dd>=Days:
                self.sim_env(mms_type=mms_cat,mms_phase="liquid",dayidx=dd-Days)
            ## daily input from housing
            if mms_cat == "MMS_indoor":
                self.daily_init(dayidx=dd,mms_info="mms_indoor_liquid")
            elif mms_cat == "MMS_open":
                self.daily_init(dayidx=dd,mms_info="mms_indoor_liquid")
            
            ## simulations at hourly timestep
            for hh in np.arange(24):
                
                ## manure pool
                self.manure_pool[hh+1] = self.manure_pool[hh] + self.manure[hh+1]
                
                ## urea/UA hydrolysis rate; orgN decomposition rate
                self.urea_hydro_rate[hh+1] = urea_hydrolysis_rate(temp=self.T_sim[hh+1],theta=1.0,delta_t=timestep)
                ## UA hydrolysis will no longer be constrained by RH (moisture) under liquid conditions/slurry management
                self.ua_conv_factor[hh+1] = ua_hydrolysis_rate(temp=self.T_sim[hh+1],rhum=100.0,ph=self.pH,
                                    delta_t=timestep)
                self.Na_decomp_rate[hh+1], self.Nr_decomp_rate[hh+1] = N_pools_decomp_rate(temp=self.T_sim[hh+1], 
                                                                                                    delta_t=timestep)                  
                
                ## urea/UA pool and orgN pools
                self.urea_pool[hh+1] = self.urea_pool[hh] + self.urea[hh+1]
                self.UA_pool[hh+1] = self.UA_pool[hh] + self.UA[hh+1]
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh] + self.avail_N[hh+1]
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh] + self.resist_N[hh+1]
                self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh] + self.unavail_N[hh+1]
                
                ## TAN production from urea/UA hydrolysis and the N decomposition rate from dung
                self.TAN_prod[hh+1] = self.urea_hydro_rate[hh+1]*self.urea_pool[hh+1]+\
                                        self.ua_conv_factor[hh+1]*self.UA_pool[hh+1]+\
                                        self.Na_decomp_rate[hh+1]*self.avail_N_pool[hh+1] +\
                                        self.Nr_decomp_rate[hh+1]*self.resist_N_pool[hh+1] 
                
                ## update urea/UA pool and orgN pools
                self.urea_pool[hh+1] = self.urea_pool[hh+1] * (1-self.urea_hydro_rate[hh+1])
                self.UA_pool[hh+1] = self.UA_pool[hh+1] * (1-self.ua_conv_factor[hh+1])
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1]*(1 - self.Na_decomp_rate[hh+1]) 
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1]* (1 - self.Nr_decomp_rate[hh+1])
                
                ## water amount range: determined by DM content (5% - 10%)
                minwater = self.manure_pool[hh+1]/(f_DM_liquid*2) - self.manure_pool[hh+1]
                maxwater = self.manure_pool[hh+1]/(f_DM_liquid/2) - self.manure_pool[hh+1]
                ## water pool
                ## Note: the water pool in [MMS indoor liquid] is "inheritated" from housing water pool
                ## water pool should be larger than the minimum water (assumed a ~10 % DM) amount of the "liquid" MMS
                stdwater = self.manure[hh+1]/f_DM_liquid - self.manure[hh+1]
                self.water[hh+1] = np.maximum(stdwater,self.water[hh+1])
                self.Total_water_pool[hh+1] = self.Total_water_pool[hh] + self.water[hh+1] - self.evap_sim[hh+1]
                
                ## open env scheme for determining the water pool, rainfall is included
                if mms_cat == "MMS_open":
                    self.Total_water_pool[hh+1] = self.Total_water_pool[hh+1] + self.rainfall[hh+1] 
                self.Total_water_pool[hh+1] = np.maximum(self.Total_water_pool[hh+1],minwater)
                self.Total_water_pool[hh+1] = np.minimum(self.Total_water_pool[hh+1],maxwater)                             
                
                ## TAN pool
                self.TAN_pool[hh+1] = self.TAN_pool[hh] + self.TAN_prod[hh+1] + self.TAN[hh+1]
                ## TAN conc
                self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/self.Total_water_pool[hh+1]
                self.TAN_amount[hh+1] = np.nan_to_num(self.TAN_amount[hh+1],posinf=0.0,neginf=0.0)
                ## TAN conc in g/m3
                self.TAN_amount[hh+1] = self.TAN_amount[hh+1] * 1e6
                ## NH3 conc in g/m3
                KNH3 = NH3_par_coeff(temp=self.T_sim[hh+1],cncH=self.cc_H)
                
                ## resistance
                self.R_star[hh+1] = (1/k_aq_NH4(temp=self.T_sim[hh+1]))+(1/(KNH3*k_gas_NH3(temp=self.T_sim[hh+1],
                                                                    u=self.u_sim[hh+1],Z=2,zo=zo_water)))
                ## emission (volatilization process similar to the pit in animal houses; liquid-base)
                self.modelled_emiss[hh+1] = NH3_volpit(pit_tanconc=self.TAN_amount[hh+1],conc_in=0.0,
                                                Rpit=self.R_star[hh+1])*timestep*3600
                ## determining the maximum emission; emission cannot exceed TAN pool
                self.modelled_emiss[hh+1] = np.minimum(self.modelled_emiss[hh+1],self.TAN_pool[hh+1])
                ## final emission flux
                self.NH3flux[hh+1] = self.modelled_emiss[hh+1]
                ## update TAN pool
                self.TAN_pool[hh+1] = self.TAN_pool[hh+1] - self.NH3flux[hh+1]
            
            if dd >= Days:
                ddidx = dd - Days
            else:
                ddidx = dd
            if mms_cat == "MMS_indoor":
                self.daily_output(dayidx=ddidx,mms_info="mms_indoor_liquid")
            elif mms_cat == "MMS_open":
                self.daily_output(dayidx=ddidx,mms_info="mms_open_liquid")
            self.mms_to_landspreading(dayidx=ddidx)
        return

    ## Simulation: Cat A.1 manure stored in barns (as solid) 
    ## manure/animal waste is processed/dried to reduce water content; 
    ## water pool is determined by:
    ## a) water transferred from housing b) referenced water content for solid manure management c) evaporation
    ## NEW) incorporate a [surface compensation point] scheme in the model:
    ##       gaseouos concentrations (x): x_indoor; x_surface; x_bulk
    ##       aqueous concrntrations: [TAN]_bulk (TAN conc in the manure); [TAN]_surf (TAN conc at the surface, in equilibrium with x_surf)
    ##       equilibrium between [TAN]_surf and x_surf: x_surf = (k_H_D/([H+]+k_NH4+))*[TAN]_surf = KNH3*[TAN]_surf; 
    ##                   Note: k_H_D = (161500/(temp+273.15))*np.exp(-10380/(temp+273.15)) 
    ##       fluxes (F): F_atm, surface to indoor; F_tosurf: manure to surface;;; 
    ##                   F_atm=(x_surf-x_indoor)/R_star; 
    ##                   F_tosurfaq=([TAN]_bulk-[TAN]_surf)/R_manureaq; F_tosurfgas=(x_bulk-x_surf)/R_manuregas; 
    ##                   F_atm = F_tosurf(aq+gas)
    ##    [TAN]_bulk is the prognostic variable and is determined by source and loss based on the mass balance approach
    ##    solve [TAN]_surf: [TAN]_surf = (x_indoor/R_ab+[TAN]_bulk*(1/R_manureaq+KNH3/R_manuregas))/
    ##                                      (1/R_manureaq+KNH3*(R_manureg+R_ab))
    ##    indoor concentration is neglected, so x_indoor is set to 0 for global simulations

    ## Simulation: Cat B.1 manure stored in open environment (as solid) 
    ## manure/animal waste is left on open land to reduce water content (evaporation); 
    ## water pool is determined by:
    ## 1) maximum water held by manure (manure pool x absorb factor) 2) water transferred from housing
    ##      excess water is going to be wash off flux 3) evaporation
    ## NEW) incorporate a [surface compensation point] scheme in the model as Cat B [MMS indoor solid]
    ## NEW) no infiltration is considered in solid phase open-env simulations, no underlying soil pools
    ##      1) paved or concrete surface, 2) solid phase, infiltration is small, so ignored
    ##       gaseouos concentrations (x): x_atm; x_surface; x_bulk
    ##       aqueous concrntrations: [TAN]_bulk; [TAN]_surf; 
    ##       equilibrium between aqueous and gaseous conc: x = KNH3 * [TAN]
    ##       fluxes (F_upwards): F_atm, surface to atmosphere; F_tosurf: manure to surface;;; 
    ##                   F_atm=(x_surf-x_atm)/R_ab;  F_runoff = qr*[TAN]_surf
    ##                   F_tosurfaq=([TAN]_bulk-[TAN]_surf)/R_manureaq; F_tosurfgas = (x_bulk - x_surf)/R_manuregas
    ##                   F_tosurf(aq+gas) = F_atm + F_runoff
    ## Note) NO3- pools in the bulk manure and soil have similar processes as TAN, which do not include 1) adsorption, 2) gaseous diffustion
    ##       NO3- aqueous diffusion is moderated by a scaling factor regarding the different diffusivity of NO3- from NH4+ 

    def MMS_solid_sim(self,mms_cat,start_day_idx,end_day_idx):
        print('current simulation is for: ['+str(mms_cat)+', solid], '+\
                                str(self.livestock)+', '+str(self.production_system))
        if mms_cat == "MMS_indoor":
            f_MMSarea = self.f_MMS_indoor_solid/self.f_MMS
            MMS_area["mms_indoor_solid_area"] = self.housingarea*(1.0-self.f_loss-self.f_sold)*f_MMSarea*MMS_area_factor["mms_indoor_solid"]
            self.mmsarea = MMS_area["mms_indoor_solid_area"].values
        elif mms_cat == "MMS_open":
            f_MMSarea = self.f_MMS_open_solid/self.f_MMS
            MMS_area["mms_open_solid_area"] = self.housingarea*(1.0-self.f_loss-self.f_sold)*f_MMSarea*MMS_area_factor["mms_open_solid"]
            self.mmsarea = MMS_area["mms_open_solid_area"].values

        if self.livestock.lower()=="poultry":
            print("Simulation for poultry.")

        for dd in np.arange(start_day_idx,end_day_idx):
            ## environemnts
            if dd < Days:
                self.sim_env(mms_type=mms_cat,mms_phase="solid",dayidx=dd)
            if dd>=Days:
                self.sim_env(mms_type=mms_cat,mms_phase="solid",dayidx=dd-Days)
            ## daily input from housing
            if mms_cat == "MMS_indoor":
                self.daily_init(dayidx=dd,mms_info="mms_indoor_solid")
            elif mms_cat == "MMS_open":
                self.daily_init(dayidx=dd,mms_info="mms_indoor_solid")
            
            ## simulations at hourly timestep
            for hh in np.arange(24):
                ## manure pool
                self.manure_pool[hh+1] = self.manure_pool[hh] + self.manure[hh+1] 

                ## urea/UA hydrolysis rate; orgN decomposition rate
                self.urea_hydro_rate[hh+1] = urea_hydrolysis_rate(temp=self.T_sim[hh+1],theta=1.0,delta_t=timestep)
                ## note the difference between the liquid and solid management; 
                ## UA hydrolysis rate is constrained by RH under solid management
                self.ua_conv_factor[hh+1] = ua_hydrolysis_rate(temp=self.T_sim[hh+1],rhum=self.RH_sim[hh+1],ph=self.pH,
                                    delta_t=timestep)
                self.Na_decomp_rate[hh+1], self.Nr_decomp_rate[hh+1] = N_pools_decomp_rate(temp=self.T_sim[hh+1], 
                                                                                                    delta_t=timestep)
                
                ## urea/UA pool and orgN pools
                self.urea_pool[hh+1] = self.urea_pool[hh] + self.urea[hh+1]
                self.UA_pool[hh+1] = self.UA_pool[hh] + self.UA[hh+1]
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh] + self.avail_N[hh+1]
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh] + self.resist_N[hh+1]
                self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh] + self.unavail_N[hh+1]

                ## TAN production from urea/UA hydrolysis and the N decomposition rate from dung
                self.TAN_prod[hh+1] = self.urea_hydro_rate[hh+1]*self.urea_pool[hh+1]+\
                                        self.ua_conv_factor[hh+1]*self.UA_pool[hh+1]+\
                                        self.Na_decomp_rate[hh+1]*self.avail_N_pool[hh+1] +\
                                        self.Nr_decomp_rate[hh+1]*self.resist_N_pool[hh+1]

                ## TAN pool
                self.TAN_pool[hh+1] = self.TAN_pool[hh] + self.TAN_prod[hh+1] + self.TAN[hh+1]

                ## update urea/UA pool and orgN pools
                self.urea_pool[hh+1] = self.urea_pool[hh+1] * (1-self.urea_hydro_rate[hh+1])
                self.UA_pool[hh+1] = self.UA_pool[hh+1] * (1-self.ua_conv_factor[hh+1])
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1]*(1 - self.Na_decomp_rate[hh+1]) 
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1]* (1 - self.Nr_decomp_rate[hh+1])

                ## indoor env scheme
                if mms_cat == "MMS_indoor":
                    ## water amount in "solid" manure
                    ## self.manure refers to the DM mass; therefore, total manure mass = DM mass/DM%
                    ## water in the "solid" manure = water% x total manure mass
                    ## maximum water that can be held by manure
                    maxwater = self.manure_pool[hh+1] * absorb_factor
                    ## referenced water content for solid manure
                    stdwater = self.manure_pool[hh+1]/(self.DM_content/100) - self.manure_pool[hh+1]
                    ## water pool
                    ## Note the difference between the water pool of [MMS indoor solid] and [MMS indoor liquid]
                    ## water pool of [MMS indoor solid] is directly determined by the amount of manure as we assumed a dry matter content 
                    self.Total_water_pool[hh+1] = self.Total_water_pool[hh] + self.water[hh+1] - self.evap_sim[hh+1]
                    ## water content is controlled to be between the maximum and standard water content
                    self.Total_water_pool[hh+1] = np.maximum(self.Total_water_pool[hh+1],stdwater)
                    self.Total_water_pool[hh+1] = np.minimum(self.Total_water_pool[hh+1],maxwater)  

                    ## water content of the manure
                    vtotal,manurewc,manure_WFPS = manure_properties(solidmass=self.manure_pool[hh+1],
                                                            watermass=self.Total_water_pool[hh+1])
                    manurewc = np.minimum(manurewc,manure_porosity)
                    manure_WFPS = np.minimum(manure_WFPS,1.0)

                    ## manure resistance
                    ## manure resistance is determined by: R = z/(2*tor_manure*D); 
                    ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                    ##        tortuosity dependence for aqueous diffusion is not considered here
                    ## layer thickness of manure (DM+water) in meter: z = Vmanure
                    z_total = vtotal
                    self.R_manurel[hh+1] = diff_resistance(distance=z_topmanure/2,phase='aqueous',theta_sat=manure_porosity,
                                                            theta=manurewc,temp=self.T_sim[hh+1])
                    self.R_manureg[hh+1] = diff_resistance(distance=z_topmanure/2,phase='gaseous',theta_sat=manure_porosity,
                                                            theta=manurewc,temp=self.T_sim[hh+1])
                    ## resistance
                    self.R_star[hh+1] = 1/k_gas_NH3(temp=self.T_sim[hh+1],u=self.u_sim[hh+1],Z=2,zo=z_manuresurf)

                    ## TAN conc; g/m3
                    ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to manure) NH4+
                    ## gaseous NH3 in air-filled pore space; epsilon(manure porosity) - theta
                    ## aqueous TAN in water-filled pore space; theta
                    ## solid phase TAN adsorbed on solid particles; 1 - epsilon
                    ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd_manure = 1.0 m^3/m^3 
                    ## Kd_manure = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                    ## manure density varies, ~ 0.3-1.9 g/cm^3, we assmume the porosity of mamnure is 40%
                    ## the volume (thickness) is determined by solid mass and bulk density (derived from porosity)
                    ## NH3(g) = KNH3*[TAN(aq)]
                    ## MTAN = Vmanure*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
                    ## so, [TAN(aq)] = MTAN/Vmanure * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
                    KNH3 = NH3_par_coeff(temp=self.T_sim[hh+1],cncH=self.cc_H)
                    # self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/(z_total*(manurewc+KNH3*(manure_porosity-manurewc)+\
                    #     (1-manure_porosity)*Kd))
                    # self.TAN_amount[hh+1] = TAN_concentration(mtan=self.TAN_pool[hh+1],zlayer=z_total,
                    #                                     theta_sat=manure_porosity,theta=manurewc,knh3=KNH3,kd=Kd_manure)
                    # self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/(self.Total_water_pool[hh+1]+self.manure_pool[hh+1])
                    self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/(self.Total_water_pool[hh+1]+self.manure_pool[hh+1]*Kd_manure/(manure_PD/1e3))
                    self.TAN_amount[hh+1] = np.nan_to_num(self.TAN_amount[hh+1],posinf=0.0,neginf=0.0)
                    self.TAN_amount[hh+1] = self.TAN_amount[hh+1] * 1e6
                    ## NH3 conc in the bulk manure
                    self.NH3_gas_bulk[hh+1] = KNH3 * self.TAN_amount[hh+1]
                    ## NH3 conc is 0 when manure water content equals to porosity
                    # self.NH3_gas_bulk[hh+1][manurewc==manure_porosity] = 0.0
                    correction_gas = np.nan_to_num((manure_porosity-manurewc)/(manure_porosity-manurewc),posinf=0.0,neginf=0.0)
                    self.NH3_gas_bulk[hh+1] = self.NH3_gas_bulk[hh+1] * correction_gas

                    ## TAN conc at the surface (g/m3)
                    # self.TAN_surf_amount[hh+1] = (self.TAN_amount[hh+1]*(1/self.R_manurel[hh+1]+KNH3/self.R_manureg[hh+1])/\
                    #         (KNH3*(1/self.R_star[hh+1]+1/self.R_manureg[hh+1])+1/self.R_manurel[hh+1]))
                    self.TAN_surf_amount[hh+1] = surf_TAN_cnc(tan_cnc=self.TAN_amount[hh+1],rliq=self.R_manurel[hh+1],
                                                            rgas=self.R_manureg[hh+1],knh3=KNH3,ratm=self.R_star[hh+1],qrunoff=0.0)
                    ## Gaseous NH3 at the surface
                    self.NH3_gas_surf[hh+1] = self.TAN_surf_amount[hh+1]*KNH3

                    ## determining the maximum emission; 
                    emiss_idx = (self.NH3_gas_surf[hh+1]/self.R_star[hh+1])*timestep*3600
                    ## fraction of [NH4+(aq)]
                    fNH4 = frac_NH4(theta=manurewc,theta_sat=manure_porosity,
                                            temp=self.T_sim[hh+1],cncH=self.cc_H,kd=Kd_manure)
                    ## nirification TAN in the bulk manure; 
                    # nitrif_idx = TAN_nitrif(tan_pool=self.TAN_pool[hh+1]*(1-f_manure_anoxic),temp=self.T_sim[hh+1],
                    #                         theta=manurewc,theta_sat=manure_porosity,pH=self.pH,
                    #                         fert_type='manure',frac_nh4=fNH4)*timestep*3600
                    nitrif_idx = (self.TAN_pool[hh+1]*fNH4)*nitrification_rate_manure(manure_temp=self.T_sim[hh+1],
                                                                        WFPS=manure_WFPS,pH=self.pH)*timestep*3600
                    ## fluxes from the bulk manure to the soil interface
                    ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                    otherflux1 = False
                    otherflux2 = False
                    emiss_idx,nitrif_idx,otherflux1,otherflux2 = N_pathways(mN=self.TAN_pool[hh+1],
                                            flux1=emiss_idx,flux2=nitrif_idx,flux3=otherflux1,flux4=otherflux2)

                    self.NH3flux[hh+1] = emiss_idx
                    self.TANnitrif_manure[hh+1] = nitrif_idx

                    ## update TAN pool
                    self.TAN_pool[hh+1] = self.TAN_pool[hh+1] - self.NH3flux[hh+1] - self.TANnitrif_manure[hh+1]
                    self.TAN_pool[hh+1] = np.maximum(self.TAN_pool[hh+1],0.0)
                    ## NO3- pool of the bulk manure
                    self.NO3_pool[hh+1] = self.NO3_pool[hh] + self.TANnitrif_manure[hh+1]

                elif mms_cat == "MMS_open":
                    ## water pool
                    ## water pool of [MMS open solid] is determined by:
                    ##   source: manure water, rain (precipitation)
                    ##   loss: evaporation, infiltration to soil/interface
                    ## minimum water amount is equivalent to the mositure equilibrium content of the manure  
                    ## Note: infiltration of manure water to soil is assumed to be 10mm/day (Vira et al., 2020GMD)
                    ##       and this should be differentiate with TAN infiltration to soil
                    ## justify the water amount;
                    maxwater = self.manure_pool[hh+1] * absorb_factor
                    stdwater = self.manure_pool[hh+1]/(self.DM_content/100) - self.manure_pool[hh+1]
                    self.Total_water_pool[hh+1] = self.Total_water_pool[hh+1] + self.water[hh+1] 
                    # water_idx = self.Total_water_pool[hh+1] - maxwater
                    # self.Total_water_pool[hh+1][water_idx>0] = maxwater[water_idx>0]
                    self.Total_water_pool[hh+1] = np.minimum(self.Total_water_pool[hh+1],maxwater)
                    self.Total_water_pool[hh+1] = self.Total_water_pool[hh+1] + self.rainfall[hh+1] - self.evap_sim[hh+1]
                    # water_idx = self.Total_water_pool[hh+1] - maxwater
                    # self.rain_avail_washoff[hh+1][water_idx>0] = water_idx[water_idx>0]/(timestep*3600)
                    self.rain_avail_washoff[hh+1] = (self.Total_water_pool[hh+1] - maxwater)/(timestep*3600)
                    self.rain_avail_washoff[hh+1] = np.maximum(self.rain_avail_washoff[hh+1],0)
                    # self.Total_water_pool[hh+1][water_idx>0] = maxwater[water_idx>0] 
                    # water_idx = self.Total_water_pool[hh+1] - stdwater
                    # self.Total_water_pool[hh+1][water_idx<0] = stdwater[water_idx<0]
                    self.Total_water_pool[hh+1] = np.maximum(self.Total_water_pool[hh+1],stdwater)

                    ## water content of the manure
                    vtotal,manurewc,manure_WFPS = manure_properties(solidmass=self.manure_pool[hh+1],
                                                            watermass=self.Total_water_pool[hh+1])
                    # manurewc[manurewc>manure_porosity] = manure_porosity
                    manurewc = np.minimum(manurewc,manure_porosity)
                    # manure_WFPS[manure_WFPS>1.0] = 1.0
                    manure_WFPS = np.minimum(manure_WFPS,1.0)
                    ## manure resistance
                    ## manure resistance is determined by: R = z/(2*tor_manure*D); 
                    ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                    ##        tortuosity dependence for aqueous diffusion is not considered here
                    ## layer thickness of manure (DM+water) in meter: z = Vmanure
                    z_total = vtotal
                    # print(dd,hh+1,"manure wc", manurewc[192,386])
                    # print(dd,hh+1,"z total", z_total[192,386])
                    
                    self.R_manurel[hh+1] = diff_resistance(distance=z_topmanure/2,phase='aqueous',theta_sat=manure_porosity,
                                                            theta=manurewc,temp=self.T_sim[hh+1])
                    self.R_manureg[hh+1] = diff_resistance(distance=z_topmanure/2,phase='gaseous',theta_sat=manure_porosity,
                                                            theta=manurewc,temp=self.T_sim[hh+1])
                    
                    ## washoff: 1) manure, 2) urea, 3) available org N, 4) reistant org N, 5) TAN
                    ## washoff coefficient (%) = washoff water (mm) * washoff (%/mm)
                    nonN_washoff_rate = self.rain_avail_washoff[hh+1]*f_washoff_nonN/1e3
                    N_washoff_rate = self.rain_avail_washoff[hh+1]*f_washoff_N/1e3
                    self.manure_washoff[hh+1] = nonN_washoff_rate*self.manure_pool[hh+1]*timestep*3600
                    self.UA_washoff[hh+1] = N_washoff_rate*self.UA_pool[hh+1]*timestep*3600
                    self.avail_N_washoff[hh+1] = N_washoff_rate*self.avail_N_pool[hh+1]*timestep*3600
                    self.resist_N_washoff[hh+1] = N_washoff_rate*self.resist_N_pool[hh+1]*timestep*3600
                    self.unavail_N_washoff[hh+1] = N_washoff_rate*self.unavail_N_pool[hh+1]*timestep*3600

                    urea_amount = N_concentration(mN=self.urea_pool[hh+1],zlayer=z_total,theta=manurewc)
                    ureasurfamount = surf_Ncnc(N_cnc=urea_amount,rliq=self.R_manurel[hh+1],qrunoff=self.rain_avail_washoff[hh+1])
                    # ureasurfamount[np.isnan(ureasurfamount)] = 0.0
                    ureasurfamount = np.nan_to_num(ureasurfamount,posinf=0.0,neginf=0.0)
                    # self.ureawashoff[hh+1] = ureasurfamount*self.rain_avail_washoff[hh+1]*timestep*3600
                    self.ureawashoff[hh+1] = surf_runoff(N_surfcnc=ureasurfamount,
                                                        qrunoff=self.rain_avail_washoff[hh+1])*timestep*3600
                    # washoff_idx = self.urea_pool[hh+1] - self.ureawashoff[hh+1]                                  
                    # self.ureawashoff[hh+1][washoff_idx<0] = self.urea_pool[hh+1][washoff_idx<0]
                    self.ureawashoff[hh+1] = np.minimum(self.ureawashoff[hh+1],self.urea_pool[hh+1])
                    # print(dd,hh+1,"urea conc", urea_amount[192,386])
                    # print(dd,hh+1,"urea surf conc", ureasurfamount[192,386])

                    ## update manure pool: subtracting washoff
                    self.manure_pool[hh+1] = self.manure_pool[hh+1] - self.manure_washoff[hh+1]
                    ## update urea pool and orgN pools
                    self.urea_pool[hh+1] = self.urea_pool[hh+1]  - self.ureawashoff[hh+1]
                    self.UA_pool[hh+1] = self.UA_pool[hh+1] - self.UA_washoff[hh+1]
                    self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1] - self.avail_N_washoff[hh+1]
                    self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1] - self.resist_N_washoff[hh+1]

                    ## TAN conc; g/m3
                    ## TAN will partitioned into gaseous NH3, aqueous and solid (adsorption to manure) NH4+
                    ## gaseous NH3 in air-filled pore space; epsilon(manure porosity) - theta
                    ## aqueous TAN in water-filled pore space; theta
                    ## solid phase TAN adsorbed on solid particles; 1 - epsilon
                    ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd_manure = 1.0 m^3/m^3 
                    ## Kd_manure = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                    ## manure density varies, ~ 0.3-1.9 g/cm^3, we assmume the porosity of mamnure is 40%
                    ## the volume (thickness) is determined by solid mass and bulk density (derived from porosity)
                    ## NH3(g) = KNH3*[TAN(aq)]
                    ## MTAN = Vmanure*(theta*[TAN(aq)]+(epsilon-theta)*NH3(g)+(1-epsilon)*[TAN(s)])
                    ## so, [TAN(aq)] = MTAN/Vmanure * (1/(theta+KNH3*(epsilon-theta)+Kd*(1-epsilon)))
                    KNH3 = NH3_par_coeff(temp=self.T_sim[hh+1],cncH=self.cc_H)
                    # self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/(z_total*(manurewc+KNH3*(manure_porosity-manurewc)+\
                    #     (1-manure_porosity)*Kd))
                    # self.TAN_amount[hh+1] = TAN_concentration(mtan=self.TAN_pool[hh+1],zlayer=z_total,
                    #                                     theta_sat=manure_porosity,theta=manurewc,knh3=KNH3,kd=Kd_manure)
                    # self.TAN_amount[hh+1][self.Total_water_pool[hh+1]==0] = 0
                    self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/(self.Total_water_pool[hh+1]+self.manure_pool[hh+1]*Kd_manure/(manure_PD/1e3))
                    self.TAN_amount[hh+1] = np.nan_to_num(self.TAN_amount[hh+1],posinf=0.0,neginf=0.0)
                    self.TAN_amount[hh+1] = self.TAN_amount[hh+1] * 1e6
                    ## NH3 conc in the bulk manure
                    self.NH3_gas_bulk[hh+1] = KNH3 * self.TAN_amount[hh+1]
                    ## NH3 conc is 0 when manure water content equals to porosity
                    # self.NH3_gas_bulk[hh+1][manurewc==manure_porosity] = 0.0
                    correction_gas = np.nan_to_num((manure_porosity-manurewc)/(manure_porosity-manurewc),posinf=0.0,neginf=0.0)
                    self.NH3_gas_bulk[hh+1] = self.NH3_gas_bulk[hh+1] * correction_gas

                    ## TAN conc at the surface (solved); atmospheric NH3 is ignored
                    # self.TAN_surf_amount[hh+1] = (self.TAN_amount[hh+1]*(1/self.R_manurel[hh+1]+KNH3/self.R_manureg[hh+1])/\
                    #         (self.rain_avail_washoff[hh+1]+KNH3*(1/self.R_star[hh+1]+1/self.R_manureg[hh+1])+1/self.R_manurel[hh+1]))
                    self.TAN_surf_amount[hh+1] = surf_TAN_cnc(tan_cnc=self.TAN_amount[hh+1],rliq=self.R_manurel[hh+1],
                                                            rgas=self.R_manureg[hh+1],knh3=KNH3,ratm=self.R_atm[hh+1],
                                                            qrunoff=self.rain_avail_washoff[hh+1])
                    ## Gaseous NH3 at the surface
                    self.NH3_gas_surf[hh+1] = self.TAN_surf_amount[hh+1]*KNH3

                    # print(dd,hh+1,"TAN pool", self.TAN_pool[hh+1,192,386])
                    # print(dd,hh+1,"water pool", self.Total_water_pool[hh+1,192,386])
                    # print(dd,hh+1,"TAN conc", self.TAN_amount[hh+1,192,386])
                    # print(dd,hh+1,"TAN surf conc", self.TAN_surf_amount[hh+1,192,386])
                    # print(dd,hh+1,"NH3 surf conc", self.NH3_gas_surf[hh+1,192,386])

                    ## determining the maximum emission; 
                    emiss_idx = (self.NH3_gas_surf[hh+1]/self.R_atm[hh+1])*timestep*3600
                    ## determining the maximum TAN runoff;
                    ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
                    ## else need to multiply by timestep*3600
                    runoff_idx = self.rain_avail_washoff[hh+1]*self.TAN_surf_amount[hh+1]*timestep*3600
                    
                    ## fraction of [NH4+(aq)]
                    fNH4 = frac_NH4(theta=manurewc,theta_sat=manure_porosity,
                                            temp=self.T_sim[hh+1],cncH=self.cc_H,kd=Kd_manure)
                    ## nirification TAN in the bulk manure; 
                    # nitrif_idx = TAN_nitrif(tan_pool=self.TAN_pool[hh+1]*(1-f_manure_anoxic),temp=self.T_sim[hh+1],
                    #                         theta=manurewc,theta_sat=manure_porosity,pH=self.pH,
                    #                         fert_type='manure',frac_nh4=fNH4)*timestep*3600
                    nitrif_idx = (self.TAN_pool[hh+1]*fNH4)*nitrification_rate_manure(manure_temp=self.T_sim[hh+1],
                                                                        WFPS=manure_WFPS,pH=self.pH)*timestep*3600
                    # print(dd,hh+1,"NH3 flux idx", emiss_idx[192,386])

                    ## fluxes from the bulk manure to the soil interface
                    ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                    otherflux = False
                    emiss_idx,runoff_idx,nitrif_idx,otherflux = N_pathways(mN=self.TAN_pool[hh+1],
                                            flux1=emiss_idx,flux2=runoff_idx,flux3=nitrif_idx,flux4=otherflux)
                    self.NH3flux[hh+1] = emiss_idx
                    self.TANwashoff[hh+1] = runoff_idx
                    self.TANnitrif_manure[hh+1] = nitrif_idx

                    # print(dd,hh+1,"NH3 flux", self.NH3flux[hh+1,192,386])

                    ## update TAN pool
                    self.TAN_pool[hh+1] = self.TAN_pool[hh+1] - self.NH3flux[hh+1] - self.TANwashoff[hh+1] - self.TANnitrif_manure[hh+1]
                    self.TAN_pool[hh+1] = np.maximum(self.TAN_pool[hh+1],0.0)
                    ## NO3- pool of the bulk manure
                    self.NO3_pool[hh+1] = self.NO3_pool[hh+1] + self.TANnitrif_manure[hh+1]
                    ## NO3- conc of the bulk manure; g/m3
                    self.NO3_amount[hh+1] = N_concentration(mN=self.NO3_pool[hh+1],zlayer=z_total,theta=manurewc)
                    NO3surfamount = surf_Ncnc(N_cnc=self.NO3_amount[hh+1],rliq=self.R_manurel[hh+1]/f_DNO3,qrunoff=self.rain_avail_washoff[hh+1])
                    self.NO3washoff[hh+1] = NO3surfamount*self.rain_avail_washoff[hh+1]*timestep*3600
                    self.NO3_pool[hh+1] = self.NO3_pool[hh+1] - self.NO3washoff[hh+1]
            if dd >= Days:
                ddidx = dd - Days
            else:
                ddidx = dd
            if mms_cat == "MMS_indoor":
                self.daily_output(dayidx=ddidx,mms_info="mms_indoor_solid")
            elif mms_cat == "MMS_open":
                self.daily_output(dayidx=ddidx,mms_info="mms_open_solid")
            self.mms_to_landspreading(dayidx=ddidx)
        return
  
    ## Simulation: Cat B.3 manure stored in open environment (as solid) //// under development 08/Sep
    def MMS_open_lagoon_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: MMS open (lagoon,liquid), '+str(self.livestock)+', '+str(self.production_system))
        f_MMSarea = self.f_MMS_open_lagoon/self.f_MMS
        MMS_area["mms_open_lagoon_area"] = self.housingarea*(1.0-self.f_loss-self.f_sold)*f_MMSarea*MMS_area_factor["mms_open_lagoon"]
        
        if self.livestock.lower()=="poultry":
            print("Simulation for poultry.")

        for dd in np.arange(start_day_idx,end_day_idx):
            self.sim_env(mms_type="MMS_open",mms_phase="liquid",dayidx=dd)
            self.daily_init(dayidx=dd,mms_info="mms_open_liquid")
            for hh in np.arange(24):
                ## NH3 conc in g/m3
                KNH3 = NH3_par_coeff(temp=self.T_sim[hh+1],cncH=self.cc_H)
                ## resistance
                self.R_star[hh+1] = (1/k_aq_NH4(temp=self.T_sim[hh+1]))+(1/(KNH3*k_gas_NH3(temp=self.T_sim[hh+1],
                                                                    u=self.u_sim[hh+1],Z=2,zo=zo_water)))
                ## 
                self.modelled_emiss[hh+1] = NH3_volpit(pit_tanconc=lagoon_TAN_conc*1e6,conc_in=0.0,
                                                Rpit=self.R_star[hh+1])*timestep*3600

                self.NH3flux[hh+1] = self.modelled_emiss[hh+1]
            if dd >= Days:
                ddidx = dd - Days
            else:
                ddidx = dd
            self.daily_output(dayidx=ddidx,mms_info="mms_open_lagoon")
        return

    def daily_init(self,dayidx,mms_info):
        self.manure_pool[0] = self.manure_pool[-1]
        self.avail_N_pool[0] = self.avail_N_pool[-1]
        self.resist_N_pool[0] = self.resist_N_pool[-1]
        self.unavail_N_pool[0] = self.unavail_N_pool[-1]
        self.urea_pool[0] = self.urea_pool[-1]
        self.UA_pool[0] = self.UA_pool[-1]
        self.TAN_pool[0] = self.TAN_pool[-1]
        self.Total_water_pool[0] = self.Total_water_pool[-1]
        if mms_info == "mms_open_solid":
            self.NO3_pool[0] = self.NO3_pool[-1]
        if mms_info == "mms_indoor_solid":
            self.NO3_pool[0] = self.NO3_pool[-1]

        if dayidx >= Days:
            dayidx = dayidx - Days
        ## inputs from housing (added in the mid-day)
        ## the equations used to represent these pools need to be explained here:
        ## each pool should be in a unit of, mass/unit area, i.e., g/m^2
        ## and each MMS corresponds to a specific MMS area, so the explicit equation would be:
        ##     self.[pool] = self.[poolmass]*(1.0-f_loss-f_sold)*f_MMS_[MMS type]/MMS_area["MMS type"]
        ## then, note that  
        ##     MMS_area["MMS type"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_[MMS type]*MMS_area_factor["MMS type"]
        ##     self.[pool] = self.[poolmass]/(self.housingarea*MMS_area_factor["MMS_type"])
        ## this is different to the HOUSING module as N excretion has been divided by the housing area before used as input
        self.manure[12] = self.manure_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        ## N input in multiple forms
        self.urea[12] = self.urea_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        self.UA[12] = self.UA_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        self.avail_N[12] = self.avail_N_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        self.resist_N[12] = self.resist_N_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        self.unavail_N[12] = self.unavail_N_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        ## TAN from housing to storage
        self.TAN[12] = self.TAN_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        ## water from housing to storage (this is unconstrained)
        self.water[12] = self.water_added[dayidx]/(self.housingarea*MMS_area_factor[mms_info])
        return

    def daily_output(self,dayidx,mms_info):
        # if dayidx >= Days:
        #     dayidx = dayidx - Days
        self.o_NH3flux[dayidx] = np.nansum(self.NH3flux[1:],axis=0)
        if mms_info == "mms_indoor_solid":
            self.o_NH4nitrif[dayidx] = np.nansum(self.TANnitrif_manure[1:],axis=0)
        if mms_info == "mms_open_solid":
            self.o_Nwashoff[dayidx] = np.nansum(self.TANwashoff[1:],axis=0)+np.nansum(self.ureawashoff[1:],axis=0)+\
                            np.nansum(self.UA_washoff[1:],axis=0)+\
                            np.nansum(self.avail_N_washoff[1:],axis=0)+np.nansum(self.resist_N_washoff[1:],axis=0)
            self.o_NH4nitrif[dayidx] = np.nansum(self.TANnitrif_manure[1:],axis=0)

    def manure_app_calendar(self):
        manureappcalpath = infile_path+crop_data_path+manure_appcalendar
        appcalds = open_ds(manureappcalpath)
        self.spring_app_cal = appcalds.spring_plant_mean.values
        self.winter_app_cal = appcalds.winter_plant_mean.values
        idx = (self.spring_app_cal!=0)&(self.winter_app_cal==0)&(np.isnan(self.winter_app_cal))
        self.winter_app_cal[idx] = self.spring_app_cal[idx] + int(180)

    def mms_to_landspreading(self,dayidx):
        ## spring planting season
        self.spring_Napp[self.spring_app_cal==dayidx] = (self.TAN_pool[-1][self.spring_app_cal==dayidx]+\
                                                        self.urea_pool[-1][self.spring_app_cal==dayidx])*\
                                                        self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_UANapp[self.spring_app_cal==dayidx] = self.UA_pool[-1][self.spring_app_cal==dayidx]*\
                                                            self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_manureapp[self.spring_app_cal==dayidx] = self.manure_pool[-1][self.spring_app_cal==dayidx]*\
                                                            self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_waterapp[self.spring_app_cal==dayidx] = self.Total_water_pool[-1][self.spring_app_cal==dayidx]*\
                                                            self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_availNapp[self.spring_app_cal==dayidx] = self.avail_N_pool[-1][self.spring_app_cal==dayidx]*\
                                                            self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_resistNapp[self.spring_app_cal==dayidx] = self.resist_N_pool[-1][self.spring_app_cal==dayidx]*\
                                                            self.mmsarea[self.spring_app_cal==dayidx]
        self.spring_unavailNapp[self.spring_app_cal==dayidx] = self.unavail_N_pool[-1][self.spring_app_cal==dayidx]*\
                                                                self.mmsarea[self.spring_app_cal==dayidx]

        self.manure_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.Total_water_pool[-1][self.spring_app_cal==dayidx]  = 0.0
        self.avail_N_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.resist_N_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.unavail_N_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.urea_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.UA_pool[-1][self.spring_app_cal==dayidx] = 0.0
        self.TAN_pool[-1][self.spring_app_cal==dayidx]  = 0.0
        
        ## winter planting season
        self.winter_Napp[self.winter_app_cal==dayidx] = (self.TAN_pool[-1][self.winter_app_cal==dayidx]+\
                                                        self.urea_pool[-1][self.winter_app_cal==dayidx])*\
                                                        self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_UANapp[self.winter_app_cal==dayidx] = self.UA_pool[-1][self.winter_app_cal==dayidx]*\
                                                        self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_manureapp[self.winter_app_cal==dayidx] = self.manure_pool[-1][self.winter_app_cal==dayidx]*\
                                                        self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_waterapp[self.winter_app_cal==dayidx] = self.Total_water_pool[-1][self.winter_app_cal==dayidx]*\
                                                        self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_availNapp[self.winter_app_cal==dayidx] = self.avail_N_pool[-1][self.winter_app_cal==dayidx]*\
                                                        self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_resistNapp[self.winter_app_cal==dayidx] = self.resist_N_pool[-1][self.winter_app_cal==dayidx]*\
                                                            self.mmsarea[self.winter_app_cal==dayidx]
        self.winter_unavailNapp[self.winter_app_cal==dayidx] = self.unavail_N_pool[-1][self.winter_app_cal==dayidx]*\
                                                            self.mmsarea[self.winter_app_cal==dayidx]

        self.manure_pool[-1][self.winter_app_cal==dayidx] = 0.0
        self.Total_water_pool[-1][self.winter_app_cal==dayidx]  = 0.0
        self.avail_N_pool[-1][self.winter_app_cal==dayidx] = 0.0
        self.resist_N_pool[-1][self.winter_app_cal==dayidx]  = 0.0
        self.unavail_N_pool[-1][self.winter_app_cal==dayidx]  = 0.0
        self.urea_pool[-1][self.winter_app_cal==dayidx]  = 0.0
        self.UA_pool[-1][self.winter_app_cal==dayidx]  = 0.0
        self.TAN_pool[-1][self.winter_app_cal==dayidx]  = 0.0
    
    def sim_out(self,mms_cat,phase):
        nlat = int(180.0/dlat)
        nlon = int(360.0/dlon)
        ntime = Days
        lats = 90 - 0.5*np.arange(nlat)
        lons = -180 + 0.5*np.arange(nlon)
        yearidx = str(sim_year)+'-01-01'
        times = pd.date_range(yearidx,periods=ntime)
        if phase == "liquid":
            MMS_NH3emiss = self.o_NH3flux*self.mmsarea
            MMS_totalN = self.total_N_MMS*self.mmsarea/self.housingarea
            outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss=(['time','lat','lon'],MMS_NH3emiss),
                    NtotalMMS=(['lat','lon'],MMS_totalN),
                    sourcearea=(['lat','lon'],self.mmsarea),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-MMS: NH3 emissions from "+\
                            self.production_system+" "+self.livestock+" "+mms_cat+" "+phase+" in " +str(sim_year),
                    info = self.production_system+" "+self.livestock+" "+mms_cat+" "+phase,
                    units="gN per grid",
                ),
            )
            outds.NH3emiss.attrs["unit"] = 'gN/day'
            outds.NH3emiss.attrs["long name"] = 'NH3 emission from '+mms_cat+" "+phase
            outds.NtotalMMS.attrs["unit"] = 'gN'
            outds.NtotalMMS.attrs["long name"] = 'Total N of '+mms_cat+" "+phase
            outds.sourcearea.attrs["unit"] = 'm2'
            outds.sourcearea.attrs["long name"] = 'Source area for NH3 emission from '+mms_cat+" "+phase

        elif phase == "solid":
            MMS_NH3emiss = self.o_NH3flux*self.mmsarea
            MMS_NH4nitrif = self.o_NH4nitrif*self.mmsarea
            MMS_totalN = self.total_N_MMS*self.mmsarea/self.housingarea
            if mms_cat == "MMS_open":
                MMS_Nwashoff = self.o_Nwashoff*self.mmsarea
                outds = xr.Dataset(
                    data_vars=dict(
                        NH3emiss=(['time','lat','lon'],MMS_NH3emiss),
                        NH4nitrif=(['time','lat','lon'],MMS_NH4nitrif),
                        Nwashoff=(['time','lat','lon'],MMS_Nwashoff),
                        NtotalMMS=(['lat','lon'],MMS_totalN),
                        sourcearea=(['lat','lon'],self.mmsarea),
                                ),
                    coords = dict(
                        time=(["time"], times),
                        lon=(["lon"], lons),
                        lat=(["lat"], lats),
                                ),
                    attrs=dict(
                        description="AMCLIM-MMS: NH3 emissions from "+\
                                self.production_system+" "+self.livestock+" "+mms_cat+" "+phase+" in " +str(sim_year),
                        info = self.production_system+" "+self.livestock+" "+mms_cat+" "+phase,
                        units="gN per grid",
                    ),
                )

                outds.NH3emiss.attrs["unit"] = 'gN/day'
                outds.NH3emiss.attrs["long name"] = 'NH3 emission from '+mms_cat+" "+phase
                outds.NH4nitrif.attrs["unit"] = 'gN/day'
                outds.NH4nitrif.attrs["long name"] = 'Nitrification from '+mms_cat+" "+phase
                outds.NH4nitrif.attrs["unit"] = 'gN/day'
                outds.NH4nitrif.attrs["long name"] = 'N washoff from '+mms_cat+" "+phase
                outds.NtotalMMS.attrs["unit"] = 'gN'
                outds.NtotalMMS.attrs["long name"] = 'Total N of '+mms_cat+" "+phase
                outds.sourcearea.attrs["unit"] = 'm2'
                outds.sourcearea.attrs["long name"] = 'Source area for NH3 emission from '+mms_cat+" "+phase

            else:
                outds = xr.Dataset(
                    data_vars=dict(
                        NH3emiss=(['time','lat','lon'],MMS_NH3emiss),
                        NH4nitrif=(['time','lat','lon'],MMS_NH4nitrif),
                        NtotalMMS=(['lat','lon'],MMS_totalN),
                        sourcearea=(['lat','lon'],self.mmsarea),
                                ),
                    coords = dict(
                        time=(["time"], times),
                        lon=(["lon"], lons),
                        lat=(["lat"], lats),
                                ),
                    attrs=dict(
                        description="AMCLIM-MMS: NH3 emissions from "+\
                                self.production_system+" "+self.livestock+" "+mms_cat+" "+phase+" in " +str(sim_year),
                        info = self.production_system+" "+self.livestock+" "+mms_cat+" "+phase,
                        units="gN per grid",
                    ),
                )

                outds.NH3emiss.attrs["unit"] = 'gN/day'
                outds.NH3emiss.attrs["long name"] = 'NH3 emission from '+mms_cat+" "+phase
                outds.NH4nitrif.attrs["unit"] = 'gN/day'
                outds.NH4nitrif.attrs["long name"] = 'Nitrification from '+mms_cat+" "+phase
                outds.NtotalMMS.attrs["unit"] = 'gN'
                outds.NtotalMMS.attrs["long name"] = 'Total N of '+mms_cat+" "+phase
                outds.sourcearea.attrs["unit"] = 'm2'
                outds.sourcearea.attrs["long name"] = 'Source area for NH3 emission from '+mms_cat+" "+phase

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in outds.data_vars}

        outds.to_netcdf(output_path+self.livestock+'.'+self.production_system+'.'+mms_cat+'.'+phase+\
                            '.'+str(sim_year)+'.nc',encoding=encoding)
        print("ncfile saved.")

        outds = xr.Dataset(
                data_vars=dict(
                    spring_TAN=(['lat','lon'],self.spring_Napp),
                    spring_UAN=(['lat','lon'],self.spring_UANapp),
                    spring_manure=(['lat','lon'],self.spring_manureapp),
                    spring_water=(['lat','lon'],self.spring_waterapp),
                    spring_availN=(['lat','lon'],self.spring_availNapp),
                    spring_resistN=(['lat','lon'],self.spring_resistNapp),
                    spring_unavailN=(['lat','lon'],self.spring_unavailNapp),
                    winter_TAN=(['lat','lon'],self.winter_Napp),
                    winter_UAN=(['lat','lon'],self.winter_UANapp),
                    winter_manure=(['lat','lon'],self.winter_manureapp),
                    winter_water=(['lat','lon'],self.winter_waterapp),
                    winter_availN=(['lat','lon'],self.winter_availNapp),
                    winter_resistN=(['lat','lon'],self.winter_resistNapp),
                    winter_unavailN=(['lat','lon'],self.winter_unavailNapp),
                            ),
                coords = dict(
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-MMS: manure application to land in spring and winter season "+\
                            self.production_system+" "+self.livestock+" "+mms_cat+" "+phase+" in " +str(sim_year),
                    info = self.production_system+" "+self.livestock+" "+mms_cat+" "+phase,
                    units="gN per grid",
                ),
            )

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in outds.data_vars}

        outds.to_netcdf(output_path+self.livestock+'.'+self.production_system+'.'+mms_cat+'.'+phase+\
                            '.'+str(sim_year)+'.manureapp.nc',encoding=encoding)
        print("manure app ncfile saved.")
        return

    def output_MMS_pathway(self):
        nlat = int(180.0/dlat)
        nlon = int(360.0/dlon)
        ntime = Days
        lats = 90 - 0.5*np.arange(nlat)
        lons = -180 + 0.5*np.arange(nlon)
        yearidx = str(sim_year)+'-01-01'
        times = pd.date_range(yearidx,periods=ntime)

        outds = xr.Dataset(
                data_vars=dict(
                    N_loss=(['lat','lon'],self.total_N_MMS*self.f_loss),
                    N_sold=(['lat','lon'],self.total_N_MMS*self.f_sold),
                    N_preserved_liquid=(['lat','lon'],self.total_N_MMS*self.f_MMS_preserve_liquid),
                    N_preserved_solid=(['lat','lon'],self.total_N_MMS*self.f_MMS_preserve_solid),
                    N_fuel=(['lat','lon'],self.total_N_MMS*self.f_MMS_fuel),
                    N_indoor_liquid=(['lat','lon'],self.total_N_MMS*self.f_MMS_indoor_liquid),
                    N_indoor_solid=(['lat','lon'],self.total_N_MMS*self.f_MMS_indoor_solid),
                    N_outdoor_liquid=(['lat','lon'],self.total_N_MMS*self.f_MMS_open_liquid),
                    N_outdoor_solid=(['lat','lon'],self.total_N_MMS*self.f_MMS_open_solid),
                    N_outdoor_lagoon=(['lat','lon'],self.total_N_MMS*self.f_MMS_open_lagoon),
                            ),
                coords = dict(
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-MMS: N pathways to each MMS of "+\
                            self.production_system+" "+self.livestock+" in " +str(sim_year),
                    info = self.production_system+" "+self.livestock,
                    units="gN per grid",
                ),
            )

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in outds.data_vars}

        outds.to_netcdf(output_path+self.livestock+'.'+self.production_system+\
                            '.'+str(sim_year)+'.MMSpathway.nc',encoding=encoding)
        print("manure pathway ncfile saved.")
        return


    def MMS_sim_main(self,mms_cat,phase,start_day_idx,end_day_idx):
        spinup = np.floor((end_day_idx - Days)/Days)
        print("Spinup year is: "+str(spinup)+" yr(s)")
        self.manure_app_calendar()
        if phase == "liquid":
            self.MMS_liquid_sim(mms_cat,start_day_idx,end_day_idx)
        elif phase == "solid":
            self.MMS_solid_sim(mms_cat,start_day_idx,end_day_idx)
        else:
            self.MMS_open_lagoon_sim(start_day_idx,end_day_idx)
        self.sim_out(mms_cat,phase)
        return

