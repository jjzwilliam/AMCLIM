
####################################
## import essential AMCLIM modules
####################################
'''file_path = os.getcwd()
module_path = str(Path(file_path).parent)
print(module_path)

if module_path not in sys.path:
    sys.path.append(module_path)'''
from INPUT.input import *
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *


######################################
## excreted N data
######################################
excretN_info = animal_file['Excreted_N'][lvl_idx]
animal_head = animal_file['Animal_head'][lvl_idx]
animal_weight = animal_file['Animal_weight'][lvl_idx]

animal_weight.values[np.where((animal_head!=0)&(animal_weight==0)&(~np.isnan(animal_head)))] = np.nanmedian(animal_weight.values[np.where(animal_weight!=0)])
animal_weight.values[np.where((animal_head!=0)&(np.isnan(animal_weight))&(~np.isnan(animal_head)))] = np.nanmedian(animal_weight.values[np.where(animal_weight!=0)])
## pig density 1: 120 kg/m^2
#animal_density = 120.0
animal_density = stocking_desity[livestock]
# print(str(livestock)+" stocking density is "+str(animal_density)+" kg/m^2")
massgrid = animal_head*animal_weight.values
housing_area = animal_head*animal_weight.values/animal_density

F_Ncorrect = 1.0
excret_N = F_Ncorrect*excretN_info/housing_area   ## unit: kg N per m^2 per year
excret_N = xr_to_np(excret_N)

###################################
## housing parameters
###################################
## fraction of N_avail, N_resist and N_unavail; 50, 45, 5 per cent
## ref: CLM_FANv1 (Riddick et al., 2016)
f_avail = 0.5
f_resist = 0.45
f_unavail = 0.05
## fraction of UAN in poultry excretion N; 60 per cent
f_uan = 0.6
## fraction of N in poultry excretion; 5 per cent
f_excretn = 0.05
## system pH is assumed to be manure pH varied by animals
pH = pH_info[livestock.upper()]
## surface roughness height of slatted floor; default 2mm
zo_house = 0.002
## surface roughness of water surface (pit storage, mostly liquid)
zo_water = 0.002
## assuming the roughness height of manure storage barn is ~ 0.5m (<ref height of 2m)
zo_barn = 0.5 

## shape in [lat,lon]
f_loss = np.zeros(mtrx[1:])
f_sold = np.zeros(mtrx[1:])
f_housing_litter = np.zeros(mtrx[1:])
f_housing_pit = np.zeros(mtrx[1:])

for mms in loss_list:
    try:
        f_mms = MMS_file[mms][lvl_idx].values
        f_mms[np.isnan(f_mms)] = 0.0
        f_loss = f_loss + f_mms   
    except:pass
for mms in sold_list:
    try:
        f_mms = MMS_file[mms][lvl_idx].values
        f_mms[np.isnan(f_mms)] = 0.0
        f_sold = f_sold + f_mms
    except:pass
for mms in MMS_house_storage_solid_list:
    try:
        f_mms = MMS_file[mms][lvl_idx].values
        f_mms[np.isnan(f_mms)] = 0.0
        f_housing_litter = f_housing_litter + f_mms
    except:pass
for mms in MMS_house_storage_liquid_list:
    try:
        f_mms = MMS_file[mms][lvl_idx].values
        f_mms[np.isnan(f_mms)] = 0.0
        f_housing_pit = f_housing_pit + f_mms
    except:pass


############################################
## HOUSING MODULE: various housing systems
############################################
class HOUSING_MODULE:
    ## Initiating class and assign values to object properties;
    ## array_shape: define dimension (model resolution) of arrays: e.g., (days,lats,lons)
    ## N_input: annual excreted N by livestock divided by farming area; kgN per unit area per year; 
    ##          shape: same as the np.zeros(array_shape), e.g., np.zeros(mtrx)
    ## housing_type: 1. 'slat/pit house' 
    ##               2. 'barn'
    ##               3. 'poultry house'
    ##               ...

    def __init__(self,N_input,livestock_name,house_env,housing_type,housing_area,fslat=0.8,fpit=1.0,idealised_sim=False,idealised_setting=None):
        if idealised_sim is True:
            print("NOTICE - This is an idealised simulation.\n Set [idealised_sim=False] for normal simulations.")
            try:
                if idealised_setting is not None:
                    self.idealised_animal_head = animal_head.copy()
                    self.idealised_animal_head.values[~np.isnan(self.idealised_animal_head.values)] = idealised_setting['head']
                    self.idealised_animal_weight = animal_weight.copy()
                    self.idealised_animal_weight.values[~np.isnan(self.idealised_animal_weight.values)] = idealised_setting['weight']
                    self.idealised_housing_area = self.idealised_animal_head*self.idealised_animal_weight/\
                                                idealised_setting['density']
                    self.idealised_N_input = F_Ncorrect*self.idealised_animal_head*idealised_setting['Nexcret_rate']/\
                                                self.idealised_housing_area
                    self.idealised_N_input = xr_to_np(self.idealised_N_input)
                    self.durine_N, self.durea, self.dmanure_N, self.durine, self.dmanure, self.manure_wc,self.pH = livestock_waste_info(livestock_type=livestock_name, 
                                                    waste_N=self.idealised_N_input)
                else:
                    raise ValueError('Invalid [idealised settings]')
            except ValueError as exp:
                print("Error occurrs! Please check [idealised settings]")
        else:
            ## animal waste info    
            self.durine_N, self.durea, self.dmanure_N, self.durine, self.dmanure, self.manure_wc,self.pH = livestock_waste_info(livestock_type=livestock_name, 
                        waste_N=N_input)

        ## show current config settings, e.g., livestock, production sys, etc.
        print('HOUSING Module - current livestock is: '+str(livestock))
        if production_system is None:
            raise ValueError("ERROR - incorrect production system for "+str(livestock))
        else:
            print('HOUSING Module - current production system is: '+str(production_system))
        if housing_system is None:
            raise ValueError("ERROR - incorrect housing system for "+str(livestock))
        else:
            print('HOUSING Module - current housing system is: '+str(housing_system))
        
        self.durine = self.durine * 1000
        self.manure_wc = self.manure_wc * 1000
        ## pH and H+ ions concentration
        self.cc_H = np.float(10**(-self.pH))

        field_shape = (lats,lons)
        array_shape = (25,lats,lons)
        ## output shape
        outarray_shape = (Days,lats,lons)

        ## environmental fields
        self.T_sim = np.zeros(array_shape)
        self.T_gnd = np.zeros(array_shape)
        self.u_sim = np.zeros(array_shape)
        self.RH_sim = np.zeros(array_shape)

        #####################################
        ## define reaction rates
        #####################################
        self.urea_hydro_rate = np.zeros(array_shape)
        self.Na_decomp_rate = np.zeros(array_shape)
        self.Nr_decomp_rate = np.zeros(array_shape)

        ##################################
        ## define prognostic variables
        ##################################
        ## for housing simulation
        self.NH3_inconc = np.zeros(array_shape)
        self.NH3_out = np.zeros(array_shape)
        print("Housing env is: "+str(house_env))
        print("Housing system is: "+str(housing_type))
        ## feces input
        self.manure = np.zeros(array_shape)
        ## urine input
        self.urine = np.zeros(array_shape)
        ## urea input
        self.urea = np.zeros(array_shape)
        ## input of available N component
        self.avail_N = np.zeros(array_shape)
        ## input of resistant N component
        self.resist_N = np.zeros(array_shape)
        ## input of unavailable N component
        self.unavail_N = np.zeros(array_shape)

        if housing_type.lower() == 'slat/pit house':
            self.evap_slat = np.zeros(array_shape)
            self.evap_pit = np.zeros(array_shape)
            self.R_slat = np.zeros(array_shape)
            self.R_pit = np.zeros(array_shape)
            ## feces pool
            self.manure_pool_slat = np.zeros(array_shape)
            self.manure_pool_pit = np.zeros(array_shape)
            ## min amount of water in manure (function of T, RH)
            ## equilibrium moisture content, water content cannot be lower than this threshold under "natural" processes, 
            ## i.e., without drying processes
            self.manure_minwc_slat = np.zeros(array_shape)
            self.manure_minwc_pit = np.zeros(array_shape)
            ## water amount in fresh feces
            self.manure_initwc_slat = np.zeros(array_shape)
            self.manure_initwc_pit = np.zeros(array_shape)
            ## urea pool
            self.urea_pool_slat = np.zeros(array_shape)
            self.urea_pool_pit = np.zeros(array_shape)
            ## Nitrogen pool that is available (easily) to form TAN
            self.avail_N_pool_slat = np.zeros(array_shape)
            self.avail_N_pool_pit = np.zeros(array_shape)
            ## Nitrogen pool that is resistant (slowly) to form TAN
            self.resist_N_pool_slat = np.zeros(array_shape)
            self.resist_N_pool_pit = np.zeros(array_shape)
            ## Nitrogen pool that is unavailable (cannot) to form TAN
            self.unavail_N_pool_slat = np.zeros(array_shape)
            self.unavail_N_pool_pit = np.zeros(array_shape)
            ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool
            self.TAN_prod_slat = np.zeros(array_shape)
            self.TAN_prod_pit = np.zeros(array_shape)
            ## TAN pool
            self.TAN_pool_slat = np.zeros(array_shape)
            self.TAN_pool_pit = np.zeros(array_shape)
            ## TAN pool conc (aqueous phase, g/m3)
            self.TAN_amount_slat = np.zeros(array_shape)
            self.TAN_amount_pit = np.zeros(array_shape)
            ## total water pool of the system (manure water; urine+manure water+[washing water])
            self.Total_water_pool_slat = np.zeros(array_shape)
            self.Total_water_pool_pit = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in g/m3)
            self.NH3_gas_slat = np.zeros(array_shape)
            self.NH3_gas_pit = np.zeros(array_shape)
            ## emission potential
            self.modelled_emiss_slat = np.zeros(array_shape)
            self.modelled_emiss_pit = np.zeros(array_shape)
            ## final emission
            ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3flux] 
            ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to 
            ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
            ## during land spreading stages, such as canopy recapture, deap injection...
            self.NH3flux_slat = np.zeros(array_shape)
            self.NH3flux_pit = np.zeros(array_shape)

            ## pools left after housing; transfer to storage
            self.manure_pool_pit_to_storage = np.zeros(outarray_shape)
            self.avail_N_pool_pit_to_storage = np.zeros(outarray_shape)
            self.resist_N_pool_pit_to_storage = np.zeros(outarray_shape)
            self.unavail_N_pool_pit_to_storage = np.zeros(outarray_shape)
            self.urea_pool_pit_to_storage = np.zeros(outarray_shape)
            self.TAN_pool_pit_to_storage = np.zeros(outarray_shape)
            self.Total_water_pool_pit_to_storage = np.zeros(outarray_shape)

            ## output of NH3 flux
            self.o_NH3flux_slat = np.zeros(outarray_shape)
            self.o_NH3flux_pit = np.zeros(outarray_shape)

            ## fraction of the floor area, assuming 80% is solid floor
            self.fslat = fslat
            ## fraction of the gap of the area, assuming 20% (1-40%) is gap
            self.fgap = 1.0 - fslat
            ## relative surface area of the underneath pit, different to fgap
            self.fpit = fpit

            ## area of slatted floor: multiplied by the fraction of slat/pit housing
            self.floor_area = housing_area*(1.0-f_loss-f_sold)*f_housing_pit
            ## area of pit: multiplied by the fraction of slat/pit housing
            self.pit_area = fpit * housing_area*(1.0-f_loss-f_sold)*f_housing_pit
        
        else:
            self.evap = np.zeros(array_shape)
            self.R_star = np.zeros(array_shape)
            ## feces pool
            self.manure_pool = np.zeros(array_shape)
            ## min amount of water in manure (function of T, RH)
            ## equilibrium moisture content, water content cannot be lower than this threshold under "natural" processes, 
            ## i.e., without drying processes
            self.manure_water = np.zeros(array_shape)
            self.manure_minwc = np.zeros(array_shape)
            ## water amount in fresh feces
            self.manure_initwc = np.zeros(array_shape)
            ## Nitrogen pool that is available (easily) to form TAN
            self.avail_N_pool = np.zeros(array_shape)
            ## Nitrogen pool that is resistant (slowly) to form TAN
            self.resist_N_pool = np.zeros(array_shape)
            ## Nitrogen pool that is unavailable (cannot) to form TAN
            self.unavail_N_pool = np.zeros(array_shape)
            ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool
            self.TAN_prod = np.zeros(array_shape)
            ## TAN pool
            self.TAN_pool = np.zeros(array_shape)
            ## TAN pool conc (aqueous phase, g/m3)
            self.TAN_amount = np.zeros(array_shape)
            ## total water pool of the system (manure water; urine+manure water+[washing water])
            self.Total_water_pool = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in g/m3)
            self.NH3_gas = np.zeros(array_shape)
            ## emission potential
            self.modelled_emiss = np.zeros(array_shape)
            ## final emission
            ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3flux] 
            ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to 
            ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
            ## during land spreading stages, such as canopy recapture, deap injection...
            self.NH3flux = np.zeros(array_shape)

            ## pools left after housing; transfer to storage
            self.manure_pool_to_storage = np.zeros(outarray_shape)
            self.avail_N_pool_to_storage = np.zeros(outarray_shape)
            self.resist_N_pool_to_storage = np.zeros(outarray_shape)
            self.unavail_N_pool_to_storage = np.zeros(outarray_shape)
            self.TAN_pool_to_storage = np.zeros(outarray_shape)
            self.Total_water_pool_to_storage = np.zeros(outarray_shape)
            self.NH3flux_from_barn = np.zeros(outarray_shape)

            ## output of NH3 flux
            self.o_NH3flux = np.zeros(outarray_shape)

            ## area of housing area
            self.floor_area = housing_area*(1.0-(1.0-f_loss-f_sold)*f_housing_pit)

            if housing_type.lower() == 'barn':
                ## urea pool
                self.urea_pool = np.zeros(array_shape)
                self.urea_pool_to_storage = np.zeros(outarray_shape)
            elif housing_type.lower() == 'poultry house':
                ## 
                self.N_birds_UAN = N_input * 1000/365 * f_uan
                self.manure = self.N_birds_UA/f_excretn
                ## uric acid pool
                self.UA_pool = np.zeros(array_shape)
                self.UA_pool_to_storage = np.zeros(array_shape)
                ## uric acid hydrolysis rate
                self.ua_conv_factor = ua_hydrolysis_rate(temp=self.T_sim,rhum=self.RH_sim,ph=self.pH)
    
    
    ## surface resistance assumed to be 100 s/m (for bedding etc.)
    R_surf = 100.0
    ## background NH3 level, ug/m3
    X_air = 0.0
    ## cleaning frequency; i.e. duration of manure in houses until collected and transfered to storage
    cleaning_cycl = 90
    cleaning_freq = int(cleaning_cycl * nhours / timestep)

    def met_input_interp(self,template):
        ##################################
        ## fill land input data
        ##################################
        self.u_sim = field_var_fill(sd_template=template,input_field=self.u_sim) ## m/s
        self.RH_sim = field_var_fill(sd_template=template,input_field=self.RH_sim)  ## s/m
        return

    ## housing environmental conditions; and diagnostic variables
    def sim_env(self,house_env,housing_type,dayidx):
        hhidx = dayidx*24
        # temp_data = temp_file.t2m[hhidx:hhidx+24] - 273.15
        # rhum_data = rhum_file.rh2m[hhidx:hhidx+24]
        # wind_data = wind_file.Wind_Speed_10m_Mean[hhidx:hhidx+24]
        temp_data = temp_file.t2m[dayidx] - 273.15
        rhum_data = rhum_file.Relative_Humidity_2m_06h[dayidx]
        wind_data = wind_file.Wind_Speed_10m_Mean[dayidx]
        ## housing environmental conditions
        if house_env.lower() == 'insulated':
            # print("HOUSING ENV: House with slatted floor")
            self.T_sim[1:], self.u_sim[1:], self.RH_sim[1:] = housing_env(temp_data,rhum_data,livestock,production_system)
            self.T_gnd = self.T_sim
        elif house_env.lower() == 'naturally ventilated':
            # print("HOUSING ENV: Naturally ventilated barn")
            self.T_sim[1:], self.T_gnd[1:], self.u_sim[1:] = barn_env(temp_data,
                        wind_profile(uref=wind_data,height_ref=wind_data_height,height_out=ref_height,zo=zo_house))
            self.RH_sim[1:] = rhum_data
        # elif house_env.lower() == 'poultry house':
        #     # print("HOUSING ENV: Poultry house")
        #     self.T_sim[1:], self.u_sim[1:], self.RH_sim[1:] = housing_env(temp_data,rhum_data,livestock,production_system)
        else:
            # print("HOUSING ENV: Other housing systems")
            self.T_sim[1:] = temp_data
            self.u_sim[1:] = wind_data
            self.RH_sim[1:] = rhum_data
        
        ## coverting xarray to numpy array
        self.T_sim = xr_to_np(self.T_sim)
        self.u_sim = xr_to_np(self.u_sim)
        self.RH_sim = xr_to_np(self.RH_sim)
        self.met_input_interp(template=animal_file['Excreted_N'][lvl_idx])

        ###################################################
        ## calculating diagnostic variables
        ###################################################
        if housing_type.lower() == 'slat/pit house':
            ## evaporation; aerodynamic method; (m/s)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            self.evap_slat = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house)
            self.evap_pit = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house)
            ## convert evap from m/s to g/m2/day
            self.evap_slat = self.evap_slat*1e6*timestep*3600*self.fslat
            self.evap_pit = self.evap_pit*1e6*timestep*3600
        else:
            ## evaporation; aerodynamic method; (m/s)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            self.evap = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house)
            self.evap = self.evap*1e6*timestep*3600

            if housing_type.lower() == 'poultry house':
                ## housing resistance; 16700 s/m for poultry houses
                self.R_star[:] = 16700
        


    ###################################
    ## define model functions
    ###################################
    ## simulation: for houses with slatted floor (two-source model: slat+pit); industrial production system
    def slat_pit_housing_sim(self,start_day_idx,end_day_idx,house_env,f_slat,f_gap):
        fpitcorrect = 1/self.fpit
        for dd in np.arange(start_day_idx,end_day_idx):
            if dd<Days:
                self.sim_env(house_env,"slat/pit house",dd)
            else:
                self.sim_env(house_env,"slat/pit house",dd-Days)
            self.daily_init(housing_type="slat/pit house")
            for hh in np.arange(24):
                self.manure[hh+1] = self.dmanure
                self.urine[hh+1] = self.durine

                ## manure pool
                ## the amount of manure left on slat depends on the solid floor ratio, here f_slat = 80%
                ## we assume manure left on slat will drop to the pit on the next day
                self.manure_pool_slat[hh+1] = self.manure_pool_slat[hh] + self.manure[hh+1]*f_slat
                self.manure_pool_pit[hh+1] = self.manure_pool_pit[hh] + self.manure[hh+1]*f_gap*fpitcorrect 

                ## N input in multiple forms
                self.urea[hh+1] = self.durea
                self.avail_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_avail
                self.resist_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_resist
                self.unavail_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_unavail

                ## urea hydrolysis rate
                self.urea_hydro_rate[hh+1] = urea_hydrolysis_rate(temp=self.T_sim[hh+1],theta=1.0,delta_t=timestep)
                ## decomposition rate of available and resistant N components
                self.Na_decomp_rate[hh+1], self.Nr_decomp_rate[hh+1] = N_pools_decomp_rate(temp=self.T_sim[hh+1], 
                                                                                                    delta_t=timestep)

                ## Urea pool
                ## slat urea pool is reset at a daily basis
                self.urea_pool_slat[hh+1] = self.urea_pool_slat[hh] + self.urea[hh+1]*f_slat
                ## urea left on slat will go to the pit urea pool
                self.urea_pool_pit[hh+1] = self.urea_pool_pit[hh] + self.urea[hh+1]*f_gap*fpitcorrect 

                ## Org N pools in various forms
                self.avail_N_pool_slat[hh+1] = self.avail_N_pool_slat[hh] + self.avail_N[hh+1]*f_slat
                self.avail_N_pool_pit[hh+1] = self.avail_N_pool_pit[hh] + self.avail_N[hh+1]*f_gap*fpitcorrect 

                self.resist_N_pool_slat[hh+1] = self.resist_N_pool_slat[hh] + self.resist_N[hh+1]*f_slat
                self.resist_N_pool_pit[hh+1] = self.resist_N_pool_pit[hh] + self.resist_N[hh+1]*f_gap*fpitcorrect

                self.unavail_N_pool_slat[hh+1] = self.unavail_N_pool_slat[hh] + self.unavail_N[hh+1]*f_slat
                self.unavail_N_pool_pit[hh+1] = self.unavail_N_pool_pit[hh] + self.unavail_N[hh+1] *f_gap*fpitcorrect

                ## TAN production from urea hydrolysis and org N decomposition from slat and pit, respectively
                ## the conditions (T, RH) of slat and pit are assume to be consistent, so is the urea hydrolysis rate
                ## and the N decomposition rate from dung
                self.TAN_prod_slat[hh+1] = self.urea_hydro_rate[hh+1]*self.urea_pool_slat[hh+1]+ \
                                            self.Na_decomp_rate[hh+1]*self.avail_N_pool_slat[hh+1] + \
                                            self.Nr_decomp_rate[hh+1]*self.resist_N_pool_slat[hh+1]

                self.TAN_prod_pit[hh+1] = self.urea_hydro_rate[hh+1]*self.urea_pool_pit[hh+1]+\
                                        self.Na_decomp_rate[hh+1]*self.avail_N_pool_pit[hh+1] +\
                                        self.Nr_decomp_rate[hh+1]*self.resist_N_pool_pit[hh+1]

                ## update urea pool and orgN pools
                self.urea_pool_slat[hh+1] = self.urea_pool_slat[hh+1]*(1-self.urea_hydro_rate[hh+1])
                ## urea left on slat will go to the pit urea pool
                self.urea_pool_pit[hh+1] = self.urea_pool_pit[hh+1]*(1-self.urea_hydro_rate[hh+1]) 

                ## Org N pools in various forms
                self.avail_N_pool_slat[hh+1] = self.avail_N_pool_slat[hh+1]*(1-self.Na_decomp_rate[hh+1])
                self.avail_N_pool_pit[hh+1] = self.avail_N_pool_pit[hh+1]*(1-self.Na_decomp_rate[hh+1])

                self.resist_N_pool_slat[hh+1] = self.resist_N_pool_slat[hh+1]*(1-self.Nr_decomp_rate[hh+1]) 
                self.resist_N_pool_pit[hh+1] = self.resist_N_pool_pit[hh+1]*(1-self.Nr_decomp_rate[hh+1])

                ## water from the fresh dung
                self.manure_initwc_slat[hh+1] = self.manure_wc*f_slat
                self.manure_initwc_pit[hh+1] = self.manure_wc*f_gap*fpitcorrect

                ## water amount when mositure content reach equilibrium
                # mois_coeff = min_manurewc(temp=self.T_sim[hh+1],rhum=self.RH_sim[hh+1])
                mois_coeff_slat = 1.0
                mois_coeff_pit = 10.0
                self.manure_minwc_slat[hh+1] = self.manure_pool_slat[hh+1]*mois_coeff_slat
                self.manure_minwc_pit[hh+1] = self.manure_pool_pit[hh+1]*mois_coeff_pit

                ## water pool of slat and pit
                ## slat water pool
                self.Total_water_pool_slat[hh+1] = self.Total_water_pool_slat[hh]+self.urine[hh+1]*f_slat +\
                                                    self.manure_initwc_slat[hh+1] - self.evap_slat[hh+1]
                water_slat_idx = self.Total_water_pool_slat[hh+1] - self.manure_minwc_slat[hh+1]
                self.Total_water_pool_slat[hh+1][water_slat_idx<0] = self.manure_minwc_slat[hh+1][water_slat_idx<0]
                ## pit water pool
                self.Total_water_pool_pit[hh+1] = self.Total_water_pool_pit[hh]+self.urine[hh+1]*f_gap*fpitcorrect +\
                                                    self.manure_initwc_pit[hh+1]-self.evap_pit[hh]
                water_pit_idx = self.Total_water_pool_pit[hh+1] - self.manure_minwc_pit[hh+1]
                self.Total_water_pool_pit[hh+1][water_pit_idx<0] = self.manure_minwc_pit[hh+1][water_pit_idx<0]


                ## TAN pool of slat and pit
                self.TAN_pool_slat[hh+1] = self.TAN_pool_slat[hh] + self.TAN_prod_slat[hh+1]
                ## TAN left on slat will go to the pit TAN pool
                self.TAN_pool_pit[hh+1] = self.TAN_pool_pit[hh] + self.TAN_prod_pit[hh+1]

                ## TAN conc on slat and in pit (currently in g/ml)
                self.TAN_amount_slat[hh+1] = self.TAN_pool_slat[hh+1]/self.Total_water_pool_slat[hh+1]
                self.TAN_amount_slat[hh+1][self.Total_water_pool_slat[hh+1]==0] = 0.0

                self.TAN_amount_pit[hh+1] = self.TAN_pool_pit[hh+1]/self.Total_water_pool_pit[hh+1]
                self.TAN_amount_pit[hh+1][self.Total_water_pool_pit[hh+1]==0] = 0.0

                ## TAN conc in g/m3
                self.TAN_amount_slat[hh+1] = self.TAN_amount_slat[hh+1] * 1e6
                self.TAN_amount_pit[hh+1] = self.TAN_amount_pit[hh+1] * 1e6

                ## NH3 conc in g/m3
                KNH3 = NH3_par_coeff(temp=self.T_sim[hh+1],cncH=self.cc_H)
                self.NH3_gas_slat[hh+1] = KNH3*self.TAN_amount_slat[hh+1]
                self.NH3_gas_pit[hh+1] = KNH3*self.TAN_amount_pit[hh+1]

                ## resistances
                self.R_slat[hh+1] = 1/k_gas_NH3(temp=self.T_sim[hh+1],u=self.u_sim[hh+1],Z=2,zo=zo_house)
                self.R_pit[hh+1] = (1/k_aq_NH4(temp=self.T_sim[hh+1]))+(1/(KNH3*k_gas_NH3(temp=self.T_sim[hh+1],
                                                                    u=self.u_sim[hh+1],Z=2,zo=zo_water)))

                ## determining the maximum emission; emission cannot exceed TAN pool
                self.modelled_emiss_slat[hh+1] = NH3_volslat(slat_conc=self.NH3_gas_slat[hh+1],conc_in=0.0,
                                                    Rslat=self.R_slat[hh+1])*timestep*3600
                emiss_slat_idx = self.TAN_pool_slat[hh+1] - self.modelled_emiss_slat[hh+1] 
                self.modelled_emiss_slat[hh+1][emiss_slat_idx<0] = self.TAN_pool_slat[hh+1][emiss_slat_idx<0]

                self.modelled_emiss_pit[hh+1] = NH3_volpit(pit_tanconc=self.TAN_amount_pit[hh+1],conc_in=0.0,
                                                Rpit=self.R_pit[hh+1])*timestep*3600
                emiss_pit_idx = self.TAN_pool_pit[hh+1] - self.modelled_emiss_pit[hh+1]
                self.modelled_emiss_pit[hh+1][emiss_pit_idx<0] = self.TAN_pool_pit[hh+1][emiss_pit_idx<0]

                ## final emission flux from slat and pit; flux are additive
                self.NH3flux_slat[hh+1] = self.modelled_emiss_slat[hh+1]
                self.NH3flux_pit[hh+1] = self.modelled_emiss_pit[hh+1]

                ## update TAN pool
                self.TAN_pool_slat[hh+1] = self.TAN_pool_slat[hh+1] - self.NH3flux_slat[hh+1]
                self.TAN_pool_pit[hh+1] = self.TAN_pool_pit[hh+1] - self.NH3flux_pit[hh+1]

            if dd >= Days:
                ddidx = dd - Days
            else:
                ddidx = dd
            self.daily_output(housing_type="slat/pit house",dayidx=ddidx)
            
        return

    ## simulation: for houses with concrete floor (singe-source model); 
    ## intermediate/backyard production system (pig/ruminants)
    def barn_housing_sim(self,start_day_idx,end_day_idx,house_env):
        for dd in np.arange(start_day_idx,end_day_idx):
            if dd<Days:
                self.sim_env(house_env,"barn",dd)
            else:
                self.sim_env(house_env,"barn",dd-Days)
            self.daily_init(housing_type="barn")
            for hh in np.arange(24):
                self.manure[hh+1] = self.dmanure
                self.urine[hh+1] = self.durine
                ## manure pool
                self.manure_pool[hh+1] = self.manure_pool[hh] + self.manure[hh+1]

                ## N input in multiple forms
                self.urea[hh+1] = self.durea
                self.avail_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_avail
                self.resist_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_resist
                self.unavail_N[hh+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_unavail

                ## urea hydrolysis rate
                self.urea_hydro_rate[hh+1] = urea_hydrolysis_rate(temp=self.T_gnd[hh+1],theta=1.0,delta_t=timestep)
                ## decomposition rate of available and resistant N components
                self.Na_decomp_rate[hh+1], self.Nr_decomp_rate[hh+1] = N_pools_decomp_rate(temp=self.T_gnd[hh+1], 
                                                                                                    delta_t=timestep)

                ## Urea pool
                self.urea_pool[hh+1] = self.urea_pool[hh] + self.urea[hh+1]

                ## Org N pools in various forms
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh] + self.avail_N[hh+1]
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh] + self.resist_N[hh+1]            
                self.unavail_N_pool[hh+1] = self.unavail_N_pool[hh] + self.unavail_N[hh+1] 

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[hh+1] = self.urea_hydro_rate[hh+1]*self.urea_pool[hh+1]+\
                                        self.Na_decomp_rate[hh+1]*self.avail_N_pool[hh+1] +\
                                        self.Nr_decomp_rate[hh+1]*self.resist_N_pool[hh+1]

                ## update urea pool and orgN pools
                self.urea_pool[hh+1] = self.urea_pool[hh+1]*(1 - self.urea_hydro_rate[hh+1]) 
                self.avail_N_pool[hh+1] = self.avail_N_pool[hh+1]*(1 - self.Na_decomp_rate[hh+1]) 
                self.resist_N_pool[hh+1] = self.resist_N_pool[hh+1]* (1 - self.Nr_decomp_rate[hh+1])  

                ## water from the fresh dung
                self.manure_initwc[hh+1] = self.manure_wc
                ## water amount when mositure content reach equilibrium
                # mois_coeff = min_manurewc(temp=self.T_sim[hh+1],rhum=self.RH_sim[hh+1])
                mois_coeff = 1.0 
                self.manure_minwc[hh+1] = self.manure_pool[hh+1]*mois_coeff
                ## water pool 
                self.Total_water_pool[hh+1] = self.Total_water_pool[hh]+self.urine[hh+1]+self.manure_initwc[hh+1]-self.evap[hh+1]
                water_idx = self.Total_water_pool[hh+1] - self.manure_minwc[hh+1] 
                self.Total_water_pool[hh+1][water_idx<0] = self.manure_minwc[hh+1][water_idx<=0] 

                ## TAN pool 
                self.TAN_pool[hh+1] = self.TAN_pool[hh]+self.TAN_prod[hh+1]

                ## TAN conc (currently in g/ml)
                self.TAN_amount[hh+1] = self.TAN_pool[hh+1]/self.Total_water_pool[hh+1]
                self.TAN_amount[hh+1][self.Total_water_pool[hh+1]==0] = 0.0
                ## TAN conc in g/m3
                self.TAN_amount[hh+1] = self.TAN_amount[hh+1] * 1e6

                ## NH3 conc in g/m3
                KNH3 = NH3_par_coeff(temp=self.T_gnd[hh+1],cncH=self.cc_H)
                self.NH3_gas[hh+1] = self.TAN_amount[hh+1]*KNH3

                ## resistance
                self.R_star[hh+1] = 1/k_gas_NH3(temp=self.T_sim[hh+1],u=self.u_sim[hh+1],Z=2,zo=zo_house)+self.R_surf

                ## determining the maximum emission; emission cannot exceed TAN pool
                self.modelled_emiss[hh+1] = NH3_volslat(slat_conc=self.NH3_gas[hh+1],conc_in=0.0,
                                                        Rslat=self.R_star[hh+1])*timestep*3600
                emiss_idx = self.TAN_pool[hh+1] - self.modelled_emiss[hh+1]
                self.modelled_emiss[hh+1][emiss_idx<0] = self.TAN_pool[hh+1][emiss_idx<0]

                ## final emission flux 
                self.NH3flux[hh+1] = self.modelled_emiss[hh+1]

                ## update TAN pool
                self.TAN_pool[hh+1] = self.TAN_pool[hh+1] - self.NH3flux[hh+1]
            if dd >= Days:
                ddidx = dd - Days
            else:
                ddidx = dd
            self.daily_output(housing_type="barn",dayidx=ddidx)
            
        return

    ## simulation: for poultry housing;
    ## broilers/layers; housing resistance is assumed to be 16700 s/m
    def poultry_housing_sim(self,start_day_idx,end_day_idx,house_env):
        for dd in np.arange(start_day_idx,end_day_idx):
            if dd<Days:
                self.sim_env(house_env,"barn",dd)
            else:
                self.sim_env(house_env,"barn",dd-Days)
            for hh in np.arange(24):
                ## excretion pool; Eq.1
                self.manure_pool[hh+1] = self.manure_pool[hh] + self.manure
                self.TAN_prod[hh+1] = self.ua_conv_factor[hh+1] * self.UA_pool[hh]
                ## UA pool; Eq.2
                UA_idx = self.UA_pool[hh] - self.TAN_prod[hh+1]
                self.UA_pool[hh+1][UA_idx>0] = self.UA_pool[hh][UA_idx>0] + self.N_birds_UAN[UA_idx>0] - self.TAN_prod[hh+1][UA_idx>0]
                self.UA_pool[hh+1][UA_idx<=0] = self.N_birds_UAN[UA_idx<=0]

                ## manure_water calculation; Eq.9  
                self.manure_water[hh+1][self.mois_coeff[hh+1]>=5] = self.manure_pool[hh+1][self.mois_coeff[hh+1]>=5] * self.mois_coeff[hh+1][self.mois_coeff[hh+1]>=5] / 100
                self.manure_water[hh+1][self.mois_coeff[hh+1]<5] = 0
                self.Total_water_pool[hh+1] = self.manure_water[hh+1]
                
                ## TAN pool; Eq.3
                TAN_idx = self.TAN_pool[hh] - self.NH3flux[hh]
                self.TAN_pool[hh+1][TAN_idx<=0] = self.TAN_prod[hh+1][TAN_idx<=0]
                self.TAN_pool[hh+1][TAN_idx>0] = TAN_idx[TAN_idx>0] + self.TAN_prod[hh+1][TAN_idx>0]
                self.TAN_pool_ug[hh+1] = self.TAN_pool[hh+1] * 1e6
                
                ## TAN concentration
                self.TAN_amount[hh+1][self.Total_water_pool[hh+1]==0] = 0
                self.TAN_amount[hh+1][self.Total_water_pool[hh+1]!=0] = self.TAN_pool[hh+1][self.Total_water_pool[hh+1]!=0] / self.Total_water_pool[hh+1][self.Total_water_pool[hh+1]!=0]
                self.TAN_amount_M[hh+1] = self.TAN_amount[hh+1] / 14 * 1000

                ## Gamma value; Eq.7
                self.Gamma_manure[hh+1] =  self.TAN_amount_M[hh+1] / (self.cc_H + self.k_NH4[hh+1])
                
                ## Gaseous concentration of NH3 at the surface; Eq.6
                self.NH3_gas_M[hh+1] = self.Henry_constant[hh+1] * self.Gamma_manure[hh+1]
                self.NH3_gas_ug[hh+1] = self.NH3_gas_M[hh+1] * 14 * 1e9

                ## Emissions; see Eq.8 in general, and Sect.2.2.2; Eq.14~17
                emiss_idx = (self.NH3_gas_ug[hh+1]*3600*timestep/self.R_star) - self.TAN_pool_ug[hh+1]
                self.modelled_emiss[hh+1][emiss_idx>0] = self.TAN_pool_ug[hh+1][emiss_idx>0]
                self.modelled_emiss[hh+1][emiss_idx<0] = self.NH3_gas_ug[hh+1][emiss_idx<0]*3600*timestep/self.R_star
                self.NH3flux[hh+1] = self.Modelled_emiss[hh+1]/1e6 
        return

    ## initialisation: model initialisation for housing simulation
    def housing_init(self,housing_type):
        self.manure[:] = 0.0
        self.urine[:] = 0.0
        self.urea[:] = 0.0
        self.avail_N[:] = 0.0
        self.resist_N[:] = 0.0
        self.unavail_N[:] = 0.0

        if housing_type.lower() == 'slat/pit house':
            self.manure_pool_slat[:] = 0.0
            self.manure_pool_pit[:] = 0.0
            self.manure_minwc_slat[:] = 0.0
            self.manure_minwc_pit[:] = 0.0
            self.manure_initwc_slat[:] = 0.0
            self.manure_initwc_pit[:] = 0.0
            self.urea_pool_slat[:] = 0.0
            self.urea_pool_pit[:] = 0.0
            self.avail_N_pool_slat[:] = 0.0
            self.avail_N_pool_pit[:] = 0.0
            self.resist_N_pool_slat[:] = 0.0
            self.resist_N_pool_pit[:] = 0.0
            self.unavail_N_pool_slat[:] = 0.0
            self.unavail_N_pool_pit[:] = 0.0
            self.TAN_prod_slat[:] = 0.0
            self.TAN_prod_pit[:] = 0.0
            self.TAN_pool_slat[:] = 0.0
            self.TAN_pool_pit[:] = 0.0
            self.TAN_amount_slat[:] = 0.0
            self.TAN_amount_pit[:] = 0.0
            self.Total_water_pool_slat[:] = 0.0
            self.Total_water_pool_pit[:] = 0.0
            self.modelled_emiss_slat[:] = 0.0
            self.modelled_emiss_pit[:] = 0.0
            self.NH3flux_slat[:] = 0.0
            self.NH3flux_pit[:] = 0.0

        else:
            self.manure_pool[:] = 0.0
            self.manure_water[:] = 0.0
            self.manure_initwc[:] = 0.0
            self.manure_minwc[:] = 0.0
            self.avail_N_pool[:] = 0.0
            self.resist_N_pool[:] = 0.0
            self.unavail_N_pool[:] = 0.0
            self.TAN_prod[:] = 0.0
            self.TAN_pool[:] = 0.0
            self.TAN_amount[:] = 0.0
            self.Total_water_pool[:] = 0.0
            self.modelled_emiss[:] = 0.0
            self.NH3flux[:] = 0.0
            if housing_type.lower() == 'barn':
                self.urea_pool[:] = 0.0
            elif housing_type.lower() == 'poultry_house':
                self.UA_pool[:] = 0.0

    ## initialisation: 2nd initialisation for housing simulation
    def housing_2nd_init(self,housing_type):
        if housing_type.lower() == 'slat/pit house':
            aa2 = self.manure_pool_slat[-1]
            aa3 = self.manure_pool_pit[-1]
            bb2 = self.urea_pool_slat[-1]
            bb3 = self.urea_pool_pit[-1]
            cc2 = self.TAN_pool_slat[-1]
            cc3 = self.TAN_pool_pit[-1]
            dd2 = self.NH3flux_slat[-1]
            dd3 = self.NH3flux_pit[-1]
            ee2 = self.Total_water_pool_slat[-1]
            ee3 = self.Total_water_pool_pit[-1]

            self.manure_pool_slat[0] = aa2
            self.manure_pool_pit[0] = aa3
            self.urea_pool_slat[0] = bb2
            self.urea_pool_pit[0] = bb3
            self.TAN_pool_slat[0] = cc2
            self.TAN_pool_pit[0] = cc3
            self.NH3flux_slat[0] = dd2
            self.NH3flux_pit[0] = dd3
            self.Total_water_pool_slat[0] = ee2
            self.Total_water_pool_pit[0] = ee3  
        else:
            aa1 = self.manure_pool[-1]
            
            cc1 = self.TAN_pool[-1]
            dd1 = self.NH3flux[-1]
            ee1 = self.Total_water_pool[-1]

            self.manure_pool[0] = aa1
            
            self.TAN_pool[0] = cc1
            self.NH3flux[0] = dd1
            self.Total_water_pool[0] = ee1

            if housing_type.lower() == 'barn':
                bb1 = self.urea_pool[-1]
                self.urea_pool[0] = bb1
            elif housing_type.lower() == 'poultry house':
                bb1 = self.UA_pool[-1]
                self.UA_pool[0] = bb1

    def daily_init(self,housing_type):
        if housing_type == 'slat/pit house':
            fpitcorrect = 1/self.fpit
            ## slat manure pool goes into pit manure pool
            self.manure_pool_pit[0] = self.manure_pool_pit[-1] + self.manure_pool_slat[-1]*fpitcorrect
            ## slat urea pool goes into pit urea pool
            self.urea_pool_pit[0] = self.urea_pool_pit[-1] + self.urea_pool_slat[-1]*fpitcorrect
            ## slat organic N pools go into pit orgN pools
            self.avail_N_pool_pit[0] = self.avail_N_pool_pit[-1] + self.avail_N_pool_slat[-1]*fpitcorrect
            self.resist_N_pool_pit[0] = self.resist_N_pool_pit[-1] + self.resist_N_pool_slat[-1]*fpitcorrect
            self.unavail_N_pool_pit[0] = self.unavail_N_pool_pit[-1] + self.unavail_N_pool_slat[-1]*fpitcorrect
            ## slat water pool goes into pit water pool
            self.Total_water_pool_pit[0] = self.Total_water_pool_pit[-1] + self.Total_water_pool_slat[-1]*fpitcorrect
            ## slat TAN pool goes into pit TAN pool
            self.TAN_pool_pit[0] = self.TAN_pool_pit[-1] + self.TAN_pool_slat[-1]*fpitcorrect
        else:
            self.manure_pool[0] = self.manure_pool[-1]
            self.manure_water[0] = self.manure_water[-1]
            self.avail_N_pool[0] = self.avail_N_pool[-1]
            self.resist_N_pool[0] = self.resist_N_pool[-1]
            self.unavail_N_pool[0] = self.unavail_N_pool[-1]
            self.urea_pool[0] = self.urea_pool[-1]
            self.TAN_pool[0] = self.TAN_pool[-1]
            self.Total_water_pool[0] = self.Total_water_pool[-1]

    ## output daily sum
    def daily_output(self,housing_type,dayidx):
        if housing_type == "slat/pit house":
            self.o_NH3flux_slat[dayidx] = np.nansum(self.NH3flux_slat[1:],axis=0)
            self.o_NH3flux_pit[dayidx] = np.nansum(self.NH3flux_pit[1:],axis=0)
        if housing_type == "barn":
            self.o_NH3flux[dayidx] = np.nansum(self.NH3flux[1:],axis=0)

    ## initialisation for house (barn) cleaning; all variables at the cleaning day are reset to zero
    def cleaning_pit(self,dayidx,housing_type):
        if housing_type.lower() == 'slat/pit house':
            ## pools: from housing to MMS
            ## Note: by multiplying the housing area, we get the total mass of each pool rather than mass/unit area
            self.manure_pool_pit_to_storage[dayidx] = self.manure_pool_pit[-1]*self.pit_area
            self.avail_N_pool_pit_to_storage[dayidx] = self.avail_N_pool_pit[-1]*self.pit_area
            self.resist_N_pool_pit_to_storage[dayidx] = self.resist_N_pool_pit[-1]*self.pit_area
            self.unavail_N_pool_pit_to_storage[dayidx] = self.unavail_N_pool_pit[-1]*self.pit_area
            self.urea_pool_pit_to_storage[dayidx] = self.urea_pool_pit[-1]*self.pit_area
            self.TAN_pool_pit_to_storage[dayidx] = self.TAN_pool_pit[-1]*self.pit_area
            self.Total_water_pool_pit_to_storage[dayidx] = self.Total_water_pool_pit[-1]*self.pit_area
        
            self.manure_pool_pit[:] = 0.0
            self.manure_initwc_pit[:] = 0.0
            self.urea_pool_pit[:] = 0.0
            self.avail_N_pool_pit[:] = 0.0
            self.resist_N_pool_pit[:] = 0.0
            self.unavail_N_pool_pit[:] = 0.0
            self.TAN_pool_pit[:] = 0.0
            self.Total_water_pool_pit[:] = 0.0
        else:
            #print(housing_type)
            ## pools: from housing to MMS
            self.manure_pool_to_storage[dayidx] = self.manure_pool[-1]*self.floor_area
            self.avail_N_pool_to_storage[dayidx] = self.avail_N_pool[-1]*self.floor_area
            self.resist_N_pool_to_storage[dayidx] = self.resist_N_pool[-1]*self.floor_area
            self.unavail_N_pool_to_storage[dayidx] = self.unavail_N_pool[-1]*self.floor_area
            self.TAN_pool_to_storage[dayidx] = self.TAN_pool[-1]*self.floor_area
            ## get rid of numpy rounding error: maximum NH3 flux is equivalent to TAN pool and should not exceed TAN pool
            ## however, when do calculations, NH3 flux is sometimes rounded to a value that is bigger than TAN pool
            ## e.g. NH3 flux = 2.7233141344, TAN pool = 2.7233141343; this leads to the [TAN to storage] a negative value,
            ## which is against reality/mass conservative principle due to the rounding error in numpy
            self.TAN_pool_to_storage[dayidx][self.TAN_pool_to_storage[dayidx]<0] = 0.0
            #print("TAN pool to storage: ",self.TAN_pool_to_storage[dayidx,41,84])
            self.Total_water_pool_to_storage[dayidx] = self.Total_water_pool[-1]*self.floor_area
            ## NH3 flux is not multiplied by housing area, do this in main script!
    #         self.NH3flux_from_barn[dayidx] = self.NH3flux[dayidx]

            self.manure_pool[:] = 0.0
            self.manure_water[:] = 0.0
            self.avail_N_pool[:] = 0.0
            self.resist_N_pool[:] = 0.0
            self.unavail_N_pool[:] = 0.0
            self.TAN_pool[:] = 0.0
            self.Total_water_pool[:] = 0.0

            if housing_type.lower() == 'barn':
                self.urea_pool_to_storage[dayidx] = self.urea_pool[-1]*self.floor_area
                self.urea_pool[:] = 0.0
            elif housing_type.lower() == 'poultry house':
                self.UA_pool_to_storage[dayidx] = self.UA_pool[-1]*self.floor_area
                self.UA_pool[:] = 0.0

    def housing_to_storage_init(self,housing_type):
        if housing_type.lower() == 'slat/pit house':
            self.manure_pool_pit_to_storage[:] = 0.0
            self.avail_N_pool_pit_to_storage[:] = 0.0
            self.resist_N_pool_pit_to_storage[:] = 0.0
            self.unavail_N_pool_pit_to_storage[:] = 0.0
            self.TAN_pool_pit_to_storage[:] = 0.0
            self.Total_water_pool_pit_to_storage[:] = 0.0
        else:
            self.manure_pool_to_storage[:] = 0.0
            self.avail_N_pool_to_storage[:] = 0.0
            self.resist_N_pool_to_storage[:] = 0.0
            self.unavail_N_pool_to_storage[:] = 0.0
            self.TAN_pool_to_storage[:] = 0.0
            self.Total_water_pool_to_storage[:] = 0.0
            self.NH3flux_from_barn[:] = 0.0
            if housing_type.lower() == 'barn':
                self.urea_pool_to_storage[:] = 0.0
            elif housing_type.lower() == 'poultry house':
                self.UA_pool_to_storage[:] = 0.0

    def sim_out(self,housing_type):
        nlat = int(180.0/dlat)
        nlon = int(360.0/dlon)
        ntime = Days
        lats = 90 - 0.5*np.arange(nlat)
        lons = -180 + 0.5*np.arange(nlon)
        yearidx = str(sim_year)+'-01-01'
        times = pd.date_range(yearidx,periods=ntime)
        if housing_type.lower() == 'slat/pit house':
            slat_emiss = self.o_NH3flux_slat*self.floor_area
            pit_emiss = self.o_NH3flux_pit*self.pit_area
            outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss_slat=(['time','lat','lon'],slat_emiss),
                    NH3emiss_pit=(['time','lat','lon'],pit_emiss),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-Housing: \
                        NH3 emissions from "+\
                            production_system+" "+livestock+" "+housing_type+"housing in " +str(sim_year),
                    info = production_system+" "+livestock+" "+housing_type,
                    units="gN per grid",
                ),
            )
            outds.NH3emiss_slat.attrs["unit"] = 'gN/day'
            outds.NH3emiss_slat.attrs["long name"] = 'NH3 emission from housing (slat)'

            outds.NH3emiss_pit.attrs["unit"] = 'gN/day'
            outds.NH3emiss_pit.attrs["long name"] = 'NH3 emission from housing (pit)'

        elif housing_type.lower() == 'barn':
            housing_NH3emiss = self.o_NH3flux*self.floor_area
            outds = xr.Dataset(
                data_vars=dict(
                    NH3emiss=(['time','lat','lon'],housing_NH3emiss),
                            ),
                coords = dict(
                    time=(["time"], times),
                    lon=(["lon"], lons),
                    lat=(["lat"], lats),
                            ),
                attrs=dict(
                    description="AMCLIM-Housing: \
                        NH3 emissions from "+\
                            production_system+" "+livestock+" "+housing_type+"housing in " +str(sim_year),
                    info = production_system+" "+livestock+" "+housing_type,
                    units="gN per grid",
                ),
            )
            outds.NH3emiss.attrs["unit"] = 'gN/day'
            outds.NH3emiss.attrs["long name"] = 'NH3 emission from housing'

        comp = dict(zlib=True, complevel=9)
        encoding = {var: comp for var in outds.data_vars}

        outds.to_netcdf(output_path+livestock+'.'+production_system+'.'+housing_type+\
                            '.'+str(sim_year)+'.nc',encoding=encoding)
        print("ncfile saved.")
        return
    
    def slat_pit_sim_main(self,house_env,housing_type,start_idx,end_idx,cleaning_frequency):
        print('HOUSING Sim - current simulation is for: '+str(house_env)+', slat/pit housing, '+str(livestock))
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.slat_pit_housing_sim(dd,dd+cleaning_frequency,house_env,self.fslat,self.fgap)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.slat_pit_housing_sim(dd,end_idx,house_env,self.fslat,self.fgap)
        self.housing_2nd_init(housing_type)
        self.slat_pit_housing_sim(0,start_idx,house_env,self.fslat,self.fgap)
        return

    def barn_sim_main(self,house_env,housing_type,start_idx,end_idx,cleaning_frequency,litter=False):
        print('HOUSING Sim - current simulation is for: '+str(house_env)+', barn, '+str(livestock))
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        if litter is True:
            self.R_surf = 100.0
            self.floor_area = self.floor_area*f_housing_litter
        else:
            self.floor_area = self.floor_area*(1.0-f_housing_litter)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.barn_housing_sim(dd,dd+cleaning_frequency,house_env)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.barn_housing_sim(dd,end_idx,house_env)
        self.housing_2nd_init(housing_type)
        self.barn_housing_sim(0,start_idx,house_env)
        return

    def poultry_house_sim_main(self,house_env,housing_type,start_idx,end_idx,cleaning_frequency):
        self.sim_env(housing_type)
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.poultry_housing_sim(dd,dd+cleaning_frequency)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.poultry_housing_sim(dd,end_idx)
        self.housing_2nd_init(housing_type)
        self.poultry_housing_sim(0,start_idx)
        return 

    
