
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
excretN_info = animal_file['Excreted_N'][0]
animal_head = animal_file['Animal_head'][0]
try:
    animal_weight = animal_file['Animal_weight'][0]
except:
    pass
#print(np.nanmedian(animal_weight.values[np.where(animal_weight!=0)]))
animal_weight[np.where((animal_head!=0)&(animal_weight==0))] = np.nanmedian(animal_weight.values[np.where(animal_weight!=0)])
animal_Nrate = animal_file['N_rate'][0]
#print(np.nanmedian(animal_Nrate.values[np.where(animal_Nrate!=0)]))
animal_Nrate[np.where((animal_head!=0)&(animal_Nrate==0))] = np.nanmedian(animal_Nrate.values[np.where(animal_Nrate!=0)])
## pig density 1: 120 kg/m^2
#animal_density = 120.0
animal_density = stocking_desity[livestock]
print(str(livestock)+" stocking density is "+str(animal_density)+" kg/m^2")
massgrid = animal_head*animal_weight
housing_area = animal_head*animal_weight/animal_density
## total animal mass per grid
# massgrid = 1000*excretN_info/(animal_Nrate*365)
# housing_area = massgrid/animal_density
## pig density 2: 1 head/m^2
# animal_density = 1.0
# housing_area = animal_head/animal_density
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
## house surface roughness height; default 2mm
zo_house = 0.002

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

    def __init__(self,array_shape,N_input,livestock_name,housing_type):
        print("Current livestock is: "+str(livestock_name))
        ## animal waste info    
        self.durine_N, self.durea, self.dmanure_N, self.durine, self.dmanure, self.manure_wc, self.pH = livestock_waste_info(livestock_type=livestock_name, waste_N=N_input)
        self.durine = self.durine * 1000
        self.manure_wc = self.manure_wc * 1000
        ## pH and H+ ions concentration
        self.cc_H = np.float(10**(-self.pH))

        ## environmental fields
        self.T_sim = np.zeros(array_shape)
        self.u_sim = np.zeros(array_shape)
        self.RH_sim = np.zeros(array_shape)

        #####################################
        ## define diagnostic fields
        #####################################
        self.mois_coeff = np.zeros(array_shape)
        self.daily_urea_hydro_rate = np.zeros(array_shape)
        self.daily_Na_decomp_rate = np.zeros(array_shape)
        self.daily_Nr_decomp_rate = np.zeros(array_shape)
        self.Henry_constant = np.zeros(array_shape)
        self.k_NH4 = np.zeros(array_shape)

        ##################################
        ## define prognostic variables
        ##################################
        ## for housing simulation
        self.NH3_inconc = np.zeros(array_shape)
        self.NH3_out = np.zeros(array_shape)
        print("Housing system is: "+str(housing_type))
        ## feces input
        self.manure = np.zeros(array_shape)
        ## daily urine input
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
            self.R_star_slat = np.zeros(array_shape)
            self.R_star_pit = np.zeros(array_shape)
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
            ## TAN pool in ug/m2
            self.TAN_pool_ug_slat = np.zeros(array_shape)
            self.TAN_pool_ug_pit = np.zeros(array_shape)
            ## TAN pool conc (aqueous phase)
            self.TAN_amount_slat = np.zeros(array_shape)
            self.TAN_amount_pit = np.zeros(array_shape)
            ## TAN pool in molar concentration
            self.TAN_amount_M_slat = np.zeros(array_shape)
            self.TAN_amount_M_pit = np.zeros(array_shape)
            ## total water pool of the system (manure water; urine+manure water+[washing water])
            self.Total_water_pool_slat = np.zeros(array_shape)
            self.Total_water_pool_pit = np.zeros(array_shape)
            ## ratio of [NH4+]/[H+] in the system
            self.Gamma_manure_slat = np.zeros(array_shape)
            self.Gamma_manure_pit = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in molar mass)
            self.NH3_gas_M_slat = np.zeros(array_shape)
            self.NH3_gas_M_pit = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in ug)
            self.NH3_gas_ug_slat = np.zeros(array_shape)
            self.NH3_gas_ug_pit = np.zeros(array_shape)
            ## emission potential
            self.modelled_emiss_slat = np.zeros(array_shape)
            self.modelled_emiss_pit = np.zeros(array_shape)
            ## final emission
            ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3_flux] 
            ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to 
            ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
            ## during land spreading stages, such as canopy recapture, deap injection...
            self.NH3_flux_slat = np.zeros(array_shape)
            self.NH3_flux_pit = np.zeros(array_shape)

            ## pools left after housing; transfer to storage
            self.manure_pool_pit_to_storage = np.zeros(array_shape)
            self.avail_N_pool_pit_to_storage = np.zeros(array_shape)
            self.resist_N_pool_pit_to_storage = np.zeros(array_shape)
            self.unavail_N_pool_pit_to_storage = np.zeros(array_shape)
            self.TAN_pool_pit_to_storage = np.zeros(array_shape)
            self.Total_water_pool_pit_to_storage = np.zeros(array_shape)
        
        else:
            self.evap_ = np.zeros(array_shape)
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
            ## TAN pool in ug/m2
            self.TAN_pool_ug = np.zeros(array_shape)
            ## TAN pool conc (aqueous phase)
            self.TAN_amount = np.zeros(array_shape)
            ## TAN pool in molar concentration
            self.TAN_amount_M = np.zeros(array_shape)
            ## total water pool of the system (manure water; urine+manure water+[washing water])
            self.Total_water_pool = np.zeros(array_shape)
            ## ratio of [NH4+]/[H+] in the system
            self.Gamma_manure = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in molar mass)
            self.NH3_gas_M = np.zeros(array_shape)
            ## surface NH3 concentrtion at equilirium (in ug)
            self.NH3_gas_ug = np.zeros(array_shape)
            ## emission potential
            self.modelled_emiss = np.zeros(array_shape)
            ## final emission
            ## the difference between [modelled_emiss] (so-called "emission potential") and [NH3_flux] 
            ## (so-called "final emission") is that whether there is other factors that limits the amount of flux to 
            ## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
            ## during land spreading stages, such as canopy recapture, deap injection...
            self.NH3_flux = np.zeros(array_shape)

            ## pools left after housing; transfer to storage
            self.manure_pool_to_storage = np.zeros(array_shape)
            self.avail_N_pool_to_storage = np.zeros(array_shape)
            self.resist_N_pool_to_storage = np.zeros(array_shape)
            self.unavail_N_pool_to_storage = np.zeros(array_shape)
            self.TAN_pool_to_storage = np.zeros(array_shape)
            self.Total_water_pool_to_storage = np.zeros(array_shape)
            self.NH3_flux_from_barn = np.zeros(array_shape)

            if housing_type.lower() == 'barn':
                ## urea pool
                self.urea_pool = np.zeros(array_shape)
                self.urea_pool_to_storage = np.zeros(array_shape)
            elif housing_type.lower() == 'poultry house':
                ## 
                self.N_birds_UAN = N_input * 1000/365 * f_uan
                self.manure = self.N_birds_UA/f_excretn
                ## uric acid pool
                self.UA_pool = np.zeros(array_shape)
                self.UA_pool_to_storage = np.zeros(array_shape)
                ## uric acid hydrolysis rate
                self.daily_ua_conv_factor = ua_hydrolysis_rate(temp=self.T_sim,rhum=self.RH_sim,ph=self.pH)
    
    ## fraction of the floor area, assuming 40% is solid floor
    fslat = 0.4
    ## fraction of the gap of the area, assuming 60% (1-40%) is gap
    fgap = 0.6
    ## background NH3 level, ug/m3
    X_air = 0.300
    ## cleaning frequency; i.e. duration of manure in houses until collected and transfered to storage
    cleaning_cycl = 90
    cleaning_freq = int(cleaning_cycl * nhours / timestep)

    ## housing environmental conditions; and diagnostic variables
    def sim_env(self,housing_type):
        ## housing environmental conditions
        if housing_type.lower() == 'slat/pit house':
            print("House with slatted floor")
            self.T_sim, self.u_sim, self.RH_sim = housing_env(temp_data,rhum_data,livestock,production_system)
        elif housing_type.lower() == 'barn':
            print("Naturally ventilated barn")
            self.T_sim, self.u_sim = barn_env(temp_data,wind_data)
            self.RH_sim = rhum_data
        elif housing_type.lower() == 'poultry house':
            print("Poultry house")
            self.T_sim, self.u_sim, self.RH_sim = housing_env(temp_data,rhum_data,livestock,production_system)
        else:
            print("Other housing systems")
            self.T_sim = temp_data
            self.u_sim = wind_data
            self.RH_sim = rhum_data
        
        ## coverting xarray to numpy array
        self.T_sim = xr_to_np(self.T_sim)
        self.u_sim = xr_to_np(self.u_sim)
        self.RH_sim = xr_to_np(self.RH_sim)

        ###################################################
        ## calculating diagnostic variables
        ###################################################
        if housing_type.lower() == 'slat/pit house':
            ## daily evaporation; aerodynamic method; (g/m^2)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            self.evap_slat = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house)*self.fslat*1000
            self.evap_pit = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=0.6,zo=zo_house)*1000
            ## housing resistance; s/m
            self.R_star_slat = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house))
            ## correction factor for indoor resistance R_star
            F_rcorret = 1.0
            self.R_star_pit = F_rcorret * resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap_pit/1000)
        else:
            ## daily evaporation; aerodynamic method; (g/m^2)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            self.evap = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_house)*1000
            ## housing resistance; s/m
            self.R_star = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap/1000)
            if housing_type.lower() == 'poultry house':
                ## housing resistance; 16700 s/m for poultry houses
                self.R_star[:] = 16700

        ## mositure equilirium, mositure content of manure
        self.mois_coeff = (-np.log(1.01-(self.RH_sim/100))/(0.0000534*(self.T_sim+273.15)))**(1/1.41)
        ## daily urea hydrolysis rate
        self.daily_urea_hydro_rate = urea_hydrolysis_rate(temp=self.T_sim,delta_t=timestep)
        ## daily decomposition rate of available and resistant N components
        self.daily_Na_decomp_rate, self.daily_Nr_decomp_rate = N_pools_decomp_rate(temp=self.T_sim, delta_t=timestep)
        ## Henry's law constant and dissociation equilibria; Eq.6
        self.Henry_constant = (161500/(self.T_sim + 273.15)) * np.exp(-10380/(self.T_sim + 273.15))
        ## dissociation constant of NH4+; 298.15 K is 25 degC (room temperature)
        self.k_NH4 = 5.67e-10*np.exp(-6286*(1/(self.T_sim + 273.15)-1/298.15))

    ###################################
    ## define model functions
    ###################################
    ## simulation: for houses with slatted floor (two-source model: slat+pit); industrial production system
    def slat_pit_housing_sim(self,start_day_idx,end_day_idx,f_slat,f_gap):
        for dd in np.arange(start_day_idx,end_day_idx-1):
            ## daily manure and urine on per unit area
            self.manure[dd+1] = self.dmanure
            self.urine[dd+1] = self.durine
            
            ## manure pool
            ## the amount of manure left on slat depends on the solid floor ratio, here f_slat = 40%
            ## we assume manure left on slat will drop to the pit on the next day
            self.manure_pool_slat[dd+1] = self.manure[dd+1]*f_slat
            self.manure_pool_pit[dd+1] = self.manure_pool_pit[dd] + self.manure[dd+1]*f_gap + self.manure_pool_slat[dd]

            ## N input in multiple forms
            self.urea[dd+1] = self.durea
            self.avail_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_avail
            self.resist_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_resist
            self.unavail_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_unavail

            ## TAN production from urea hydrolysis from slat and pit, respectively
            ## the conditions (T, RH) of slat and pit are assume to be consistent, so is the urea hydrolysis rate
            ## and the N decomposition rate from dung
            self.TAN_prod_pit[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool_pit[dd]+\
                                    self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool_pit[dd] +\
                                    self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool_pit[dd]

            ## Urea pool
            ## slat urea pool is reset at a daily basis
            self.urea_pool_slat[dd+1] = self.urea[dd+1]*f_slat
            ## urea left on slat will go to the pit urea pool
            self.urea_pool_pit[dd+1] = self.urea_pool_pit[dd] * (1 - self.daily_urea_hydro_rate[dd+1]) + \
                                    self.urea[dd+1]*f_gap + self.urea_pool_slat[dd]*(1-self.daily_urea_hydro_rate[dd])

            ## Org N pools in various forms
            self.avail_N_pool_slat[dd+1] = self.avail_N[dd+1]*f_slat
            self.avail_N_pool_pit[dd+1] = self.avail_N_pool_pit[dd]*(1 - self.daily_Na_decomp_rate[dd+1])+ \
                                        self.avail_N[dd+1]*f_gap + \
                                        self.avail_N_pool_slat[dd]*(1-self.daily_Na_decomp_rate[dd])

            self.resist_N_pool_slat[dd+1] = self.resist_N[dd+1]*f_slat
            self.resist_N_pool_pit[dd+1] = self.resist_N_pool_pit[dd]* (1 - self.daily_Nr_decomp_rate[dd+1])+ \
                                        self.resist_N[dd+1]*f_gap+ \
                                        self.resist_N_pool_slat[dd]*(1 - self.daily_Nr_decomp_rate[dd])

            self.unavail_N_pool_slat[dd+1] = self.unavail_N[dd+1]*f_slat
            self.unavail_N_pool_pit[dd+1] = self.unavail_N_pool_pit[dd] + self.unavail_N[dd+1] *f_slat + self.unavail_N_pool_slat[dd]

            ## TAN production from urea hydrolysis and org N decomposition
            self.TAN_prod_slat[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool_slat[dd+1]+ \
                            self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool_slat[dd+1] + \
                            self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool_slat[dd+1]

            ## water from the fresh dung
            self.manure_initwc_slat[dd+1] = self.manure_wc*f_slat
            self.manure_initwc_pit[dd+1] = self.manure_wc*f_gap

            ## water amount when mositure content reach equilibrium
            self.manure_minwc_slat[dd+1] = self.manure_pool_slat[dd+1]*self.mois_coeff[dd+1] / 100
            self.manure_minwc_pit[dd+1] = self.manure_pool_pit[dd+1]*self.mois_coeff[dd+1]/100

            ## water pool of slat and pit
            ## slat water pool
            water_slat_idx = self.urine[dd+1]*f_slat + self.manure_initwc_slat[dd+1] - self.evap_slat[dd+1] -self. manure_minwc_slat[dd+1]
            self.Total_water_pool_slat[dd+1][water_slat_idx>=0] = self.urine[dd+1][water_slat_idx>=0]*f_slat +\
                                                            self.manure_initwc_slat[dd+1][water_slat_idx>=0] -\
                                                            self.evap_slat[dd+1][water_slat_idx>=0]
            self.Total_water_pool_slat[dd+1][water_slat_idx<0] = self.manure_minwc_slat[dd+1][water_slat_idx<0]
            ## pit water pool
            water_pit_idx = self.Total_water_pool_pit[dd]-self.evap_pit[dd]-self.manure_minwc_pit[dd+1]
            self.Total_water_pool_pit[dd+1][water_pit_idx>0] = self.Total_water_pool_pit[dd][water_pit_idx>0] +\
                                                            self.urine[dd+1][water_pit_idx>0]*f_gap +\
                                                            self.manure_initwc_pit[dd+1][water_pit_idx>0] -\
                                                            self.evap_pit[dd][water_pit_idx>0] +\
                                                            self.Total_water_pool_slat[dd][water_pit_idx>0]
            self.Total_water_pool_pit[dd+1][water_pit_idx<=0] = self.manure_minwc_pit[dd+1][water_pit_idx<=0]+\
                                                            self.urine[dd+1][water_pit_idx<=0]*f_gap +\
                                                            self.manure_initwc_pit[dd+1][water_pit_idx<=0] +\
                                                            self.Total_water_pool_slat[dd][water_pit_idx<=0]

            ## TAN pool of slat and pit
            self.TAN_pool_slat[dd+1] = self.TAN_prod_slat[dd+1]
            ## TAN left on slat will go to the pit TAN pool
            TAN_pit_idx = self.TAN_pool_pit[dd] - self.NH3_flux_pit[dd]
            self.TAN_pool_pit[dd+1][TAN_pit_idx>0] = TAN_pit_idx[TAN_pit_idx>0]+self.TAN_prod_pit[dd+1][TAN_pit_idx>0]+\
                                                self.TAN_pool_slat[dd][TAN_pit_idx>0]-self.NH3_flux_slat[dd][TAN_pit_idx>0]
            self.TAN_pool_pit[dd+1][TAN_pit_idx<=0] = self.TAN_prod_pit[dd+1][TAN_pit_idx<=0]+self.TAN_pool_slat[dd][TAN_pit_idx<=0]-\
                                                self.NH3_flux_slat[dd][TAN_pit_idx<=0]

            ## TAN pool in ug
            self.TAN_pool_ug_slat[dd+1] = self.TAN_pool_slat[dd+1] * 1e6
            self.TAN_pool_ug_pit[dd+1] = self.TAN_pool_pit[dd+1] * 1e6

            ## TAN conc on slat and in pit
            self.TAN_amount_slat[dd+1][self.Total_water_pool_slat[dd+1]==0] = 0
            self.TAN_amount_slat[dd+1][self.Total_water_pool_slat[dd+1]!=0] = self.TAN_pool_slat[dd+1][self.Total_water_pool_slat[dd+1]!=0]/\
                                                        self.Total_water_pool_slat[dd+1][self.Total_water_pool_slat[dd+1]!=0]
            self.TAN_amount_pit[dd+1][self.Total_water_pool_pit[dd+1]==0] = 0
            self.TAN_amount_pit[dd+1][self.Total_water_pool_pit[dd+1]!=0] = self.TAN_pool_pit[dd+1][self.Total_water_pool_pit[dd+1]!=0]/\
                                                        self.Total_water_pool_pit[dd+1][self.Total_water_pool_pit[dd+1]!=0]

            ## TAN molar conc
            self.TAN_amount_M_slat[dd+1] = self.TAN_amount_slat[dd+1]/14*1000
            self.TAN_amount_M_pit[dd+1] = self.TAN_amount_pit[dd+1]/14*1000

            ## Gamma value
            self.Gamma_manure_slat[dd+1] =  self.TAN_amount_M_slat[dd+1]/(self.cc_H + self.k_NH4[dd+1])
            self.Gamma_manure_pit[dd+1] =  self.TAN_amount_M_pit[dd+1]/(self.cc_H + self.k_NH4[dd+1])

            ## Gaseous NH3 at the surface of slat TAN pool and pit TAN pool
            self.NH3_gas_M_slat[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure_slat[dd+1]
            self.NH3_gas_M_pit[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure_pit[dd+1]

            ## in ug
            self.NH3_gas_ug_slat[dd+1] = self.NH3_gas_M_slat[dd+1]*14*1e9
            self.NH3_gas_ug_pit[dd+1] = self.NH3_gas_M_pit[dd+1]*14*1e9

            ## determining the maximum emission; emission cannot exceed TAN pool
            emiss_slat_idx = (self.NH3_gas_ug_slat[dd+1]*3600*timestep/self.R_star_slat[dd+1]) - self.TAN_pool_ug_slat[dd+1]
            self.modelled_emiss_slat[dd+1][emiss_slat_idx>=0] = self.TAN_pool_ug_slat[dd+1][emiss_slat_idx>=0]
            self.modelled_emiss_slat[dd+1][emiss_slat_idx<0] = self.NH3_gas_ug_slat[dd+1][emiss_slat_idx<0]*3600*timestep/\
                                                        self.R_star_slat[dd+1][emiss_slat_idx<0]

            emiss_pit_idx = (self.NH3_gas_ug_pit[dd+1]*3600*timestep/self.R_star_pit[dd+1]) - self.TAN_pool_ug_pit[dd+1]
            self.modelled_emiss_pit[dd+1][emiss_pit_idx>=0] = self.TAN_pool_ug_pit[dd+1][emiss_pit_idx>=0]
            self.modelled_emiss_pit[dd+1][emiss_pit_idx<0] = self.NH3_gas_ug_pit[dd+1][emiss_pit_idx<0]*3600*timestep/\
                                                        self.R_star_pit[dd+1][emiss_pit_idx<0]

            ## final emission flux from slat and pit; flux are additive
            self.NH3_flux_slat[dd+1] = self.modelled_emiss_slat[dd+1]/1e6
            self.NH3_flux_pit[dd+1] = self.modelled_emiss_pit[dd+1]/1e6
        return

    ## simulation: for houses with concrete floor (singe-source model); 
    ## intermediate/backyard production system (pig/ruminants)
    def barn_housing_sim(self,start_day_idx,end_day_idx):
        for dd in np.arange(start_day_idx,end_day_idx-1):
            ## daily manure and urine on per unit area
            self.manure[dd+1] = self.dmanure
            self.urine[dd+1] = self.durine
            ## manure pool
            self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]
    
            ## N input in multiple forms
            self.urea[dd+1] = self.durea
            self.avail_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_avail
            self.resist_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_resist
            self.unavail_N[dd+1] = (self.dmanure_N + (self.durine_N-self.durea))*f_unavail
            
            ## Urea pool
            self.urea_pool[dd+1] = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd]) + self.urea[dd+1]
            ## Org N pools in various forms
            self.avail_N_pool[dd+1] = self.avail_N_pool[dd]*(1 - self.daily_Na_decomp_rate[dd]) + self.avail_N[dd+1]
            self.resist_N_pool[dd+1] = self.resist_N_pool[dd]* (1 - self.daily_Nr_decomp_rate[dd]) + self.resist_N[dd+1]            
            self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1] 
            
            ## TAN production from urea hydrolysis and the N decomposition rate from dung
            self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd+1]+\
                                    self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd+1] +\
                                    self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd+1]
            
            ## water from the fresh dung
            self.manure_initwc[dd+1] = self.manure_wc
            ## water amount when mositure content reach equilibrium
            self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100
            ## water pool 
            water_idx = self.Total_water_pool[dd] + self.urine[dd+1] + self.manure_initwc[dd+1] - \
                                self.evap[dd+1] - self.manure_minwc[dd+1] 
            self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] +\
                                                            self.urine[dd+1][water_idx>0] +\
                                                            self.manure_initwc[dd+1][water_idx>0] -\
                                                            self.evap[dd+1][water_idx>0] 
            self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd+1][water_idx<=0] 

            ## TAN pool 
            self.TAN_pool[dd+1] = self.TAN_pool[dd]-self.NH3_flux[dd]+self.TAN_prod[dd+1]
            ## TAN pool in ug
            self.TAN_pool_ug[dd+1] = self.TAN_pool[dd+1] * 1e6

            ## TAN conc 
            self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                                                        self.Total_water_pool[dd+1][self.Total_water_pool[dd+1]!=0]
            ## TAN molar conc
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/14*1000
            ## Gamma value 
            self.Gamma_manure[dd+1] =  self.TAN_amount_M[dd+1]/(self.cc_H + self.k_NH4[dd+1])
            ## Gaseous NH3 at the surface 
            self.NH3_gas_M[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure[dd+1]
            ## in ug
            self.NH3_gas_ug[dd+1] = self.NH3_gas_M[dd+1]*14*1e9
            
            ## determining the maximum emission; emission cannot exceed TAN pool
            emiss_idx = (self.NH3_gas_ug[dd+1]*3600*timestep/self.R_star[dd+1]) - self.TAN_pool_ug[dd+1]
            self.modelled_emiss[dd+1][emiss_idx>=0] = self.TAN_pool_ug[dd+1][emiss_idx>=0]
            self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*timestep/self.R_star[dd+1][emiss_idx<0]
            ## final emission flux 
            self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return

    ## simulation: for poultry housing;
    ## broilers/layers; housing resistance is assumed to be 16700 s/m
    def poultry_housing_sim(self,start_day_idx,end_day_idx):
        for dd in np.arange(start_day_idx,end_day_idx-1):
            ## excretion pool; Eq.1
            self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure
            self.TAN_prod[dd+1] = self.daily_ua_conv_factor[dd+1] * self.UA_pool[dd]
            ## UA pool; Eq.2
            UA_idx = self.UA_pool[dd] - self.TAN_prod[dd+1]
            self.UA_pool[dd+1][UA_idx>0] = self.UA_pool[dd][UA_idx>0] + self.N_birds_UAN[UA_idx>0] - self.TAN_prod[dd+1][UA_idx>0]
            self.UA_pool[dd+1][UA_idx<=0] = self.N_birds_UAN[UA_idx<=0]

            ## manure_water calculation; Eq.9  
            self.manure_water[dd+1][self.mois_coeff[dd+1]>=5] = self.manure_pool[dd+1][self.mois_coeff[dd+1]>=5] * self.mois_coeff[dd+1][self.mois_coeff[dd+1]>=5] / 100
            self.manure_water[dd+1][self.mois_coeff[dd+1]<5] = 0
            self.Total_water_pool[dd+1] = self.manure_water[dd+1]
            
            ## TAN pool; Eq.3
            TAN_idx = self.TAN_pool[dd] - self.NH3_Flux[dd]
            self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]
            self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0] + self.TAN_prod[dd+1][TAN_idx>0]
            self.TAN_pool_ug[dd+1] = self.TAN_pool[dd+1] * 1e6
            
            ## TAN concentration
            self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
            self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0] / self.Total_water_pool[dd+1][self.Total_water_pool[dd+1]!=0]
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1] / 14 * 1000

            ## Gamma value; Eq.7
            self.Gamma_manure[dd+1] =  self.TAN_amount_M[dd+1] / (self.cc_H + self.k_NH4[dd+1])
            
            ## Gaseous concentration of NH3 at the surface; Eq.6
            self.NH3_gas_M[dd+1] = self.Henry_constant[dd+1] * self.Gamma_manure[dd+1]
            self.NH3_gas_ug[dd+1] = self.NH3_gas_M[dd+1] * 14 * 1e9

            ## Emissions; see Eq.8 in general, and Sect.2.2.2; Eq.14~17
            emiss_idx = (self.NH3_gas_ug[dd+1]*3600*timestep/self.R_star) - self.TAN_pool_ug[dd+1]
            self.modelled_emiss[dd+1][emiss_idx>0] = self.TAN_pool_ug[dd+1][emiss_idx>0]
            self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*timestep/self.R_star
            self.NH3_Flux[dd+1] = self.Modelled_emiss[dd+1]/1e6 
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
            self.TAN_pool_ug_slat[:] = 0.0
            self.TAN_pool_ug_pit[:] = 0.0
            self.TAN_amount_slat[:] = 0.0
            self.TAN_amount_pit[:] = 0.0
            self.TAN_amount_M_slat[:] = 0.0
            self.TAN_amount_M_pit[:] = 0.0
            self.Total_water_pool_slat[:] = 0.0
            self.Total_water_pool_pit[:] = 0.0
            self.Gamma_manure_slat[:] = 0.0
            self.Gamma_manure_pit[:] = 0.0
            self.NH3_gas_M_slat[:] = 0.0
            self.NH3_gas_M_pit[:] = 0.0
            self.NH3_gas_ug_slat[:] = 0.0
            self.NH3_gas_ug_pit[:] = 0.0
            self.modelled_emiss_slat[:] = 0.0
            self.modelled_emiss_pit[:] = 0.0
            self.NH3_flux_slat[:] = 0.0
            self.NH3_flux_pit[:] = 0.0

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
            self.TAN_pool_ug[:] = 0.0
            self.TAN_amount[:] = 0.0
            self.TAN_amount_M[:] = 0.0
            self.Total_water_pool[:] = 0.0
            self.Gamma_manure[:] = 0.0
            self.NH3_gas_M[:] = 0.0
            self.NH3_gas_ug[:] = 0.0
            self.modelled_emiss[:] = 0.0
            self.NH3_flux[:] = 0.0
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
            dd2 = self.NH3_flux_slat[-1]
            dd3 = self.NH3_flux_pit[-1]
            ee2 = self.Total_water_pool_slat[-1]
            ee3 = self.Total_water_pool_pit[-1]

            self.manure_pool_slat[0] = aa2
            self.manure_pool_pit[0] = aa3
            self.urea_pool_slat[0] = bb2
            self.urea_pool_pit[0] = bb3
            self.TAN_pool_slat[0] = cc2
            self.TAN_pool_pit[0] = cc3
            self.NH3_flux_slat[0] = dd2
            self.NH3_flux_pit[0] = dd3
            self.Total_water_pool_slat[0] = ee2
            self.Total_water_pool_pit[0] = ee3  
        else:
            aa1 = self.manure_pool[-1]
            
            cc1 = self.TAN_pool[-1]
            dd1 = self.NH3_flux[-1]
            ee1 = self.Total_water_pool[-1]

            self.manure_pool[0] = aa1
            
            self.TAN_pool[0] = cc1
            self.NH3_flux[0] = dd1
            self.Total_water_pool[0] = ee1

            if housing_type.lower() == 'barn':
                bb1 = self.urea_pool[-1]
                self.urea_pool[0] = bb1
            elif housing_type.lower() == 'poultry house':
                bb1 = self.UA_pool[-1]
                self.UA_pool[0] = bb1

    ## initialisation for house (barn) cleaning; all variables at the cleaning day are reset to zero
    def cleaning_pit(self,day_idx,housing_type):
        if housing_type.lower() == 'slat/pit house':
            ## pools: from housing to MMS
            self.manure_pool_pit_to_storage[day_idx] = self.manure_pool_pit[day_idx]+self.manure_pool_slat[day_idx]
            self.avail_N_pool_pit_to_storage[day_idx] = self.avail_N_pool_pit[day_idx]
            self.resist_N_pool_pit_to_storage[day_idx] = self.resist_N_pool_pit[day_idx]
            self.unavail_N_pool_pit_to_storage[day_idx] = self.unavail_N_pool_pit[day_idx]
            self.TAN_pool_pit_to_storage[day_idx] = self.TAN_pool_pit[day_idx]+self.urea_pool_pit[day_idx]
            self.Total_water_pool_pit_to_storage[day_idx] = self.Total_water_pool_pit[day_idx]
            
            self.manure_pool_slat[day_idx] = 0.0
            self.manure_pool_pit[day_idx] = 0.0
            self.manure_initwc_slat[day_idx] = 0.0
            self.manure_initwc_pit[day_idx] = 0.0
            self.urea_pool_slat[day_idx] = 0.0
            self.urea_pool_pit[day_idx] = 0.0
            self.avail_N_pool_slat[day_idx] = 0.0
            self.avail_N_pool_pit[day_idx] = 0.0
            self.resist_N_pool_slat[day_idx] = 0.0
            self.resist_N_pool_pit[day_idx] = 0.0
            self.unavail_N_pool_slat[day_idx] = 0.0
            self.unavail_N_pool_pit[day_idx] = 0.0
            self.TAN_pool_slat[day_idx] = 0.0
            self.TAN_pool_pit[day_idx] = 0.0
            self.TAN_pool_ug_slat[day_idx] = 0.0
            self.TAN_pool_ug_pit[day_idx] = 0.0
            self.Total_water_pool_slat[day_idx] = 0.0
            self.Total_water_pool_pit[day_idx] = 0.0
            self.modelled_emiss_slat[day_idx] = 0.0
            self.modelled_emiss_pit[day_idx] = 0.0
            self.NH3_flux_slat[day_idx] = 0.0
            self.NH3_flux_pit[day_idx] = 0.0
        else:
            ## pools: from housing to MMS
            self.manure_pool_to_storage[day_idx] = self.manure_pool[day_idx]
            self.avail_N_pool_to_storage[day_idx] = (1 - self.daily_Na_decomp_rate[day_idx])*self.avail_N_pool[day_idx]
            self.resist_N_pool_to_storage[day_idx] = (1 - self.daily_Nr_decomp_rate[day_idx])*self.resist_N_pool[day_idx]
            self.unavail_N_pool_to_storage[day_idx] = self.unavail_N_pool[day_idx]
            self.TAN_pool_to_storage[day_idx] = self.TAN_pool[[day_idx]]-self.NH3_flux[day_idx]
            self.Total_water_pool_to_storage[day_idx] = self.Total_water_pool[day_idx]
            self.NH3_flux_from_barn[day_idx] = self.NH3_flux[day_idx]

            self.manure_pool[day_idx] = 0.0
            self.manure_water[day_idx] = 0.0
            self.avail_N_pool[day_idx] = 0.0
            self.resist_N_pool[day_idx] = 0.0
            self.unavail_N_pool[day_idx] = 0.0
            self.TAN_pool[day_idx] = 0.0
            self.Total_water_pool[day_idx] = 0.0
            self.modelled_emiss[day_idx] = 0.0
            self.NH3_flux[day_idx] = 0.0

            if housing_type.lower() == 'barn':
                self.urea_pool_to_storage[day_idx] = (1 - self.daily_urea_hydro_rate[day_idx])*self.urea_pool[day_idx]
                self.urea_pool[day_idx] = 0.0
            elif housing_type.lower() == 'poultry house':
                self.UA_pool_to_storage[day_idx] = self.UA_pool[day_idx]
                self.UA_pool[day_idx] = 0.0

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
            self.NH3_flux_from_barn[:] = 0.0
            if housing_type.lower() == 'barn':
                self.urea_pool_to_storage[:] = 0.0
            elif housing_type.lower() == 'poultry house':
                self.UA_pool_to_storage[:] = 0.0
    
    def slat_pit_sim_main(self,housing_type,start_idx,end_idx,cleaning_frequency):
        self.sim_env(housing_type)
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.slat_pit_housing_sim(dd,dd+cleaning_frequency+1,self.fslat,self.fgap)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.slat_pit_housing_sim(dd,end_idx,self.fslat,self.fgap)
        self.housing_2nd_init(housing_type)
        self.slat_pit_housing_sim(0,start_idx,self.fslat,self.fgap)
        return

    def barn_sim_main(self,housing_type,start_idx,end_idx,cleaning_frequency):
        self.sim_env(housing_type)
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.barn_housing_sim(dd,dd+cleaning_frequency+1)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.barn_housing_sim(dd,end_idx)
        self.housing_2nd_init(housing_type)
        self.barn_housing_sim(0,start_idx)
        return

    def poultry_house_sim_main(self,housing_type,start_idx,end_idx,cleaning_frequency):
        self.sim_env(housing_type)
        self.housing_init(housing_type)
        self.housing_to_storage_init(housing_type)
        for dd in np.arange(start_idx,end_idx,cleaning_frequency):
            if dd + cleaning_frequency < end_idx:
                self.poultry_housing_sim(dd,dd+cleaning_frequency+1)
                self.cleaning_pit(dd+cleaning_frequency,housing_type)
            else:
                self.poultry_housing_sim(dd,end_idx)
        self.housing_2nd_init(housing_type)
        self.poultry_housing_sim(0,start_idx)
        return 

    
