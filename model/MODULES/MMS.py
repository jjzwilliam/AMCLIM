####################################
## import python packages
####################################
import sys
import os
from pathlib import Path

####################################
## import essential AMCLIM modules
####################################
'''file_path = os.getcwd()
module_path = str(Path(file_path).parent)
print(module_path)

if module_path not in sys.path:
    sys.path.append(module_path)'''
from CONFIG.config import *
from MODULES.PARAMETERS import *
from MODULES.FUNC import *
#sys.exit()

####################################
## calculating diagnostic variables
####################################
T_barn, u_barn = barn_env(temp_data,wind_data)
RH_barn = rhum_data

T_out = temp_data
u_out = wind_data
RH_out = rhum_data
evap_out = evap_data

## coverting xarray to numpy array
T_barn = xr_to_np(T_barn)
T_out = xr_to_np(T_out)
u_barn = xr_to_np(u_barn)
u_out = xr_to_np(u_barn)
RH_barn = xr_to_np(RH_barn)
RH_out = xr_to_np(RH_out)
evap_out = xr_to_np(evap_out)
## daily evaporation; aerodynamic method; (g/m^2)
## Q_vent above pit is set to be 0.6 m/s (full efficiency)
evap_barn = water_evap_a(temp=T_barn,rhum=RH_barn,u=u_barn)*1000
## storage resistance; s/m
R_star_storage = resistance_water_air(temp=T_in,rhum=RH,evap_flux=water_evap_a(temp=T_in,rhum=RH,u=u_in))
## correction factor for indoor resistance R_star
F_rcorret = 1.0


## pH and H+ ions concentration
cc_H = np.float(10**(-pH))
## mositure equilirium, mositure content of manure
mois_coeff_barn = (-np.log(1.01-(RH_barn/100))/(0.0000534*(T_barn+273.15)))**(1/1.41)
mois_coeff_out = (-np.log(1.01-(RH_out/100))/(0.0000534*(T_out+273.15)))**(1/1.41)
if livestock.lower()=="broiler":
    ## daily UA hydrolysis rate
    daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH,ph=pH)
elif livestock.lower()=="layer":
    daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH,ph=pH)
else:
    ## daily urea hydrolysis rate
    daily_urea_hydro_rate = urea_hydrolysis_rate(temp=T_in,delta_t=24)
    ## daily decomposition rate of available and resistant N components
    daily_Na_decomp_rate, daily_Nr_decomp_rate = N_pools_decomp_rate(temp=T_in, delta_t=24)

## Henry's law constant and dissociation equilibria; Eq.6
Henry_constant = (161500/(Total_Temp_housing + 273.15)) * np.exp(-10380/(Total_Temp_housing + 273.15))
## dissociation constant of NH4+
k_NH4 = 5.67e-10*np.exp(-6286*(1/(Total_Temp_housing + 273.15)-1/298.15))
## background NH3 level, ug/m3
X_air = 0.300
## indoor NH3 level, read from datasets
# NH3_inconc_ug = NH3_inconc * 1e6

##################################
## define 
##################################
class MMS_module:
    def __init__(self,array_shape):
        ## feces input from housing
        self.manure_added = np.zeros(array_shape)
        ## feces pool
        self.manure_pool = np.zeros(array_shape)
        ## min amount of water in manure (function of T, RH)
        ## equilibrium moisture content, water content cannot be lower than this threshold under "natural" processes,
        ## i.e., without drying processes
        self.manure_water = np.zeros(array_shape)
        self.manure_minwc = np.zeros(array_shape)
        ## water amount in fresh feces
        self.manure_initwc = np.zeros(array_shape)
        ## urine input from housing
        self.urine_added = np.zeros(array_shape)
        ## urea input from housing
        self.urea_added = np.zeros(array_shape)
        ## urea pool
        self.urea_pool = np.zeros(array_shape)
        ## uric acid input from housing
        self.UA_added = np.zeros(array_shape)
        ## uric acid pool
        self.UA_pool = np.zeros(array_shape)
        ## input of available N component from housing
        self.avail_N_added = np.zeros(array_shape)
        ## Nitrogen pool that is available (easily) to form TAN
        self.avail_N_pool = np.zeros(array_shape)
        ## input of resistant N component from housing
        self.resist_N_added = np.zeros(array_shape)
        ## Nitrogen pool that is resistant (slowly) to form TAN
        self.resist_N_pool = np.zeros(array_shape)
        ## input of unavailable N component from housing
        self.unavail_N_added = np.zeros(array_shape)
        ## Nitrogen pool that is unavailable (cannot) to form TAN
        self.unavail_N_pool = np.zeros(array_shape)
        ## TAN added from housing
        self.TAN_added = np.zeros(array_shape)
        ## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool
        self.TAN_prod = np.zeros(array_shape)
        ## TAN pool
        self.TAN_pool = np.zeros(array_shape)
        ## TAN pool in ug/m2
        self.TAN_pool_ug = np.zeros(array_shape)
        ## TAN pool conc (aqueous phase)
        self.TAN_amount = np.zeros(array_shape)
        ## TAN pool in molar concentration
        self.TAN_amout_M = np.zeros(array_shape)
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

    def sim_env(self):
        if MMS_type == '':
            T_barn, u_barn = barn_env(temp_data,wind_data)
            RH_barn = rhum_data
            ## daily evaporation; aerodynamic method; (g/m^2)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            evap_barn = water_evap_a(temp=T_barn,rhum=RH_barn,u=u_barn)*1000
        else:
            T_out = temp_data
            u_out = wind_data
            RH_out = rhum_data
            evap_out = evap_data
        return

    ## mositure equilirium, mositure content of manure
    mois_coeff_barn = (-np.log(1.01-(RH_barn/100))/(0.0000534*(T_barn+273.15)))**(1/1.41)
    if livestock.lower()=="broiler":
    ## daily UA hydrolysis rate
        daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH,ph=pH)
    elif livestock.lower()=="layer":
        daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH,ph=pH)
    else:
        ## daily urea hydrolysis rate
        daily_urea_hydro_rate = urea_hydrolysis_rate(temp=T_in,delta_t=24)
        ## daily decomposition rate of available and resistant N components
        daily_Na_decomp_rate, daily_Nr_decomp_rate = N_pools_decomp_rate(temp=T_in, delta_t=24)

    ## Henry's law constant and dissociation equilibria; Eq.6
    Henry_constant = (161500/(Total_Temp_housing + 273.15)) * np.exp(-10380/(Total_Temp_housing + 273.15))
    ## dissociation constant of NH4+
    k_NH4 = 5.67e-10*np.exp(-6286*(1/(Total_Temp_housing + 273.15)-1/298.15))

    def MMS_barn_sim(self,start_day_idx,end_day_idx):
        for dd in np.arange(start_day_idx,end_day_idx-1):
            ## daily manure and urine on per unit area
            self.manure_added[dd+1] = 
            self.urine_added[dd+1] = 

            ## manure pool
            self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]

            ## N input in multiple forms
            self.urea_added[dd+1] = 
            self.avail_N_added[dd+1] = 
            self.resist_N_added[dd+1] = 
            self.unavail_N_added[dd+1] = 

            ## TAN production from urea hydrolysis and the N decomposition rate from dung
            self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                    self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                    self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
            ## TAN from housing to storage
            self.TAN_added[dd+1] = 

            ## Urea pool
            urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1])
            self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
            self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

            ## Org N pools in various forms
            self.avail_N_pool[dd+1] = self.avail_N_pool[dd]*(1 - self.daily_Na_decomp_rate[dd+1]) + self.avail_N_added[dd+1]
            self.resist_N_pool[dd+1] = self.resist_N_pool[dd]* (1 - self.daily_Nr_decomp_rate[dd+1]) + self.resist_N_added[dd+1]
            self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1]


            ## water from the fresh dung
            # self.manure_initwc[dd+1] = 

            ## water amount when mositure content reach equilibrium
            self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

            ## water pool
            water_idx = Total_water_pool[dd]-evap[dd]-self.manure_minwc[dd+1]
            self.Total_water_pool[dd+1][water_idx>=0] = self.Total_water_pool[dd][water_idx>0] +\
                                                            self.urine_added[dd+1][water_idx>0] +\
                                                            self.manure_initwc[dd+1][water_idx>0] -\
                                                            self.evap[dd][water_idx>0]
            self.Total_water_pool[dd+1][water_idx<0] = self.manure_minwc[dd+1][water_idx<=0] +\
                                                            self.urine_added[dd+1][water_idx<=0] +\
                                                            self.manure_initwc[dd+1][water_idx<=0]

            ## TAN pool
            self.TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd]
            self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]
            self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]

            ## TAN pool in ug
            self.TAN_pool_ug[dd+1] = self.TAN_pool[dd+1] * 1e6

            ## TAN conc
            self.TAN_amount[dd+1][Total_water_pool[dd+1]==0] = 0
            self.TAN_amount[dd+1][Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][Total_water_pool[dd+1]!=0]/\
                                                        self.Total_water_pool[dd+1][Total_water_pool[dd+1]!=0]

            ## TAN molar conc
            self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/14*1000

            ## Gamma value
            self.Gamma_manure[dd+1] =  self.TAN_amount_M[dd+1]/(self.cc_H + self.k_NH4[dd+1])

            ## Gaseous NH3 at the surface
            self.NH3_gas_M[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure[dd+1]

            ## in ug
            self.NH3_gas_ug[dd+1] = self.NH3_gas_M[dd+1]*14*1e9

            ## determining the maximum emission; emission cannot exceed TAN pool
            emiss_idx = (self.NH3_gas_ug[dd+1]*3600*24/self.R_star[dd+1]) - self.TAN_pool_ug[dd+1]
            self.modelled_emiss[dd+1][emiss_idx>=0] = self.TAN_pool_ug[dd+1][emiss_idx>=0]
            self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*24/\
                                                        self.R_star[dd+1][emiss_idx<0]

            ## final emission flux
            self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return
    
    def MMS_land_sim(self,start_day_idx,end_day_idx):
        # for dd in np.arange(start_day_idx,end_day_idx-1):
        
        return
    
    def MMS_liquid_sim(self,start_day_idx,end_day_idx):
        # for dd in np.arange(start_day_idx,end_day_idx):

        return