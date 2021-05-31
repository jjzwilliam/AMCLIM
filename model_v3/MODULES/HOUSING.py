####################################
## import python packages
####################################
import numpy as np
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
animal_density = 120.0 
massgrid = animal_head*animal_weight
housing_area = animal_head*animal_weight/animal_density
## total animal mass per grid
# massgrid = 1000*excretN_info/(animal_Nrate*365)
# housing_area = massgrid/animal_density
## pig density 2: 1 head/m^2
# animal_density = 1.0
# housing_area = animal_head/animal_density

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
## fraction of the floor area, assuming 40% is solid floor
fslat = 0.4
## fraction of the gap of the area, assuming 60% (1-40%) is gap
fgap = 0.6

## test:
fslat = 1.0
fgap = 0.0

## house surface roughness height; default 2mm
zo_house = 0.002

###################################
## calculating diagnostic fields
####################################
if production_system == 'industrial':
    T_in, u_in, RH_in = housing_env(temp_data,rhum_data,livestock,production_system)
elif production_system == 'intermediate':
    T_in, u_in = barn_env(temp_data,wind_data)
    RH_in = rhum_data
else:
    T_in = temp_data
    u_in = wind_data
    RH_in = rhum_data

## coverting xarray to numpy array
T_in = xr_to_np(T_in)
u_in = xr_to_np(u_in)
RH_in = xr_to_np(RH_in)

## daily evaporation; aerodynamic method; (g/m^2)
## Q_vent above pit is set to be 0.6 m/s (full efficiency)
evap = water_evap_a(temp=T_in,rhum=RH_in,u=u_in,zo=zo_house)*1000
evap_slat = water_evap_a(temp=T_in,rhum=RH_in,u=u_in,zo=zo_house)*fslat*1000
evap_pit = water_evap_a(temp=T_in,rhum=RH_in,u=0.6,zo=zo_house)*1000
## housing resistance; s/m
R_star_slat = resistance_water_air(temp=T_in,rhum=RH_in,evap_flux=water_evap_a(temp=T_in,rhum=RH_in,u=u_in,zo=zo_house))
R_star = R_star_slat

## correction factor for indoor resistance R_star
F_rcorret = 1.0
R_star_pit = F_rcorret * resistance_water_air(temp=T_in,rhum=RH_in,evap_flux=water_evap_a(temp=T_in,rhum=RH_in,u=0.6,zo=zo_house))

## animal waste N info
# excret_N = livestock_N_info(livestock_type='PIG', region='NA')
F_Ncorrect = 1.0
excret_N = F_Ncorrect*excretN_info/housing_area   ## unit: kg N per m^2 per year
excret_N = xr_to_np(excret_N)
# excret_N = F_Ncorrect*excretN_info*(1000/365)/housing_area         ## g N per m^2 per day
## animal waste info
durine_N, durea, dmanure_N, durine, dmanure, manure_wc, pH = livestock_waste_info(livestock_type=livestock, waste_N=excret_N)
durine = durine * 1000
manure_wc = manure_wc * 1000

## pH and H+ ions concentration
cc_H = np.float(10**(-pH))
## mositure equilirium, mositure content of manure
mois_coeff = (-np.log(1.01-(RH_in/100))/(0.0000534*(T_in+273.15)))**(1/1.41)
if livestock.lower()=="broiler":
    ## daily UA hydrolysis rate
    daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH_in,ph=pH)
elif livestock.lower()=="layer":
    daily_ua_conv_factor = ua_hydrolysis_rate(temp=T_in,rhum=RH_in,ph=pH)
else:
    ## daily urea hydrolysis rate
    daily_urea_hydro_rate = urea_hydrolysis_rate(temp=T_in,delta_t=24)
    ## daily decomposition rate of available and resistant N components
    daily_Na_decomp_rate, daily_Nr_decomp_rate = N_pools_decomp_rate(temp=T_in, delta_t=24)

## Henry's law constant and dissociation equilibria; Eq.6
Henry_constant = (161500/(T_in + 273.15)) * np.exp(-10380/(T_in + 273.15))
## dissociation constant of NH4+
k_NH4 = 5.67e-10*np.exp(-6286*(1/(T_in + 273.15)-1/298.15))
## background NH3 level, ug/m3
X_air = 0.300
## indoor NH3 level, read from datasets
# NH3_inconc_ug = NH3_inconc * 1e6



####################################
## define prognostic fields
####################################
## animal number
animal_n = np.zeros(mtrx)
## feces input
manure = np.zeros(mtrx)
## feces pool
manure_pool = np.zeros(mtrx)
manure_pool_slat = np.zeros(mtrx)
manure_pool_pit = np.zeros(mtrx)
## min amount of water in manure (function of T, RH)
## equilibrium moisture content, water content cannot be lower than this threshold under "natural" processes, 
## i.e., without drying processes
manure_water = np.zeros(mtrx)
manure_minwc = np.zeros(mtrx)
manure_minwc_slat = np.zeros(mtrx)
manure_minwc_pit = np.zeros(mtrx)
## water amount in fresh feces
manure_initwc = np.zeros(mtrx)
manure_initwc_slat = np.zeros(mtrx)
manure_initwc_pit = np.zeros(mtrx)
## daily urine input
urine = np.zeros(mtrx)
## urea input
urea = np.zeros(mtrx)
## urea pool
urea_pool = np.zeros(mtrx)
urea_pool_slat = np.zeros(mtrx)
urea_pool_pit = np.zeros(mtrx)
## uric acid pool
UA_pool = np.zeros(mtrx)
## input of available N component
avail_N = np.zeros(mtrx)
## Nitrogen pool that is available (easily) to form TAN
avail_N_pool = np.zeros(mtrx)
avail_N_pool_slat = np.zeros(mtrx)
avail_N_pool_pit = np.zeros(mtrx)
## input of resistant N component
resist_N = np.zeros(mtrx)
## Nitrogen pool that is resistant (slowly) to form TAN
resist_N_pool = np.zeros(mtrx)
resist_N_pool_slat = np.zeros(mtrx)
resist_N_pool_pit = np.zeros(mtrx)
## input of unavailable N component
unavail_N = np.zeros(mtrx)
## Nitrogen pool that is unavailable (cannot) to form TAN
unavail_N_pool = np.zeros(mtrx)
unavail_N_pool_slat = np.zeros(mtrx)
unavail_N_pool_pit = np.zeros(mtrx)
## TAN production from urea hydrolysis, conversion from N_avail and N_resist pool
TAN_prod = np.zeros(mtrx)
TAN_prod_slat = np.zeros(mtrx)
TAN_prod_pit = np.zeros(mtrx)
## TAN pool
TAN_pool = np.zeros(mtrx)
TAN_pool_slat = np.zeros(mtrx)
TAN_pool_pit = np.zeros(mtrx)
## TAN pool in ug/m2
TAN_pool_ug = np.zeros(mtrx)
TAN_pool_ug_slat = np.zeros(mtrx)
TAN_pool_ug_pit = np.zeros(mtrx)
## TAN pool conc (aqueous phase)
TAN_amount = np.zeros(mtrx)
TAN_amount_slat = np.zeros(mtrx)
TAN_amount_pit = np.zeros(mtrx)
## TAN pool in molar concentration
TAN_amount_M = np.zeros(mtrx)
TAN_amount_M_slat = np.zeros(mtrx)
TAN_amount_M_pit = np.zeros(mtrx)
## total water pool of the system (manure water; urine+manure water+[washing water])
Total_water_pool = np.zeros(mtrx)
Total_water_pool_slat = np.zeros(mtrx)
Total_water_pool_pit = np.zeros(mtrx)
## ratio of [NH4+]/[H+] in the system
Gamma_manure = np.zeros(mtrx)
Gamma_manure_slat = np.zeros(mtrx)
Gamma_manure_pit = np.zeros(mtrx)
## surface NH3 concentrtion at equilirium (in molar mass)
NH3_gas_M = np.zeros(mtrx)
NH3_gas_M_slat = np.zeros(mtrx)
NH3_gas_M_pit = np.zeros(mtrx)
## surface NH3 concentrtion at equilirium (in ug)
NH3_gas_ug = np.zeros(mtrx)
NH3_gas_ug_slat = np.zeros(mtrx)
NH3_gas_ug_pit = np.zeros(mtrx)
## emission potential
modelled_emiss = np.zeros(mtrx)
modelled_emiss_slat = np.zeros(mtrx)
modelled_emiss_pit = np.zeros(mtrx)
## final emission
## the difference between [modelled_emiss] (so-called "emission potential") and [NH3_flux] 
## (so-called "final emission") is that whether there is other factors that limits the amount of flux to 
## the atmosphere, i.e. measures mitigating emissions during housing stage such as coverings, straw ...;
## during land spreading stages, such as canopy recapture, deap injection...
NH3_flux = np.zeros(mtrx)
NH3_flux_slat = np.zeros(mtrx)
NH3_flux_pit = np.zeros(mtrx)

## pools left after housing; transfer to storage
manure_pool_pit_to_storage = np.zeros(mtrx)
avail_N_pool_pit_to_storage = np.zeros(mtrx)
resist_N_pool_pit_to_storage = np.zeros(mtrx)
unavail_N_pool_pit_to_storage = np.zeros(mtrx)
TAN_pool_pit_to_storage = np.zeros(mtrx)
Total_water_pool_pit_to_storage = np.zeros(mtrx)

manure_pool_to_storage = np.zeros(mtrx)
urea_pool_to_storage = np.zeros(mtrx)
avail_N_pool_to_storage = np.zeros(mtrx)
resist_N_pool_to_storage = np.zeros(mtrx)
unavail_N_pool_to_storage = np.zeros(mtrx)
TAN_pool_to_storage = np.zeros(mtrx)
Total_water_pool_to_storage = np.zeros(mtrx)
NH3_flux_from_barn = np.zeros(mtrx)

## for housing simulation
NH3_inconc = np.zeros(mtrx)
NH3_out = np.zeros(mtrx)


###################################
## define model functions
###################################

## simulation: for houses with slatted floor (two-source model: slat+pit); industrial production system
def slat_pit_housing_sim(start_day_idx,end_day_idx,f_slat,f_gap):
    for dd in np.arange(start_day_idx,end_day_idx-1):
        
        ## daily manure and urine on per unit area
        manure[dd+1] = dmanure
        urine[dd+1] = durine
        
        ## manure pool
        ## the amount of manure left on slat depends on the solid floor ratio, here f_slat = 40%
        ## we assume manure left on slat will drop to the pit on the next day
        manure_pool_slat[dd+1] = manure[dd+1]*f_slat
        manure_pool_pit[dd+1] = manure_pool_pit[dd] + manure[dd+1]*f_gap + manure_pool_slat[dd]

        ## N input in multiple forms
        urea[dd+1] = durea
        avail_N[dd+1] = (dmanure_N + (durine_N-durea))*f_avail
        resist_N[dd+1] = (dmanure_N + (durine_N-durea))*f_resist
        unavail_N[dd+1] = (dmanure_N + (durine_N-durea))*f_unavail

        ## TAN production from urea hydrolysis from slat and pit, respectively
        ## the conditions (T, RH) of slat and pit are assume to be consistent, so is the urea hydrolysis rate
        ## and the N decomposition rate from dung
        TAN_prod_pit[dd+1] = daily_urea_hydro_rate[dd+1]*urea_pool_pit[dd]+\
                                daily_Na_decomp_rate[dd+1]*avail_N_pool_pit[dd] +\
                                daily_Nr_decomp_rate[dd+1]*resist_N_pool_pit[dd]

        ## Urea pool
        ## slat urea pool is reset at a daily basis
        urea_pool_slat[dd+1] = urea[dd+1]*f_slat
        ## urea left on slat will go to the pit urea pool
        urea_pool_pit[dd+1] = urea_pool_pit[dd]* (1 - daily_urea_hydro_rate[dd+1]) + urea[dd+1]*f_gap + urea_pool_slat[dd]*(1-daily_urea_hydro_rate[dd])

        ## Org N pools in various forms
        avail_N_pool_slat[dd+1] = avail_N[dd+1]*f_slat
        avail_N_pool_pit[dd+1] = avail_N_pool_pit[dd]*(1 - daily_Na_decomp_rate[dd+1])+ \
                                    avail_N[dd+1]*f_gap + \
                                    avail_N_pool_slat[dd]*(1-daily_Na_decomp_rate[dd])

        resist_N_pool_slat[dd+1] = resist_N[dd+1]*f_slat
        resist_N_pool_pit[dd+1] = resist_N_pool_pit[dd]* (1 - daily_Nr_decomp_rate[dd+1])+ \
                                    resist_N[dd+1]*f_gap+ \
                                    resist_N_pool_slat[dd]*(1 - daily_Nr_decomp_rate[dd])

        unavail_N_pool_slat[dd+1] = unavail_N[dd+1]*f_slat
        unavail_N_pool_pit[dd+1] = unavail_N_pool_pit[dd] + unavail_N[dd+1] *f_slat + unavail_N_pool_slat[dd]

        ## TAN production from urea hydrolysis and org N decomposition
        TAN_prod_slat[dd+1] = daily_urea_hydro_rate[dd+1]*urea_pool_slat[dd+1]+ \
                           daily_Na_decomp_rate[dd+1]*avail_N_pool_slat[dd+1] + \
                           daily_Nr_decomp_rate[dd+1]*resist_N_pool_slat[dd+1]

        ## water from the fresh dung
        manure_initwc_slat[dd+1] = manure_wc*f_slat
        manure_initwc_pit[dd+1] = manure_wc*f_gap

        ## water amount when mositure content reach equilibrium
        manure_minwc_slat[dd+1] = manure_pool_slat[dd+1]*mois_coeff[dd+1] / 100
        manure_minwc_pit[dd+1] = manure_pool_pit[dd+1]*mois_coeff[dd+1]/100

        ## water pool of slat and pit
        ## slat water pool
        water_slat_idx = urine[dd+1]*f_slat + manure_initwc_slat[dd+1] - evap_slat[dd+1] - manure_minwc_slat[dd+1]
        Total_water_pool_slat[dd+1][water_slat_idx>=0] = urine[dd+1][water_slat_idx>=0]*f_slat +\
                                                        manure_initwc_slat[dd+1][water_slat_idx>=0] -\
                                                        evap_slat[dd+1][water_slat_idx>=0]
        Total_water_pool_slat[dd+1][water_slat_idx<0] = manure_minwc_slat[dd+1][water_slat_idx<0]
        ## pit water pool
        water_pit_idx = Total_water_pool_pit[dd]-evap_pit[dd]-manure_minwc_pit[dd+1]
        Total_water_pool_pit[dd+1][water_pit_idx>=0] = Total_water_pool_pit[dd][water_pit_idx>0] +\
                                                        urine[dd+1][water_pit_idx>0]*f_gap +\
                                                        manure_initwc_pit[dd+1][water_pit_idx>0] -\
                                                        evap_pit[dd][water_pit_idx>0] +\
                                                        Total_water_pool_slat[dd][water_pit_idx>0]
        Total_water_pool_pit[dd+1][water_pit_idx<0] = manure_minwc_pit[dd+1][water_pit_idx<=0]+\
                                                        urine[dd+1][water_pit_idx<=0]*f_gap +\
                                                        manure_initwc_pit[dd+1][water_pit_idx<=0] +\
                                                        Total_water_pool_slat[dd][water_pit_idx<=0]

        ## TAN pool of slat and pit
        TAN_pool_slat[dd+1] = TAN_prod_slat[dd+1]
        ## TAN left on slat will go to the pit TAN pool
        TAN_pit_idx = TAN_pool_pit[dd] - NH3_flux_pit[dd]
        TAN_pool_pit[dd+1][TAN_pit_idx>0] = TAN_pit_idx[TAN_pit_idx>0]+TAN_prod_pit[dd+1][TAN_pit_idx>0]+\
                                            TAN_pool_slat[dd][TAN_pit_idx>0]-NH3_flux_slat[dd][TAN_pit_idx>0]
        TAN_pool_pit[dd+1][TAN_pit_idx<=0] = TAN_prod_pit[dd+1][TAN_pit_idx<=0]+TAN_pool_slat[dd][TAN_pit_idx<=0]-\
                                            NH3_flux_slat[dd][TAN_pit_idx<=0]

        ## TAN pool in ug
        TAN_pool_ug_slat[dd+1] = TAN_pool_slat[dd+1] * 1e6
        TAN_pool_ug_pit[dd+1] = TAN_pool_pit[dd+1] * 1e6

        ## TAN conc on slat and in pit
        TAN_amount_slat[dd+1][Total_water_pool_slat[dd+1]==0] = 0
        TAN_amount_slat[dd+1][Total_water_pool_slat[dd+1]!=0] = TAN_pool_slat[dd+1][Total_water_pool_slat[dd+1]!=0]/\
                                                    Total_water_pool_slat[dd+1][Total_water_pool_slat[dd+1]!=0]
        TAN_amount_pit[dd+1][Total_water_pool_pit[dd+1]==0] = 0
        TAN_amount_pit[dd+1][Total_water_pool_pit[dd+1]!=0] = TAN_pool_pit[dd+1][Total_water_pool_pit[dd+1]!=0]/\
                                                    Total_water_pool_pit[dd+1][Total_water_pool_pit[dd+1]!=0]

        ## TAN molar conc
        TAN_amount_M_slat[dd+1] = TAN_amount_slat[dd+1]/14*1000
        TAN_amount_M_pit[dd+1] = TAN_amount_pit[dd+1]/14*1000

        ## Gamma value
        Gamma_manure_slat[dd+1] =  TAN_amount_M_slat[dd+1]/(cc_H + k_NH4[dd+1])
        Gamma_manure_pit[dd+1] =  TAN_amount_M_pit[dd+1]/(cc_H + k_NH4[dd+1])

        ## Gaseous NH3 at the surface of slat TAN pool and pit TAN pool
        NH3_gas_M_slat[dd+1] = Henry_constant[dd+1]*Gamma_manure_slat[dd+1]
        NH3_gas_M_pit[dd+1] = Henry_constant[dd+1]*Gamma_manure_pit[dd+1]

        ## in ug
        NH3_gas_ug_slat[dd+1] = NH3_gas_M_slat[dd+1]*14*1e9
        NH3_gas_ug_pit[dd+1] = NH3_gas_M_pit[dd+1]*14*1e9

        ## determining the maximum emission; emission cannot exceed TAN pool
        emiss_slat_idx = (NH3_gas_ug_slat[dd+1]*3600*24/R_star_slat[dd+1]) - TAN_pool_ug_slat[dd+1]
        modelled_emiss_slat[dd+1][emiss_slat_idx>=0] = TAN_pool_ug_slat[dd+1][emiss_slat_idx>=0]
        modelled_emiss_slat[dd+1][emiss_slat_idx<0] = NH3_gas_ug_slat[dd+1][emiss_slat_idx<0]*3600*24/\
                                                    R_star_slat[dd+1][emiss_slat_idx<0]

        emiss_pit_idx = (NH3_gas_ug_pit[dd+1]*3600*24/R_star_pit[dd+1]) - TAN_pool_ug_pit[dd+1]
        modelled_emiss_pit[dd+1][emiss_pit_idx>=0] = TAN_pool_ug_pit[dd+1][emiss_pit_idx>=0]
        modelled_emiss_pit[dd+1][emiss_pit_idx<0] = NH3_gas_ug_pit[dd+1][emiss_pit_idx<0]*3600*24/\
                                                    R_star_pit[dd+1][emiss_pit_idx<0]

        ## final emission flux from slat and pit; flux are additive
        NH3_flux_slat[dd+1] = modelled_emiss_slat[dd+1]/1e6
        NH3_flux_pit[dd+1] = modelled_emiss_pit[dd+1]/1e6
    return

## simulation: for houses with concrete floor (singe-source model); 
## intermediate/backyard production system (pig/ruminants)
def barn_sim(start_day_idx,end_day_idx):
    for dd in np.arange(start_day_idx,end_day_idx-1):
        ## daily manure and urine on per unit area
        manure[dd+1] = dmanure
        urine[dd+1] = durine
        
        ## manure pool
        manure_pool[dd+1] = manure_pool[dd] + manure[dd+1]
        
        ## N input in multiple forms
        urea[dd+1] = durea
        avail_N[dd+1] = (dmanure_N + (durine_N-durea))*f_avail
        resist_N[dd+1] = (dmanure_N + (durine_N-durea))*f_resist
        unavail_N[dd+1] = (dmanure_N + (durine_N-durea))*f_unavail
        
         ## Urea pool
        urea_pool[dd+1] = urea_pool[dd]*(1 - daily_urea_hydro_rate[dd]) + urea[dd+1]

        ## Org N pools in various forms
        avail_N_pool[dd+1] = avail_N_pool[dd]*(1 - daily_Na_decomp_rate[dd]) + avail_N[dd+1]
        resist_N_pool[dd+1] = resist_N_pool[dd]* (1 - daily_Nr_decomp_rate[dd]) + resist_N[dd+1]            
        unavail_N_pool[dd+1] = unavail_N_pool[dd] + unavail_N[dd+1] 
        
        ## TAN production from urea hydrolysis and the N decomposition rate from dung
        TAN_prod[dd+1] = daily_urea_hydro_rate[dd+1]*urea_pool[dd+1]+\
                                daily_Na_decomp_rate[dd+1]*avail_N_pool[dd+1] +\
                                daily_Nr_decomp_rate[dd+1]*resist_N_pool[dd+1]
        
        ## water from the fresh dung
        manure_initwc[dd+1] = manure_wc
        
        ## water amount when mositure content reach equilibrium
        manure_minwc[dd+1] = manure_pool[dd+1]*mois_coeff[dd+1]/100
        
        ## water pool 
        water_idx = Total_water_pool[dd] + urine[dd+1] + manure_initwc[dd+1] - evap[dd+1]-manure_minwc[dd+1] 
        Total_water_pool[dd+1][water_idx>0] = Total_water_pool[dd][water_idx>0] +\
                                                        urine[dd+1][water_idx>0] +\
                                                        manure_initwc[dd+1][water_idx>0] -\
                                                        evap[dd+1][water_idx>0] 
        Total_water_pool[dd+1][water_idx<=0] = manure_minwc[dd+1][water_idx<=0] 

        ## TAN pool 
        TAN_pool[dd+1] = TAN_pool[dd]-NH3_flux[dd]+TAN_prod[dd+1]


        ## TAN pool in ug
        TAN_pool_ug[dd+1] = TAN_pool[dd+1] * 1e6

        ## TAN conc 
        TAN_amount[dd+1][Total_water_pool[dd+1]==0] = 0
        TAN_amount[dd+1][Total_water_pool[dd+1]!=0] = TAN_pool[dd+1][Total_water_pool[dd+1]!=0]/\
                                                    Total_water_pool[dd+1][Total_water_pool[dd+1]!=0]

        ## TAN molar conc
        TAN_amount_M[dd+1] = TAN_amount[dd+1]/14*1000
        
        ## Gamma value 
        Gamma_manure[dd+1] =  TAN_amount_M[dd+1]/(cc_H + k_NH4[dd+1])

        ## Gaseous NH3 at the surface 
        NH3_gas_M[dd+1] = Henry_constant[dd+1]*Gamma_manure[dd+1]
        
        ## in ug
        NH3_gas_ug[dd+1] = NH3_gas_M[dd+1]*14*1e9
        
        ## determining the maximum emission; emission cannot exceed TAN pool
        emiss_idx = (NH3_gas_ug[dd+1]*3600*24/R_star[dd+1]) - TAN_pool_ug[dd+1]
        modelled_emiss[dd+1][emiss_idx>=0] = TAN_pool_ug[dd+1][emiss_idx>=0]
        modelled_emiss[dd+1][emiss_idx<0] = NH3_gas_ug[dd+1][emiss_idx<0]*3600*24/\
                                                    R_star[dd+1][emiss_idx<0]

        ## final emission flux 
        NH3_flux[dd+1] = modelled_emiss[dd+1]/1e6
    return

## simulation: for poultry housing;
## broilers/layers; housing resistance is assumed to be 16700 s/m
def poultry_housing_sim(start_day_idx,end_day_idx):
    
    for dd in np.arange(start_day_idx,end_day_idx-1):
        
        ## excretion pool; Eq.1
        manure_pool[dd+1] = manure_pool[dd] + manure

        TAN_prod[dd+1] = daily_ua_conv_factor[dd+1] * UA_pool[dd]
        
        ## UA pool; Eq.2
        UA_idx = UA_pool[dd] - TAN_prod[dd+1]
        UA_pool[dd+1][UA_idx>0] = UA_pool[dd][UA_idx>0] + N_birds_UAN[UA_idx>0] - \
                                                TAN_prod[dd+1][UA_idx>0]
        UA_pool[dd+1][UA_idx<=0] = N_birds_UAN[UA_idx<=0]


        ## manure_water calculation; Eq.9  
        manure_water[dd+1][mois_coeff[dd+1]>=5] = manure_pool[dd+1][mois_coeff[dd+1]>=5] * mois_coeff[dd+1][mois_coeff[dd+1]>=5] / 100
        manure_water[dd+1][mois_coeff[dd+1]<5] = 0
        Total_water_pool[dd+1] = manur_water[dd+1]
        
        ## TAN pool; Eq.3
        TAN_idx = TAN_pool[dd] - NH3_Flux[dd]
        TAN_pool[dd+1][TAN_idx<=0] = TAN_prod[dd+1][TAN_idx<=0]
        TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0] + TAN_prod[dd+1][TAN_idx>0]
        TAN_pool_ug[dd+1] = TAN_pool[dd+1] * 1e6
        
        ## TAN concentration
        TAN_amount[dd+1][Total_water_pool[dd+1]==0] = 0
        TAN_amount[dd+1][Total_water_pool[dd+1]!=0] = TAN_pool[dd+1][Total_water_pool[dd+1]!=0] \
                                                         / Total_water_pool[dd+1][Total_water_pool[dd+1]!=0]
        TAN_amount_M[dd+1] = TAN_amount[dd+1] / 14 * 1000

        ## Gamma value; Eq.7
        Gamma_manure[dd+1] =  TAN_amount_M[dd+1] / (cc_H+k_NH4[dd+1])
        
        ## Gaseous concentration of NH3 at the surface; Eq.6
        NH3_gas_M[dd+1] =Henry_constant[dd+1] * Gamma_manure[dd+1]
        NH3_gas_ug[dd+1] = NH3_gas_M[dd+1] * 14 * 1e9

        ## Emissions; see Eq.8 in general, and Sect.2.2.2; Eq.14~17
        emiss_idx = (NH3_gas_ug[dd+1]*3600*24/R_ab) - TAN_pool_ug[dd+1]
        modelled_emiss[dd+1][emiss_idx>0] = TAN_pool_ug[dd+1][emiss_idx>0]
        modelled_emiss[dd+1][emiss_idx<0] = NH3_gas_ug[dd+1][emiss_idx<0]*3600*24/R_star_p

        NH3_Flux[dd+1] = Modelled_emiss[dd+1]/1e6 
    return

## initialisation: model initialisation for housing simulation
def housing_init():
    manure[:] = 0.0
    manure_pool[:] = 0.0
    manure_pool_slat[:] = 0.0
    manure_pool_pit[:] = 0.0
    manure_water[:] = 0.0
    manure_minwc_slat[:] = 0.0
    manure_minwc_pit[:] = 0.0
    manure_initwc_slat[:] = 0.0
    manure_initwc_pit[:] = 0.0
    UA_pool[:] = 0.0
    urine[:] = 0.0
    urea[:] = 0.0
    urea_pool[:] = 0.0
    urea_pool_slat[:] = 0.0
    urea_pool_pit[:] = 0.0
    avail_N[:] = 0.0
    avail_N_pool[:] = 0.0
    avail_N_pool_slat[:] = 0.0
    avail_N_pool_pit[:] = 0.0
    resist_N[:] = 0.0
    resist_N_pool[:] = 0.0
    resist_N_pool_slat[:] = 0.0
    resist_N_pool_pit[:] = 0.0
    unavail_N[:] = 0.0
    unavail_N_pool[:] = 0.0
    unavail_N_pool_slat[:] = 0.0
    unavail_N_pool_pit[:] = 0.0
    TAN_prod[:] = 0.0
    TAN_prod_slat[:] = 0.0
    TAN_prod_pit[:] = 0.0
    TAN_pool[:] = 0.0
    TAN_pool_slat[:] = 0.0
    TAN_pool_pit[:] = 0.0
    TAN_pool_ug[:] = 0.0
    TAN_pool_ug_slat[:] = 0.0
    TAN_pool_ug_pit[:] = 0.0
    TAN_amount[:] = 0.0
    TAN_amount_slat[:] = 0.0
    TAN_amount_pit[:] = 0.0
    TAN_amount_M[:] = 0.0
    TAN_amount_M_slat[:] = 0.0
    TAN_amount_M_pit[:] = 0.0
    Total_water_pool[:] = 0.0
    Total_water_pool_slat[:] = 0.0
    Total_water_pool_pit[:] = 0.0
    Gamma_manure[:] = 0.0
    Gamma_manure_slat[:] = 0.0
    Gamma_manure_pit[:] = 0.0
    NH3_gas_M[:] = 0.0
    NH3_gas_M_slat[:] = 0.0
    NH3_gas_M_pit[:] = 0.0
    NH3_gas_ug[:] = 0.0
    NH3_gas_ug_slat[:] = 0.0
    NH3_gas_ug_pit[:] = 0.0
    modelled_emiss[:] = 0.0
    modelled_emiss_slat[:] = 0.0
    modelled_emiss_pit[:] = 0.0
    NH3_flux[:] = 0.0
    NH3_flux_slat[:] = 0.0
    NH3_flux_pit[:] = 0.0
    return

## initialisation: 2nd initialisation for housing simulation
def housing_2nd_init():
    aa1 = manure_pool[-1]
    aa2 = manure_pool_slat[-1]
    aa3 = manure_pool_pit[-1]
    bb1 = urea_pool[-1]
    bb2 = urea_pool_slat[-1]
    bb3 = urea_pool_pit[-1]
    cc1 = TAN_pool[-1]
    cc2 = TAN_pool_slat[-1]
    cc3 = TAN_pool_pit[-1]
    dd1 = NH3_flux[-1]
    dd2 = NH3_flux_slat[-1]
    dd3 = NH3_flux_pit[-1]
    ee1 = Total_water_pool[-1]
    ee2 = Total_water_pool_slat[-1]
    ee3 = Total_water_pool_pit[-1]
    ff1 = UA_pool[-1]

    manure_pool[0] = aa1
    urea_pool[0] = bb1
    TAN_pool[0] = cc1
    NH3_flux[0] = dd1
    Total_water_pool[0] = ee1
    UA_pool[0] = ff1

    manure_pool_slat[0] = aa2
    manure_pool_pit[0] = aa3
    urea_pool_slat[0] = bb2
    urea_pool_pit[0] = bb3
    TAN_pool_slat[0] = cc2
    TAN_pool_pit[0] = cc3
    NH3_flux_slat[0] = dd2
    NH3_flux_pit[0] = dd3
    Total_water_pool_slat[0] = ee2
    Total_water_pool_pit[0] = ee3   
    return

## initialisation for house (barn) cleaning; all variables at the cleaning day are reset to zero
def cleaning_pit(day_idx):
    ## pools: from housing to MMS
    manure_pool_pit_to_storage[day_idx] = manure_pool_pit[day_idx]+manure_pool_slat[day_idx]
    avail_N_pool_pit_to_storage[day_idx] = avail_N_pool_pit[day_idx]
    resist_N_pool_pit_to_storage[day_idx] = resist_N_pool_pit[day_idx]
    unavail_N_pool_pit_to_storage[day_idx] = unavail_N_pool_pit[day_idx]
    TAN_pool_pit_to_storage[day_idx] = TAN_pool_pit[day_idx]+urea_pool_pit[day_idx]
    Total_water_pool_pit_to_storage[day_idx] = Total_water_pool_pit[day_idx]
    
    # manure[day_idx] = 0.0
    manure_pool[day_idx] = 0.0
    manure_pool_slat[day_idx] = 0.0
    manure_pool_pit[day_idx] = 0.0
    # manure_water[day_idx] = 0.0
    # manure_minwc_slat[day_idx] = 0.0
    # manure_minwc_pit[day_idx] = 0.0
    manure_initwc_slat[day_idx] = 0.0
    manure_initwc_pit[day_idx] = 0.0
    UA_pool[day_idx] = 0.0
    # urine[day_idx] = 0.0
    # urea[day_idx] = 0.0
    urea_pool[day_idx] = 0.0
    urea_pool_slat[day_idx] = 0.0
    urea_pool_pit[day_idx] = 0.0
    # avail_N[day_idx] = 0.0
    avail_N_pool[day_idx] = 0.0
    avail_N_pool_slat[day_idx] = 0.0
    avail_N_pool_pit[day_idx] = 0.0
    # resist_N[day_idx] = 0.0
    resist_N_pool[day_idx] = 0.0
    resist_N_pool_slat[day_idx] = 0.0
    resist_N_pool_pit[day_idx] = 0.0
    # unavail_N[day_idx] = 0.0
    unavail_N_pool[day_idx] = 0.0
    unavail_N_pool_slat[day_idx] = 0.0
    unavail_N_pool_pit[day_idx] = 0.0
    # TAN_prod[day_idx] = 0.0
    # TAN_prod_slat[day_idx] = 0.0
    # TAN_prod_pit[day_idx] = 0.0
    TAN_pool[day_idx] = 0.0
    TAN_pool_slat[day_idx] = 0.0
    TAN_pool_pit[day_idx] = 0.0
    TAN_pool_ug[day_idx] = 0.0
    TAN_pool_ug_slat[day_idx] = 0.0
    TAN_pool_ug_pit[day_idx] = 0.0
    # TAN_amount[day_idx] = 0.0
    # TAN_amount_slat[day_idx] = 0.0
    # TAN_amount_pit[day_idx] = 0.0
    # TAN_amount_M[day_idx] = 0.0
    # TAN_amount_M_slat[day_idx] = 0.0
    # TAN_amount_M_pit[day_idx] = 0.0
    Total_water_pool[day_idx] = 0.0
    Total_water_pool_slat[day_idx] = 0.0
    Total_water_pool_pit[day_idx] = 0.0
    # Gamma_manure[day_idx] = 0.0
    # Gamma_manure_slat[day_idx] = 0.0
    # Gamma_manure_pit[day_idx] = 0.0
    # NH3_gas_M[day_idx] = 0.0
    # NH3_gas_M_slat[day_idx] = 0.0
    # NH3_gas_M_pit[day_idx] = 0.0
    # NH3_gas_ug[day_idx] = 0.0
    # NH3_gas_ug_slat[day_idx] = 0.0
    # NH3_gas_ug_pit[day_idx] = 0.0
    modelled_emiss[day_idx] = 0.0
    modelled_emiss_slat[day_idx] = 0.0
    modelled_emiss_pit[day_idx] = 0.0
    NH3_flux[day_idx] = 0.0
    NH3_flux_slat[day_idx] = 0.0
    NH3_flux_pit[day_idx] = 0.0
    return

def cleaning_barn(day_idx):
    ## pools: from housing to MMS
    manure_pool_to_storage[day_idx] = manure_pool[day_idx]
    urea_pool_to_storage[day_idx] = (1 - daily_urea_hydro_rate[day_idx])*urea_pool[day_idx]
    avail_N_pool_to_storage[day_idx] = (1 - daily_Na_decomp_rate[day_idx])*avail_N_pool[day_idx]
    resist_N_pool_to_storage[day_idx] = (1 - daily_Nr_decomp_rate[day_idx])*resist_N_pool[day_idx]
    unavail_N_pool_to_storage[day_idx] = unavail_N_pool[day_idx]
    TAN_pool_to_storage[day_idx] = TAN_pool[[day_idx]]-NH3_flux[day_idx]
    Total_water_pool_to_storage[day_idx] = Total_water_pool[day_idx]
    NH3_flux_from_barn[day_idx] = NH3_flux[day_idx]

    manure_pool[day_idx] = 0.0
    manure_water[day_idx] = 0.0
    UA_pool[day_idx] = 0.0
    urea_pool[day_idx] = 0.0
    avail_N_pool[day_idx] = 0.0
    resist_N_pool[day_idx] = 0.0
    unavail_N_pool[day_idx] = 0.0
    TAN_pool[day_idx] = 0.0
    Total_water_pool[day_idx] = 0.0
    modelled_emiss[day_idx] = 0.0
    NH3_flux[day_idx] = 0.0
    return

def housing_to_storage_init():
    manure_pool_pit_to_storage[:] = 0.0
    avail_N_pool_pit_to_storage[:] = 0.0
    resist_N_pool_pit_to_storage[:] = 0.0
    unavail_N_pool_pit_to_storage[:] = 0.0
    TAN_pool_pit_to_storage[:] = 0.0
    Total_water_pool_pit_to_storage[:] = 0.0
    
    manure_pool_to_storage[:] = 0.0
    urea_pool_to_storage[:] = 0.0
    avail_N_pool_to_storage[:] = 0.0
    resist_N_pool_to_storage[:] = 0.0
    unavail_N_pool_to_storage[:] = 0.0
    TAN_pool_to_storage[:] = 0.0
    Total_water_pool_to_storage[:] = 0.0
    NH3_flux_from_barn[:] = 0.0
    return
