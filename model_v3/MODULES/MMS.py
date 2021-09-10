
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
#sys.exit()

##################################
## get MMS data
##################################
## MMS_file has been read in INPUT.input.py
MMS_type_list = []
for MMS in MMS_file.data_vars:
    MMS_type_list.append(str(MMS))

## provisional MMS categories
## loss (untraceable): fishpond, discahrge, public sewage
loss_list = ['fishpond','discharge','publsewage','dumping']
## sold (untraceable): sold
sold_list = ['sold']
## Cat D: no significant emission; used as fuel: biogas(digester,liquid), burned (solid)
MMS_fuel_list = ['mmsbiogas','mmsburned']
## Cat C: N mostly preserved for further use (as solid): thermal drying 
MMS_preserve_solid_list = ['mmsthermal']
## Cat C: N mostly preserved for further use (as liquid), e.g., manure stored in liquid phase with emission mitigation measures: lagoon (typically with cover), liquid crust
MMS_preserve_liquid_list = ['mmslagoon','mmsliqcrust']
## Cat A: manure stored in barns (as solid): composting, deep litter, litter (poultry), no litter (poultry), pit (layer),solid storage
MMS_barn_solid_list = ['mmscompost','mmsdeeplitt','mmslitter','mmsnolitt','mmssolid']
##  Cat A: manure stored in barns (as liquid): aerobic processing, liquid, pit1, pit2
MMS_barn_liquid_list = ['mmsaerproc','mmsliquid','mmsliqoth','mmspit1','mmspit2']
## Cat B: manure in open environment; left on land (as solid): aerobic processing, daily spreading, dry lot, pasture, pasture+paddock
MMS_open_solid_list = ['mmsconfin','mmsdaily','mmsdrylot','mmspasture','mmspastpad']
## Cat B: manure in open environment (as liquid): aerobic lagoon, liquid
MMS_open_liquid_list = ['mmsaerobic']

## shape in [lat,lon]
f_loss = np.zeros(mtrx[1:])
f_sold = np.zeros(mtrx[1:])
f_MMS_fuel = np.zeros(mtrx[1:])
f_MMS_preserve_solid = np.zeros(mtrx[1:])
f_MMS_preserve_liquid = np.zeros(mtrx[1:])
f_MMS_barn_solid = np.zeros(mtrx[1:])
f_MMS_barn_liquid = np.zeros(mtrx[1:])
f_MMS_open_solid = np.zeros(mtrx[1:])
f_MMS_open_liquid = np.zeros(mtrx[1:]) 

for mms in loss_list:
    try:f_loss = f_loss + MMS_file[mms].values   
    except:pass
for mms in sold_list:
    try:f_sold = f_sold + MMS_file[mms].values
    except:pass
for mms in MMS_fuel_list:
    try:f_MMS_fuel = f_MMS_fuel + MMS_file[mms].values
    except:pass
for mms in MMS_preserve_solid_list:
    try:f_MMS_preserve_solid = f_MMS_preserve_solid + MMS_file[mms].values
    except:pass
for mms in MMS_preserve_liquid_list:
    try:f_MMS_preserve_liquid = f_MMS_preserve_liquid + MMS_file[mms].values
    except:pass
for mms in MMS_barn_solid_list:
    try:f_MMS_barn_solid = f_MMS_barn_solid + MMS_file[mms].values
    except:pass
for mms in MMS_barn_liquid_list:
    try:f_MMS_barn_liquid = f_MMS_barn_liquid + MMS_file[mms].values
    except:pass
for mms in MMS_open_solid_list:
    try:f_MMS_open_solid = f_MMS_open_solid + MMS_file[mms].values
    except:pass
for mms in MMS_open_liquid_list:
    try:f_MMS_open_liquid = f_MMS_open_liquid + MMS_file[mms].values
    except:pass

## this refers to the storage surface area in farms relative to the housing area
## e.g., mms_barn_solid = 0.1, 
##       which means the area of a barn that stores solid manure has 1/10 of housing area in the local farm;
##       housing area = 10 km^2 in a grid; barn_solid area = 1 km^2 (1/10) * f_MMS_barn_solid
## these values need to be reviewed (?)
MMS_area_factor = {
    "mms_barn_solid":0.1,
    "mms_barn_liquid":1.0,
    "mms_open_solid":0.1,
    "mms_open_liquid":0.5}

## areas of each MMS, initial values are None
MMS_area = {
    "mms_barn_solid_area":None,
    "mms_barn_liquid_area":None,
    "mms_open_solid_area":None,
    "mms_open_liquid_area":None}

###################################
## MMS parameters
###################################
## assuming the roughness height of manure storage barn is ~ 0.5m (<ref height of 2m)
zo_barn = 0.5  
## assuming the roughness height of manure pile (open land) is ~ 1.0m (<ref height of 2m)
zo_manureland = 1.0
## assuming the dry matter (DM) content of solid manure is 20%
DM_content = solid_m_DM[livestock]
## assuming the density of manure; 1t kg/m^3 or 1g/cm^3
manure_density = pho_m[livestock]
## assuming source layer of the top soil is 2 cm (0.02 m) thick
z_soil = 0.02
## assuming infiltration of manure water to the soil is 3mm/day (3000 g/m^2/day) ref: Sommer&Jacobsen,1999 and Vira et al.,2020 GMD
ki = 3000.0
## assuming soil characteristics: 1) sand (%), 2) clay (%), 3) bulk density (g/cm^3), 4) particle density (g/cm^3)
soil_sand = 80
soil_clay = 8
soil_bd = 1.5
soil_pd = 2.66
##################################
## define 
##################################
class MMS_module:
    def __init__(self,array_shape,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,water_added,pH_value,area_housing):
        ## feces input from housing
        self.manure_added = manure_added
        self.manure = np.zeros(array_shape)
        ## feces pool
        self.manure_pool = np.zeros(array_shape)
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
        self.urea_added = urea_added
        self.urea = np.zeros(array_shape)
        ## urea pool
        self.urea_pool = np.zeros(array_shape)
        ## uric acid input from housing
        self.UA_added = np.zeros(array_shape)
        ## uric acid pool
        self.UA_pool = UA_added
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
        ## TAN added from housing
        self.TAN_added = TAN_added
        self.TAN = np.zeros(array_shape)
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
        ## TAN conc at the surface; for [MMS (barn,open) solid]
        self.TAN_surf_amount_M = np.zeros(array_shape)
        ## TAN conc at the soil surface/interface between manure and land (soil); for [MMS open solid]
        self.TAN_soil_amount_M = np.zeros(array_shape)
        ## water added from housing
        self.water_added = water_added
        self.water = np.zeros(array_shape)
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
        ## diffusion of aqueous TAN to soil
        self.diffusiveflux = np.zeros(array_shape)
        ## infiltration of aqueous TAN to subsurface (not diffusive)
        self.infilflux = np.zeros(array_shape)

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
        ## barn/buidling resistance
        self.R_star = np.zeros(array_shape)
        ## manure resistance
        self.R_manure = np.zeros(array_shape)
        ## soil resistance
        self.R_soil = np.zeros(array_shape)
        ## mositure equilirium, mositure content of manure
        self.mois_coeff = np.zeros(array_shape)
        ## daily UA hydrolysis rate
        self.daily_ua_conv_factor = np.zeros(array_shape)
        ## daily urea hydrolysis rate
        self.daily_urea_hydro_rate = np.zeros(array_shape)
        ## daily decomposition rate of available and resistant N components
        self.daily_Na_decomp_rate = np.zeros(array_shape)
        self.daily_Nr_decomp_rate = np.zeros(array_shape)
        ## Henry's law constant and dissociation equilibria; Eq.6
        self.Henry_constant = np.zeros(array_shape)
        ## dissociation constant of NH4+; 298.15 K is 25 degC (room temperature)
        self.k_NH4 = np.zeros(array_shape)
        ## molecular diffusivity of NH4+ in the water
        self.D_aq_NH4 = np.zeros(array_shape)
        ## tortuosity for diffusion
        self.tor_soil = np.zeros(array_shape)
        ## infiltration rate; m/s
        self.qinfil = np.zeros(array_shape)
        ## pH and H+ ions concentration
        self.pH = pH_value
        self.cc_H = np.float(10**(-pH_value))
        ## housing area that is used to determine MMS area
        self.housingarea = area_housing

    def sim_env(self,mms_type):
        if mms_type == 'MMS_barn':
            self.T_sim,self.u_sim = barn_env(temp_data,\
                wind_profile(uref=wind_data,height_ref=wind_data_height,height_out=ref_height,zo=zo_barn))
            self.RH_sim = rhum_data
            self.T_sim = xr_to_np(self.T_sim)
            self.RH_sim = xr_to_np(self.RH_sim)
            self.u_sim = xr_to_np(self.u_sim)
            ## daily evaporation; aerodynamic method; (g/m^2)
            ## Q_vent above pit is set to be 0.6 m/s (full efficiency)
            self.evap_sim = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=self.u_sim,zo=zo_barn)*1000
            self.R_star = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap_sim/1000)
        else:
            self.T_sim = temp_data
            self.u_sim = wind_data
            self.RH_sim = rhum_data
            self.evap_sim = evap_data
            self.soilmoist = soilmoist_data
            self.persm = persm_data
            self.R_star = resistance_aero_boundary(temp=self.T_sim,rhum=self.RH_sim,\
                u=wind_profile(uref=self.u_sim,height_ref=wind_data_height,height_out=ref_height,zo=zo_manureland),\
                    H=sshf_data,Z=ref_height,zo=zo_manureland)
            self.T_sim = xr_to_np(self.T_sim)
            self.RH_sim = xr_to_np(self.RH_sim)
            self.u_sim = xr_to_np(self.u_sim)
            self.evap_sim = xr_to_np(self.evap_sim)
            self.tor_soil = soil_tuotorsity(theta_sat=0,theta=0,phase="aqueous")
            ## 
            self.qinfil = infiltration_rate(soilmoist_percent=self.soilmoist,sand_percent=soil_sand,\
            clay_percent=soil_clay,bulk_density=soil_bd,particle_density=soil_pd)

        ## mositure equilirium, mositure content of manure
        self.mois_coeff = (-np.log(1.01-(self.RH_sim/100))/(0.0000534*(self.T_sim+273.15)))**(1/1.41)
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
        self.D_aq_NH4 = diffusivity_NH4(temp=self.T_sim)
        return
    
    ## Simulation: Cat A manure stored in barns (as liquid)
    ## water pool is transfered from housing to MMS barn 
    def MMS_barn_liquid_sim(self,start_day_idx,end_day_idx):
        ## area of MMS_barn_liquid
        MMS_area["mms_barn_liquid_area"] = self.housingarea*(1.0 - f_loss - f_sold)*f_MMS_barn_liquid*MMS_area_factor['mms_barn_liquid']
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx-1):

                ## the equations used to represent these pools need to be explained here:
                ## each pool should be in a unit of, mass/unit area, i.e., g/m^2
                ## and each MMS corresponds to a specific MMS area, so the explicit equation would be:
                ##     self.[pool] = self.[poolmass]*(1.0-f_loss-f_sold)*f_MMS_[MMS type]/MMS_area["MMS type"]
                ## then, note that  
                ##     MMS_area["MMS type"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_[MMS type]*MMS_area_factor["MMS type"]
                ##     self.[pool] = self.[poolmass]/(self.housingarea*MMS_area_factor["MMS_type"])
                ## this is different to the HOUSING module as N excretion has been divided by the housing area before used as input

                ## daily manure and urine on per unit area
                self.manure[dd+1] = self.manure_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])

                ## manure pool
                self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]

                ## N input in multiple forms
                self.urea[dd+1] = self.urea_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])
                self.avail_N[dd+1] = self.avail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])
                self.resist_N[dd+1] = self.resist_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])
                self.unavail_N[dd+1] = self.unavail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                        self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                        self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
                ## TAN from housing to storage
                self.TAN[dd+1] = self.TAN_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1])
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

                ## Org N pools in various forms
                self.avail_N_pool[dd+1] = self.avail_N_pool[dd]*(1 - self.daily_Na_decomp_rate[dd+1]) + self.avail_N[dd+1]
                self.resist_N_pool[dd+1] = self.resist_N_pool[dd]* (1 - self.daily_Nr_decomp_rate[dd+1]) + self.resist_N[dd+1]
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1]

                ## water from the fresh dung
                # self.manure_initwc[dd+1] = * f_MMS_barn

                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

                ## water pool
                ## Note: the water pool in [MMS barn liquid] is "inheritated" from housing water pool
                water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd+1]
                self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] +\
                                                                self.water_added[dd+1][water_idx>0]/(self.housingarea*MMS_area_factor["mms_barn_liquid"]) -\
                                                                self.evap_sim[dd][water_idx>0]
                self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd+1][water_idx<=0] +\
                                                                self.water_added[dd+1][water_idx<=0]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])
                                                                
                ## TAN pool
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

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
                self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*timestep/\
                                                            self.R_star[dd+1][emiss_idx<0]

                ## final emission flux
                self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return

    ## Simulation: Cat B manure stored in barns (as solid) //// under development 31/Aug
    ## manure/animal waste is processed/dried to reduce water content; water pool is determined by:
    ## a) mositure equilibrium (minimum amount of manure moisture content) b) assuming 20 % (?) DM content (Sommer&Hutchings,2001)
    ## c) incorporate a [surface compensation point] scheme in the model:
    ##       gaseouos concentrations (chi): chi_indoor; chi_surface
    ##       aqueous concrntrations: [TAN]_bulk (TAN conc in the manure); [TAN]_surf (TAN conc at the surface, in equilibrium with chi_surf)
    ##       equilibrium between [TAN]_surf and chi_surf: chi_surf = (k_H_D/([H+]+k_NH4+))*[TAN]_surf; 
    ##                   Note: k_H_D = (161500/(temp+273.15))*np.exp(-10380/(temp+273.15)) 
    ##       fluxes (F): F_atm, surface to indoor; F_tosurf: manure to surface;;; 
    ##                   F_atm=(chi_surf-chi_indoor)/R_star; F_tosurf=([TAN]_bulk-[TAN]_surf)/R_manure; F_atm = F_tosurf
    ##    [TAN]_bulk is the prognostic variable and is determined by source and loss based on the mass balance approach
    ##    solve chi_surf: chi_surf = ([TAN]_bulk*R_star+chi_indoor*R_manure)/(R_manure*(k_H_D/([H+]+k_NH4+))+R_star)

    def MMS_barn_solid_sim(self,start_day_idx,end_day_idx):
        MMS_area["mms_barn_solid_area"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_barn_solid*MMS_area_factor["mms_barn_solid"]
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx-1):
                ## note the equations for the pools are similar to "barn liquid sim"; see comments above
                ## daily manure and urine on per unit area
                self.manure[dd+1] = self.manure_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])

                ## manure pool
                self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]

                ## manure resistance
                ## manure resistance is determined by: R = z/(2*D); 
                ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not applicable here (?)
                ## layer thickness of manure: z = manure pool/manure density
                ##        manure density is given in g/cm^3, and needs to be converted to g/m^3 by multiplying by 10^6
                self.R_manure[dd+1] = ((self.manure_pool[dd+1]/(DM_content/100))/(manure_density*1e6))/(2*self.D_aq_NH4[dd+1]) 

                ## N input in multiple forms
                self.urea[dd+1] = self.urea_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])
                self.avail_N[dd+1] = self.avail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])
                self.resist_N[dd+1] = self.resist_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])
                self.unavail_N[dd+1] = self.unavail_N_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])

                ## TAN production from urea hydrolysis and the N decomposition rate from dung
                self.TAN_prod[dd+1] = self.daily_urea_hydro_rate[dd+1]*self.urea_pool[dd]+\
                                        self.daily_Na_decomp_rate[dd+1]*self.avail_N_pool[dd] +\
                                        self.daily_Nr_decomp_rate[dd+1]*self.resist_N_pool[dd]
                ## TAN from housing to storage
                self.TAN[dd+1] = self.TAN_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])

                ## Urea pool
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1])
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

                ## Org N pools in various forms
                self.avail_N_pool[dd+1] = self.avail_N_pool[dd]*(1 - self.daily_Na_decomp_rate[dd+1]) + self.avail_N[dd+1]
                self.resist_N_pool[dd+1] = self.resist_N_pool[dd]* (1 - self.daily_Nr_decomp_rate[dd+1]) + self.resist_N[dd+1]
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1]

                ## water amount in "solid" manure
                ## self.manure refers to the DM mass; therefore, total manure mass = DM mass/DM%
                ## water in the "solid" manure = water% x total manure mass
                self.manure_water[dd+1] = (self.manure[dd+1]/(DM_content/100))*(1-DM_content/100)

                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

                ## water pool
                ## Note the difference between the water pool of [MMS barn solid] and [MMS barn liquid]
                ## water pool of [MMS barn solid] is directly determined by the amount of manure as we assumed a dry matter content
                water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd+1]
                self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] +\
                                                                self.manure_water[dd+1][water_idx>0] -\
                                                                self.evap_sim[dd][water_idx>0]
                self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd+1][water_idx<=0] +\
                                                                self.manure_water[dd+1][water_idx<=0] 
                                                                
                ## TAN pool
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

                ## TAN pool in ug
                self.TAN_pool_ug[dd+1] = self.TAN_pool[dd+1] * 1e6

                ## TAN conc
                ## TAN will partitioned into aqueous and solid (adsorption to manure) phase
                ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd = 1.0 m^3/m^3 
                ## Kd = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                ## [TAN(s)] is concentration of sorbed NH4+ with respect to the volume of manure
                ## manure density varies, ~ 0.3-1.9 g/cm^3, we assume 1.0 g/cm^3 (note the unit!)
                ## [TAN(aq)] = TAN_mass(total)/(manure_mass/manure_density)
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                                                            (self.manure_pool[dd+1][self.Total_water_pool[dd+1]!=0]/manure_density)

                ## TAN molar conc
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/14*1000

                ## TAN conc at the surface
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount_M[dd+1]*self.R_star[dd+1])/\
                                            (self.R_manure[dd+1]*(self.Henry_constant[dd+1]/(self.cc_H + self.k_NH4[dd+1]))+self.R_star[dd+1])

                ## Gamma value
                self.Gamma_manure[dd+1] =  self.TAN_surf_amount_M[dd+1]/(self.cc_H + self.k_NH4[dd+1])

                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure[dd+1]

                ## in ug
                self.NH3_gas_ug[dd+1] = self.NH3_gas_M[dd+1]*14*1e9

                ## determining the maximum emission; emission cannot exceed TAN pool
                emiss_idx = (self.NH3_gas_ug[dd+1]*3600*timestep/self.R_star[dd+1]) - self.TAN_pool_ug[dd+1]
                self.modelled_emiss[dd+1][emiss_idx>=0] = self.TAN_pool_ug[dd+1][emiss_idx>=0]
                self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*timestep/\
                                                            self.R_star[dd+1][emiss_idx<0]

                ## final emission flux
                self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return

    ## Simulation: Cat C manure stored in open environment (as solid) //// under development 07/Sep
    ## manure/animal waste is processed/dried to reduce water content; water pool is determined by:
    ## a) mositure equilibrium (minimum amount of manure moisture content) b) assuming 20 % (?) DM content (Sommer&Hutchings,2001)
    ## c) incorporate a [surface compensation point] scheme in the model as Cat B [MMS barn solid]
    ## d) incorporate a [soil surface compensation point] in the model; with additional TAN loss through aqueous diffusion to soil
    ##    and infiltration (subsurface leaching) to soil
    ##       gaseouos concentrations (chi): chi_atm; chi_surface
    ##       aqueous concrntrations: [TAN]_bulk; [TAN]_surf; [TAN]_soil (TAN conc at the interface between manure and land/soil)
    ##       fluxes (F_upwards): F_atm, surface to atmosphere; F_tosurf: manure to surface;;; 
    ##                   F_atm=(chi_surf-chi_atm)/R_ab; F_tosurf=([TAN]_bulk-[TAN]_surf)/R_manure; F_atm = F_tosurf
    ##       fluxes (F_downwards): F_tosoil: manure to interface layer between manure and soil; F_diffusion: aqueous diffusion to soil
    ##                             F_infiltration: infiltration/subsurface leaching of TAN
    ##                   F_diffusion=[TAN]_soil/R_soil; F_infil=[TAN]_soil*qinfil; F_tosoil=([TAN]_bulk-[TAN]_soil)/R_manure;
    ##                   F_tosoil = F_diffusion + F_infiltration
    ##    ([TAN]_bulk is the prognostic variable and is determined by source and loss based on the mass balance approach)
    ##    (solve chi_surf: chi_surf = ([TAN]_bulk*R_star+chi_atm*R_manure)/(R_manure*(k_H_D/([H+]+k_NH4+))+R_ab; as Cat B [MMS barn solid])
    ##    solve [TAN]_soil: [TAN]_soil = [TAN]_bulk/(R_manure/R_soil+R_manure*qinfil+1)
    def MMS_land_sim(self,start_day_idx,end_day_idx):
        MMS_area[:] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_open_solid*MMS_area_factor["mms_open_solid"]
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx-1):
                ## note the equations for the pools are similar to "barn liquid sim"; see comments above
                ## daily manure and urine on per unit area
                self.manure[dd+1] = self.manure_added[dd+1]/(self.housingarea*MMS_area_factor["mms_open_solid"])

                ## manure pool
                self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1]

                ## manure resistance
                ## manure resistance is determined by: R = z/(2*D); 
                ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not applicable here (?)
                ## layer thickness of manure: z = manure pool/manure density
                ##        manure density is given in g/cm^3, and needs to be converted to g/m^3 by multiplying by 10^6
                self.R_manure[dd+1] = ((self.manure_pool[dd+1]/(DM_content/100))/(manure_density*1e6))/(2*self.D_aq_NH4[dd+1]) 
                
                ## soil resistance
                ## soil resistance = thickness of the source layer (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                self.R_soil[dd+1] = z_soil/(self.tor_soil[dd+1]*self.D_aq_NH4[dd+1])

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
                urea_idx = self.urea_pool[dd]*(1 - self.daily_urea_hydro_rate[dd+1])
                self.urea_pool[dd+1][urea_idx>0] = urea_idx[urea_idx>0] + self.urea[dd+1][urea_idx>0]
                self.urea_pool[dd+1][urea_idx<=0] = self.urea[dd+1][urea_idx<=0]

                ## Org N pools in various forms
                self.avail_N_pool[dd+1] = self.avail_N_pool[dd]*(1 - self.daily_Na_decomp_rate[dd+1]) + self.avail_N[dd+1]
                self.resist_N_pool[dd+1] = self.resist_N_pool[dd]* (1 - self.daily_Nr_decomp_rate[dd+1]) + self.resist_N[dd+1]
                self.unavail_N_pool[dd+1] = self.unavail_N_pool[dd] + self.unavail_N[dd+1]

                ## water amount in "solid" manure
                ## self.manure refers to the DM mass; therefore, total manure mass = DM mass/DM%
                ## water in the "solid" manure = water% x total manure mass
                self.manure_water[dd+1] = (self.manure[dd+1]/(DM_content/100))*(1-DM_content/100)

                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

                ## water pool
                ## water pool of [MMS open solid] is determined by:
                ##   source: manure water, rain (precipitation)
                ##   loss: evaporation, infiltration to soil 
                ## Note: infiltration of manure water to soil is assumed to be 3mm/day (Sommer&Jacobsen,1999; Vira et al., 2020GMD)
                ##       and this should be differentiate with TAN infiltration to soil
                water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-ki-self.manure_minwc[dd+1]
                self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] +\
                                                                self.manure_water[dd+1][water_idx>0] -\
                                                                self.evap_sim[dd][water_idx>0] - ki
                self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd+1][water_idx<=0] +\
                                                                self.manure_water[dd+1][water_idx<=0] 
                                                                
                ## TAN pool (different from [MMS barn (liquid,solid)])
                ## Note: source of TAN pool: TAN production from 1) urea hydrolysis and 2) decomposition of org N
                ##       loss of TAN pool: 1) NH3 volatilization to the atmospere, 2) diffusion to soil, and 3) infiltration (subsurface leaching) to soil
                ##       only aqueous phase diffusion is considered, and gaseous diffusion is not considered in this study
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusiveflux[dd] - self.infilflux[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

                ## TAN pool in ug
                self.TAN_pool_ug[dd+1] = self.TAN_pool[dd+1] * 1e6

                ## TAN conc
                ## TAN will partitioned into aqueous and solid (adsorption to manure) phase
                ## we applied: [TAN(s)] = Kd[TAN(aq)], Kd = 1.0 m^3/m^3 
                ## Kd = 1.0 m^3/m^3 (this is probably for the convenience of calculation...); (Vira et al, 2020 GMD)
                ## [TAN(s)] is concentration of sorbed NH4+ with respect to the volume of manure
                ## manure density varies, ~ 0.3-1.9 g/cm^3, we assume 1.0 g/cm^3 (note the unit!)
                ## [TAN(aq)] = TAN_mass(total)/(manure_mass/manure_density)
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                                                            (self.manure_pool[dd+1][self.Total_water_pool[dd+1]!=0]/manure_density)

                ## TAN molar conc
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/14*1000

                ## TAN conc at the surface
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount_M[dd+1]*self.R_star[dd+1])/\
                                            (self.R_manure[dd+1]*(self.Henry_constant[dd+1]/(self.cc_H + self.k_NH4[dd+1]))+self.R_star[dd+1])

                ## TAN conc at the soil surface/interface between manure and land (soil)

                self.TAN_soil_amount_M[dd+1] = self.TAN_amount_M[dd+1]/\
                    (self.R_manure[dd+1]/self.R_soil[dd+1]+self.R_manure[dd+1]*self.qinfil[dd+1]+1)

                ## TAN loss through aqueous diffusion to soil
                self.diffusiveflux[dd+1] = self.TAN_soil_amount_M[dd+1]/self.R_soil[dd+1]*24*3600

                ## TAN loss through subsurface leaching (infiltration to soil)
                self.infilflux[dd+1] = self.TAN_soil_amount_M[dd+1]*self.qinfil[dd+1]*24*3600

                ## Gamma value
                self.Gamma_manure[dd+1] =  self.TAN_surf_amount_M[dd+1]/(self.cc_H + self.k_NH4[dd+1])

                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.Henry_constant[dd+1]*self.Gamma_manure[dd+1]

                ## in ug
                self.NH3_gas_ug[dd+1] = self.NH3_gas_M[dd+1]*14*1e9

                ## determining the maximum emission; emission cannot exceed TAN pool
                emiss_idx = (self.NH3_gas_ug[dd+1]*3600*timestep/self.R_star[dd+1]) - self.TAN_pool_ug[dd+1]
                self.modelled_emiss[dd+1][emiss_idx>=0] = self.TAN_pool_ug[dd+1][emiss_idx>=0]
                self.modelled_emiss[dd+1][emiss_idx<0] = self.NH3_gas_ug[dd+1][emiss_idx<0]*3600*timestep/\
                                                            self.R_star[dd+1][emiss_idx<0]

                ## final emission flux
                self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return
    
    ## Simulation: Cat C manure stored in open environment (as solid) //// under development 08/Sep
    def MMS_liquid_sim(self,start_day_idx,end_day_idx):
        MMS_area[:] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_open_liquid*MMS_area_factor["mms_open_liquid"]
        # if livestock.lower()=="poultry":
        #     # for dd in np.arange(start_day_idx,end_day_idx-1):
        #     print(livestock)

        # else:
        #     for dd in np.arange(start_day_idx,end_day_idx-1):


        return