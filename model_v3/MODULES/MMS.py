
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
## e.g., mms_barn_solid = 0.2, 
##       which means the area of a barn that stores solid manure has 1/5 of housing area in the local farm;
##       housing area = 10 km^2 in a grid; barn_solid area = 1 km^2 (1/5) * f_MMS_barn_solid
## these values need to be reviewed (?)
MMS_area_factor = {
    "mms_barn_solid":0.2,
    "mms_barn_liquid":0.4,
    "mms_open_solid":1.0,
    "mms_open_liquid":2.5}

## areas of each MMS, initial values are None
MMS_area = {
    "mms_barn_solid_area":None,
    "mms_barn_liquid_area":None,
    "mms_open_solid_area":None,
    "mms_open_liquid_area":None}

###################################
## MMS parameters
###################################
## adsorption constant for manure; m3/m3
Kd = 1.0
## assuming the roughness height of manure storage barn is ~ 0.5m (<ref height of 2m)
zo_barn = 0.5  
## assuming the roughness height of manure pile (open land) is ~ 1.0m (<ref height of 2m)
zo_manureland = 1.0
## manure surface thickness is set to be 2 cm
z_manuresurf = 0.02
## dry matter (DM) content of solid manure 
DM_content = solid_m_DM[livestock]
## dry matter (DM) content of liquid manure is assumed to be 5%
f_DM_liquid = 0.05
## maximum water content of manure
f_wcmax = 1 - (DM_content/100)/2
## assuming the density of manure; 1t kg/m^3 or 1g/cm^3
manure_density = rho_m[livestock]
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

# groundtemp_datalvl1 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=groundtemp_datalvl1)  ## degC
# groundtemp_datalvl2 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=groundtemp_datalvl2)  ## degC
# rhum_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=rhum_data)  ## per cent
# wind_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=wind_data)  ## m/s
# evap_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=evap_data) ## g/day
# # soilmoist_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_data)  ## m3/m3
# soilmoist_datalvl1 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_datalvl1)  ## m3/m3
# # soilmoist_datalvl2 = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilmoist_datalvl2)  ## m3/m3
# soilbd_ds = open_ds(file_path+soil_data_path+soilbdfile)
# soilbd = soilbd_ds.T_BULK_DEN.values
# soilporosity = 1 - (soilbd/(rho_soil/1000))
# persm_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=soilporosity)  ## per cent
# ram1_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=ram1_data)  ## s/m
# rb1_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=rb1_data)  ## s/m
# # runoff_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=runoff_data)  ## m/day
# # subrunoff_data = field_var_fill(sd_template=animal_file['Excreted_N'][lvl_idx],input_field=subrunoff_data)  ## m/day

##################################
## MMS module
##################################
class MMS_module:
    def __init__(self,array_shape,manure_added,urea_added,UA_added,avail_N_added,resist_N_added,unavail_N_added,\
    TAN_added,water_added,pH_value,area_housing):
        ## show current config settings, e.g., MMS type
        print('MMS Module - current MMS is for production system: '+str(MMS_type))
        ## feces input from housing
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
        ## urea input from housing
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
        ## TAN pool (the bulk manure)
        self.TAN_pool = np.zeros(array_shape)
        ## TAN pool of the soil interface
        self.TAN_pool_soil = np.zeros(array_shape)
        ## washoff flux of TAN 
        self.TAN_washoff = np.zeros(array_shape)
        ## TAN pool conc (aqueous phase)
        self.TAN_amount = np.zeros(array_shape)
        ## TAN pool in molar concentration
        self.TAN_amount_M = np.zeros(array_shape)
        ## TAN conc at the surface; for [MMS (barn,open) solid]
        self.TAN_surf_amount_M = np.zeros(array_shape)
        ## TAN conc at the soil surface/interface between manure and land (soil); for [MMS open solid]
        self.TAN_soil_amount = np.zeros(array_shape)
        ## NO3- from nitrification in bulk manrue 
        self.nitrif_NO3_manure = np.zeros(array_shape)
        ## NO3- washoff 
        self.NO3_washoff = np.zeros(array_shape)
        ## NO3- pool in the bulk manure
        self.NO3_pool = np.zeros(array_shape)
        ## NO3- conc in the bulk manure
        self.NO3_amount = np.zeros(array_shape)
        ## nitrification in the soil interface
        self.nitrif_NO3_soil = np.zeros(array_shape)
        ## NO3- pool in the soil interface
        self.NO3_pool_soil = np.zeros(array_shape)
        ## NO3- conc in the soil interface
        self.NO3_soil_amount = np.zeros(array_shape)
        ## infiltration of NO3 from manure to the soil interface
        self.NO3_infilmanure = np.zeros(array_shape)
        ## diffusive aqueous NO3- from manure to the soil interface
        self.NO3_diffusivemanure = np.zeros(array_shape)
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
        ## NH3 concentration in the bulk manure
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
        ## diffusion of aqueous TAN from manure to the soil interface
        self.diffusivefluxmanure_aq = np.zeros(array_shape)
        ## diffusion of gaseous NH3 from manure to the soil interface
        self.diffusivefluxmanure_gas = np.zeros(array_shape)
        ## infiltration of aqueous TAN to the soil interface
        self.infilflux = np.zeros(array_shape)
        ## diffusion of aqueous TAN to deeper soil
        self.diffusivefluxsoil_aq = np.zeros(array_shape)
        ## diffusion of gaseous NH3 to deeper soil
        self.diffusivefluxsoil_gas = np.zeros(array_shape)
        ## leaching of aqueous TAN to deeper soil (not diffusive)
        self.leachingflux = np.zeros(array_shape)

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
        ## barn/buidling resistance
        self.R_star = np.zeros(array_shape)
        ## manure resistance for aqueous diffusion
        self.R_manurel = np.zeros(array_shape)
        ## manure resistance for gaseous diffusion
        self.R_manureg = np.zeros(array_shape)
        ## soil resistance for aqeous diffusion
        self.R_soilaq = np.zeros(array_shape)
        ## soil resistance for gaseous diffusion
        self.R_soilg = np.zeros(array_shape)
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
        ## infiltration rate; m/s
        self.qinfil = np.zeros(array_shape)
        ## percolation flux; m/s
        self.qpsoil = np.zeros(array_shape)
        ## pH and H+ ions concentration
        self.pH = pH_value
        self.cc_H = np.float(10**(-pH_value))
        ## housing area that is used to determine MMS area
        self.housingarea = area_housing

    def met_input_interp(self,template):
        ##################################
        ## fill land input data
        ##################################
        self.T_sim = field_var_fill(sd_template=template,input_field=self.T_sim)  ## degC
        self.evap_sim = field_var_fill(sd_template=template,input_field=self.evap_sim) ## g/day
        # soilmoist_data = field_var_fill(sd_template=template,input_field=soilmoist_data)  ## m3/m3
        self.soilmoist= field_var_fill(sd_template=template,input_field=self.soilmoist)  ## m3/m3
        self.persm = field_var_fill(sd_template=template,input_field=self.persm)  ## m3/m3
        self.R_star = field_var_fill(sd_template=template,input_field=self.R_star)  ## s/m
        self.soilmoist[self.soilmoist>self.persm] = self.persm[self.soilmoist>self.persm]
        return

    def sim_env(self,mms_type,mms_phase):
        if mms_type == 'MMS_barn':
            self.T_sim,self.u_sim = barn_env(temp_data,\
                wind_profile(uref=wind_data,height_ref=wind_data_height,height_out=ref_height,zo=zo_barn))
            self.RH_sim = rhum_data
            self.T_sim = xr_to_np(self.T_sim)
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
                self.R_star = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap_sim)
                ## convert evap from m/s to g/m2/day
                self.evap_sim = self.evap_sim*1e6*timestep*3600
                print('MMS ENV: mms barn (liquid)')
            elif mms_phase == "solid":
                # self.R_star, self.evap_sim = resistance_manure(temp=self.T_sim,u=u2m,rhum=self.RH_sim)
                # self.evap_sim = self.evap_sim*1e6*timestep*3600
                self.evap_sim = water_evap_a(temp=self.T_sim,rhum=self.RH_sim,u=u2m,zo=zo_barn)
                self.R_star = resistance_water_air(temp=self.T_sim,rhum=self.RH_sim,evap_flux=self.evap_sim)
                ## convert evap from m/s to g/m2/day
                f_evap_correct = 1.0
                self.evap_sim = self.evap_sim*1e6*timestep*3600*f_evap_correct
                print('MMS ENV: mms barn (solid)')

        else:
            self.T_sim = groundtemp_datalvl1
            self.u_sim = wind_data
            self.RH_sim = rhum_data
            self.evap_sim = evap_data
            self.rainfall = rain_data
            self.R_star = ram1_data+rb1_data
            self.soilmoist = soilmoist_datalvl1
            # self.persm = persm_data
            soilbd_ds = open_ds(file_path+soil_data_path+soilbdfile)
            ## buld density unit: kg/dm3
            soilbd = soilbd_ds.T_BULK_DEN.values
            soilporosity = 1 - (soilbd/(rho_soil/1000))
            self.persm[:] = soilporosity
            self.persm[self.persm>1.0] = 0.99

            self.met_input_interp(template=animal_file['Excreted_N'][lvl_idx])

            # self.T_sim = xr_to_np(self.T_sim)
            # self.RH_sim = xr_to_np(self.RH_sim)
            # self.u_sim = xr_to_np(self.u_sim)
            # self.evap_sim = xr_to_np(self.evap_sim)
            # self.soilmoist = xr_to_np(self.soilmoist)
            # self.persm = xr_to_np(self.persm)
            # self.persm[self.persm>1.0] = 0.9999
            # self.rainfall = xr_to_np(self.rainfall)
            # self.R_star = xr_to_np(self.R_star)

            self.daily_KNO3 = nitrification_rate_soil(ground_temp=xr_to_np(groundtemp_datalvl1),theta=self.soilmoist,theta_sat=self.persm,
                                                        pH = self.pH,fer_type="manure")*timestep*3600
            self.daily_KNO3[self.daily_KNO3<0] = 0.0
            self.daily_KNO3[np.isnan(self.daily_KNO3)] = 0.0
            print('MMS ENV: mms open env')

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
    
    ## Simulation: Cat A manure stored in barns (as liquid)
    ## water pool is transfered from housing to MMS barn 
    def MMS_barn_liquid_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: MMS barn (liquid), '+str(livestock))
        ## area of MMS_barn_liquid
        MMS_area["mms_barn_liquid_area"] = self.housingarea*(1.0 - f_loss - f_sold)*f_MMS_barn_liquid*MMS_area_factor['mms_barn_liquid']
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):

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

                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

                ## water pool
                ## Note: the water pool in [MMS barn liquid] is "inheritated" from housing water pool
                # water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd+1]
                # self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] - self.evap_sim[dd][water_idx>0]
                # self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd+1][water_idx<=0]
                # self.Total_water_pool[dd+1] = self.Total_water_pool[dd+1] + \
                #     self.water_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_liquid"])
                self.Total_water_pool[dd+1] = self.manure_pool[dd+1]/f_DM_liquid                                             

                ## TAN pool
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

                ## TAN pool in ug
                TAN_pool_ug = self.TAN_pool[dd+1] * 1e6

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
                NH3_gas_ug = self.NH3_gas_M[dd+1]*14*1e9

                ## determining the maximum emission; emission cannot exceed TAN pool
                emiss_idx = (NH3_gas_ug*3600*timestep/self.R_star[dd+1]) - TAN_pool_ug
                self.modelled_emiss[dd+1][emiss_idx>=0] = TAN_pool_ug[emiss_idx>=0]
                self.modelled_emiss[dd+1][emiss_idx<0] = NH3_gas_ug[emiss_idx<0]*3600*timestep/\
                                                            self.R_star[dd+1][emiss_idx<0]

                ## final emission flux
                self.NH3_flux[dd+1] = self.modelled_emiss[dd+1]/1e6
        return

    ## Simulation: Cat B manure stored in barns (as solid) //// under development 31/Aug
    ## manure/animal waste is processed/dried to reduce water content; 
    ## water pool is determined by:
    ## a) mositure equilibrium (minimum amount of manure moisture content) b) water transferred from housing
    ## NEW) incorporate a [surface compensation point] scheme in the model:
    ##       gaseouos concentrations (chi): chi_indoor; chi_surface; chi_bulk
    ##       aqueous concrntrations: [TAN]_bulk (TAN conc in the manure); [TAN]_surf (TAN conc at the surface, in equilibrium with chi_surf)
    ##       equilibrium between [TAN]_surf and chi_surf: chi_surf = (k_H_D/([H+]+k_NH4+))*[TAN]_surf = KNH3*[TAN]_surf; 
    ##                   Note: k_H_D = (161500/(temp+273.15))*np.exp(-10380/(temp+273.15)) 
    ##       fluxes (F): F_atm, surface to indoor; F_tosurf: manure to surface;;; 
    ##                   F_atm=(chi_surf-chi_indoor)/R_star; 
    ##                   F_tosurfaq=([TAN]_bulk-[TAN]_surf)/R_manureaq; F_tosurfgas=(chi_bulk-chi_surf)/R_manuregas; 
    ##                   F_atm = F_tosurf(aq+gas)
    ##    [TAN]_bulk is the prognostic variable and is determined by source and loss based on the mass balance approach
    ##    solve [TAN]_surf: [TAN]_surf = (chi_indoor/R_ab+[TAN]_bulk*(1/R_manureaq+KNH3/R_manuregas))/
    ##                                      (1/R_manureaq+KNH3*(R_manureg+R_ab))
    def MMS_barn_solid_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: MMS barn (solid), '+str(livestock))
        MMS_area["mms_barn_solid_area"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_barn_solid*MMS_area_factor["mms_barn_solid"]
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):
                ## note the equations for the pools are similar to "barn liquid sim"; see comments above
                ## daily manure and urine on per unit area
                self.manure[dd+1] = self.manure_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])

                ## manure pool
                self.manure_pool[dd+1] = self.manure_pool[dd] + self.manure[dd+1] 

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
                # self.manure_water[dd+1] = self.water_added[dd+1]/(self.housingarea*MMS_area_factor["mms_barn_solid"])
                self.manure_water[dd+1] = (self.manure[dd+1]/(DM_content/100))*(1-DM_content/100)

                ## water amount when mositure content reach equilibrium
                self.manure_minwc[dd+1] = self.manure_pool[dd+1]*self.mois_coeff[dd+1]/100

                ## water pool
                ## Note the difference between the water pool of [MMS barn solid] and [MMS barn liquid]
                ## water pool of [MMS barn solid] is directly determined by the amount of manure as we assumed a dry matter content 
                water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd]
                self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] +\
                                                                self.manure_water[dd+1][water_idx>0] -\
                                                                self.evap_sim[dd][water_idx>0]
                self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd][water_idx<=0] +\
                                                                self.manure_water[dd+1][water_idx<=0] 
                # self.Total_water_pool[dd+1] = self.Total_water_pool[dd] + self.manure_water[dd+1]
                
                ## water content of the manure
                vtotal,manurewc,manure_WFPS = manure_properties(solidmass=self.manure_pool[dd+1],
                                                        watermass=self.Total_water_pool[dd+1])
                manurewc[manurewc>manure_porosity] = manure_porosity
                manure_WFPS[manure_WFPS>1.0] = 1.0
                manure_torl = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='aqueous')
                manure_torg = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='gaseous')

                ## manure resistance
                ## manure resistance is determined by: R = z/(2*tor_manure*D); 
                ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not considered here
                ## layer thickness of manure (DM+water) in meter: z = Vmanure
                z_total = vtotal
                self.R_manurel[dd+1] = z_total/(2*self.D_aq_NH4[dd+1]*manure_torl)
                self.R_manureg[dd+1] = z_total/(2*self.D_air_NH3[dd+1]*manure_torg)
                                                                
                ## TAN pool
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.nitrif_NO3_manure[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

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
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                    (z_total[self.Total_water_pool[dd+1]!=0]*(manurewc[self.Total_water_pool[dd+1]!=0]+\
                    KNH3[self.Total_water_pool[dd+1]!=0]*(manure_porosity-manurewc[self.Total_water_pool[dd+1]!=0])+\
                    (1-manure_porosity)*Kd))
                ## TAN molar conc; mol/L
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
                ## NH3 conc in the bulk manure
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when manure water content equals to porosity
                self.NH3_gas_bulk[dd+1][manurewc==manure_porosity] = 0.0

                ## TAN conc at the surface
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_manurel[dd+1]+KNH3/self.R_manureg[dd+1])/\
                        (KNH3*(1/self.R_star[dd+1]+1/self.R_manureg[dd+1])+1/self.R_manurel[dd+1]))/\
                            (14*1000)

                ## TAN surface conc in g/m3
                TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
                NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

                ## determining the maximum emission; 
                emiss_idx = (NH3_gas_g*3600*timestep/self.R_star[dd+1])

                ## nirification rate of TAN in the bulk manure; daily maximum nitrification rate is 0.1 per day
                KNO3_manure = nitrification_rate_manure(manure_temp=self.T_sim[dd+1],WFPS=manure_WFPS)*timestep*3600
                KNO3_manure[KNO3_manure>0.1] = 0.1
                ## correction for WPFS response
                KNO3_manure[KNO3_manure<0.0] = 0.0
                KNO3_manure[np.isnan(KNO3_manure)] = 0.0
                ## determining the maximum nitrification of TAN
                # nitrif_idx = KNO3_manure*self.TAN_pool[dd]
                ## fraction of [NH4+(aq)]
                f_NH4 = manurewc/(manurewc+\
                    KNH3*(manure_porosity-manurewc)+\
                    (1-manure_porosity)*Kd)*(self.cc_H/(self.cc_H+self.k_NH4[dd+1]))
                nitrif_idx = KNO3_manure*self.TAN_pool[dd]*f_NH4*(1-f_manure_anoxic)

                ## fluxes from the bulk manure to the soil interface
                ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                manureall_loss = emiss_idx + nitrif_idx
                manureloss_idx = self.TAN_pool[dd+1] - manureall_loss
                # if dd+1 == 143:
                #     print("manureloss_idx: ",manureloss_idx[130,363]) 
                self.NH3_flux[dd+1][manureloss_idx>=0] = emiss_idx[manureloss_idx>=0]
                self.nitrif_NO3_manure[dd+1][manureloss_idx>=0] = nitrif_idx[manureloss_idx>=0]
                ## otherwise, we use a weighted distribution
                self.NH3_flux[dd+1][manureloss_idx<0] = self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (emiss_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.nitrif_NO3_manure[dd+1][manureloss_idx<0]= self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (nitrif_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])

                ## NO3- pool of the bulk manure
                self.NO3_pool[dd+1] = self.NO3_pool[dd] + self.nitrif_NO3_manure[dd+1]
                
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
    ## Note) NO3- pools in the bulk manure and soil have similar processes as TAN, which do not include 1) adsorption, 2) gaseous diffustion
    ##       NO3- aqueous diffusion is moderated by a scaling factor regarding the different diffusivity of NO3- from NH4+    
    def MMS_land_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: MMS open (land,solid), '+str(livestock))
        MMS_area["mms_open_solid_area"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_open_solid*MMS_area_factor["mms_open_solid"]
        soilclayds = open_ds(file_path+soil_data_path+soilclayfile)
        soilclay = soilclayds.T_CLAY.values
        Kdsoil = ammonium_adsorption(clay_content=soilclay)
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
                self.qinfil[dd+1] = (self.Total_water_pool[dd]-absorb_factor*self.manure_pool[dd])/1e6
                self.qinfil[dd+1][self.qinfil[dd+1]>dailymaxinfil]=dailymaxinfil
                self.qinfil[dd+1][self.qinfil[dd+1]<0] = 0.0
                ## covert m/day to m/s
                self.qinfil[dd+1] = self.qinfil[dd+1]/(timestep*3600)
                ## soil infiltration flux is given in m/s
                self.qpsoil[dd+1] = infiltration_rate_method(dailyinfil=self.qinfil[dd+1]*timestep*3600,
                                                                theta_sat=self.persm[dd+1],theta=self.soilmoist[dd+1])
                self.qpsoil[dd+1][self.qpsoil[dd+1]<0] = 0.0
                ## justify the water amount; infil flux is the infiltration within the manure (between 0-10mm/day)
                water_idx = self.Total_water_pool[dd]-self.evap_sim[dd]-self.manure_minwc[dd]-self.qinfil[dd+1]*timestep*3600*1e6
                self.Total_water_pool[dd+1][water_idx>0] = self.Total_water_pool[dd][water_idx>0] + self.rainfall[dd+1][water_idx>0] + \
                                                                self.manure_water[dd+1][water_idx>0] - \
                                                                    self.evap_sim[dd][water_idx>0] - \
                                                                        self.qinfil[dd+1][water_idx>0]*timestep*3600*1e6
                self.Total_water_pool[dd+1][water_idx<=0] = self.manure_minwc[dd][water_idx<=0] + self.rainfall[dd+1][water_idx<=0] + \
                                                                self.manure_water[dd+1][water_idx<=0] 
                ## maximum water capcity of the manure pool (g/m2); assuming manure water holding capcity + maximum infiltration
                max_wc = absorb_factor*self.manure_pool[dd] + dailymaxinfil*1e6
                ## justipy whether current water pool exceeds the maximum water holidng capacity
                water_idx2 =  max_wc - self.Total_water_pool[dd+1]
                ## if exceeds: the surplus amount of water acts of "washoff" water, and the water pool equals the maximum wc
                ## rain available for washoff has the unit of m (as an accumulation of a day's rainfall in m)
                self.rain_avail_washoff[dd+1][water_idx2<0] = (-1)*water_idx2[water_idx2<0]/1e6
                self.Total_water_pool[dd+1][water_idx2<0] = max_wc[water_idx2<0]
                ## if not exceeds: "washoff" water is 0 as water is absorbed by the manure
                # self.rain_avail_washoff[dd+1] = 0.0
                self.rain_avail_washoff[dd+1][water_idx2>=0] = 0.0
                self.rain_avail_washoff[dd+1] = self.rain_avail_washoff[dd+1]/(timestep*3600)

                ## water pool of the soil interface
                water_soil_idx = self.qinfil[dd+1]*timestep*3600+z_soil*self.soilmoist[dd+1] - z_soil*self.persm[dd+1]
                self.water_pool_soil[dd+1][water_soil_idx>0] = (z_soil*self.persm[dd+1][water_soil_idx>0])*1e6
                self.water_pool_soil[dd+1][water_soil_idx<=0] = (self.qinfil[dd+1][water_soil_idx<=0]*timestep*3600 + \
                                                                z_soil*self.soilmoist[dd+1][water_soil_idx<=0])*1e6
                soilwc = (self.water_pool_soil[dd+1]/1e6)/z_soil
                soilwc[soilwc>self.persm[dd+1]] = self.persm[dd+1][soilwc>self.persm[dd+1]]

                ## washoff: 1) manure, 2) urea, 3) available org N, 4) reistant org N, 5) TAN
                ## washoff coefficient (m) = washoff water (mm) * washoff (%/mm)
                nonN_washoff_rate = self.rain_avail_washoff[dd+1]*f_washoff_nonN*1e3
                N_washoff_rate = self.rain_avail_washoff[dd+1]*f_washoff_N*1e3
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
                vtotal,manurewc,manure_WFPS = manure_properties(solidmass=self.manure_pool[dd+1],
                                                        watermass=self.Total_water_pool[dd+1])
                manurewc[manurewc>manure_porosity] = manure_porosity
                manure_WFPS[manure_WFPS>1.0] = 1.0
                manure_torl = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='aqueous')
                manure_torg = soil_tuotorsity(theta_sat=manure_porosity,theta=manurewc,phase='gaseous')

                ## manure resistance
                ## manure resistance is determined by: R = z/(2*tor_manure*D); 
                ##        z is the layer thickness of manure; D is molecular diffusivity of NH4+ in water
                ##        tortuosity dependence for aqueous diffusion is not considered here
                ## layer thickness of manure (DM+water) in meter: z = Vmanure
                z_total = vtotal
                self.R_manurel[dd+1] = z_total/(2*self.D_aq_NH4[dd+1]*manure_torl)
                self.R_manureg[dd+1] = z_total/(2*self.D_air_NH3[dd+1]*manure_torg)
                # self.R_manureg[dd+1][manure_torg!=0] = z_total[manure_torg!=0]/\
                #                                     (2*self.D_air_NH3[dd+1][manure_torg!=0]*manure_torg[manure_torg!=0])  
                # ## when water content is zero, gaseous diffusion is ceased by infinite resistance
                # self.R_manureg[dd+1][manure_torg==0] = np.inf
                
                ## soil resistance
                ## soil resistance = distance (z) / (tortuosity for diffusion x diffusivity of the species)
                ## Note the differences between Vira et al., 2020 GMD and Moring et al., 2016 BG; we used Moring et al.
                tor_soil_aq = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=soilwc,phase="aqueous")
                tor_soil_gas = soil_tuotorsity(theta_sat=self.persm[dd+1],theta=soilwc,phase="gaseous")
                self.R_soilaq[dd+1] = d_deepsoil/(tor_soil_aq*self.D_aq_NH4[dd+1])
                self.R_soilg[dd+1] = d_deepsoil/(tor_soil_gas*self.D_air_NH3[dd+1])
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
                TAN_idx = self.TAN_pool[dd] - self.NH3_flux[dd] - self.diffusivefluxmanure_aq[dd] - self.diffusivefluxmanure_gas[dd]-\
                    self.infilflux[dd] - self.nitrif_NO3_manure[dd] - self.TAN_washoff[dd]
                self.TAN_pool[dd+1][TAN_idx>0] = TAN_idx[TAN_idx>0]+self.TAN_prod[dd+1][TAN_idx>0]+self.TAN[dd+1][TAN_idx>0]
                self.TAN_pool[dd+1][TAN_idx<=0] = self.TAN_prod[dd+1][TAN_idx<=0]+self.TAN[dd+1][TAN_idx<=0]

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
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0
                self.TAN_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.TAN_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                    (z_total[self.Total_water_pool[dd+1]!=0]*(manurewc[self.Total_water_pool[dd+1]!=0]+\
                    KNH3[self.Total_water_pool[dd+1]!=0]*(manure_porosity-manurewc[self.Total_water_pool[dd+1]!=0])+\
                    (1-manure_porosity)*Kd))
                ## TAN molar conc; mol/L
                self.TAN_amount_M[dd+1] = self.TAN_amount[dd+1]/(14*1000)
                ## NH3 conc in the bulk manure
                self.NH3_gas_bulk[dd+1] = KNH3 * self.TAN_amount[dd+1]
                ## NH3 conc is 0 when manure water content equals to porosity
                self.NH3_gas_bulk[dd+1][manurewc==manure_porosity] = 0.0

                ## TAN conc at the surface (solved); atmospheric NH3 is ignored
                self.TAN_surf_amount_M[dd+1] = (self.TAN_amount[dd+1]*(1/self.R_manurel[dd+1]+KNH3/self.R_manureg[dd+1])/\
                        (self.rain_avail_washoff[dd+1]+KNH3*(1/self.R_star[dd+1]+1/self.R_manureg[dd+1])+1/self.R_manurel[dd+1]))/\
                            (14*1000)
                ## TAN surface conc in g/m3
                TAN_surf_amount = self.TAN_surf_amount_M[dd+1]*(14*1000)
                ## Gaseous NH3 at the surface
                self.NH3_gas_M[dd+1] = self.TAN_surf_amount_M[dd+1]*KNH3
                ## in ug
                NH3_gas_g = self.NH3_gas_M[dd+1]*14*1000

                ## determining the maximum emission; 
                emiss_idx = (NH3_gas_g*3600*timestep/self.R_star[dd+1])
                ## determining the maximum TAN runoff;
                ## if rain available for runoff already in m/day; don't need to multiply be timestep*3600;
                ## else need to multiply by timestep*3600
                runoff_idx = self.rain_avail_washoff[dd+1]*TAN_surf_amount*timestep*3600
                ## determining the maximum TAN aqueous diffusion to the soil interface;
                ## diffusion is considered to be unidirectional - fro the bulk manure to the soil interface
                diffaq_idx = (self.TAN_amount[dd+1]-self.TAN_soil_amount[dd])/self.R_manurel[dd+1]*timestep*3600
                diffaq_idx[diffaq_idx<=0] = 0.0
                ## when soil mositure is 0, aqueous diffusion stops
                diffaq_idx[soilwc==0] = 0.0
                ## determining the maximum TAN gaseous diffusion to the soil interface
                diffgas_idx = (self.NH3_gas_bulk[dd+1]-self.NH3_gas_soil[dd])/self.R_manureg[dd+1]*timestep*3600
                diffgas_idx[diffgas_idx<0] = 0.0
                ## when the soil moisture reaches the saturation, gaseous diffusion stops
                diffgas_idx[soilwc==self.persm[dd+1]] = 0.0
                ## determining the maximum infiltration of TAN
                infil_idx = self.qinfil[dd+1]*self.TAN_amount[dd+1]*timestep*3600

                ## nirification rate of TAN in the bulk manure; daily maximum nitrification rate is 0.1 per day
                KNO3_manure = nitrification_rate_manure(manure_temp=self.T_sim[dd+1],WFPS=manure_WFPS)*timestep*3600
                KNO3_manure[KNO3_manure>0.1] = 0.1
                ## correction for WPFS response
                KNO3_manure[KNO3_manure<0.0] = 0.0
                KNO3_manure[np.isnan(KNO3_manure)] = 0.0
                ## determining the maximum nitrification of TAN
                # nitrif_idx = KNO3_manure*self.TAN_pool[dd]
                ## fraction of [NH4(aq)]
                f_NH4_manure = manurewc/(manurewc+\
                    KNH3*(manure_porosity-manurewc)+\
                    (1-manure_porosity)*Kd)*(self.cc_H/(self.cc_H+self.k_NH4[dd+1]))
                nitrif_idx = KNO3_manure*self.TAN_pool[dd]*f_NH4_manure*(1-f_manure_anoxic)


                ## fluxes from the bulk manure to the soil interface
                ## if TAN pool > sum of all losses, then each loss equals to its corresponding maximum value
                manureall_loss = emiss_idx + runoff_idx + diffaq_idx + diffgas_idx + infil_idx + nitrif_idx
                # if dd+1 == 143:
                #     print("manureall_loss: ",manureall_loss[130,363]) 
                manureloss_idx = self.TAN_pool[dd+1] - manureall_loss
                # if dd+1 == 143:
                #     print("manureloss_idx: ",manureloss_idx[130,363]) 
                self.NH3_flux[dd+1][manureloss_idx>=0] = emiss_idx[manureloss_idx>=0]
                self.TAN_washoff[dd+1][manureloss_idx>=0] = runoff_idx[manureloss_idx>=0]
                self.diffusivefluxmanure_aq[dd+1][manureloss_idx>=0] = diffaq_idx[manureloss_idx>=0]
                self.diffusivefluxmanure_gas[dd+1][manureloss_idx>=0] = diffgas_idx[manureloss_idx>=0]
                self.infilflux[dd+1][manureloss_idx>=0] = infil_idx[manureloss_idx>=0]
                self.nitrif_NO3_manure[dd+1][manureloss_idx>=0] = nitrif_idx[manureloss_idx>=0]
                ## otherwise, we use a weighted distribution
                self.NH3_flux[dd+1][manureloss_idx<0] = self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (emiss_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.TAN_washoff[dd+1][manureloss_idx<0] = self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (runoff_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.diffusivefluxmanure_aq[dd+1][manureloss_idx<0]= self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (diffaq_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.diffusivefluxmanure_gas[dd+1][manureloss_idx<0] = self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (diffgas_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.infilflux[dd+1][manureloss_idx<0] = self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (infil_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])
                self.nitrif_NO3_manure[dd+1][manureloss_idx<0]= self.TAN_pool[dd+1][manureloss_idx<0]*\
                                                            (nitrif_idx[manureloss_idx<0]/manureall_loss[manureloss_idx<0])

                ## NO3- pool of the bulk manure
                NO3_idx = self.NO3_pool[dd] - self.NO3_infilmanure[dd] - self.NO3_diffusivemanure[dd] - self.NO3_washoff[dd+1]
                self.NO3_pool[dd+1][NO3_idx>0] = NO3_idx[NO3_idx>0] + self.nitrif_NO3_manure[dd+1][NO3_idx>0]
                self.NO3_pool[dd+1][NO3_idx<=0] = self.nitrif_NO3_manure[dd+1][NO3_idx<=0]

                ## NO3- conc of the bulk manure; g/mL
                self.NO3_amount[dd+1][self.Total_water_pool[dd+1]!=0] = self.NO3_pool[dd+1][self.Total_water_pool[dd+1]!=0]/\
                                                                            self.Total_water_pool[dd+1][self.Total_water_pool[dd+1]!=0]
                self.NO3_amount[dd+1][self.Total_water_pool[dd+1]==0] = 0.0
                ## NO3- conc of the bulk manure in g/m3
                self.NO3_amount[dd+1] = self.NO3_amount[dd+1]*1e6
                
                ## diffusive aquous NO3 and infiltration of NO3 from the bulk manure to the soil interface
                ## the diffusion distance for NO3 is considered to be from the surface (oxygen available zone) to the bottom;
                ## therefore, 0.5 is applied (double the resistance by doubling the distance)
                NO3_diffidx = (self.NO3_amount[dd+1] - self.NO3_soil_amount[dd])*(0.5*f_DNO3/self.R_manurel[dd+1])*timestep*3600
                NO3_diffidx[NO3_diffidx<0] = 0.0
                NO3_diffidx[soilwc==0] = 0.0
                NO3_infilidx = self.qinfil[dd+1]*self.NO3_amount[dd+1]*timestep*3600
                NO3_lossall = NO3_diffidx + NO3_infilidx 
                manureloss_idx = self.NO3_pool[dd+1] - NO3_lossall
                self.NO3_diffusivemanure[dd+1][manureloss_idx>=0] = NO3_diffidx[manureloss_idx>=0]
                self.NO3_infilmanure[dd+1][manureloss_idx>=0] = NO3_infilidx[manureloss_idx>=0]
                self.NO3_diffusivemanure[dd+1][manureloss_idx<0] = self.NO3_pool[dd+1][manureloss_idx<0]*\
                                                                        NO3_diffidx[manureloss_idx<0]/NO3_lossall[manureloss_idx<0]                                           
                self.NO3_infilmanure[dd+1][manureloss_idx<0] = self.NO3_pool[dd+1][manureloss_idx<0]*\
                                                                        NO3_infilidx[manureloss_idx<0]/NO3_lossall[manureloss_idx<0]                                         

                ## TAN pool in the soil interface
                TAN_soil_idx = self.TAN_pool_soil[dd] - self.diffusivefluxsoil_aq[dd] - self.diffusivefluxsoil_gas[dd] - \
                    self.leachingflux[dd] - self.nitrif_NO3_soil[dd] 
                self.TAN_pool_soil[dd+1][TAN_soil_idx>0] = TAN_soil_idx[TAN_soil_idx>0] + self.infilflux[dd+1][TAN_soil_idx>0] + \
                                                            self.diffusivefluxmanure_aq[dd+1][TAN_soil_idx>0] + \
                                                            self.diffusivefluxmanure_gas[dd+1][TAN_soil_idx>0]
                self.TAN_pool_soil[dd+1][TAN_soil_idx<=0] =  self.infilflux[dd+1][TAN_soil_idx<=0]+ \
                                                                self.diffusivefluxmanure_aq[dd+1][TAN_soil_idx<=0] +\
                                                                 self.diffusivefluxmanure_gas[dd+1][TAN_soil_idx<=0]

                ## TAN conc at the soil surface/interface between manure and land (soil); g/m3
                self.TAN_soil_amount[dd+1][self.water_pool_soil[dd+1]==0] = 0.0
                self.TAN_soil_amount[dd+1][self.water_pool_soil[dd+1]!=0] = self.TAN_pool_soil[dd+1][self.water_pool_soil[dd+1]!=0]/\
                    (z_soil*(soilwc[self.water_pool_soil[dd+1]!=0]+\
                    KNH3[self.water_pool_soil[dd+1]!=0]*(self.persm[dd+1][self.water_pool_soil[dd+1]!=0]-\
                    soilwc[self.water_pool_soil[dd+1]!=0])+\
                    (1-self.persm[dd+1][self.water_pool_soil[dd+1]!=0])*Kdsoil[self.water_pool_soil[dd+1]!=0]))
                ## NH3 concentration in the soil pore space
                self.NH3_gas_soil[dd+1] = KNH3 * self.TAN_soil_amount[dd+1]
                ## NH3 conc is zero when soil moisture content reaches the saturation
                self.NH3_gas_soil[dd+1][soilwc==self.persm[dd+1]] = 0.0 

                ## fluxes to deeper soil
                ## TAN loss through aqueous diffusion and leaching to deeper soil
                soildiffaq_idx = self.TAN_soil_amount[dd+1]/self.R_soilaq[dd+1]*timestep*3600
                soildiffaq_idx[self.soilmoist[dd+1]==0.0] = 0.0
                soildiffgas_idx = self.NH3_gas_soil[dd+1]/self.R_soilg[dd+1]*timestep*3600
                soildiffgas_idx[self.soilmoist[dd+1]==self.persm[dd+1]] = 0.0
                soilleaching_idx = self.qpsoil[dd+1]*self.TAN_soil_amount[dd+1]*timestep*3600
                ## fraction of [NH4+(aq)] in soil
                f_NH4_soil = soilwc/(soilwc+\
                    KNH3*(self.persm[dd+1]-soilwc)+\
                    (1-self.persm[dd+1])*Kdsoil)*(self.cc_H/(self.cc_H+self.k_NH4[dd+1]))
                soilnitrif_idx =  self.daily_KNO3[dd+1]*self.TAN_pool_soil[dd]*f_NH4_soil
                soilall_loss = soildiffaq_idx + soildiffgas_idx + soilleaching_idx + soilnitrif_idx
                soilloss_idx = self.TAN_pool_soil[dd+1] - soilall_loss
                self.diffusivefluxsoil_aq[dd+1][soilloss_idx>=0] = soildiffaq_idx[soilloss_idx>=0]
                self.diffusivefluxsoil_gas[dd+1][soilloss_idx>=0] = soildiffgas_idx[soilloss_idx>=0]
                self.leachingflux[dd+1][soilloss_idx>=0] = soilleaching_idx[soilloss_idx>=0]
                self.nitrif_NO3_soil[dd+1][soilloss_idx>=0] = soilnitrif_idx[soilloss_idx>=0]
                self.diffusivefluxsoil_aq[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                    soildiffaq_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
                self.diffusivefluxsoil_gas[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                    soildiffgas_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
                self.leachingflux[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                    soilleaching_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]
                self.nitrif_NO3_soil[dd+1][soilloss_idx<0] = self.TAN_pool_soil[dd+1][soilloss_idx<0]*\
                                                                    soilnitrif_idx[soilloss_idx<0]/soilall_loss[soilloss_idx<0]

                ## NO3- pool of the soil interface
                NO3_soil_idx = self.NO3_pool_soil[dd] - self.NO3_leaching[dd] - self.NO3_diffusivesoil[dd] 
                self.NO3_pool_soil[dd+1][NO3_soil_idx>0] = NO3_soil_idx[NO3_soil_idx>0] + self.nitrif_NO3_soil[dd+1][NO3_soil_idx>0] + \
                                self.NO3_infilmanure[dd+1][NO3_soil_idx>0] + self.NO3_diffusivemanure[dd+1][NO3_soil_idx>0]
                self.NO3_pool_soil[dd+1][NO3_soil_idx<=0] = self.nitrif_NO3_soil[dd+1][NO3_soil_idx<=0] + \
                                self.NO3_infilmanure[dd+1][NO3_soil_idx<=0] + self.NO3_diffusivemanure[dd+1][NO3_soil_idx<=0]

                ## NO3- conc of the soil interface; g/mL
                self.NO3_soil_amount[dd+1][self.water_pool_soil[dd+1]!=0] = self.NO3_pool_soil[dd+1][self.water_pool_soil[dd+1]!=0]/\
                                                                            self.water_pool_soil[dd+1][self.water_pool_soil[dd+1]!=0]   
                self.NO3_soil_amount[dd+1][self.water_pool_soil[dd+1]==0] = 0.0
                ## NO3- conc of the soil interface in g/m3
                self.NO3_soil_amount[dd+1] = self.NO3_soil_amount[dd+1]*1e6
                                                      
                ## NO3 loss through aqueous diffusion and leaching to deeper soil
                NO3_soildiffidx = self.NO3_soil_amount[dd+1]*(f_DNO3/self.R_soilaq[dd+1])*timestep*3600
                NO3_soilleachingidx = self.qpsoil[dd+1]*self.NO3_soil_amount[dd+1]*timestep*3600
                NO3_soilall_loss = NO3_soildiffidx + NO3_soilleachingidx
                soilloss_idx = self.NO3_pool_soil[dd+1] - (self.NO3_soil_amount[dd+1]*(f_DNO3/self.R_soilaq[dd+1] + self.qpsoil[dd+1])*timestep*3600)
                self.NO3_diffusivesoil[dd+1][soilloss_idx>=0] = NO3_soildiffidx[soilloss_idx>=0]
                self.NO3_leaching[dd+1][soilloss_idx>=0] = NO3_soilleachingidx[soilloss_idx>=0]
                self.NO3_diffusivesoil[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                                 NO3_soildiffidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
                self.NO3_leaching[dd+1][soilloss_idx<0] = self.NO3_pool_soil[dd+1][soilloss_idx<0]*\
                                                            NO3_soilleachingidx[soilloss_idx<0]/NO3_soilall_loss[soilloss_idx<0]
                
        return
    
    ## Simulation: Cat D manure stored in open environment (as solid) //// under development 08/Sep
    def MMS_liquid_sim(self,start_day_idx,end_day_idx):
        print('current simulation is for: MMS open (lagoon,liquid), '+str(livestock))
        MMS_area["mms_open_liquid_area"] = self.housingarea*(1.0-f_loss-f_sold)*f_MMS_open_liquid*MMS_area_factor["mms_open_liquid"]
        if livestock.lower()=="poultry":
            # for dd in np.arange(start_day_idx,end_day_idx-1):
            print(livestock)

        else:
            for dd in np.arange(start_day_idx,end_day_idx):
                KNH3 = self.Henry_constant[dd+1]/(self.cc_H + self.k_NH4[dd+1])
                ## wind speed at 8m
                u8 = wind_profile(uref=self.u_sim[dd+1],height_ref=10,height_out=8,zo=8e-5)
                ## mass transfer coefficient in the liquid boundary layer
                k_aq = k_aq_NH4(wind_8m=u8,temp=self.T_sim[dd+1])
                ## mass tranfer coefficient in the gas boundary layer
                k_gas = k_gas_NH3(wind_8m=u8,temp=self.T_sim[dd+1])
                ## surface TAN conc in mol/L
                self.TAN_surf_amount_M[dd+1] = ((k_aq*lagoon_TAN_conc)/(KNH3*k_gas+k_aq))/14*1000
                ## gaseous NH3 conc at surface in mol/L
                self.NH3_gas_M[dd+1] = KNH3*self.TAN_surf_amount_M[dd+1]
                ## NH3 conc in g/m3
                NH3_gas = self.NH3_gas_M[dd+1] *14*1000
                self.NH3_flux[dd+1] = NH3_gas*k_gas*timestep*3600
        return

    def MMS_sim_main(self,mms_cat,phase,start_day_idx,end_day_idx):
        if mms_cat == "MMS_barn":
            self.sim_env(mms_type=mms_cat,mms_phase=phase)
            if phase == 'solid':
                self.MMS_barn_solid_sim(start_day_idx,end_day_idx)
            elif phase == 'liquid':
                self.MMS_barn_liquid_sim(start_day_idx,end_day_idx)
        elif mms_cat == "MMS_open":
            self.sim_env(mms_type=mms_cat,mms_phase=phase)
            if phase == 'solid':
                self.MMS_land_sim(start_day_idx,end_day_idx)
            elif phase == 'liquid':
                self.MMS_liquid_sim(start_day_idx,end_day_idx)
        return
