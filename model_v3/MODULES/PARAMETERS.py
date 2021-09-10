####################################
## import python packages
####################################
import numpy as np
import xarray as xr
import pandas as pd
from collections import defaultdict
from pathlib import Path
import scipy.stats
import os
import sys

###################################
## import essential AMCLIM modules
###################################
file_path = os.getcwd()
module_path = str(Path(file_path).parent)
if module_path not in sys.path:
    sys.path.append(module_path)

from CONFIG.config import *
from INPUT.input import *
from MODULES.FUNC import *

###################################
## define parameter functions
###################################
## housing indoor conditions, e.g. temperature, ventilation, humidity as functions of natural environment
## ref: Gyldenkærne et al., 2005 for general livestock housing env; Jiang et al, 2021 for poultry
def housing_env(temp,rhum,livestock_type,production_system_type):
    temp = np.array(temp)
    rhum_in = np.zeros(temp.shape)
    u_in = np.zeros(temp.shape)
    if livestock_type.lower() == 'cattle':
        t_in = temp
        u_in[temp<0] = 0.2
        u_in[(temp>0)&(temp<12.5)] = 0.2+temp[(temp>0)&(temp<12.5)]*(0.18/12.5)
        u_in[temp>=12.5] = 0.38
        rhum_in = rhum
    if livestock_type.lower() == 'dairy':
        t_in = temp
        u_in[temp<0] = 0.2
        u_in[(temp>0)&(temp<12.5)] = 0.2+temp[(temp>0)&(temp<12.5)]*(0.18/12.5)
        u_in[temp>=12.5] = 0.38
        rhum_in = rhum
    elif livestock_type.lower() == 'pig':
        t_in = np.zeros(temp.shape)
        t_in[temp<0] = 20+0.25*temp[temp<0]
        t_in[(temp>0)&(temp<12.5)] = 20
        t_in[temp>=12.5] = 20+(temp[temp>=12.5]-12.5)
        u_in[temp<0] = 0.2
        u_in[(temp>0)&(temp<12.5)] = 0.2+temp[(temp>0)&(temp<12.5)]*(0.18/12.5)
        u_in[temp>=12.5] = 0.38
        rhum_in = rhum - 10
    elif livestock_type.lower() == 'poultry':
        if production_system_type.lower() == 'broiler':
            t_in = 0.00020*temp**3 + 0.0010*temp**2 + 0.024*temp + 22.1        ## Jiang et a., 2021
            u_in[temp<0] = 0.2  
            u_in[(temp>0)&(temp<12.5)] = 0.2+temp[(temp>0)&(temp<12.5)]*(0.23/12.5)
            u_in[temp>=12.5] = 0.43 
            rhum_in = rhum
        elif production_system_type.lower() == 'layer':
            t_in = 0.00014*temp**3 + 0.0023*temp**2 + 0.011*temp + 23.8        ## Jiang et a., 2021
            u_in[temp<0] = 0.2
            u_in[(temp>0)&(temp<12.5)] = 0.2+temp[(temp>0)&(temp<12.5)]*(0.23/12.5)
            u_in[temp>=12.5] = 0.43
            rhum_in = rhum
    return t_in, u_in, rhum_in

## barn conditions, including open barns for livestock; ref: Gyldenkærne et al., 2005
## barn temperature is distinct to natural temperature
def barn_env(temp,wind):
    temp = np.array(temp)
    wind = np.array(wind)
    t_barn = np.zeros(temp.shape)
    u_barn = np.zeros(wind.shape)
    
    Gnd_temp = 4.0
    tempdiff = 4.0
    temp_idx = Gnd_temp - tempdiff  
    t_barn[temp<=temp_idx] = Gnd_temp
    t_barn[temp>temp_idx] = temp[temp>temp_idx] + tempdiff
    
    blocking_factor = 0.8
    u_barn[temp<=temp_idx] = wind[temp<=temp_idx]
    u_barn[temp>temp_idx] = (1-blocking_factor)*wind[temp>temp_idx]
    return t_barn, u_barn

## rate: uric acid hydrolysis to TAN; temp in degC, rhum in per cent
def ua_hydrolysis_rate(temp,rhum,ph):
    ## maximum daily hydrolysis rate is 20%
    dmax_rate = 0.2
    ft = np.array(np.exp(0.165*(temp-10) + 1.8)/np.exp(0.165*(35-10) + 1.8))
    ft[ft>1.0] = 1.0
    erh = np.array((-np.log(1.01 - (rhum/100))/(0.0000534*(temp+273.15)))**(1/1.41))
    frh = np.array(np.round(0.0025*np.exp(0.1676*erh),3))
    frh[erh>35.5] = 1.0
    fph = (1.34*ph - 7.2)/(1.34*9 - 7.2)
    ## daily conversion rate
    drate = dmax_rate*ft*frh*fph
    return drate
    
## rate: urea hydrolysis to TAN; temp in degC, delta_t is the time step, i.e 1 hour, daily res: delta_t=24
# def urea_hydrolysis_rate(temp,WFPS,delta_t):
#     ## T is the soil temperature (or air temperature for housing)
#     ## WFPS is water-filled porosity, i.e. in this simulation, is the moisture content of manure for housing
#     ## it remains unclear of k_h and WFPS
#     ## k_h = 1.2 * WFPS
def urea_hydrolysis_rate(temp,delta_t): 
    ## k_h equals 0.23 per hour at 20 degC; Goh and Sherlock, 1985; Muck, 1981;
    k_h = 0.23
    ## Ah_t is a temeprature scaling factor for k_h; temperature dependence Q10 is ~2. 
    Ah_t = 0.25 * np.exp(0.0693*temp)
    hydrolysis_rate = 1 - np.exp((-k_h*delta_t)*Ah_t)  
    return hydrolysis_rate
 
## rate: TAN production from the decompostion of N_avail and N_resist; temp in degC
## ref: Vigil and Kissel (1995)
## other ref: CLM_FANv1 (Riddick et al., 2016) and FAN_v2 (Julius et al., 2020?)
def N_pools_decomp_rate(temp,delta_t):    
    ## temp coefficient and dependence
    t_r1 = 0.0106
    t_r2 = 0.12979
    f_T = t_r1*np.exp(t_r2*temp)    
    ## coefficients
    ## unit: hr^-1
    B_a = 8.94 * 1e-7 * 3600
    B_r = 6.38 * 1e-8 * 3600    
    ka = B_a * f_T 
    kr = B_r * f_T 
    k_a = 1 - np.exp(-ka*delta_t) 
    k_r = 1 - np.exp(-kr*delta_t) 
    return k_a, k_r    

## physical variables: specific humidity; temp in degC, rhum in per cent
def humidity_measures(temp, rhum):
    M_w = 18.0                ## 18.0 g/mol, molar mass of water vapor;
    R_w = 461.52              ## 461.52 j/(kg K), specifica gas constant for water vapor; R_w = 1000R/M_w;
    ## e_sat, saturation vapor pressure METHOD Two
    T = temp + 273.15
    ## saturated water vapor pressure, use a different method compare to function (humidity_measures)
    e_sat = 0.6108*np.exp(17.27*temp/(temp+237.3))      # T is K, e_sat in kPa
    e_sat = e_sat * 1000
    ## water vapor pressure
    e_vp = e_sat*rhum/100          # vapor pressure = RH x saturation vapor pressure; RH is in percentage
    ## specific humidity
    p = 101325                ## 101,325 Pa, standard atmospheric pressure;
    M_d = 28.96               ## 28.96 g/mol, molar mass of dry air;
    q_atm = ((M_w/M_d)*e_vp)/(p-(1-(M_w/M_d))*e_vp)
    q_sat = ((M_w/M_d)*e_sat)/(p-(1-(M_w/M_d))*e_sat)
    return q_sat, q_atm, e_sat, e_vp

## physical variables: water evaporation in houses; temp in degC, rhum in per cent, Rn in J/m^2/s, u in m/s
## aerodynamic method: similar to Penman's method based on the site simulation
def water_evap_a(temp,rhum,u,zo):
    ## pressure in the atmoshpere in pa
    Pressure = 101.3 * 1e3
    ## saturated water vapor pressure, use a different method compare to function (humidity_measures)
    es = 0.6108*np.exp(17.27*temp/(temp+237.3))      # temp in degC
    ## water vapor pressure
    ea = es*rhum/100

    # zo = 0.002               ## roughness height 2mm
    Z = 2                     ## reference height 2m
    k = 0.41                  ## van Karman constant
    rho_air = 1.2754          ## 1.2754 kg/m^3, density of air;
    rho_water = 997           ## 997kg/m^3, density of water;
    B = (rho_air/rho_water)*(0.622*k**2*u)/(Pressure*(np.log(Z/zo))**2)
    evap = B*(es-ea)*1000
    evap = evap*1000*3600*24
    return evap

## physical variables: molecular diffusivity of NH4+ in water; D_aq_nh4, m^2/s (Van Der Molen et al.,1990; Vira et al., GMD2020)
## temp in degC
def diffusivity_NH4(temp):
    d_aq_nh4 = 9.8e-10*1.03**temp
    return d_aq_nh4

## soil characteristics: tortuosity for diffusion
## theta is the volumetric soil water content, and theta_sat is the volumetric soil water content at saturation (equivalent as porosity)
## theta in per cent
def soil_tuotorsity(theta_sat,theta,phase):
    ## convert per cent to float
    theta_sat = theta_sat/100
    theta = theta/100
    ## soil tuotorsity in aqueous phase and gaeous phase (Millington and Quirk, 1961)
    if phase == 'aqueous':
        soil_tor = ((theta)**(10/3))/(theta_sat**2)
    elif phase == 'gaseous':
        soil_tor = ((theta_sat-theta)**(10/3))/(theta_sat**2)
    return soil_tor

## soil characteristics: infiltration rate (m/s)
## this is an empirically-derived expression for vertical/percolation/infiltration/subsurface leaching/ of animal slurry
## data source: Patle et al., Geo.Eco.Landscale 2019
## loamy sand and sandy loam in 14 sites (quality control; remove bad quality data)
## multilinear regression: R^2 = 0.75; 
## soil para: sand (%), clay (%), bulk density (g/cm^3), particle density (g/cm^3)
## percentage of saturation soil moisture 
def infiltration_rate(soilmoist_percent,sand_percent,clay_percent,bulk_density,particle_density):
    ## soil parameters matrix
    soil_para = 0.62*sand_percent - 0.17*clay_percent - 26.42*bulk_density + 3.14*particle_density - 10.58
    ## soil moisture dependence
    infil_func = soil_para + 6.31*(soilmoist_percent/100)
    ## Ks is the permeability coefficient at saturation for loamy sand, Ks=0.714cm/h (0.0119cm/min;171.36mm/day) ref: Hu et al., 2017 J.Arid.Land
    Ks = 0.714
    ## infiltration rate Ki m/s
    Ki = (infil_func/Ks)/(100*3600)
    return Ki

## resistance: resistance for water-air exchange; temp in degC, rhum in per cent
def resistance_water_air(temp,rhum,evap_flux):
    T = temp + 273.15
    ## to avoid dividing zero
    rhum = rhum-0.001
    ## evap flux in m/s (kg/m^2/s)
    evap_flux = evap_flux/(1000*24*3600)
    Q_sat, Q_atm, e_sat, e_vp = humidity_measures(temp, rhum)
    ## molecular weight of air (28.96 g/mol), NH3 (17 g/mol) and water (18 g/mol).
    M_air = 28.96
    M_NH3 = 17
    M_H2O = 18
    ## sigma_v is the value of atomic diffusion volume in [lagoon system] (from Liley et al., 1984;
    ## Physical and chemical data)
    sigma_v_air = 20.1
    sigma_v_NH3 = 14.9
    sigma_v_H2O = 12.7
    ## pressure in the atmoshpere in pa
    Pressure = 101.3 * 1e3
    D_air_NH3 = (1e-7*(T)**1.75*((M_air+M_NH3)/M_air*M_NH3)**
                 0.5)/(Pressure*(sigma_v_air**(1/3)+sigma_v_NH3**(1/3))**2)
    D_air_H2O = (1e-7*(T)**1.75*((M_air+M_H2O)/M_air*M_H2O)**
                 0.5)/(Pressure*(sigma_v_air**(1/3)+sigma_v_H2O**(1/3))**2)
    rho_air = 1.2754          ## 1.2754 kg/m^3, density of air;
    rho_water = 997           ## 997kg/m^3, density of water;
    Rc = (rho_air/rho_water)*((Q_sat-Q_atm)/evap_flux)
    return Rc

## resistance: resistance for manure in houses and in storage; temp in degC, u in m/s
def resistance_manure(temp, u):
    ## converting degC to Kelvin
    T = temp + 273.15
    ## avoid dividing zero
    u = np.array(u)
    u[u<=0] = 1e-5
    u[np.isnan(u)] = 1e-5
    # resistance with no crust; m/day; Pinder et al., 2004; probability > 0.3
    # we assumed that this corresponds to T = 25, u = 0.1 m/s (low ventilaiton), R = 0.04 day/m ~ 3400 s/m
    r_standard = 0.04
    ## mass transfer coefficient = A*V^0.8*T^(-1.4); A is a parameterised constant; Muck & Steenhius, 1982;
    a = 0.8
    b = -1.4
    T_correction = T**b/(298.15**b)
    u_correction = u**a/(0.1**a)
    r_correct = r_standard/(T_correction*u_correction)
    r_ab_star = r_correct*24*3600
    return r_ab_star

## resistance: resistance for aerodynamic and boundary layer resistance; 
## temp in degC; rhum in %; u in m/s; H is sensible heat flux in J/(m^2 s); Z is reference height in m; zo is surface roughness in m 
def resistance_aero_boundary(temp,rhum,u,H,Z,zo):
    ## temp in K
    T = temp + 273.15
    ## latent heat of vaporization; MJ per kg
    lambda_v = 2.45
    ## atmospheric pressuure; Pa
    p = 101325
    ## saturated water vapor pressure, use a different method compare to function (humidity_measures)
    es = 0.6108*np.exp(17.27*temp/(temp+237.3))      # T is K, es in kPa here
    es = es * 1000
    ## water vapor pressure
    ea = es*rhum/100                           # RH is in percentage
    ## water vapor deficit
    #delta_e = es - ea
   
    ## specific humidity
    q = (0.622*ea)/(p-0.378*ea)
    ## virtual temperature
    T_v = T/(1+0.608*q)
    ## specific gas constant of dry air: 287 J/kg/K
    R = 287.0
    ## acceleration of gravity
    g = 9.81
    ## Karman constant
    k = 0.41
    ## heat capcity of air: 1005 J/kg/K
    cpair = 1005.0
    ## density of air
    pho = p/(R*T_v)
    ## friction velocity
    ustar = k*u/(np.log(Z/zo))
    ## Monin-Obukhov length
    L = -T*(ustar**3)*pho*cpair/(k*g*H)
    ## stability correction function: psi
    u = np.array(u)
    H = np.array(H)
    psi = np.zeros(u.shape)
    ## stable condition
    psi[H<=0] = -5*Z/L[H<=0]  
    ## unstable condition
    X = np.zeros(u.shape)
    X[H>0] = (1-16*Z/L[H>0])**0.25
    psi[H>0] = np.log(((1+X[H>0])/2)**2)+np.log((1+X[H>0]**2)/2)-2*np.arctan(X[H>0])+np.pi/2
    ## aerodynamic resistance
    Ra = (np.log(Z/zo)-psi)**2/(k**2*u)

    ## coefficient B that is used to determine boundary layer resistance 
    B = 5
    Rb = 1/(B*ustar)
    return Ra, Rb

## aniaml info: waste N should be consistent to livestock_N
## unit in kg N per head per year; returning daily values
def livestock_waste_info(livestock_type, waste_N):
    ## daily N excretion from urine and feces
    dN = 1000*waste_N/365
    N_durine = dN * frac_N[livestock_type]['urine_N']
    N_durea = N_durine * frac_urea[livestock_type]
    N_ddung = dN * frac_N[livestock_type]['dung_N']
    ## daily urination and excretion
    durine = N_durine/conc_N[livestock_type]['conc_N_urine']
    ddung = m_DM[livestock_type] * N_ddung/conc_N[livestock_type]['conc_N_dung']
    ## convert g water per kg SM to L water per kg SM
    dmanure_water = (1000-m_DM[livestock_type])* N_ddung/conc_N[livestock_type]['conc_N_dung']/1000
    pH = pH_info[livestock_type]
    return N_durine, N_durea, N_ddung, durine, ddung, dmanure_water, pH

############################################
## livestock waste info database
############################################
livestock_N = defaultdict(dict)
livestock_Nrate = defaultdict(dict)
livestock_weight = defaultdict(dict)
stocking_desity = defaultdict(dict)
frac_N = defaultdict(dict)
conc_N = defaultdict(dict)
frac_urea = defaultdict(dict)
m_DM = defaultdict(dict)
solid_m_DM = defaultdict(dict)
pho_m = defaultdict(dict)
pH_info = defaultdict(dict)
name = ['CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
## regions include: 1)North America, 2)Western Europe, 3) Eastern Europe, 4)Oceania, 5)Latin America, 6)Africa
##                  7) Middle East, 8)Asia, 9) India
region = ['NA','WE','EE','OC','LA','AF','ME','AS','IN']
N_type = ['urine_N','dung_N']
conc_type = ['conc_N_urine','conc_N_dung']
## ref: Vira, J., Hess, P., Melkonian, J., and Wieder, W. R.: An improved mechanistic model for ammonia volatilization 
##      in Earth system models: Flow of Agricultural Nitrogen, version 2 (FANv2), 
##      Geosci. Model Dev. Discuss., https://doi.org/10.5194/gmd-2019-233, in review, 2019.
##      in supplementary materials, sect.2.2, Table 1.
N_values = [[58.0,65.0,55.3,65.5,44.2,42.6,52.7,42.4,25.4],
          [97.0,105.1,70.3,80.3,70.1,60.2,70.3,60.0,47.2],
          [44.0,50.6,50.0,60.2,40.1,39.8,49.9,39.6,13.7], 
          [11.2,16.1,17.0,15.6,16.8,16.8,16.8,5.1,5.1],
          [7.1,9.3,10.0,8.7,16.0,16.0,16.0,4.3,4.3],
          [17.3,30.4,30.2,30.2,5.6,5.6,5.6,2.5,2.5],
          [7.4,15.0,15.9,20.0,12.0,12.0,12.0,12.0,12.0],
          [6.3,18.0,18.0,20.0,15.0,15.0,15.0,15.0,15.0],
          [0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5,0.5],
          [44.4,44.4,44.4,44.4,44.4,44.4,44.4,44.4,34.5]]
## Nitrogen excretion rates: kg N per 1000 kg animal mass per day; (equivalent to: g N per kg animal mass per day)
N_rates =[[0.29,0.31,0.32,0.44,0.35,0.59,0.74,0.34],
          [0.44,0.48,0.35,0.44,0.48,0.60,0.70,0.47],
          [0.31,0.33,0.35,0.50,0.36,0.63,0.79,0.34],
          [0.40,0.50,0.54,0.52,1.47,1.47,1.47,0.40],
          [0.42,0.51,0.55,0.53,1.57,1.57,1.57,0.42],
          [0.24,0.42,0.46,0.46,0.55,0.55,0.55,0.24],
          [0.42,0.85,0.90,1.13,1.17,1.17,1.17,1.17],
          [0.45,1.28,1.28,1.42,1.37,1.37,1.37,1.37],
          [0.83,0.83,0.82,0.82,0.82,0.82,0.82,0.82],
          [0.32,0.32,0.32,0.32,0.32,0.32,0.32,0.32]]
## livestock stocking density in animal houses (needs to be updated)
name = ['CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
den_stock = [200.0, 200.0, 200.0, 120.0, 120.0, 60.0, 100.0, 100.0, 30.0, 200.0]
## fraction of urine N and dung N; proportion
f_N = [[1.0/2, 1.0/2], [8.8/13.8, 5/13.8], [1.0/2, 1.0/2], [2.0/3, 1.0/3],[2.0/3, 1.0/3],[2.0/3, 1.0/3],
       [1.0/2, 1.0/2], [1.0/2, 1.0/2], [0.0, 1.0],[1.0/2, 1.0/2]]
## N concentration in urine and dung; g N per L urine,  g N per kg SM
## for poultry, urine N is 0, and c_N_urine is set to 1 for convenience
c_N = [[4.40, 4.85], [9.00, 4.85], [4.40, 4.85], [4.90, 10.45], [4.90, 10.45], [4.90, 10.45],
       [12.60, 6.40], [12.60, 6.40], [1.00, 29.60], [4.40, 4.85]]
## fraction of dry matter in solid manure; g DM per kg SM
## ref: 1. Sommer and Hutchings, Ammonia emission from field applied manure and its reduction -- invited paper,
##      Europ. J. Agronomy 15 (2001) 1-15; (for cattle, pig and poultry)
##      2. Zhao et al., Nitrogen utilization efficiency and prediction of nitrogen excretion in sheep
##      offered fresh perennial ryegrass, J. of. Animal Science, 2016; (for sheep/goat)
f_DM = [181.5, 181.5, 181.5, 222.0, 222.0, 222.0, 155.0, 155.0, 574.0, 181.5]
## the dry matter (DM) content of solid manure, assuming 20% for cattle, pigs etc, 50% for poultry
f_solid_DM = [20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 50.0, 20.0]
## assuming the density of manure; 1t/m^3 or 1g/cm^3 for cattle, pigs etc, 0.4 for poultry
pho_manure = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 1.0]
## fraction of urinal N in the form of urea
f_urea = [0.75, 0.75, 0.75, 0.75, 0.75, 0.75, 0.80, 0.80, 0, 0.75]
## pH value of livestock slurry
## ref: 1. Sommer and Hutchings, Ammonia emission from field applied manure and its reduction -- invited paper,
##      Europ. J. Agronomy 15 (2001) 1-15; (for cattle, pig and poultry)
pH_val = [7.8, 7.8, 7.8, 7.7, 7.7, 7.7, 8.0, 8.0, 8.5, 7.8]
for ii in np.arange(10):
    for jj in np.arange(8):
        livestock_N[name[ii]][region[jj]] = N_values[ii][jj]
        livestock_Nrate[name[ii]][region[jj]] = N_rates[ii][jj]
        livestock_weight[name[ii]][region[jj]] = 1000*N_values[ii][jj]/(N_rates[ii][jj]*365)
        stocking_desity[name[ii]] = den_stock[ii]
for ii in np.arange(10):
    for jj in np.arange(2):
        frac_N[name[ii]][N_type[jj]] = f_N[ii][jj]
        conc_N[name[ii]][conc_type[jj]] = c_N[ii][jj]
        frac_urea[name[ii]] = f_urea[ii]
        m_DM[name[ii]] = f_DM[ii]
        solid_m_DM[name[ii]] = f_solid_DM[ii]
        pho_m[name[ii]] = pho_manure[ii]
        pH_info[name[ii]] = pH_val[ii] 




