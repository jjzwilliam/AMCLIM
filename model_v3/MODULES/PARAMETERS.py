####################################
## import python packages
####################################
import numpy as np
from numpy.lib.function_base import piecewise
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


###################
## ENV/Management
###################
## EMPIRICAL - housing indoor conditions, e.g. temperature, ventilation, humidity as functions of natural environment
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
## EMPIRICAL - barn conditions, including open barns for livestock; ref: Gyldenkærne et al., 2005
## barn temperature is distinct to natural temperature
def barn_env(temp,wind):
    temp = np.array(temp)
    wind = np.array(wind)
    t_barn = np.zeros(temp.shape)
    t_gnd = np.zeros(temp.shape)
    u_barn = np.zeros(wind.shape)
    
    Gnd_temp = 4.0
    tempdiff = 4.0
    temp_idx = Gnd_temp - tempdiff  
    t_gnd[temp<=temp_idx] = Gnd_temp
    t_gnd[temp>temp_idx] = temp[temp>temp_idx] + tempdiff
    t_barn = temp + 3.0
    # t_gnd = temp + 0.6*(20-temp)
    
    blocking_factor1 = 0.2
    blocking_factor2 = 0.8
    u_barn[temp<=temp_idx] = (1-blocking_factor1)*wind[temp<=temp_idx]
    u_barn[temp>temp_idx] = (1-blocking_factor2)*wind[temp>temp_idx]
    return t_barn, t_gnd, u_barn

####################
## Reaction rates
####################
## EMPIRICAL - rate: uric acid hydrolysis to TAN; temp in degC, rhum in per cent
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
## EMPIRICAL - rate: urea hydrolysis to TAN; temp in degC, delta_t is the time step, i.e 1 hour, daily res: delta_t=24
# def urea_hydrolysis_rate(temp,WFPS,delta_t):
#     ## T is the soil temperature (or air temperature for housing)
#     ## WFPS is water-filled porosity, i.e. in this simulation, is the moisture content of manure for housing
#     ## it remains unclear of k_h and WFPS
#     ## k_h = 1.2 * WFPS
def urea_hydrolysis_rate(temp,theta,delta_t,k_h=0.23): 
    ## DEFAULT:k_h equals 0.23 per hour at 20 degC; Goh and Sherlock, 1985; Muck, 1981;
    ## Ah_t is a temeprature scaling factor for k_h; temperature dependence Q10 is ~2. 
    Ah_t = 0.25 * np.exp(0.0693*temp)
    hydrolysis_rate = 1 - np.exp((-k_h*theta*delta_t)*Ah_t)  
    hydrolysis_rate[np.isnan(hydrolysis_rate)] = 0.0
    return hydrolysis_rate
## EMPIRICAL - rate: TAN production from the decompostion of N_avail and N_resist; temp in degC
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
## EMPIRICAL - rate: nitrification rate in soils
## ref: Stange&Neue, 2009; Riddick at al, FANv1, 2016BG; Vira et al., FANv2, 2020GMD
## ground temp in degC; theta in m3/m3, theta_sat in m3/m3
def nitrification_rate_soil(ground_temp,theta,theta_sat,pH,fer_type):
    ## convert degC to K
    ground_temp = ground_temp + 273.15
    ## maximum  temp for microbial activity
    tmax = 313   
    if fer_type == "manure":
        ## optimum temp for microbial activity
        topt = 301
        ## empirical factor in this parameterization
        asigma = 2.4
        ## critical water content of soil (g/g)
        #mcrit = 0.12
        ## sharp paramter of moisture response function
        ##b = 2
    elif fer_type == "mineral":
        ## optimum temp for microbial activity
        topt = 303
        ## empirical factor in this parameterization
        asigma = 1.8
        ## critical water content of soil (g/g)
        #mcrit = 0.10
        ## sharp paramter of moisture response function
        #b = 16
    ## temperature responce funtion
    func_tg = ((tmax-ground_temp)/(tmax-topt))**asigma*np.exp(asigma*(ground_temp-topt)/(tmax-topt))
    ## water in soil
    # soilwc = theta*rho_water/((1-theta_sat)*rho_soil)
    ## moisture response function
    # func_soilwc = 1 - np.exp(-(soilwc/mcrit)**b)
    ## soil WFPS
    WFPS = theta/theta_sat
    ## WFPS response function; 
    a = 0.55
    b = 1.7
    c = -0.007
    d = 3.22
    func_wfps = (((WFPS-b)/(a-b))**(d*(b-a)/(a-c)))*(((WFPS-c)/(a-c))**d)
    ## soil pH responce function
    func_pH = 0.56+(np.arctan(np.pi*0.45*(-5+pH))/np.pi)
    ## maximum rate of nitrification; s^-1
    rmax = 1.16e-6
    # rmax = 1.16e-6/2   ## test rmax
    ## nitrification rate; per second
    # nitrif_rate = (2*rmax*func_tg*func_wfps)/(func_tg+func_wfps)
    nitrif_rate = rmax*func_tg*func_wfps*func_pH
    return nitrif_rate
## EMPIRICAL - rate: nitrification rate of manure (analogy to nitrification rate in soils)
## ref: Stange&Neue, 2009; Riddick at al, FANv1, 2016BG; Vira et al., FANv2, 2020GMD
## manure temp in degC; WFPS in fraction
def nitrification_rate_manure(manure_temp,WFPS):
    ## convert degC to K
    manure_temp = manure_temp + 273.15
    ## maximum and optimum temp for microbial activity
    tmax = 313   
    topt = 301
    ## empirical factor in this parameterization
    asigma = 2.4
    ## temperature response funtion
    func_tg = ((tmax-manure_temp)/(tmax-topt))**asigma*np.exp(asigma*(manure_temp-topt)/(tmax-topt))
    ## WFPS response function; a 4th order polynomial regression
    # a1 = 26.34
    # a2 = -54.51
    # a3 = 32.61
    # a4 = -4.31
    # a5 = 0.15
    # func_wfps = a1*(WFPS**4)+a2*(WFPS**3)+a3*(WFPS**2)+a4*WFPS+a5
    a = 0.55
    b = 1.7
    c = -0.007
    d = 3.22
    func_wfps = (((WFPS-b)/(a-b))**(d*(b-a)/(a-c)))*(((WFPS-c)/(a-c))**d)
    ## maximum rate of nitrification; s^-1
    rmax = 1.16e-6
    # rmax = 1.16e-6/2   ## test rmax
    ## nitrification rate; per second
    # nitrif_rate = (2*rmax*func_tg*func_wfps)/(func_tg+func_wfps)
    nitrif_rate = rmax*func_tg*func_wfps
    return nitrif_rate

#####################
## Concentration
#####################
## calculate TAN contentration
def TAN_concentration(mtan,zlayer,theta_sat,theta,knh3,kd):
    cnc = mtan/(zlayer*(theta+knh3*(theta_sat-theta)+(1-theta_sat)*kd))
    cnc[theta==0] = 0.0
    return cnc
## calculate NH3 concentration
def NH3_concentration(tan_cnc,knh3,theta_sat,theta):
    cnc = tan_cnc * knh3
    cnc[theta==theta_sat] = 0.0
    return cnc
## calculate N species concentration; urea, NO3-
def N_concentration(mN,zlayer,theta):
    cnc = mN/(zlayer*theta)
    cnc[theta==0.0] = 0.0
    return cnc
## calculate surface compensation point for TAN
def surf_TAN_cnc(tan_cnc,rliq,rgas,knh3,ratm,qrunoff):
    TAN_surf_cnc = (tan_cnc*(1/rliq+knh3/rgas)/(qrunoff+knh3*(1/ratm+1/rgas)+1/rliq))
    return TAN_surf_cnc
## calculate surface compensation point for N species, i.e., urea, NO3-
def surf_Ncnc(N_cnc,rliq,qrunoff):
    N_surf_cnc = N_cnc/(rliq*qrunoff+1)
    return N_surf_cnc

######################################
## physical/chemical coefficients
######################################
## calculate dissociation constant of NH4+, kNH4 (to be put in PARAMETERS.py)
## temp in degC
def NH4_dissoc_coeff(temp):
    knh4 = 5.67e-10*np.exp(-6286*(1/(temp + 273.15)-1/298.15))
    return knh4 
## calculated NH3 partitioning coefficient, KNH3 (to be put in PARAMETERS.py)
## temp in degC; H+ ions concentration
def NH3_par_coeff(temp,cncH):
    ## dimensonessHenry's Law number
    henry_constant = (161500/(temp + 273.15)) * np.exp(-10380/(temp + 273.15))
    ## dissociation constant of NH4+
    k_NH4 = NH4_dissoc_coeff(temp)
    knh3 = henry_constant/(cncH + k_NH4)
    return knh3
## calculate NH4+ fraction
def frac_NH4(theta,theta_sat,temp,cncH,kd):
    k_NH4 = NH4_dissoc_coeff(temp)
    knh3 = NH3_par_coeff(temp,cncH)
    f_NH4 = theta/(theta+knh3*(theta_sat-theta)+(1-theta_sat)*kd)*(cncH/(cncH+k_NH4))
    return f_NH4
## calculate molecular diffusivity of NH4+ in water (inc. urea); D_aq_nh4, m^2/s 
## ref : Van Der Molen et al.,1990; Vira et al., GMD2020
## temp in degC
def diffusivity_NH4(temp,phase):
    if phase == 'aqueous':
        d_nh4 = 9.8e-10*1.03**temp
    elif phase == 'gaseous':
        temp = temp + 273.15
        d_nh4 = (1e-7*(temp)**1.75*(1/M_air+1/M_NH3)**
                 0.5)/(1.0*(sigma_v_air**(1/3)+sigma_v_NH3**(1/3))**2)
    return d_nh4

#######################
## Fluxes/Pathways
#######################
## NH3 release from the slat
def NH3_volslat(slat_conc,conc_in,Rslat):
    NH3flux = (slat_conc-conc_in)/Rslat
    return NH3flux
## NH3 release from the pit
def NH3_volpit(pit_tanconc,conc_in,Rpit):
    NH3flux = (pit_tanconc-conc_in)/Rpit
    return NH3flux
## calculate flux: NH3 volatilization (g/m2/s)
## tan_surfcnc in g/m3, ratm in s/m
def NH3_vol(nh3_surfcnc,ratm,nh3_atmcnc=0.0):
    nh3vol = (nh3_surfcnc-nh3_atmcnc)/ratm
    return nh3vol
## calculate flux: surface runoffs (g/m2/s)
## N_surfcnc in g/m3, qrunoff in m/s
def surf_runoff(N_surfcnc,qrunoff):
    surfrunoff = N_surfcnc*qrunoff
    return surfrunoff
## calculate flux: infiltrtion/leaching (g/m2/s)
## N_cnc in g/m3, qsubrunoff in m/s
def subsurf_leaching(N_cnc,qsubrunoff):
    subsrfleaching = N_cnc*qsubrunoff
    return subsrfleaching
## calculate flux: TAN aq diffusion (g/m2/s)
## cnc in g/m3, resist in s/m  
def N_diffusion(cnc1,cnc2,resist):
    diffusion = (cnc1-cnc2)/resist
    diffusion[diffusion<0.0] = 0.0
    return diffusion
## calculate chemical transformation: nitrification (g/m2/s)
def TAN_nitrif(tan_pool,temp,theta,theta_sat,pH,fert_type,frac_nh4):
    nitrif_rate = nitrification_rate_soil(temp,theta,theta_sat,pH,fert_type)
    nitrif_rate[nitrif_rate>0.1] = 0.1
    ## correction for WPFS response
    nitrif_rate[nitrif_rate<0.0] = 0.0
    # nitrif_rate[np.isnan(nitrif_rate)] = 0.0
    tan_nitrif = tan_pool*nitrif_rate*frac_nh4
    tan_nitrif[np.isnan(tan_nitrif)] = 0.0
    return tan_nitrif
## calculate plant N uptake rate (to be removed from PARAMETER.py and to be put in LAND.py)
## Ammonium and Nitrate N in g/m2; soil C in gC
def plant_N_uptake(mNH4,mNO3,temp,uptake,substrateC=0.04,substrateN=0.004):
    ## root activity weighting parameters
    v1,v2,v3,v4 = 1.0,0.5,0.25,0.1
    ## root structural dry matter components; g/m2
    Wr1,Wr2,Wr3,Wr4 = 20,40,60,80
    ## root activity paramter; gN/(g root day)
    sigmaN20 = 0.05
    ## root activity parameters; [C], [N], gN/m2
    K_C, J_N, K_Neff = 0.05, 0.005, 5
    ## temperature response
    func_temp = 0.25*np.exp(0.0693*temp)
    func_temp[func_temp>1.0] = 1.0
    ## soil mineral concentration; 
    Neff = mNH4 + func_temp*mNO3
    ## non-linear relationship between N uptake and effective soil mineral concentration
    NC_factor = 1/(1+K_C/substrateC*(1+substrateN/J_N))
    ## gN/m2/s
    if uptake == 'nh4':
        uptake = (sigmaN20*func_temp*(v1*Wr1+v2*Wr2+v3*Wr3+v4*Wr4)*(mNH4/(Neff+K_Neff))*NC_factor)/(24*3600)
    elif uptake == 'no3':
        uptake = (sigmaN20*func_temp*(v1*Wr1+v2*Wr2+v3*Wr3+v4*Wr4)*(func_temp*mNO3/(Neff+K_Neff))*NC_factor)/(24*3600)
    return uptake

## flux1: [NH3 volatalization] or [NH3 upwards diffusion]
## flux2: [TAN washoff] or [TAN upwards diffusion]
## flux3: [TAN infiltration/leaching]
## flux4: [TAN downwards diffusion]
## flux5: [NH3 downwards diffusion]
## flux6: [ammonium uptake]
def TAN_pathways(mN,flux1=False,flux2=False,flux3=False,flux4=False,flux5=False,flux6=False):
    if flux1 is not False:
        flux1[np.isnan(flux1)] = 0.0
        flux1[np.isinf(flux1)] = 0.0
    if flux2 is not False:
        flux2[np.isnan(flux2)] = 0.0
        flux2[np.isinf(flux2)] = 0.0
    if flux3 is not False:
        flux3[np.isnan(flux3)] = 0.0
        flux3[np.isinf(flux3)] = 0.0
    if flux4 is not False:
        flux4[np.isnan(flux4)] = 0.0
        flux4[np.isinf(flux4)] = 0.0
    if flux5 is not False:
        flux5[np.isnan(flux5)] = 0.0
        flux5[np.isinf(flux5)] = 0.0
    if flux6 is not False:
        flux6[np.isnan(flux6)] = 0.0
        flux6[np.isinf(flux6)] = 0.0
    
    totalidx = flux1+flux2+flux3+flux4+flux5+flux6
    massidx = mN-totalidx
    if flux1 is not False:
        flux1[massidx<0] = (flux1[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux2 is not False:
        flux2[massidx<0] = (flux2[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux3 is not False:
        flux3[massidx<0] = (flux3[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux4 is not False:
        flux4[massidx<0] = (flux4[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux5 is not False:
        flux5[massidx<0] = (flux5[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux6 is not False:
        flux6[massidx<0] = (flux6[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    return flux1,flux2,flux3,flux4,flux5,flux6

## flux1: [urea/NO3 washoff] or [NO3 uptake]
## flux2: [urea/NO3 upwards diffusion]
## flux3: [urea/NO3 infiltration/leaching]
## flux4: [urea/NO3 downwards diffusion]
def N_pathways(mN,flux1=False,flux2=False,flux3=False,flux4=False):
    if flux1 is not False:
        flux1[np.isnan(flux1)] = 0.0
        flux1[np.isinf(flux1)] = 0.0
    if flux2 is not False:
        flux2[np.isnan(flux2)] = 0.0
        flux2[np.isinf(flux2)] = 0.0
    if flux3 is not False:
        flux3[np.isnan(flux3)] = 0.0
        flux3[np.isinf(flux3)] = 0.0
    if flux4 is not False:
        flux4[np.isnan(flux4)] = 0.0
        flux4[np.isinf(flux4)] = 0.0

    totalidx = flux1+flux2+flux3+flux4
    massidx = mN-totalidx
    if flux1 is not False:
        flux1[massidx<0] = (flux1[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux2 is not False:
        flux2[massidx<0] = (flux2[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux3 is not False:
        flux3[massidx<0] = (flux3[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    if flux4 is not False:
        flux4[massidx<0] = (flux4[massidx<0]/totalidx[massidx<0])*mN[massidx<0]
    return flux1,flux2,flux3,flux4

############################
## Other physical/chemical 
############################
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
    B = (rho_air/rho_water)*(0.622*k**2*u)/(Pressure*(np.log(Z/zo))**2)
    evap = B*(es-ea)*1000
    # evap = evap*1000*3600*24
    return evap
## EMPIRICAL - mass tranfer coefficient for NH4 in liquid boundary layer
## temp in deg C
def k_aq_NH4(temp):
    temp = temp + 273.15
    kl = 1.417*1e-12*(temp**4)
    return kl
## EMPIRICAL - mass transfer coefficient for NH3 in gas boudnary layer
def k_gas_NH3(temp,u,Z,zo):
    k = 0.41
    Sc = Sc_num(temp)
    ustar = k*u/(np.log(Z/zo))
    kg = 0.001+0.0462*ustar*(Sc**(-0.67))
    return kg
## calculate resistance for diffusion of aq,gas TAN
## for NO3-, a correction must be applied (after calling this function)
def diff_resistance(distance,phase,theta_sat,theta,temp):
    f_soil = soil_tortuosity(theta_sat,theta,phase)
    diff_val = diffusivity_NH4(temp,phase)
    rdiff = distance/(f_soil*diff_val)
    return rdiff
## resistance: resistance for water-air exchange; temp in degC, rhum in per cent; evap flux in m/s
def resistance_water_air(temp,rhum,evap_flux):
    T = temp + 273.15
    ## to avoid dividing zero
    rhum = rhum-0.00001
    ## evap flux in m/s (1000 kg/m^2/s; m/s)
    # evap_flux = evap_flux/(1000*24*3600)
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
    # D_air_NH3 = (1e-7*(T)**1.75*(1/M_air+1/M_NH3)**
    #              0.5)/(Pressure*(sigma_v_air**(1/3)+sigma_v_NH3**(1/3))**2)
    # D_air_H2O = (1e-7*(T)**1.75*(1/M_air+1/M_NH3)**
    #              0.5)/(Pressure*(sigma_v_air**(1/3)+sigma_v_H2O**(1/3))**2)
    Rc = (rho_air/rho_water)*((Q_sat-Q_atm)/evap_flux)
    return Rc
## resistance: resistance for aerodynamic and boundary layer resistance; 
## temp in degC; rhum in %; u in m/s; H is sensible heat flux in J/(m^2 s); Z is reference height in m; zo is surface roughness in m 
def resistance_aero_boundary(temp,rhum,u,H,Z,zo):
    temp = np.array(temp)
    rhum = np.array(rhum)
    u = np.array(u)
    H = np.array(H)
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
    rho_air = p/(R*T_v)
    ## friction velocity
    ustar = k*u/(np.log(Z/zo))
    ## Monin-Obukhov length
    L = -T*(ustar**3)*rho_air*cpair/(k*g*H)
    ## stability correction function: psi
    psi = np.zeros(u.shape)
    ## stable condition
    psi[H<=0] = -5*Z/L[H<=0]  
    ## unstable condition
    X = np.zeros(u.shape)
    X[H>0] = (1-16*Z/L[H>0])**0.25
    psi[H>0] = np.log(((1+X[H>0])/2)**2)+np.log((1+X[H>0]**2)/2)-2*np.arctan(X[H>0])+np.pi/2
    ## aerodynamic resistance
    Ra = (np.log(Z/zo)-psi)**2/(k**2*u)
    
    B = 5
    Rb = 1/(B*ustar)
    return Ra+Rb
## physical: wind profile
## calculating mean wind speed at a specific height by knowing the wind speed at a reference height
## uref in m/s; height_ref and height_out in m; zo in m
def wind_profile(uref,height_ref,height_out,zo):
    uout = uref*(np.log(height_out/zo)/np.log(height_ref/zo))
    return uout   
## Reynold's number
def Re_num(temp,u,Z,zo):
    k=0.41
    ## friction velocity m/s
    ustar = k*u/(np.log(Z/zo))
    ## kinematic viscosity of air; m2/s
    v = 1.56e-5*((temp+273.15)/298.15)**(3/2)
    Re = ustar*zo/v
    return Re
## Schmidt number
def Sc_num(temp):
    ## diffusivity of NH3 in the air
    DNH3 = diffusivity_NH4(temp,phase="gaseous")
    ## kinematic viscosity of air; m2/s
    v = 1.56e-5*((temp+273.15)/298.15)**(3/2)
    Sc = v/DNH3
    return Sc

#######################################
## Soil properties/characteristics
#######################################
## soil characteristics: tortuosity for diffusion
## theta is the volumetric soil water content, and theta_sat is the volumetric soil water content at saturation (equivalent as porosity)
## theta in m3/m3
def soil_tortuosity(theta_sat,theta,phase):
    ## soil tuotorsity in aqueous phase and gaeous phase (Millington and Quirk, 1961)
    ## implementing parameter tuning 
    F_correct = 0.85
    if phase == 'aqueous':
        soil_tor = (((theta)**(10/3))/(theta_sat**2))**F_correct
    elif phase == 'gaseous':
        soil_tor = (((theta_sat-theta)**(10/3))/(theta_sat**2))*F_correct
    return soil_tor
## soil characteristics: infiltration rate (m/s) - empirical method 
## this is an empirically-derived expression for vertical/percolation/infiltration/subsurface leaching/ of animal slurry
## data source: Patle et al., Geo.Eco.Landscale 2019
## loamy sand and sandy loam in 14 sites (quality control; remove bad quality data)
## multilinear regression: R^2 = 0.75; 
## soil para: sand (%), clay (%), bulk density (g/cm^3), particle density (g/cm^3)
## percentage of saturation soil moisture 
def infiltration_rate_empirical(soilmoist_percent,sand_percent,clay_percent,bulk_density,particle_density):
    ## soil parameters matrix
    soil_para = 0.62*sand_percent - 0.17*clay_percent - 26.42*bulk_density + 3.14*particle_density - 10.58
    ## soil moisture dependence
    infil_func = soil_para + 6.31*(soilmoist_percent/100)
    ## Ks is the permeability coefficient at saturation for loamy sand, Ks=0.714cm/h (0.0119cm/min;171.36mm/day) ref: Hu et al., 2017 J.Arid.Land
    Ks = 0.714
    ## infiltration rate Ki m/s
    Ki = (infil_func/Ks)/(100*3600)
    return Ki
## soil characteristics: infiltration rate (m/s) - conceptual method
## In Sommer&Jacobsen 1999, 3kg/m2 slurry was applied to loamdy sand (9.5%clay;11%silt;77%sand;BD~1.5g/cm3)
## and infiltrate within 24h, which is equivalent to 3mm/24h;
## dailyevap in m/day
def infiltration_rate_method(dailyinfil,theta_sat,theta):
    ## 10 mm of water/slurry; ref 5mm for 12h infiltration in Vira et al., FANv2 2020 GMD
    ## infiltration to a depth of 4mm; ref: source layer thickness used in Moring et al., GAG model, 2016 BG
    z_layer = 0.004
    ## infiltration rate; m/s
    Ki = (dailyinfil - z_layer*(theta_sat-theta))/(24*3600)
    return Ki
## EMPIRICAL - soil characteristics: soil pH after manure/fertilizer application
## soil pH will increase due to application of urea, peak within 24-48 h, then decrease
## urea decomposition consumes H+ ion, which leads to pH increase
## [app_timing_app] is an indexing map that remarks the timing of fertilizer application with  a shape of [time,lat,lon]
def soil_pH_postapp(base_pH,app_timing_map,fert_pH):
    pH_postapp = np.zeros(mtrx2)
    for tt in np.arange(base_pH.shape[0]-6):
        pH_postapp[tt+1][app_timing_map[tt]==1] = fert_pH
        pH_postapp[tt+2][app_timing_map[tt]==1] = fert_pH
        pH_postapp[tt+3][app_timing_map[tt]==1] = fert_pH - (1/5)*(fert_pH - base_pH[tt][app_timing_map[tt]==1])
        pH_postapp[tt+4][app_timing_map[tt]==1] = fert_pH - (2/5)*(fert_pH - base_pH[tt][app_timing_map[tt]==1])
        pH_postapp[tt+5][app_timing_map[tt]==1] = fert_pH - (3/5)*(fert_pH - base_pH[tt][app_timing_map[tt]==1])
        pH_postapp[tt+6][app_timing_map[tt]==1] = fert_pH - (4/5)*(fert_pH - base_pH[tt][app_timing_map[tt]==1])
    pH_postapp[pH_postapp<base_pH] = base_pH[pH_postapp<base_pH]
    return pH_postapp
## EMPIRICAL - soil characteristics: adsorption coefficient of NH4+ on soil particles
## clay_content in %
def ammonium_adsorption(clay_content):
    clay = clay_content/100
    kd = 0.5*(7.2733*clay**3-11.22*clay**2) + 5.7198*clay + 0.0263
    return kd

####################################################
## Livestock waste properties/characteristics
####################################################
## manure characteristics: manure 1) volumetric water content 2) porosity 3) water-filled pore space (WFPS)
## solid mass in g/m2; water mass in g/m2; rho_manure in g/cm3
def manure_properties(solidmass,watermass):
    ## total volume of manure; m3/m2
    vtotal = solidmass/(manure_BD*1e3)
    ## gravimetric water content of manure
    theta_g = watermass/solidmass
    theta_g[np.isnan(theta_g)] = 0.0
    ## volumetric water content
    theta_v = theta_g*manure_BD/rho_water
    ## WFPS: volumetric water content/porosity
    WFPS = theta_v/manure_porosity
    return vtotal,theta_v,WFPS
## manure minimum moisture
def min_manurewc(temp,rhum):
    mois_coeff = (-np.log(1.01-(rhum/100))/(0.0000534*(temp+273)))**(1/1.41)
    return mois_coeff/100
## aniaml info: waste N should be consistent to livestock_N
## unit in kg N per head per year; returning daily values
def livestock_waste_info(livestock_type, waste_N):
    ## daily N excretion from urine and feces
    dN = 1000*waste_N/(365*24)
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
rho_m = defaultdict(dict)
pH_info = defaultdict(dict)
name = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
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
name = ['BEEF_CATTLE','DAIRY_CATTLE','OTHER_CATTLE','PIG','MARKET_SWINE','BREEDING_SWINE','SHEEP','GOAT','POULTRY','BUFFALO']
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
## the dry matter (DM) content of solid manure (assumed to be equivalent to manure fertilizer)； ref (Boyd, CAB reviews, 2018)
f_solid_DM = [31.4, 24.1, 31.4, 30.8, 30.8, 30.8, 32.2, 32.2, 60.6, 31.4]
## assuming the density of manure; 1t/m^3 or 1g/cm^3 for cattle, pigs etc, 0.4 for poultry
rho_manure = [1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.4, 1.0]
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
        rho_m[name[ii]] = rho_manure[ii]
        pH_info[name[ii]] = pH_val[ii]


##############################
## Manure Management Systems
##############################
## provisional MMS categories
## loss (untraceable): fishpond, discahrge, public sewage
loss_list = ['fishpond','discharge','publsewage','dumping']
## sold (untraceable): sold
sold_list = ['sold']
## Cat I.D: no significant emission; used as fuel: biogas(digester,liquid), burned (solid)
MMS_fuel_list = ['mmsbiogas','mmsburned']
## Cat I.C.1: N mostly preserved for further use (as solid): thermal drying 
MMS_preserve_solid_list = ['mmsthermal']
## Cat I.C.2: N mostly preserved for further use (as liquid), e.g., manure stored in liquid phase with emission mitigation measures: lagoon (typically with cover), liquid crust
MMS_preserve_liquid_list = ['mmsliqcrust']
## Cat I.A.1: manure stored in barns (as solid): composting, deep litter, litter (poultry), no litter (poultry), pit (layer),solid storage (inc. ot)
MMS_indoor_solid_list = ['mmscompost','mmssolid','mmsolidot']
##  Cat I.A.2: manure stored in barns (as liquid): aerobic processing, liquid, pit1, pit2
MMS_indoor_liquid_list = ['mmslagoon']
## Cat I.B.1: manure in open/outdoor environment; left on land (as solid): aerobic processing, daily spreading, dry lot, pasture, pasture+paddock
MMS_outdoor_solid_list = ['mmsconfin','mmsdrylot','mmspasture','mmspastpad']
## Cat I.B.2: manure in open/ooutdoor environment (as liquid): aerobic lagoon, liquid
MMS_outdoor_liquid_list = ['mmsaerproc','mmsliquid','mmsliqoth']
## Cat I.B.3: manure in open environment (as liquid): aerobic lagoon, liquid
MMS_outdoor_lagoon_list = ['mmsaerobic']

## MMS for housing module (in-situ storage of manure in animal houses)
## Cat II.A.1: in-situ storage of manure (as solid)
MMS_house_storage_solid_list = ['mmsdeeplitt','mmslitter','mmsnolitt'] 
## Cat II.A.2: in-situ storgate of manure (as liquid)
MMS_house_storage_liquid_list = ['mmspit1','mmspit2']

## MMS for land module
## Cat III
MMS_land_spread_list = ['mmsdaily']

############################
## model variables
############################
## wind speed at 10m
wind_data_height = 10
## reference height is 2m
ref_height = 2
## density of air: 1.2754 kg/m3
rho_air = 1.2754
## density of water: 997 kg/m3
rho_water = 997
## density of soil particle: 2660 kg/m3; 2.66 g/cm3
rho_soil = 2660
## density of manure solids: 1460 kg/m3; 1.46g/cm3
manure_PD = 1460
## manure porosity; ref: Khater et al, 2015 cattle manure porosity
manure_porosity = 0.42
## manure bulk density
manure_BD = manure_PD*(1-manure_porosity)
## NO3- diffusivity scaling factor
f_DNO3 = 1.3e-8/9.8e-10
## molar mass of air
M_air = 28.96
## molar mass of NH3
M_NH3 = 17
## moloar mass of H2O
M_H2O = 18
## sigma_v is the value of atomic diffusion volume in [lagoon system] (from Liley et al., 1984;
## Physical and chemical data)
sigma_v_air = 20.1
sigma_v_NH3 = 14.9
sigma_v_H2O = 12.7

