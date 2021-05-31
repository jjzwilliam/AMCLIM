## import AMCLIM modules
#from CONFIG.config import *
#from INPUT.input import *
#from MODULES.PARAMETERS import *
from MODULES.HOUSING import *
#from POSTPROC.postprocessing import *

#sys.exit()

## test MODULE.HOUSING
month_start=7
print(Months_name[month_start-1])
housing_init()
slat_pit_housing_sim(Months_idx[month_start-1],Months_idx[month_start],f_slat=0.4,f_gap=0.6)

slat_emiss = np.nansum(NH3_flux_slat,axis=0)*farming_area
pit_emiss = np.nansum(NH3_flux_pit,axis=0)*farming_area
total_emiss = slat_emiss+pit_emiss


# total_emiss = np.nansum(NH3_flux_from_barn,axis=0)*farming_area

total_N_Jan = excretN_info*1000*31/365
Pv = total_emiss/total_N_Jan*100

print('NH3 emission from slat is '+str(np.nansum(slat_emiss)/1e9)+ 'Gg')
print('NH3 emission from pit is '+str(np.nansum(pit_emiss)/1e9)+' Gg')
print('Total NH3 emission is '+str(np.nansum(total_emiss)/1e9)+' Gg')
print('Monthly excreted N is '+str(np.nansum(total_N_Jan)/1e9)+' Gg')
print('Average Pv is '+str(np.nansum(total_emiss)/np.nansum(total_N_Jan)*100)+' %')

print('farming area: '+str(np.nansum(farming_area)))
print('R slat: '+str(np.nanmean(R_star_slat)))
print('R pit: '+str(np.nanmean(R_star_pit)))
print('urea hydro rate: '+str(np.nanmean(daily_urea_hydro_rate)))

print('urea: '+str(np.nanmean(durea)))
print('urine: '+str(np.nanmean(durine)))

# plot_monthlysim(total_emiss,Pv,month_start)
