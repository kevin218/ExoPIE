
import warnings
warnings.filterwarnings('ignore')
import pandexo.engine.justdoit as jdi
import numpy as np
import os,sys
from exoMAST import ObsTable as mast

exo_dict = jdi.load_exo_dict('WASP-43 b')

exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
exo_dict['observation']['sat_unit'] = '%'
exo_dict['observation']['noccultations'] = 1 #number of transits
exo_dict['observation']['R'] = 100          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                             #without having to redo the calcualtion
#exo_dict['observation']['baseline_unit'] = 'total'  #Defines how you specify out of transit observing time
                                                    #'frac' : fraction of time in transit versus out = in/out
                                                    #'total' : total observing time (seconds)
#exo_dict['observation']['baseline'] = 4.0*60.0*60.0 #in accordance with what was specified above (total observing time)
exo_dict['observation']['baseline_unit'] = 'frac'
exo_dict['observation']['baseline'] = 1.0 #in accordance with what was specified above (total observing time)
exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                             #to a wavelength dependent noise floor solution (units are ppm)

#exo_dict['planet']['type'] ='user'                       #tells pandexo you are uploading your own spectrum
#exo_dict['planet']['exopath'] = 'wasp12b.txt'
#exo_dict['planet']['w_unit'] = 'cm'                      #other options include "um","nm" ,"Angs", "sec" (for phase curves)
#exo_dict['planet']['f_unit'] = 'rp^2/r*^2'               #other options are 'fp/f*'

exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
#exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit

#ORRRRR....
#if you wanted to instead to secondary transit at constant temperature
exo_dict['planet']['f_unit'] = 'fp/f*'
exo_dict['planet']['temp'] = 1000
#exo_dict['planet']['transit_duration'] = 3600
#exo_dict['planet']['td_unit'] = 'second'
#exo_dict['planet']['type'] = 'grid'                #tells pandexo you want to pull from the grid
#exo_dict['planet']['temp'] = 1000                 #grid: 500, 750, 1000, 1250, 1500, 1750, 2000, 2250, 2500
#exo_dict['planet']['chem'] = 'noTiO'              #options: 'noTiO' and 'eqchem', noTiO is chemical eq. without TiO
#exo_dict['planet']['cloud'] = 'ray10'               #options: nothing: '0',

#jdi.print_instruments()
'''
Choose from the following:
dict_keys(['WFC3 G141', 'MIRI LRS', 'NIRISS SOSS', 'NIRSpec G140M', 'NIRSpec G140H', 'NIRSpec G235M', 'NIRSpec G235H', 'NIRSpec G395M', 'NIRSpec G395H', 'NIRSpec Prism', 'NIRCam F322W2', 'NIRCam F444W'])
'''
#result = jdi.run_pandexo(exo_dict,['NIRSpec Prism'],output_file='nirspec_prism.p')
result = jdi.run_pandexo(exo_dict,['NIRISS SOSS'],output_file='niriss_soss-W43-r100.p')
result = jdi.run_pandexo(exo_dict,['NIRSpec G395H'],output_file='nirspec_g395h-W43-r100.p')
result = jdi.run_pandexo(exo_dict,['MIRI LRS'],output_file='miri_lrs-W43-r100.p')

wave = result['FinalSpectrum']['wave']
spectrum = result['FinalSpectrum']['spectrum']
error = result['FinalSpectrum']['error_w_floor']
randspec = result['FinalSpectrum']['spectrum_w_rand']

plt.figure(1)
plt.clf()
plt.errorbar(wave, spectrum, error, fmt='o')


# Prox Cen b

exo_dict = jdi.load_exo_dict('Proxima Cen b')
exo_dict['observation']['sat_level'] = 80    #saturation level in percent of full well
exo_dict['observation']['sat_unit'] = '%'
exo_dict['observation']['noccultations'] = 1 #number of transits
exo_dict['observation']['R'] = 100          #fixed binning. I usually suggest ZERO binning.. you can always bin later
                                             #without having to redo the calcualtion
exo_dict['observation']['baseline_unit'] = 'frac'
exo_dict['observation']['baseline'] = 1.0/3 #in accordance with what was specified above (total observing time)
exo_dict['observation']['noise_floor'] = 0   #this can be a fixed level or it can be a filepath
                                             #to a wavelength dependent noise floor solution (units are ppm)
exo_dict['planet']['type'] = 'constant'                  #tells pandexo you want a fixed transit depth
#exo_dict['planet']['f_unit'] = 'rp^2/r*^2'        #this is what you would do for primary transit
exo_dict['planet']['f_unit'] = 'fp/f*'
exo_dict['planet']['temp'] = 250
exo_dict['planet']['transit_duration'] = 3600
exo_dict['planet']['td_unit'] = 'second'
#exo_dict['planet']['transit_duration'] = 1/24.
#exo_dict['planet']['td_unit'] = 'd'
exo_dict['star']['metal'] = 0
exo_dict['star']['logg'] = mast.planetLogg(0.12*1048, 0.141*9.96)

result = jdi.run_pandexo(exo_dict,['MIRI LRS'],output_file='miri_lrs-ProxCen-r100-1hrs.p')


wave = result['FinalSpectrum']['wave']
spectrum = result['FinalSpectrum']['spectrum']
error = result['FinalSpectrum']['error_w_floor']
randspec = result['FinalSpectrum']['spectrum_w_rand']

plt.figure(1)
plt.clf()
plt.errorbar(wave, spectrum*1e6, error*1e6, fmt='o')
