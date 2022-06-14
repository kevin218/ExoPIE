#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pickle, sys, os
import dynesty
from exopie import transforms as tr

Rsun    = 69.634e7  #meters
Rjup    = 6.9911e7  #meters

#Load PandExo data
handle  = open('pandexo/niriss_soss-W43-R100.p', 'rb')
model   = pickle.load(handle)
obstime = model['RawData']['electrons_out'][0]/model['RawData']['e_rate_out'][0]/3600
print(obstime)
#1.0455 hours

wave = model['FinalSpectrum']['wave']
spectrum = model['FinalSpectrum']['spectrum']
error = model['FinalSpectrum']['error_w_floor']
randspec = model['FinalSpectrum']['spectrum_w_rand']
#snr = model['RawData']['electrons_out']/np.sqrt(model['RawData']['var_out'])
snr = 1./error

'''
# WASP-43
radS    = mast.planck(wave*1e-6, 4400)
areaS   = np.pi*(0.667*Rsun)**2
intS    = radS*areaS
# WASP-43b emission
Tprng   = np.arange(400,2501,100)
radP    = mast.planck(wave*1e-6, Tprng[:,np.newaxis])
areaP   = np.pi*(1.036*Rjup)**2
intP    = radP*areaP
trdepth = areaP/areaS
'''
#ecldepth    = fluxP/fluxS*trdepth
#ecldepth    = intP/intS
#bberror     = intS/snr
#bberror     = error*fluxS/ecldepth
#ndim        = 3
#unc         = bberror
Tprng   = np.arange(400,3001,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667

#foo = tr.run_transiting(wave, snr, Tprng, Ts, Rp, Rs, sigma=[300,10,0.01], savefig='figs-tr-soss')
foo = tr.run_transiting(wave, snr, Tprng, Ts, Rp, Rs, pmin=[100,4300,0.1], pmax=[3500,4500,2.0], savefig='figs-tr-soss')
niriss_means_tr, niriss_errors_tr, niriss_quantiles_tr = foo
np.savez("tr-soss.npz", means=niriss_means_tr, errors=niriss_errors_tr, quantiles=niriss_quantiles_tr)

"""
plt.figure(1)
plt.clf()
plt.plot(Tprng, errors[:,0], 'o')

plt.figure(2)
plt.clf()
plt.plot(Tprng, errors[:,1], 'o')

plt.figure(3)
plt.clf()
plt.plot(Tprng, errors[:,2], 'o')
"""


Tprng   = np.arange(400,3001,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
foo = tr.run_nontransiting(wave, snr, Tprng, Ts, Rp, Rs, pmin=[100,4300,0.1,0.66], pmax=[3500,4500,2.0,0.675], savefig='figs-non-soss')
niriss_means_non, niriss_errors_non, niriss_quantiles_non = foo
np.savez("non-soss.npz", means=niriss_means_non, errors=niriss_errors_non, quantiles=niriss_quantiles_non)
plt.close('all')



handle  = open('pandexo/nirspec_g395h-W43-r100.p', 'rb')
model_g395  = pickle.load(handle)
print(model['RawData']['electrons_out'][0]/model['RawData']['e_rate_out'][0]/3600)

wave_g395 = model_g395['FinalSpectrum']['wave']
spectrum_g395 = model_g395['FinalSpectrum']['spectrum']
error_g395 = model_g395['FinalSpectrum']['error_w_floor']
randspec_g395 = model_g395['FinalSpectrum']['spectrum_w_rand']
#snr_g395 = model_g395['RawData']['electrons_out']/np.sqrt(model_g395['RawData']['var_out'])
snr_g395 = 1./error_g395

wave2   = np.concatenate((wave, wave_g395))
snr2    = np.concatenate(( snr,  snr_g395))

Tprng   = np.arange(400,3001,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
foo = tr.run_nontransiting(wave2, snr2, Tprng, Ts, Rp, Rs, pmin=[100,4300,0,0.66], pmax=[3500,4500,2.0,0.67], savefig='figs-non-soss+g395')
niriss_g395_means_non, niriss_g395_errors_non, niriss_g395_quantiles_non = foo
np.savez("non-soss+g395.npz", means=niriss_g395_means_non, errors=niriss_g395_errors_non, quantiles=niriss_g395_quantiles_non)
plt.close('all')


Tprng   = np.arange(400,3001,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
foo = tr.run_transiting(wave2, snr2, Tprng, Ts, Rp, Rs, pmin=[100,4300,0.9], pmax=[3500,4500,1.1], savefig='figs-tr-soss+g395')
niriss_g395_means_tr, niriss_g395_errors_tr, niriss_g395_quantiles_tr = foo
np.savez("tr-soss+g395.npz", means=niriss_g395_means_tr, errors=niriss_g395_errors_tr, quantiles=niriss_g395_quantiles_tr)
plt.close('all')

'''
#Load data
'''

tupper_tr_niriss        = niriss_quantiles_tr[:,0,2] - niriss_means_tr[:,0]
tlower_tr_niriss        = niriss_means_tr[:,0] - niriss_quantiles_tr[:,0,0]
terr_tr_niriss          = niriss_errors_tr[:,0]
Rperr_tr_niriss         = niriss_errors_tr[:,2]
tupper_tr_niriss_g395   = niriss_g395_quantiles_tr[:,0,2] - niriss_g395_means_tr[:,0]
tlower_tr_niriss_g395   = niriss_g395_means_tr[:,0] - niriss_g395_quantiles_tr[:,0,0]
terr_tr_niriss_g395     = niriss_g395_errors_tr[:,0]
Rperr_tr_niriss_g395    = niriss_g395_errors_tr[:,2]
tupper_non_niriss       = niriss_quantiles_non[:,0,2] - niriss_means_non[:,0]
tlower_non_niriss       = niriss_means_non[:,0] - niriss_quantiles_non[:,0,0]
terr_non_niriss         = niriss_errors_non[:,0]
Rperr_non_niriss        = niriss_errors_non[:,2]
tupper_non_niriss_g395  = niriss_g395_quantiles_non[:,0,2] - niriss_g395_means_non[:,0]
tlower_non_niriss_g395  = niriss_g395_means_non[:,0] - niriss_g395_quantiles_non[:,0,0]
terr_non_niriss_g395    = niriss_g395_errors_non[:,0]
Rperr_non_niriss_g395   = niriss_g395_errors_non[:,2]

plt.figure(1)
plt.clf()
lw  = 2
#plt.fill_between(Tprng, tlower_tr_niriss, tupper_tr_niriss, alpha=0.6, label='NIRISS, Transiting')
#plt.fill_between(Tprng, tlower_tr_niriss_g395, tupper_tr_niriss_g395, alpha=0.6, label='NIRISS+G395H, Transiting')
#plt.fill_between(Tprng, tlower_non_niriss, tupper_non_niriss, alpha=0.6, label='NIRISS, Non-Transiting')
#plt.fill_between(Tprng, tlower_non_niriss_g395 , tupper_non_niriss_g395 , alpha=0.6, label='NIRISS+G395H, Non-Transiting')
plt.plot(Tprng, terr_tr_niriss, '-', color='green', lw=lw, label='Transiting, 0.8-2.8 $\mu m$')
plt.plot(Tprng, terr_tr_niriss_g395, '-', color='purple', lw=lw, label='Transiting, 0.8-5.0 $\mu m$')
plt.plot(Tprng, terr_non_niriss, '--', color='green', lw=lw, label='Non-Transiting, 0.8-2.8 $\mu m$')
plt.plot(Tprng, terr_non_niriss_g395, '--', color='purple', lw=lw, label='Non-Transiting, 0.8-5.0 $\mu m$')
plt.grid(True, alpha=0.3, ls=':')
plt.xlim(500,2500)
plt.ylim(0,100)
plt.xlabel('Simulated Planet Temperature (K)', size=11)
plt.ylabel("Retrieved Planet Temperature Uncertainty (K)", size=11)
plt.legend(loc="upper right", fontsize=11)
plt.tight_layout()

plt.savefig("Tp-Constraint.png", dpi=200)



plt.figure(2)
plt.clf()
lw  = 2
#plt.plot(Tprng, Rperr_tr_niriss, '-', color='C0', lw=lw, label='NIRISS, Transiting')
#plt.plot(Tprng, Rperr_tr_niriss_g395, '-', color='C1', lw=lw, label='NIRISS+G395H, Transiting')
plt.plot(Tprng, Rperr_non_niriss, '--', color='green', lw=lw, label='Non-Transiting, 0.8-2.8 $\mu m$')
plt.plot(Tprng, Rperr_non_niriss_g395, '--', color='purple', lw=lw, label='Non-Transiting, 0.8-5.0 $\mu m$')
plt.grid(True, alpha=0.3, ls=':')
plt.xlim(500,2500)
plt.ylim(0,0.2)
plt.yticks(np.linspace(0,0.2,11),np.linspace(0,0.2,11))
plt.xlabel('Simulated Planet Temperature (K)', size=11)
plt.ylabel("Retrieved Planet Radius Uncertainty ($R_J$)", size=11)
plt.legend(loc="upper right", fontsize=11)
plt.tight_layout()

plt.savefig("Rp-Constraint.png", dpi=200)
