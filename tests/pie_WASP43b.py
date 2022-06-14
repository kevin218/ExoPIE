#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import pickle, sys
import dynesty
from exopie import transforms as tr

Rsun    = 69.634e7  #meters
Rjup    = 6.9911e7  #meters


handle  = open('pandexo/nirspec_g395h-W43-r100.p', 'rb')
model_g395  = pickle.load(handle)
# print(model['RawData']['electrons_out'][0]/model['RawData']['e_rate_out'][0]/3600)

wave_g395 = model_g395['FinalSpectrum']['wave']
spectrum_g395 = model_g395['FinalSpectrum']['spectrum']
error_g395 = model_g395['FinalSpectrum']['error_w_floor']
randspec_g395 = model_g395['FinalSpectrum']['spectrum_w_rand']
#snr_g395 = model_g395['RawData']['electrons_out']/np.sqrt(model_g395['RawData']['var_out'])
snr_g395 = 1./error_g395

Tprng   = 800
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
foo = tr.run_nontransiting(wave_g395, snr_g395, Tprng, Ts, Rp, Rs, pmin=[400,4300,0,0.66], pmax=[1200,4500,2.0,0.67], savefig='test-non')
g395_means_non, g395_errors_non, g395_quantiles_non = foo
# np.savez("non-soss+g395.npz", means=g395_means_non, errors=g395_errors_non, quantiles=g395_quantiles_non)
plt.close('all')


Tprng   = np.arange(700,801,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
foo = tr.run_transiting(wave_g395, snr_g395, Tprng, Ts, Rp, Rs, pmin=[400,4300,0.9], pmax=[1200,4500,1.1], savefig='test-tr')
g395_means_tr, g395_errors_tr, g395_quantiles_tr = foo
# np.savez("tr-soss+g395.npz", means=g395_means_tr, errors=g395_errors_tr, quantiles=g395_quantiles_tr)
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
