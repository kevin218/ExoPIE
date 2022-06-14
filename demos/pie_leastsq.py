#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as spo
import scipy.stats as sps
import pickle, sys
from exoMAST import ObsTable as mast
import dynesty
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc

Rsun    = 69.634e7  #meters
Rjup    = 6.9911e7  #meters

#Load PandExo data
handle  = open('niriss_soss-W43-R100.p', 'rb')
#handle  = open('niriss_soss-W43.p', 'rb')
model   = pickle.load(handle)

wave = model['FinalSpectrum']['wave']
spectrum = model['FinalSpectrum']['spectrum']
error = model['FinalSpectrum']['error_w_floor']
randspec = model['FinalSpectrum']['spectrum_w_rand']
snr = model['RawData']['electrons_out']/np.sqrt(model['RawData']['var_out'])

# WASP-43
#waveset = np.logspace(-1, 2, num=1000) * u.um
#fluxS   = np.array(blackbody_lambda(wave*u.um, 4400*u.K))
#fluxS   = (blackbody_lambda(wave*u.um, 4400*u.K)*np.pi*(0.667*Rsun)**2).to(u.W / u.sr / u.m)
radS    = mast.planck(wave*1e-6, 4400)
areaS   = np.pi*(0.667*Rsun)**2
intS    = radS*areaS
# WASP-43b emission
Tprng   = np.arange(400,2501,100)
radP    = mast.planck(wave*1e-6, Tprng[:,np.newaxis])
areaP   = np.pi*(1.036*Rjup)**2
intP    = radP*areaP
#for i in range(len(Tprng)):
    #fluxP.append((blackbody_lambda(wave*u.um, Tprng[i]*u.K)*np.pi*(1.036*Rjup)**2).to(u.W / u.sr / u.m))
#fluxP   = np.array(fluxP)
trdepth = areaP/areaS

#ecldepth    = fluxP/fluxS*trdepth
ecldepth    = intP/intS
bberror     = intS/snr
#bberror     = error*fluxS/ecldepth
'''
plt.figure(2)
plt.clf()
plt.errorbar(wave, ecldepth[2]*1e6, yerr=error*1e6, fmt='o', zorder=0)
plt.plot(wave, ecldepth[2]*1e6, '-', zorder=1)
'''

plt.figure(1, figsize=(6.4, 3.2))
plt.clf()
ax = plt.subplot()
#ax.plot(wave, fluxS+fluxP[5], '--', color='C1')
ax.errorbar(wave, intS+intP[5], bberror, fmt='.', zorder=0)
ax.plot(wave, intS, color='C1', zorder=1)
ax.set_xlabel(r'Wavelength ($\mu m$)', size=11, labelpad=0)
ax.set_ylabel(r'Spectral Intensity (W/(sr m))', size=11, labelpad=0)
#ax.set_xscale('log')
ax.set_yscale('log')
plt.tight_layout()

#plt.savefig('Emission-NIR.png', dpi=200)

# Generate simulated data
ii      = 3
unc     = bberror
data    = intS+intP[ii]
#data    = np.random.normal(intS+intP[ii], unc)

# Fit star + planet to simulated data
#planetStarEmission
def model_func(params, wave, trdepth, data, unc):
    Tp, Ts, Rp = params
    Bp      = mast.planck(wave*1e-6, Tp)
    Bs      = mast.planck(wave*1e-6, Ts)
    areaP   = np.pi*(Rp*Rjup)**2
    #Rs      = np.sqrt(areaP/trdepth/np.pi)
    areaS   = areaP/trdepth
    model   = Bs*areaS + Bp*areaP
    return (data - model)/unc

params = (Tprng[ii], 4500, 1.0)
output, cov_x, infodict, mesg, err = spo.leastsq(model_func, params, args=(wave,trdepth,data,unc), full_output=True)
print(mesg)
print(output)
print(np.sqrt(np.diagonal(cov_x)))

#allparams, bestop, numaccept, numit = demc.demc()

'''
#arr=[wlgrid, y_meas, err, y_mod_best, wlgrid_hires,y_median, y_low_2sig, y_low_1sig, y_high_1sig, y_high_2sig]
n_spec  = 15
model   = []
for i in range(n_spec):
    handle  = open('./Retrievals3/PHASE'+str(i)+'_WFC3_IRAC_spectra.pic', 'rb')
    model.append(pickle.load(handle, encoding="latin1"))
nmodel  = len(model)

#Combine into lists
wlgrid  = model[0][4]
hstwave = model[0][0]
smodel  = []
hstflux = []
hstferr = []
binmodel= []
modm1   = []
modp1   = []
for i in range(nmodel):
    smodel.append(smooth.smooth(model[i][5],15)*1e3)    #Smooth model spectra in ppt
    modm1.append(smooth.smooth(model[i][7],15))
    modp1.append(smooth.smooth(model[i][8],15))
    hstflux.append(model[i][1]*1e3)
    hstferr.append(model[i][2]*1e3)
    binmodel.append(model[i][3])
hstflux = np.array(hstflux)
hstferr = np.array(hstferr)
binmodel= np.array(binmodel)
'''

"""
binsize     = 0.0625  #In orbits
binphase    = np.arange(binsize/2.,1-binsize/2.,binsize)+binsize/2.

xmin    = 0.6
xmax    = 5.0

plt.figure(1, figsize=(6.4, 3.2))
plt.clf()
ax = plt.subplot()
#ax.plot(waveset.value, fluxP/fluxS*depth*1e6, color='k')
for i in [0,1,13,14]:
    ax.plot(wlgrid, smodel[i], '-', label=binphase[i])
ax.set_xlabel(r'Wavelength ($\mu m$)', size=10, labelpad=0)
ax.set_ylabel(r'Planet Emission (ppm)', size=10, labelpad=0)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(1, 1e3)
ax.set_xticks([0.6,0.8,1,1.3,2,3,4,5])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.legend(loc='upper left')
plt.subplots_adjust(top=0.97, bottom=0.14, left=0.09, right=0.98)

plt.savefig("Emisison-WASP43-4Phases.png", dpi=200)




# WASP-43
waveset = np.logspace(-1, 2, num=1000) * u.um
fluxS   = blackbody_lambda(waveset, 4400*u.K)
# WASP-43b emission
Tprng   = np.arange(500,1301,100)
fluxP   = []
for i in range(len(Tprng)):
    fluxP.append(blackbody_lambda(waveset, Tprng[i]*u.K))
fluxP   = np.array(fluxP)
depth   = 0.025371

ymin    = 10
ymax    = 1e4
xmin    = 0.6
xmax    = 5.0

plt.figure(2, figsize=(6.4, 3.2))
plt.clf()
ax = plt.subplot()
#ax.plot(waveset.value, fluxP/fluxS*depth*1e6, color='k')
for i in np.arange(len(Tprng)-1,-1,-1):
    ax.plot(waveset.value, fluxP[i]/fluxS*depth*1e6, '-', label=str(Tprng[i])+' K')
ax.set_xlabel(r'Wavelength ($\mu m$)', size=11, labelpad=0)
ax.set_ylabel(r'Planet-Star Flux Ratio (ppm)', size=11, labelpad=0)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks([0.6,0.8,1,1.3,2,3,4,5])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.axhline(20, color='k', linestyle='dashed')
plt.legend(loc='upper left')
plt.subplots_adjust(top=0.97, bottom=0.14, left=0.09, right=0.98)

plt.savefig('Emission-NIR.png', dpi=200)





# Mid-IR range

# TRAPPIST-1
waveset = np.logspace(-1, 2, num=1000) * u.um
fluxS   = blackbody_lambda(waveset, 2500*u.K)
Tprng   = np.arange(250,451,50)
fluxP   = []
for i in range(len(Tprng)):
    fluxP.append(blackbody_lambda(waveset, Tprng[i]*u.K))
fluxP   = np.array(fluxP)
depth   = 0.004931  #T1-e


ymin    = 10
ymax    = 1e3
xmin    = 5.0
xmax    = 12.0

plt.figure(3, figsize=(6.4, 3.2))
plt.clf()
ax = plt.subplot()
#ax.plot(waveset.value, fluxP/fluxS*depth*1e6, color='k')
for i in np.arange(len(Tprng)-1,-1,-1):
    ax.plot(waveset.value, fluxP[i]/fluxS*depth*1e6, '-', label=str(Tprng[i])+' K')
ax.set_xlabel(r'Wavelength ($\mu m$)', size=11, labelpad=0)
ax.set_ylabel(r'Planet-Star Flux Ratio (ppm)', size=11, labelpad=0)
ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks([5,6,7,8,9,10,12])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.axhline(20, color='k', linestyle='dashed')
plt.legend(loc='upper left')
plt.subplots_adjust(top=0.97, bottom=0.14, left=0.09, right=0.98)

plt.savefig('Emission-MIR.png', dpi=200)
"""
