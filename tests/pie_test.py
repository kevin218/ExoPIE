#!/usr/bin/python

import pytest
import numpy as np
import matplotlib.pyplot as plt
import pickle, sys, os
import dynesty
from exopie import transforms as tr

# Load PandExo data
try:
    handle  = open('./demos/pandexo/niriss_soss-W43-R100.p', 'rb')
except:
    handle  = open('../demos/pandexo/niriss_soss-W43-R100.p', 'rb')
model   = pickle.load(handle)

# Define spectrum and SNR
wave = model['FinalSpectrum']['wave']
spectrum = model['FinalSpectrum']['spectrum']
error = model['FinalSpectrum']['error_w_floor']
randspec = model['FinalSpectrum']['spectrum_w_rand']
snr = 1./error

# Run non-transiting case
Tprng   = 1000
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
savedir = 'test-non'
foo = tr.run_nontransiting(wave, snr, Tprng, Ts, Rp, Rs, pmin=[400,4350,0,0.664], pmax=[1600,4450,2.0,0.67], savefig=savedir)
niriss_means_tr, niriss_errors_tr, niriss_quantiles_tr = foo
plt.close('all')

good = np.abs(niriss_means_tr-np.array([Tprng, Ts, Rp, Rs])) < 3*niriss_errors_tr
for isgood in good.squeeze():
    assert isgood
os.system(f"rm -r " + savedir)


# Load PandExo data
try:
    handle  = open('./demos/pandexo/nirspec_g395h-W43-r100.p', 'rb')
except:
    handle  = open('../demos/pandexo/nirspec_g395h-W43-r100.p', 'rb')
model_g395  = pickle.load(handle)

# Define spectrum and SNR
wave_g395 = model_g395['FinalSpectrum']['wave']
spectrum_g395 = model_g395['FinalSpectrum']['spectrum']
error_g395 = model_g395['FinalSpectrum']['error_w_floor']
randspec_g395 = model_g395['FinalSpectrum']['spectrum_w_rand']
snr_g395 = 1./error_g395

# Run transiting case
Tprng   = np.arange(900,1000,100)
Ts      = 4400
Rp      = 1.036
Rs      = 0.667
savedir = 'test-tr'
foo = tr.run_transiting(wave_g395, snr_g395, Tprng, Ts, Rp, Rs, pmin=[400,4300,0.9], pmax=[1600,4500,1.1], savefig=savedir)
g395_means_tr, g395_errors_tr, g395_quantiles_tr = foo
plt.close('all')

good = np.abs(g395_means_tr-np.array([Tprng[0], Ts, Rp])) < 3*g395_errors_tr
for isgood in good.squeeze():
    assert isgood
os.system(f"rm -r " + savedir)
