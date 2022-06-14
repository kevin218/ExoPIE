
import numpy as np
import sys
sys.path.append('/Users/stevekb1/Documents/code/exoMAST/exoMAST_Obs')
import exoMAST_TableSNR as mast

Rsun    = 69.634e7  #meters
Rjup    = 6.9911e7  #meters

Ts      = 4400
Rs      = 0.667
Tp      = 1000
Rp      = 1.036

wave    = np.logspace(-1,2, num=1000)
radS    = mast.planck(wave*1e-6, Ts)
areaS   = np.pi*(Rs*Rsun)**2
intS    = radS*areaS
radP    = mast.planck(wave*1e-6, Tp)
areaP   = np.pi*(Rp*Rjup)**2
intP    = radP*areaP
trdepth = areaP/areaS
#unc     = intS/snr

ymin    = 1e28
ymax    = 1e31
xmin    = 0.6
xmax    = 5

plt.figure(3, figsize=(5, 3.2))
plt.clf()
ax = plt.subplot()
ax.plot(wave, intS, color='C0', label='Star Only')
#ax.plot(waveset.value, fluxP.to(u.W /u.m**2/u.um/u.sr).value, color='C0')
ax.plot(wave, intS+1000*intP, '--', color='C1', label=r'Star + Planet Nightside$\times$1000')
ax.set_xlabel(r'Wavelength ($\mu m$)', size=12, labelpad=0)
ax.set_ylabel(r'Spectral Intensity (W/sr/m)', size=12, labelpad=0)
ax.loglog()
ax.set_xlim(xmin, xmax)
ax.set_ylim(ymin, ymax)
ax.set_xticks([0.6,0.8,1,1.3,2,3,4,5])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.get_xaxis().set_minor_formatter(matplotlib.ticker.NullFormatter())
plt.yticks(visible=False)
#plt.xticks([0.6,0.8,1,1.3,2,3,4,5],['0.6','0.8',1,1.3,2,3,4,5])
plt.legend(loc='lower left', fontsize=12)
plt.text(0.85,5e29, "Reference\nWavelength\nRange", color='k', fontsize=12, verticalalignment='center', horizontalalignment='center')
plt.fill_betweenx([ymin,ymax], xmin, 1.3, facecolor='C2', alpha=0.1)
plt.axvline(1.3, ls='dotted', color='C2')
plt.text(3.5,4e30, "Planetary\nInfrared\nExcess", color='k', fontsize=12, verticalalignment='center', horizontalalignment='center')
plt.arrow(3.5, 1.8e30, 0, -1.4e30, width= 0.15, head_length=1.5e29, fc='C1', ec='k')
plt.subplots_adjust(top=0.97, bottom=0.14, left=0.06, right=0.98)
#plt.tight_layout()
plt.savefig("Blackbody-Nightside-1000x.png", dpi=200)
#plt.savefig("Blackbody-WASP43b-Nightside-1000x.png", dpi=200)
