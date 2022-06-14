
import numpy as np
import scipy.optimize as spo
import scipy.stats as sps
import dynesty, sys, os
from dynesty import plotting as dyplot
from dynesty import utils as dyfunc
from exoMAST_Obs import exoMAST_TableSNR as mast

Rsun    = 69.634e7  #meters
Rjup    = 6.9911e7  #meters


def loglike_transiting(params, wave, data, unc, trdepth):
    """The log-likelihood function for transiting planets with known depth."""
    Tp, Ts, Rp = params
    Bp      = mast.planck(wave*1e-6, Tp)
    Bs      = mast.planck(wave*1e-6, Ts)
    areaP   = np.pi*(Rp*Rjup)**2
    areaS   = areaP/trdepth
    model   = Bs*areaS + Bp*areaP

    return -0.5 * (np.sum(((data-model)/unc)**2 + np.log(2*np.pi*unc**2)))

def loglike_nontransiting(params, wave, data, unc):
    """The log-likelihood function for non-transiting planets."""
    Tp, Ts, Rp, Rs = params
    Bp      = mast.planck(wave*1e-6, Tp)
    Bs      = mast.planck(wave*1e-6, Ts)
    areaP   = np.pi*(Rp*Rjup)**2
    areaS   = np.pi*(Rs*Rsun)**2
    model   = Bs*areaS + Bp*areaP

    return -0.5 * (np.sum(((data-model)/unc)**2 + np.log(2*np.pi*unc**2)))

# Define our uniform prior.
def ptform_uniform(u, pmin, pmax):
    """Transforms the uniform random variables `u ~ Unif[0., 1.)`
    to the parameters of interest."""

    #pmin    = np.array([ 800,4390,0.9])
    #pmax    = np.array([1200,4410,1.1])
    pmin    = np.array(pmin)
    pmax    = np.array(pmax)
    diff    = pmax - pmin
    return u*diff + pmin

# Define our Gaussian prior.
def ptform_gaussian(u, mu, s):
    """Transforms the uniform random variables `u ~ Unif[0., 1.)`
    to the parameters of interest."""

    x   = np.copy(u)
    # Normal
    mu  = np.array(mu)
    s   = np.array(s)
    #mu  = np.array([ Tp, 4400, 1.04])
    #s   = np.array([300,   10, 0.01])
    x   = sps.norm.ppf(u, loc=mu, scale=s)

    return x

def run_transiting(wave, snr, Tprng, Ts, Rp, Rs, sigma=None, pmin=None, pmax=None, savefig='.'):
    """

    """
    if sigma != None:
        ndim    = len(sigma)
        isgauss = True
    elif pmin != None:
        ndim    = len(pmin)
        isgauss = False
    else:
        print("Input error")
        return
    radS    = mast.planck(wave*1e-6, Ts)
    areaS   = np.pi*(Rs*Rsun)**2
    intS    = radS*areaS
    radP    = mast.planck(wave*1e-6, Tprng[:,np.newaxis])
    areaP   = np.pi*(Rp*Rjup)**2
    intP    = radP*areaP
    trdepth = areaP/areaS
    unc     = intS/snr

    ntemp       = len(Tprng)
    means       = np.zeros((ntemp,ndim))
    errors      = np.zeros((ntemp,ndim))
    quantiles   = np.zeros((ntemp,ndim,3))
    for ii in range(ntemp):
        data    = intS+intP[ii]
        mu      = [Tprng[ii],Ts,Rp]
        if isgauss:
            sampler = dynesty.NestedSampler(loglike_transiting, ptform_gaussian, ndim,
                                        nlive=500, bound='multi', sample='unif',
                                        logl_args=(wave, data, unc, trdepth),
                                        ptform_args=(mu,sigma))
        else:
            sampler = dynesty.NestedSampler(loglike_transiting, ptform_uniform, ndim,
                                        nlive=500, bound='multi', sample='unif',
                                        logl_args=(wave, data, unc, trdepth),
                                        ptform_args=(pmin,pmax))
        #sampler = dynesty.DynamicNestedSampler(loglike, ptform_gaussian, ndim,
        #                                bound='multi', sample='unif')
        #sampler = dynesty.DynamicNestedSampler(loglike, ptform_uniform, ndim,
        #                                       bound='balls', sample='rwalk')
        #print(sampler.citations)
        sampler.run_nested()
        results = sampler.results
        results.summary()

        #quantiles=[0.159, 0.5, 0.841]
        cfig, caxes = dyplot.cornerplot(results, color='b', show_titles=True, truths=mu, title_fmt='.4f')
        cfig.savefig(savefig+"/Corner-"+str(Tprng[ii])+"K.png",dpi=200)
        tfig, taxes = dyplot.traceplot(results)
        tfig.savefig(savefig+"/Trace-"+str(Tprng[ii])+"K.png",dpi=200)
        rfig, raxes = dyplot.runplot(results)
        rfig.savefig(savefig+"/Run-"+str(Tprng[ii])+"K.png",dpi=200)

        # Extract sampling results.
        samples = results.samples  # samples
        weights = np.exp(results.logwt - results.logz[-1])  # normalized weights
        # Compute 16%, 84%, 50% quantiles.
        quantiles[ii] = np.array([dyfunc.quantile(samps, [0.159, 0.5, 0.841], weights=weights)
                                 for samps in samples.T])
        # Compute weighted mean and covariance.
        means[ii], cov = dyfunc.mean_and_cov(samples, weights)
        errors[ii]   = np.sqrt(np.diagonal(cov))
    return means, errors, quantiles

def run_nontransiting(wave, snr, Tprng, Ts, Rp, Rs, sigma=None, pmin=None, pmax=None, savefig='.'):
    """

    """
    if sigma != None:
        ndim    = len(sigma)
        isgauss = True
    elif pmin != None:
        ndim    = len(pmin)
        isgauss = False
    else:
        print("Input error")
        return
    radS    = mast.planck(wave*1e-6, Ts)
    areaS   = np.pi*(Rs*Rsun)**2
    intS    = radS*areaS
    radP    = mast.planck(wave*1e-6, Tprng[:,np.newaxis])
    areaP   = np.pi*(Rp*Rjup)**2
    intP    = radP*areaP
    trdepth = areaP/areaS
    unc     = intS/snr

    ntemp       = len(Tprng)
    means       = np.zeros((ntemp,ndim))
    errors      = np.zeros((ntemp,ndim))
    quantiles   = np.zeros((ntemp,ndim,3))
    for ii in range(ntemp):
        data    = intS+intP[ii]
        mu      = [Tprng[ii],Ts,Rp,Rs]
        if isgauss:
            sampler = dynesty.NestedSampler(loglike_nontransiting, ptform_gaussian, ndim,
                                        nlive=500, bound='multi', sample='unif',
                                        logl_args=(wave, data, unc),
                                        ptform_args=(mu,sigma))
        else:
            sampler = dynesty.NestedSampler(loglike_nontransiting, ptform_uniform, ndim,
                                        nlive=500, bound='multi', sample='unif',
                                        logl_args=(wave, data, unc),
                                        ptform_args=(pmin,pmax))
        #print(sampler.citations)
        sampler.run_nested()
        results = sampler.results
        results.summary()

        # Create directory
        if os.path.isdir(savefig) == False:
            os.mkdir(savefig, 511)

        #quantiles=[0.159, 0.5, 0.841]
        cfig, caxes = dyplot.cornerplot(results, color='b', show_titles=True, truths=mu, title_fmt='.4f')
        cfig.savefig(savefig+"/Corner-"+str(Tprng[ii])+"K.png",dpi=200)
        tfig, taxes = dyplot.traceplot(results)
        tfig.savefig(savefig+"/Trace-"+str(Tprng[ii])+"K.png",dpi=200)
        rfig, raxes = dyplot.runplot(results)
        rfig.savefig(savefig+"/Run-"+str(Tprng[ii])+"K.png",dpi=200)

        # Extract sampling results.
        samples = results.samples  # samples
        weights = np.exp(results.logwt - results.logz[-1])  # normalized weights
        # Compute 10%-90% quantiles.
        quantiles[ii] = np.array([dyfunc.quantile(samps, [0.159, 0.5, 0.841], weights=weights)
                                 for samps in samples.T])
        # Compute weighted mean and covariance.
        means[ii], cov = dyfunc.mean_and_cov(samples, weights)
        errors[ii]   = np.sqrt(np.diagonal(cov))
    return means, errors, quantiles
