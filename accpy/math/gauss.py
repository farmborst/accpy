# -*- coding: utf-8 -*-
''' accpy.math.gauss
author:     felix.kramer(at)physik.hu-berlin.de
'''
from bisect import bisect
from numpy import (exp, array, ndarray, pi, sqrt, ones_like, mean, argmax, log,
                   trapz, sum, shape, sort)


def Gauss1D(pars, args):
    '''
    x      - independant variable (argument)
    amp    - amplitude (parameter)
    mu     - expected value (parameter)
    sigma  - standard deviation (parameter)
    offset - background height (parameter)
    '''
    x = args
    if isinstance(pars, dict):
        amp, mu, sigma, offset = pars['amp'], pars['mu'], pars['sig'], pars['off']
    elif isinstance(pars, (list, ndarray)):
        amp, mu, sigma, offset = pars
    return amp * exp( -(((x-mu)/sigma)**2)/2 ) / sigma/sqrt(2*pi) + offset


def NGauss1D(Npars, args):
    l = [Gauss1D(pars, args) for pars in Npars]
    y = sum(l, axis=0)
    return y


def DGauss1D(pars, args):
    x, _ = args
    amp, mu, sig, off = pars
    x_sub_mu = x - mu
    df_by_damp = exp( -((x_sub_mu / sig)**2) / 2 ) / sqrt(2*pi) / sig
    tmp = amp * df_by_damp / sig**2
    df_by_dmu = tmp * x_sub_mu
    df_by_dsig = tmp * (x_sub_mu - sig) * (x_sub_mu + sig) / sig
    df_by_doff = ones_like(x)
    return array([df_by_damp, df_by_dmu, df_by_dsig, df_by_doff])


def GuessGauss1D(args, ydata):
    '''
    infos:
         full width at half maximum
         FWHM = 2*sqrt(2*ln(2))*sigma
         full width at tenth of maximum
         FWTM = 2*sqrt(2*ln(10))*sigma
    '''
    x = args
    i_max = argmax(ydata)
    mu = x[i_max]
#    tmp = sort(ydata[:i_max])
#    offset = (tmp[:int(i_max/100)] + tmp[-int(i_max/100):])/2
    offset = mean(ydata[:int(i_max/100)])
    print(offset)
    i_max = argmax(ydata)
    hm = (max(ydata[:i_max]) + offset)/2
    i_hm = bisect(ydata[:], hm)
    fwhm = 2*abs(x[i_hm] - mu) # 1st val > hm
    sigma = fwhm/(2*sqrt(2*log(2)))
    amp = trapz(ydata-offset, x)
    print([amp, mu, sigma, offset])
    return [amp, mu, sigma, offset]


def Gauss2D(pars, args):
    '''
    x, y       - input arguments
    parameters as dict or list
        amp        - amplitude
        mux, muy   - expected value
        sigx, sigy - standard deviation
        offset     - constant background

    infos:
         full width at half maximum
         FWHM = 2*sqrt(2*ln(2))*sigma
         full width at tenth of maximum
         FWTM = 2*sqrt(2*ln(10))*sigma
    '''
    x, y = args
    if isinstance(pars, dict):
        amp = pars['amp']
        mux, muy = pars['mux'], pars['muy']
        sigx, sigy = pars['sigx'], pars['sigy']
        offset = pars['off']
    elif isinstance(pars, (list, ndarray)):
        amp, mux, muy, sigx, sigy, offset = pars
    return amp * exp( -(((x-mux)/sigx)**2 + ((y-muy)/sigy)**2)/2 ) /sigx/sigy/(2*pi) + offset


def DGauss2D(pars, args):
    x, _ = args
    amp, mu, sig, off = pars
    x_sub_mu = x - mu
    df_by_damp = exp( -((x_sub_mu / sig)**2) / 2 ) / sqrt(2*pi) / sig
    tmp = amp * df_by_damp / sig**2
    df_by_dmu = tmp * x_sub_mu
    df_by_dsig = tmp * (x_sub_mu - sig) * (x_sub_mu + sig) / sig
    df_by_doff = ones_like(x)
    return array([df_by_damp, df_by_dmu, df_by_dsig, df_by_doff])


def NGauss2D(Npars, args):
    if isinstance(Npars, dict):
        subparams = [{key : Npars[key][i] for key in Npars} for i in range(len(Npars['amp']))]
        l = [Gauss2D(pars, args) for pars in subparams]
        z = sum(l, axis=0)
    elif isinstance(Npars, (list, ndarray)):
        l = [Gauss2D(pars, args) for pars in Npars]
        z = sum(l, axis=0)
    return z