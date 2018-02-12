# -*- coding: utf-8 -*-
''' accpy.math.specialfuns
author:     felix.kramer(at)physik.hu-berlin.de
'''
from numpy import exp, sqrt, pi, sin, cos, linspace, logspace, log2


def gauss(x, amp, mu, sigma, offset):
    ''' gauss(x, a, mu, sigma, offset) =
        a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset
    x      - input variable
    amp    - amplitude (parameter)
    mu     - expected value (parameter)
    sigma  - standard deviation (parameter)
    offset - (parameter)

    infos:
         full width at half maximum
         FWHM = 2*sqrt(2*ln(2))*sigma
         full width at tenth of maximum
         FWTM = 2*sqrt(2*ln(10))*sigma
    '''
    return amp * exp( -(((x-mu)/sigma)**2)/2 ) / sigma/sqrt(2*pi) + offset


def gauss2D((x, y), amp, (mux, muy), (sigx, sigy), offset, theta=0):
    ''' gauss2D(x, a, mu, sigma, offset) =
        amp*exp(-(((x-x0)/sigma)**2)/2)/sigma/sqrt(2*pi) + offset
    x, y       - input variables
    amp        - amplitude (parameter)
    mux, muy   - expected value (parameter)
    sigx, sigy - standard deviation (parameter)
    theta      - rotation (parameter)
    offset     - (parameter)

    infos:
         full width at half maximum
         FWHM = 2*sqrt(2*ln(2))*sigma
         full width at tenth of maximum
         FWTM = 2*sqrt(2*ln(10))*sigma
    '''
    if theta != 0:
        sigx2, sigy2 = sigx**2, sigy**2
        sinth2, costh2 = sin(theta)**2, cos(theta)**2
        a = (costh2/sigx2 + sinth2/sigy2)/2
        b = sin(2*theta)/4 * (1/sigy2 - 1/sigx2)
        c = (sinth2/sigx2 + costh2/sigy2)/2
        return amp * exp( -a*(x-mux)**2 - 2*b*(x-mux)*(y-muy) - c*(y-muy)**2 ) /sigx/sigy/(2*pi) + offset
    else:
        return amp * exp( -(((x-mux)/sigx)**2 + ((y-muy)/sigy)**2)/2 ) /sigx/sigy/(2*pi) + offset


def hanning(N):
    return (cos(linspace(-pi, pi, N))+1)/2


def squarespace(xmin, xmax, n):
    xmax *= 1 + 1e-1
    out = xmax - logspace(log2(xmin + 1e-1*xmax), log2(xmax), n, base=2)
    return out[::-1]
