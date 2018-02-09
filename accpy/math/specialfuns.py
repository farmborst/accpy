# -*- coding: utf-8 -*-
''' accpy.math.specialfuns
author:     felix.kramer(at)physik.hu-berlin.de
'''
from numpy import exp, sqrt, pi, cos, linspace, logspace, log2


def gauss(x, a, mu, sigma, offset):
    ''' gauss(x, a, mu, sigma, offset) = a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset
    x      - input variable
    a      - (parameter)
    mu     - expected value (parameter)
    sigma  - standard deviation (parameter)
    offset - (parameter)

    infos:
         full width at half maximum
         FWHM = 2*sqrt(2*ln(2))*sigma
         full width at tenth of maximum
         FWTM = 2*sqrt(2*ln(10))*sigma
    '''
    return a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset


def hanning(N):
    return (cos(linspace(-pi, pi, N))+1)/2


def squarespace(xmin, xmax, n):
    xmax *= 1 + 1e-1
    out = xmax - logspace(log2(xmin + 1e-1*xmax), log2(xmax), n, base=2)
    return out[::-1]
