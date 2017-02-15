# -*- coding: utf-8 -*-
''' accpy.math.specialfuns
author:     felix.kramer(at)physik.hu-berlin.de
'''
from numpy import exp, sqrt, pi, cos, linspace


def gauss(x, a, mu, sigma, offset):
    return a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset


def hanning(N):
    return (cos(linspace(-pi, pi, N))+1)/2