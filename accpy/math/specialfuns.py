# -*- coding: utf-8 -*-
''' accpy.math.specialfuns
author:     felix.kramer(at)physik.hu-berlin.de
'''
from numpy import exp, sqrt, pi


def gauss(x, a, mu, sigma, offset):
    return a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset
