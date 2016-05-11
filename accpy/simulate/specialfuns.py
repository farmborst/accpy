# -*- coding: utf-8 -*-
"""
Felix Kramer
accpy.simulate.specialfuns
"""
from numpy import exp, sqrt, pi


def gauss(x, a, mu, sigma, offset):
    return a*exp(-(((x-mu)/sigma)**2)/2)/sigma/sqrt(2*pi)+offset
