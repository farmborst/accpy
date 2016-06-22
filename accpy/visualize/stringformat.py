# -*- coding: utf-8 -*-
"""
Created on Thu Jun 16 19:25:00 2016

@author: user
"""
from __future__ import division
from math import log10


def SI(d, significantfigures=4):
    # Unit prefixes
    prefixes = ['a', 'f', 'p', 'n', r'$\mu$', 'm', '',
                'k', 'M', 'G', 'T', 'P', 'E']
    exp10 = log10(abs(d))
    i = int((exp10) // 3)
    prefix = prefixes[6+i]
    d = round(d / 1e3**i, significantfigures)
    string = '{d} {prefix}'.format(d=d, prefix=prefix)
    return string


def SId(d, significantfigures=4):
    # Unit prefixes
    prefixes = ['a', 'f', 'p', 'n', r'$\mu$', 'm', '',
                'k', 'M', 'G', 'T', 'P', 'E']
    exp10 = log10(abs(d))
    i = int((exp10) // 3)
    prefix = prefixes[6+i]
    metric = 1e3**i
    d = round(d / metric, significantfigures)
    return prefix, metric
