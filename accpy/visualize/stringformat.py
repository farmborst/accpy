# -*- coding: utf-8 -*-
"""accpy.visualize.stringformat
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division
from math import log10


class uc:
    class greek:
        Alpha   = u'\u0391'
        Beta    = u'\u0392'
        Gamma   = u'\u0393'
        Delta   = u'\u0394'
        Epsilon = u'\u0395'
        Zeta    = u'\u0396'
        Eta     = u'\u0397'
        Theta   = u'\u0398'
        Iota    = u'\u0399'
        Kappa   = u'\u039a'
        Lambda  = u'\u039b'
        Mu      = u'\u039c'
        Nu      = u'\u039d'
        Xi      = u'\u039e'
        Omicron = u'\u039f'
        Pi      = u'\u03a0'
        Rho     = u'\u03a1'
        Sigma   = u'\u03a3'
        Tau     = u'\u03a4'
        Upsilon = u'\u03a5'
        Phi     = u'\u03a6'
        Chi     = u'\u03a7'      
        Psi     = u'\u03a8'
        Omega   = u'\u03a9'
        alpha   = u'\u03b1'
        beta    = u'\u03b2'
        gamma   = u'\u03b3'
        delta   = u'\u03b4'
        epsilon = u'\u03b5'
        zeta    = u'\u03b6'
        eta     = u'\u03b7'
        theta   = u'\u03b8'
        iota    = u'\u03b9'
        kappa   = u'\u03ba'
        lamda   = u'\u03bb'
        mu      = u'\u03bc'
        nu      = u'\u03bd'
        xi      = u'\u03be'
        omicron = u'\u03bf'
        pi      = u'\u03c0'
        rho     = u'\u03c1'
        sigma   = u'\u03c3'
        tau     = u'\u03c4'
        upsilon = u'\u03c5'
        phi     = u'\u03c6'
        chi     = u'\u03c7'
        psi     = u'\u03c8'
        omega   = u'\u03c9'
    class math:
        permille = u'\u2030'
        squared  = u'\u00B2'
    pi = u'\u03c0'
    ppt = u'\u2030'
    squared = u'\u00B2'
    alpha = u'\u03B1'
    beta = u'\u03B2'
    delta =	u"\u03B4"
    epsilon = u'\u03B5'


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


def time2str(t):
    dd = t//86400
    t %= 86400
    hh = t//3600
    t %= 3600
    mm = t//60
    t %= 60
    ss = t//1
    t %= 1
    ms = (t*1e3)//1
    if dd > 0:
        timestring = '%g days %g hours %g minutes %g.%g seconds' % (dd, hh, mm, ss, ms)
    elif hh > 0:
        timestring = '%g hours %g minutes %g.%g seconds' % (hh, mm, ss, ms)
    elif mm > 0:
        timestring = '%g minutes %g.%g seconds' % (mm, ss, ms)
    elif ss > 0:
        timestring = '%g.%g seconds' % (ss, ms)
    else:
        timestring = '%g milliseconds' % (ms)
    return timestring
