# -*- coding: utf-8 -*-
"""accpy.visualize.stringformat
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division
from math import log10


class uc:
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
