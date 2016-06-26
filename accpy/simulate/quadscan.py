#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.quadscan
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
import numpy as np

def simulate_quadscan(twiss4, krange, driftlength):
    for i in range(100):
        a = np.random.rand(1000, 1000)
        b = np.random.rand(1000, 1000)
        c = np.dot(a, b)
    figs = []
    return figs