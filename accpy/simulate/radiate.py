# -*- coding: utf-8 -*-
""" accpy.simulate.radiate
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division
from . import const


def radiationloss(energy, q, E0, rho):
    '''synchrotronradiation(energy, rho)
    input:
        particle energy
        particle charge
        particle rest energy
        bending radius
    return:
        energyloss per round due to synchrotron radiation at given energy and
        bending radius / eV
    notice:
        for electrons
    '''
    # energyloss per turn per elektron
    e = const.qe
    e0 = const.e0
    restenergy = const.Ee/e
    energy = q*energy**4/rho/3/e0/(E0**4)
    return energy