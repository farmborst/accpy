# -*- coding: utf-8 -*-
''' accpy.simulate.particles
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import size
from . import const


def part2mq(particle):
    if particle == 'electron':
        m = const.me
        q = const.qe
    elif particle == 'proton':
        m = const.mp
        q = const.qe
    return m, q


def part2mqey(E, UC, particle):
    m, q = part2mq(particle)
    E0 = m*const.cl**2/q
    gamma = E/E0+1
    P_UC = size(UC, 1)        # nr of elements in unit cell
    return m, q, E0, gamma, P_UC


def part2E0(particle):
    m, q = part2mq(particle)
    E0 = m*const.cl**2/q
    return E0


def part2mqE0(particle):
    m, q = part2mq(particle)
    E0 = m*const.cl**2/q
    return m, q, E0