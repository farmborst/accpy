#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.ramp
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import pi, sin, linspace, sqrt
from . import const
from .radiate import radiationloss
from .particles import part2mqE0
from ..visualize.plot import plotramp


def energy(T, t_inj, t_ext, t_ext2, E_inj, E_ext, particle, ND, LD):
    m, q, E0 = part2mqE0(particle)
    f = 1/T                     # frequency of booster ramp
    w = 2*pi*f                  # angular frequency of booster ramp
    t_max = (t_ext+t_ext2)/2    # peak time

    # analytic calculations show for y = A*sin(2*pi*f(t+B))+C = A*sin(w(t+B))+C
    t_0 = 1/(4*f)-t_max         # max(y) = A*sin(w(tmax+B))+C = A + C -> B
    amp = (E_ext-E_inj)/(sin(w*(t_ext+t_0))-sin(w*(t_inj+t_0)))
    off = E_inj-amp*sin(w*(t_inj+t_0))
    energy = lambda t: amp*sin(w*(t+t_0))+off

    t = linspace(0, T, 1e3)
    E = energy(t)

    lorentzgamma = abs(E)/E0 + 1     # 1/sqrt(1-(v/c)^2)
    lorentzbeta = sqrt(1-1/lorentzgamma**2)     # v/c
    UD = ND*LD                       # total orbit length of all dipoles
    rho = UD/2/const.pi              # bending radius
    B = E/const.cl/lorentzbeta/rho
    loss = radiationloss(E, q, E0, rho)
    figs = plotramp(T, t, E, B, t_inj, t_ext, t_ext2, loss)
    return figs
