#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.ramp
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import pi, sin, linspace, sqrt, mean
from . import const
from .particles import part2mqE0
from ..visualize.plot import plotramp


def energy(amp, w, t_0, off):
    ''' Energy over ramp as function of time
    '''
    E = lambda t: amp*sin(w*(t+t_0))+off
    return E


def lorentz(E, E0):
    ''' Relativistic Lorentz factors as functions of time
    '''
    gamma = lambda t: abs(E(t))/E0+1         # 1/sqrt(1-(v/c)^2)
    beta = lambda t: sqrt(1-1/gamma(t)**2)   # v/c
    return gamma, beta


def Bfluxdensity(E, lorentzbeta, rho):
    ''' Magnetic flux density in T
    '''
    B = lambda t: E(t)/const.cl/lorentzbeta(t)/rho
    return B


def radiationloss(E, q, rho, E0):
    ''' per turn and electron / eV
    '''
    loss = lambda t: q*E(t)**4/rho/3/const.e0/(E0**4)
    return loss


def requiredvoltage(E, loss, Trev):
    # required voltage in cavity
    ''' required acceleration voltage in cavity to compensate radiation loss
        and reach desired energy / eV
    '''
    voltage = lambda t: E(t+Trev) - E(t) + loss(t)
    return voltage


def simulate_ramp(T, t_inj, t_ext, t_ext2, E_inj, E_ext, particle, ND, LD, U, points):
    m, q, E0 = part2mqE0(particle)
    UD = ND*LD                       # total orbit length of all dipoles
    rho = UD/2/const.pi              # bending radius

    # energy ramp
    f = 1/T                     # frequency of booster ramp
    w = 2*pi*f                  # angular frequency of booster ramp
    t_max = (t_ext+t_ext2)/2    # peak time
    # analytic calculations show for y = A*sin(2*pi*f(t+B))+C = A*sin(w(t+B))+C
    t_0 = 1/(4*f)-t_max         # max(y) = A*sin(w(tmax+B))+C = A + C -> B
    amp = (E_ext-E_inj)/(sin(w*(t_ext+t_0))-sin(w*(t_inj+t_0)))
    off = E_inj-amp*sin(w*(t_inj+t_0))

    # funtions of time
    E = energy(amp, w, t_0, off)
    lorentzgamma, lorentzbeta = lorentz(E, E0)
    B = Bfluxdensity(E, lorentzbeta, rho)
    loss = radiationloss(E, q, rho, E0)
    Trev = lambda t: U/(lorentzbeta(t)*const.cl)

    t = linspace(0, T, points)
    volt = requiredvoltage(E, loss, mean(Trev(t)))

    figs = plotramp(T, t, E(t), B(t), t_inj, t_ext, t_ext2, loss(t), volt(t))
    return figs
