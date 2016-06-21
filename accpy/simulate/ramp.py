#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.ramp
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import sin, cos, arcsin, linspace, sqrt, mean
from .const import pi, cl, e0
from .particles import part2mqE0
from ..visualize.plot import plotramp


def energy(amp, w, t_0, off):
    ''' Energy over ramp as function of time
    '''
    E = lambda t: amp*sin(w*(t+t_0))+off
    Edot = lambda t: w*amp*cos(w*(t+t_0))
    return E, Edot


def lorentz(E, E0):
    ''' Relativistic Lorentz factors as functions of time
    '''
    gamma = lambda t: abs(E(t))/E0+1         # 1/sqrt(1-(v/c)^2)
    beta = lambda t: sqrt(1-1/gamma(t)**2)   # v/c
    return gamma, beta


def Bfluxdensity(E, lorentzbeta, rho):
    ''' Magnetic flux density in T
    '''
    B = lambda t: E(t)/cl/lorentzbeta(t)/rho
    return B


def radiationloss(E, q, rho, E0):
    ''' per turn and electron / eV
    '''
    loss = lambda t: q*E(t)**4/rho/3/e0/(E0**4)
    return loss


def requiredvoltage(E, loss, Trev):
    ''' required acceleration voltage in cavity to compensate radiation loss
        and reach desired energy / eV
    '''
    voltage = lambda t: E(t+Trev) - E(t) + loss(t)
    return voltage


def synchrophase(requiredvoltage, cavityvoltage):
    ''' Synchronous phase at which electrons see the correct voltage in the
    cavity with cavity field:
        seenvoltage = cavityvoltage*sin(phase)'''
    phase = lambda t: arcsin(requiredvoltage(t)/cavityvoltage(t))/2/pi
    return phase


def synchrofrequency(energy, cavityvoltage, phase, Trev, E0, f_hf, slipfactor):
    ''' Synchronous phase at which electrons see the correct voltage in the
    cavity'''
    h = lambda t: Trev(t)*f_hf  # harmonic number
    pc = lambda t: sqrt(energy(t)**2-E0**2)
    # peaks up at energy zerocrossing and down at voltage zerocrossing
    # maybe abs(cavityvoltage) ?
    fsyn = lambda t: sqrt(h(t)*slipfactor*cavityvoltage(t)*cos(phase(t))/2/pi/pc(t))/Trev(t)
    return fsyn


def bunchlength(fsyn, energyspread, lorentzbeta, slipfactor):
    bdur = lambda t: energyspread(t)*abs(slipfactor)/(fsyn(t)*2*pi)
    blen = lambda t: bdur(t)*lorentzbeta(t)*cl
    return bdur, blen

def overvoltage(OV, voltage):
    cavityvoltage = lambda t: OV*voltage(t)
    return cavityvoltage


def simulate_ramp(T, t_inj, t_ext, t_ext2, E_inj, E_ext, particle, ND, LD, U,
                  points, f_hf, slipfactor):
    t = linspace(0, T, points)
    m, q, E0 = part2mqE0(particle)
    UD = ND*LD                      # total orbit length of all dipoles
    rho = UD/2/pi                   # bending radius

    # energy ramp
    f = 1/T                     # frequency of booster ramp
    w = 2*pi*f            # angular frequency of booster ramp
    t_max = (t_ext+t_ext2)/2    # peak time
    # analytic calculations show for y = A*sin(2*pi*f(t+B))+C = A*sin(w(t+B))+C
    t_0 = 1/(4*f)-t_max         # max(y) = A*sin(w(tmax+B))+C = A + C -> B
    amp = (E_ext-E_inj)/(sin(w*(t_ext+t_0))-sin(w*(t_inj+t_0)))
    off = E_inj-amp*sin(w*(t_inj+t_0))

    # funtions of time
    E, Edot = energy(amp, w, t_0, off)
    lorentzgamma, lorentzbeta = lorentz(E, E0)
    B = Bfluxdensity(E, lorentzbeta, rho)
    loss = radiationloss(E, q, rho, E0)
    Trev = lambda t: U/(lorentzbeta(t)*cl)
    volt = requiredvoltage(E, loss, mean(Trev(t)))
    cavityvoltages = [overvoltage(OV, volt) for OV in [1, 2, 5, 10]]
    particlephases = [synchrophase(volt, cav) for cav in cavityvoltages]
    fsyns = [synchrofrequency(E, cav, phase, Trev, E0, f_hf, slipfactor) for phase, cav in zip(particlephases, cavityvoltages)]

    phases = [phase(t) for phase in particlephases]
    freqs = [fsyn(t) for fsyn in fsyns]

    figs = plotramp(T, t, E(t), B(t), t_inj, t_ext, t_ext2, loss(t), volt(t),
                    phases, freqs)
    return figs
