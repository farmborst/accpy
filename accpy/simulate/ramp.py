#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.ramp
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import sin, cos, arcsin, linspace, sqrt, mean, trapz
from .const import pi, cl, e0, re, hb, qe
from .particles import part2mqey
from ..visualize.plot import plotramp
from ..lattices.bessy2 import lattice
from ..simulate.tracking import tracktwiss4
from ..simulate.lsd import oneturn
from ..simulate.radiate import dipolering
from ..simulate.slicing import cellslice
from ..simulate.rmatrices import UCS2R


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
    fsyn = lambda t: sqrt(h(t)*slipfactor(t)*cavityvoltage(t)*cos(phase(t))/2/pi/pc(t))/Trev(t)
    return fsyn


def overvoltage(OV, voltage):
    cavityvoltage = lambda t: OV*voltage(t)
    return cavityvoltage


def synchroints(sdip, xtwissdip, disperdip, rho, E0):
    # SOME SYNCHRTRON CONSTANTS
    hbar = hb/qe
    Cq = 55*hbar*cl/32/sqrt(3)/E0            # Chao: 3.8319e-13 m
    Ca = re*cl/3/E0**3                       # Chao: 2.113e-24 m²/(eV)³/s
    Cy = 4*pi/3*re/E0**3                    # Chao: 8.846e-32 m/(eV)³

    # H FUNCTION: H(s) = beta*D'^2 + 2*alpha*D*D' + gamma*D^2
    Hsx = (xtwissdip[1, 1, :]*disperdip[0, :]**2
           - 2*xtwissdip[0, 1, :]*disperdip[0, :]*disperdip[1, :]
           + xtwissdip[0, 0, :]*disperdip[1, :]**2)

    # SYNCHROTRIN INTEGRALS
    SYNIN1 = trapz(disperdip[0, :], sdip)/rho
    SYNIN2 = 2*pi/rho
    SYNIN4x = SYNIN1/rho/rho
    SYNIN4y = 0
    SYNIN5x = trapz(Hsx, sdip)/rho/rho/rho

    # DAMPING PARTITION NUMBERS
    Jx = 1 - SYNIN4x/SYNIN2
    Jy = 1 - SYNIN4y/SYNIN2
    Js = 2 + (SYNIN4x+SYNIN4y)/SYNIN2
    return SYNIN1, SYNIN2, SYNIN4x, SYNIN4y, SYNIN5x, Jx, Jy, Js, Cq


def momentumdynamics(sdip, disperdip, lorentzgamma, C, rho):
    alpha_mc = 1/rho/C*trapz(disperdip[0, :], sdip)     # moment compaction
    # alpha_mc = M1T[4, 5]/U                            # moment compaction
    gamma_tr = 1/sqrt(alpha_mc)                         # transition energy
    eta_mc = lambda t: 1/gamma_tr**2-1/lorentzgamma(t)**2   # slip factor
    return eta_mc, gamma_tr, alpha_mc

def dampingdecrements(Ca, C, E, SYNIN2, Jx, Jy, Js):
    alphax = lambda t: Ca/C*E(t)**3*SYNIN2*Jx
    alphay = lambda t: Ca/C*E(t)**3*SYNIN2*Jy
    alphas = lambda t: Ca/C*E(t)**3*SYNIN2*Js
    return alphax, alphay, alphas


def bunchlength(fsyn, energyspread, lorentzbeta, slipfactor):
    bdur = lambda t: energyspread(t)*abs(slipfactor(t))/(fsyn(t)*2*pi)
    blen = lambda t: bdur(t)*lorentzbeta(t)*cl
    return bdur, blen


def quantumexcitation(Cq, Ca, lorentzgamma, E, SYNIN5x, SYNIN2, rho):
    ''' QUANTUM EXCITATION (per round)
    quantex = lambda t: 55*re*hbar*c*c/24/sqrt(3)/E0*gamma(t)**5*SYNIN5x/C
    following formulas give same results!!!
    Chao (3.1.4.3):   quantex = 55*re*hbar*c/48/sqrt(3)/E0**6*SYNIN5x
    my quantum excitation formulas are based on Wiedemann (13.15)
    '''
    x = lambda t: Cq*lorentzgamma(t)**2*Ca*2*E(t)**3*SYNIN5x/SYNIN2/rho**2
    s = lambda t: Cq*re*cl*2/3*lorentzgamma(t)**5*rho
    return x, s


# RADIAL EMITTANCE
def Xequilibriumemittance(Cq, lorentzgamma, SYNIN2, SYNIN5x, Jx):
    # from emit_dot != 0 we get:
    emit = lambda t: Cq*lorentzgamma(t)**2/Jx*SYNIN5x/SYNIN2
    return emit


def Xemittancedot(E, Edot, emitx, quantex, alphax):
    def emitdot(t, emitx):
        # adiabatic damping: -emitx*Edot(t)/E(t)
        # radiation damping: -2*emitx*alphax(t)
        # quantumexcitation: quantex(t)/2  ("radiation anti-damping")
        y = quantex(t)-(Edot(t)/E(t)+2*alphax(t))*emitx
        return y
    return emitdot


# AXIAL EMITTANCE
def Yequilibriumemittance(Cq, ytwiss, rho, Jy):
    # from emit_dot != 0 we get:
    emit = Cq*mean(ytwiss[0, 0, :])/(2*Jy*rho)
    return emit


def Yemittancedot(E, Edot, emity, alphay):
    def emitdot(t, emity):
        # adiabatic damping: -emity*Edot(t)/E(t)
        # radiation damping: -2*emity*alphay(t)
        # quantumexcitation: none
        # See Hemmie 'Reduction of antidamping' or Borland SRFEL-003
        y = -(Edot(t)/E(t)+2*alphay(t))*emity
        return y
    return emitdot


# LONGITUDINAL EMITTANCE
def Sequilibriumemittance(Cq, lorentzgamma, Js, rho):
    # from emit_dot != 0 we get:
    emit = lambda t: sqrt(Cq*lorentzgamma(t)**2/(Js*rho))
    return emit


def Semittancedot(E, Edot, emitz, Cq, Js, alphas, lorentzgamma, rho):
    def emitdot(t, emitz):
        # adiabatic damping: none
        # radiation damping: -2*emit*alphax(t)
        # quantumexcitation: quantex(t)/2
        # y = quantes(t)-(Edot(t)/E(t)+2*alphas(t))*emitz
        y = (Edot(t)/E(t)+2*alphas(t))*sqrt(Cq*lorentzgamma(t)**2/(Js*rho))-(Edot(t)/E(t)+2*alphas(t))*emitz
        return y
    return emitdot


def simulate_ramp(T, t_inj, t_ext, t_ext2, E_inj, E_ext, latt, points, f_hf):
    # get parameters and unit cell of lattice
    (closed, particle, E, I, UC, diagnostics, N_UC,     # always
     HF_f, HF_V,                                        # closed lattice
     xtwiss0, ytwiss0, xdisp0) = lattice(latt)          # open lattice
    m, q, E0, gamma, P_UC = part2mqey(E, UC, particle)
    xtwiss0, ytwiss0, xdisp0, rho, D_UC, UD, LD = oneturn(UC, P_UC, N_UC, gamma)
    C = sum(UC[1, :])*N_UC  # circumference

    # energy ramp
    t = linspace(0, T, points)
    f = 1/T                     # frequency of booster ramp
    w = 2*pi*f            # angular frequency of booster ramp
    t_max = (t_ext+t_ext2)/2    # peak time
    # analytic calculations show for y = A*sin(2*pi*f(t+B))+C = A*sin(w(t+B))+C
    t_0 = 1/(4*f)-t_max         # max(y) = A*sin(w(tmax+B))+C = A + C -> B
    amp = (E_ext-E_inj)/(sin(w*(t_ext+t_0))-sin(w*(t_inj+t_0)))
    off = E_inj-amp*sin(w*(t_inj+t_0))
    # twisstracking to get required observables
    s, UCS, P_UCS = cellslice(UC, P_UC, points)
    # calculate according sliced R matrix
    R = UCS2R(P_UCS, UCS, gamma)
    # track twiss and dispersion for 1 unit cell
    xtwiss, ytwiss, xdisp, xytwiss = tracktwiss4(R, P_UCS, closed, xtwiss0, ytwiss0, xdisp0)
    # calculate according ring of dipoles
    sdip, disperdip, xtwissdip, ytwissdip = dipolering(s, N_UC, UD, P_UCS, UCS, xdisp, xtwiss, ytwiss, points, D_UC)
    # synchrotron integrals
    SYNIN1, SYNIN2, SYNIN4x, SYNIN4y, SYNIN5x, Jx, Jy, Js, Cq = synchroints(sdip, xtwissdip, disperdip, rho, E0)

    # funtions of time
    E, Edot = energy(amp, w, t_0, off)
    lorentzgamma, lorentzbeta = lorentz(E, E0)
    lorentzbetagamma = lambda t: lorentzbeta(t)*lorentzgamma(t)
    slipfactor, gamma_tr, alpha_mc = momentumdynamics(sdip, disperdip, lorentzgamma, C, rho)
    B = Bfluxdensity(E, lorentzbeta, rho)
    loss = radiationloss(E, q, rho, E0)
    Trev = lambda t: C/(lorentzbeta(t)*cl)
    volt = requiredvoltage(E, loss, mean(Trev(t)))
    cavityvoltages = [overvoltage(OV, volt) for OV in [1, 2, 5, 10]]
    particlephases = [synchrophase(volt, cav) for cav in cavityvoltages]
    fsyns = [synchrofrequency(E, cav, phase, Trev, E0, f_hf, slipfactor) for phase, cav in zip(particlephases, cavityvoltages)]

    Xemitequi = Xequilibriumemittance(Cq, lorentzgamma, SYNIN2, SYNIN5x, Jx)
    Yemitequi = Yequilibriumemittance(Cq, ytwiss, rho, Jy)
    Semitequi = Sequilibriumemittance(Cq, lorentzgamma, Js, rho)

    bdurequi, blenequi = bunchlength(fsyns[0], Semitequi, lorentzbeta, slipfactor)


    phases = [phase(t) for phase in particlephases]
    freqs = [fsyn(t) for fsyn in fsyns]

    figs = plotramp(T, t, E(t), B(t), t_inj, t_ext, t_ext2, loss(t), volt(t),
                    phases, freqs, Xemitequi, Yemitequi, Semitequi,
                    bdurequi, blenequi, lorentzbetagamma)
    return figs
