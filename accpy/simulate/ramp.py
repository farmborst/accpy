#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.ramp
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (sin, cos, arcsin, linspace, sqrt, mean, trapz, array, where,
                   repeat)
from itertools import product
from .const import pi, cl, e0, re, hb, qe
from .particles import part2mqey
from .tracking import tracktwiss4
from .lsd import oneturn
from .radiate import dipolering
from .slicing import cellslice
from .rmatrices import UCS2R
from ..visualize.plot import plotramp
from ..lattices.bessy2 import lattice
from ..math.ode import odeint


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
    return SYNIN1, SYNIN2, SYNIN4x, SYNIN4y, SYNIN5x, Jx, Jy, Js, Cq, Ca


def momentumdynamics(sdip, disperdip, lorentzgamma, C, rho):
    alpha_mc = 1/rho/C*trapz(disperdip[0, :], sdip)     # moment compaction
    # alpha_mc = M1T[4, 5]/U                            # moment compaction
    gamma_tr = 1/sqrt(alpha_mc)                         # transition energy
    eta_mc = lambda t: 1/gamma_tr**2-1/lorentzgamma(t)**2   # slip factor
    return eta_mc, gamma_tr, alpha_mc


def dampingdecrements(Ca, Cdip, E, SYNIN2, Jx, Jy, Js):
    alphax = lambda t: Ca/Cdip*E(t)**3*SYNIN2*Jx
    alphay = lambda t: Ca/Cdip*E(t)**3*SYNIN2*Jy
    alphas = lambda t: Ca/Cdip*E(t)**3*SYNIN2*Js
    return alphax, alphay, alphas


def bunchlength(lorentzbeta, slipfactor, bdur):
    blen = lambda t: bdur(t)*lorentzbeta(t)*cl
    return blen


def bunchduration(fsyn, energyspread, slipfactor):
    bdur = lambda t: energyspread(t)*abs(slipfactor(t))/(fsyn(t)*2*pi)
    return bdur


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


def Xemittancedot(E, Edot, quantex, alphax):
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


def Yemittancedot(E, Edot, alphay):
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


def Semittancedot(E, Edot, Cq, Js, alphas, lorentzgamma, rho):
    def emitdot(t, emits):
        # adiabatic damping: none
        # radiation damping: -2*emit*alphax(t)
        # quantumexcitation: quantex(t)/2
        # y = quantes(t)-(Edot(t)/E(t)+2*alphas(t))*emits
        y = (Edot(t)/E(t)+2*alphas(t))*sqrt(Cq*lorentzgamma(t)**2/(Js*rho))-(Edot(t)/E(t)+2*alphas(t))*emits
        return y
    return emitdot


def simulate_ramp(T, t_inj, t_ext, t_ext2, E_inj, E_ext, latt, points, f_hf,
                  V_HFs, emitxs, emitys, emitss):
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
    w = 2*pi*f                  # angular frequency of booster ramp
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
    Cdip = sdip[-1]
    # synchrotron integrals
    SYNIN1, SYNIN2, SYNIN4x, SYNIN4y, SYNIN5x, Jx, Jy, Js, Cq, Ca = synchroints(sdip, xtwissdip, disperdip, rho, E0)

    # funtions of time
    f_E, f_Edot = energy(amp, w, t_0, off)
    f_lorentzgamma, f_lorentzbeta = lorentz(f_E, E0)
    slipfactor, gamma_tr, alpha_mc = momentumdynamics(sdip, disperdip, f_lorentzgamma, C, rho)
    f_B = Bfluxdensity(f_E, f_lorentzbeta, rho)
    f_loss = radiationloss(f_E, q, rho, E0)
    Trev = lambda t: C/(f_lorentzbeta(t)*cl)
    f_volt = requiredvoltage(f_E, f_loss, mean(Trev(t)))
    volt = f_volt(t)
    overvoltagefactors = array(V_HFs)/max(volt)
    cavityvoltages = [overvoltage(OV, f_volt) for OV in overvoltagefactors]
    particlephases = [synchrophase(f_volt, cav) for cav in cavityvoltages]
    fsyns = [synchrofrequency(f_E, cav, phase, Trev, E0, f_hf, slipfactor) for phase, cav in zip(particlephases, cavityvoltages)]

    f_Xemitequi = Xequilibriumemittance(Cq, f_lorentzgamma, SYNIN2, SYNIN5x, Jx)
    Yemitequi = Yequilibriumemittance(Cq, ytwiss, rho, Jy)
    f_Semitequi = Sequilibriumemittance(Cq, f_lorentzgamma, Js, rho)

    # funtions of time and initial value (odes)
    alphax, alphay, alphas = dampingdecrements(Ca, Cdip, f_E, SYNIN2, Jx, Jy, Js)
    f_quantex, f_quantes = quantumexcitation(Cq, Ca, f_lorentzgamma, f_E, SYNIN5x, SYNIN2, rho)
    f_Xemitdot = Xemittancedot(f_E, f_Edot, f_quantex, alphax)
    f_Yemitdot = Yemittancedot(f_E, f_Edot, alphay)
    f_Semitdot = Semittancedot(f_E, f_Edot, Cq, Js, alphas, f_lorentzgamma, rho)
    f_bdurequis = [bunchduration(x, f_Semitequi, slipfactor) for x in fsyns]
    f_blenequis = [bunchlength(f_lorentzbeta, slipfactor, x) for x in f_bdurequis]

    # prepare data for plots
    E = f_E(t)
    B = f_B(t)
    loss = f_loss(t)
    i1 = where(t > t_inj)[0][0]
    i2 = where(t > t_ext)[0][0]
    i4 = where(t < t_ext2)[0][-1:]
    t_max = (t_ext+t_ext2)/2
    i3 = where(t > t_max)[0][0]
    i5 = where(volt < 0)[0][0]
    i6 = where(volt == max(volt))[0][0]
    i7 = where(E > 0)[0][0]
    i8 = where(E > 0)[0][-1]
    i9 = where(E > E[i1])[0][-1]
    tt = array([t[i1], t[i2], t[i3], t[i4]])
    tt2 = array([t[i1], t[i2], t[i3], t[i4], t[i5], t[i6]])
    EE = array([E[i1], E[i2], E[i3], E[i4]])
    BB = array([B[i1], B[i2], B[i3], B[i4]])
    LL = array([loss[i1], loss[i2], loss[i3], loss[i4]])
    VV = array([volt[i1], volt[i2], volt[i3], volt[i4], volt[i5], volt[i6]])
    # time where energy > 0
    tEgZ = t[i7:i8]
    EEgZ = E[i7:i8]
    # time after injection where energy > E_inj
    tAI = t[i1:i9]
    EAI = E[i1:i9]
    # time after injection where energy and acceleration voltage > 0
    tVgZ = t[i1:i5]
    EVgZ = E[i1:i5]

    # get data from functions of t for plots
    phases = [phase(t) for phase in particlephases]
    freqs = [fsyn(tVgZ) for fsyn in fsyns]

    # get data from functions of tEgZ for plots
    bdurequis = [f_bdurequi(tVgZ) for f_bdurequi in f_bdurequis]
    blenequis = [f_blenequi(tVgZ) for f_blenequi in f_blenequis]
    Xemitequi = f_Xemitequi(tAI)
    Yemitequi = repeat(Yemitequi, len(tAI))
    Semitequi = f_Semitequi(tAI)*1e3

    # get data from functions of tAI for plots
    Xemits = [odeint(f_Xemitdot, tAI, ex) for ex in emitxs]
    Yemits = [odeint(f_Yemitdot, tAI, ey, atol=1e-18)+Yemitequi for ey in emitys]
    # careful "late binding" !!!!!!!!!!!!
    f_Semits = [lambda t, es=es: odeint(f_Semitdot, t, es) for es in emitss]
    Semits = [f_Semit(tAI)*1e3 for f_Semit in f_Semits]
    f_bdurs = [bunchduration(fsyn, f_Semit, slipfactor) for fsyn, f_Semit in product(fsyns, f_Semits)]
    f_blens = [bunchlength(f_lorentzbeta, slipfactor, x) for x in f_bdurs]
    bdurs = [f_bdur(tVgZ) for f_bdur in f_bdurs]
    blens = [f_blen(tVgZ) for f_blen in f_blens]

    # calculate normalized emittances
    lorentzbetagammaAI = f_lorentzgamma(tAI)*f_lorentzbeta(tAI)
    NXemitequi = lorentzbetagammaAI*Xemitequi
    NYemitequi = lorentzbetagammaAI*Yemitequi
    NXemits = [lorentzbetagammaAI*ex for ex in Xemits]
    NYemits = [lorentzbetagammaAI*ey for ey in Yemits]


    figs = plotramp(T, t, tt, tt2, tEgZ, tAI, tVgZ, E, EE, EEgZ, EAI, EVgZ, B,
                    BB, loss, LL, volt, VV, phases, freqs, Xemitequi,
                    Yemitequi, Semitequi, bdurequis, blenequis, V_HFs,
                    Xemits, Yemits, Semits, NXemitequi, NYemitequi, NXemits,
                    NYemits, bdurs, blens)

    return figs
