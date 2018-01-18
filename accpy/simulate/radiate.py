# -*- coding: utf-8 -*-
""" accpy.simulate.radiate
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division
from . import const
from numpy import (hstack, dstack, zeros, linspace, sqrt, trapz, pi, sin,
                   nanmean)


def dipolering(s, N_UC, UD, P_UCS, UCS, xdisp, xtwiss, ytwiss, slic, D_UC):
    UCS1T = hstack(UCS for i in range(N_UC))
    xdisp_t = hstack(xdisp for i in range(N_UC))
    xtwiss_t = dstack(xtwiss for i in range(N_UC))
    ytwiss_t = dstack(ytwiss for i in range(N_UC))
    PD_TOT = N_UC*slic*D_UC
    disperdip = zeros([2, PD_TOT])
    xtwissdip = zeros([2, 2, PD_TOT])
    ytwissdip = zeros([2, 2, PD_TOT])
    j = 0
    for i in range(N_UC*P_UCS):
        if UCS1T[0, i] == 1:
            disperdip[:, j] = xdisp_t[:2, i]
            xtwissdip[:, :, j] = xtwiss_t[:, :, i]
            ytwissdip[:, :, j] = ytwiss_t[:, :, i]
            j += 1
    sdip = linspace(0, UD, PD_TOT)
    return sdip, disperdip, xtwissdip, ytwissdip


def synchroints(N_UC, s, gamma, xtwissdip, disperdip, sdip, rho, E, E0, I,
                q, m, ytwiss):
    U = N_UC*s[-1:]
    beta = sqrt(1-1/gamma**2)
    Hsx = (xtwissdip[1, 1, :]*disperdip[0, :]**2
           - 2*xtwissdip[0, 1, :]*disperdip[0, :]*disperdip[1, :]
           + xtwissdip[0, 0, :]*disperdip[1, :]**2)
    SYNIN1 = trapz(disperdip[0, :], sdip)/rho
    SYNIN2 = 2*pi/rho
    SYNIN4x = SYNIN1/rho/rho
    SYNIN4y = 0
    SYNIN5x = trapz(Hsx, sdip)/rho/rho/rho
    Jx = 1 - SYNIN4x/SYNIN2
    Jy = 1 - SYNIN4y/SYNIN2
    Js = 2 + (SYNIN4x+SYNIN4y)/SYNIN2
    # momentum dynamics
    alpha_mc = 1/rho/U*trapz(disperdip[0, :], sdip)     # moment compaction
    # alpha_mc = M1T[4, 5]/U                            # moment compaction
    gamma_tr = 1/sqrt(alpha_mc)                         # transition energy
    eta_mc = 1/gamma_tr**2-1/gamma**2                   # slip factor
    # synchrotron radiation
    # constants
    c = const.cl                    # speed of light / (m/s)
    e_0 = const.e0                  # vacuum permittivity / (As/Vm)
    r_e = const.re                  # classical radius of electron / (m)
    mu_0 = const.u0                 # vac permeability / (N/A^2)=(Henry/m)
    h_bar = const.hb                # reduced Planck constant / (Js)
    # spectrum
    T_rev = U/beta/const.cl         # total revolution time for particles
    Tr = 2*pi*rho/beta/const.cl     # radiating rev time for particles
    t_ph = 2*rho*sin(1/gamma)/c     # time for photons on straight
    t_e = (2/gamma)*rho/(beta*c)    # time for electron on arc
    dt = t_e-t_ph                   # time until electron and photon meet
    omega_c = 2/dt                  # critical wavenumber / (1/s)
    lambda_c = 2*pi*c/omega_c       # critical wavelength / (m)
    E_c = h_bar*omega_c/const.qe    # critical energy / (eV)
    # power per electron from Hinterberger / (W)
    P_rad = (c*q**2*E**4)/(2*pi*rho**2*3*e_0*E0**4)
    # energyloss per circumnavigation and electron / (Ws)
    U_rad = P_rad*Tr/const.qe
    # total radiated power
    N_e = I*Tr/const.qe             # number of electrons in accelerator
    U_ges = U_rad*const.qe*N_e      # radiated energy / circumnavigation / (Ws)
    P_ges = U_ges/Tr                # total power / (W)
    # equilibrium emittances and damping times (!amplitude damping *2)
    Cq = (55*const.hb)/(32*sqrt(3)*m*const.cl)
    emiteqx = Cq*gamma*gamma*SYNIN5x/Jx/SYNIN2
    tau_x = T_rev*E/(U_rad*Jx)
    emiteqy = Cq*nanmean(ytwiss[0, 0, :])/(2*Jy*rho)
    tau_y = T_rev*E/(U_rad*Jy)
    # longitudinal (synchrotron motion)
    sigma_E = sqrt(Cq*gamma**2/(Js*rho))  # equilibrium energy spread
    tau_s = T_rev*E/(U_rad*Js)

    v_s = 50e3 # synchrotron frequency
    omega_s = v_s*2*pi
    sigma_tau = sigma_E*abs(eta_mc)/(omega_s)
    sigma_s = sigma_tau*beta*c
    Q_s = omega_s*T_rev/2/pi
    return (Cq, Jx, emiteqx, tau_x, Jy, E, emiteqy, tau_y, alpha_mc, eta_mc,
            gamma_tr, Q_s, Js, sigma_E, sigma_tau, sigma_s, tau_s, U_rad,
            P_ges, E_c, lambda_c)
