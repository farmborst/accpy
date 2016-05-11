#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.lsd
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (eye, dot, hstack, dstack, trapz, pi, zeros,
                   sqrt, linspace, sin, nanmean)
from matplotlib.pyplot import show
from . import const
from .slicing import cellslice
from .rmatrices import rmatrix, UCS2R
from .tracking import (initialtwiss, tracktwiss4)
from ..lattices import bessy2
from ..visualize.figures import plotstandards, figsave
from ..visualize.plot import (plotdisptraj, plotopticpars_closed,
                              plotopticpars_open, plotoptic)



def lsd(latt, slic=int(1e4), save=False, ft='pdf',
        plotstandard='presentation_1920x1080', scale=[1, 1]):
    ''' LSD - Lattice Simulation and Design
    '''
    plotstandards(plotstandard, scale)
    lattice = bessy2.lattice
    # get parameters and unit cell of lattice
    m, q, E0, E, I, gamma, UC, P_UC, diagnostics, closed, \
        rho, UD, D_UC, N_UC, \
        xtwiss, ytwiss, xdisp, HF_f, HF_V = lattice(latt)

    if closed:
        # one turn R matrix of unit cell
        M = eye(6)
        for i in range(P_UC):
            # R matrices of unsliced unit cell
            M = dot(rmatrix(UC[:, i], gamma), M)
        xtwiss, ytwiss, xdisp = initialtwiss(M)
        # points in one turn
        P_TOT = P_UC*N_UC
        # one turn R matrix of ring
        M1T = eye(6)
        for i in range(N_UC):
            M1T = dot(M, M1T)

    # get sliced unit cell for finer tracking
    s, UCS, P_UCS = cellslice(UC, P_UC, slic)

    # calculate according sliced R matrix
    R = UCS2R(P_UCS, UCS, gamma)

    # track twiss and dispersion
    xtwiss, ytwiss, xdisp, xytwiss = tracktwiss4(R, P_UCS, closed, xtwiss, ytwiss, xdisp)

    if closed:
        # tune Q_u:=1/2pi*int(ds/beta_u(s))
        Qx = N_UC*trapz(1./xtwiss[0, 0, :], s)/2/pi
        Qy = N_UC*trapz(1./ytwiss[0, 0, :], s)/2/pi
        # nat chromaticity xi_u:=1/4pi*int(k_u(s)*beta_u(s)) with k_y = - k_x
        Xx = N_UC*trapz(-UCS[4, :]*xtwiss[0, 0, :], s)/4/pi
        Xy = N_UC*trapz(UCS[4, :]*ytwiss[0, 0, :], s)/4/pi
        # calculate according ring of dipoles
        U = N_UC*s[-1:]
        beta = sqrt(1-1/gamma**2)
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
        # synchrotron integrals
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
        alpha_mc = 1/rho/U*trapz(disperdip[0, :], sdip)  # moment compaction
        # alpha_mc = M1T[4, 5]/U                            # moment compaction
        gamma_tr = 1/sqrt(alpha_mc)                      # transition energy
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
        Tr = 2*pi*rho/beta/const.cl  # radiating rev time for particles
        t_ph = 2*rho*sin(1/gamma)/c  # time for photons on straight
        t_e = (2/gamma)*rho/(beta*c)    # time for electron on arc
        dt = t_e-t_ph                   # time until electron and photon meet
        omega_c = 2/dt                  # critical wavenumber / (1/s)
        lambda_c = 2*pi*c/omega_c    # critical wavelength / (m)
        E_c = h_bar*omega_c/const.qe    # critical energy / (eV)
        # power per electron from Hinterberger / (W)
        P_rad = (c*q**2*E**4)/(2*pi*rho**2*3*e_0*E0**4)
        # energyloss per circumnavigation and electron / (Ws)
        U_rad = P_rad*Tr/const.qe
        # total radiated power
        N_e = I*Tr/const.qe         # number of electrons in accelerator
        U_ges = U_rad*const.qe*N_e  # radiated energy / circumnavigation / (Ws)
        P_ges = U_ges/Tr            # total power / (W)
        # equilibrium emittances and damping times (!amplitude damping *1/2)
        Cq = (55*const.hb)/(32*sqrt(3)*m*const.cl)
        Cgamma = 4*pi*const.re/(3*E0**3)
        Calpha = 2113.1e-27
        emiteqx = Cq*gamma*gamma*SYNIN5x/Jx/SYNIN2
        tau_x = T_rev*E/(U_rad*Jx)
        emiteqy = Cq*nanmean(ytwiss[0, 0, :])/(2*Jy*rho)
        tau_y = T_rev*E/(U_rad*Jy)
        # longitudinal (synchrotron motion)
        sigma_E = sqrt(Cq*gamma**2/(Js*rho))  # equilibrium energy spread
        tau_s = T_rev*E/(U_rad*Js)
    #    vs = HF_f
    #    ws = HF_f*2*pi
    #    h = T_rev*HF_f
    #    ns = eta_mc
    #    ps = E
    #    phi_s
    #    omega_syn = ws*sqrt(h*ns/2/pi/ps/vs*q*HF_V*np.cos(phi_s))
    #    Q_syn = omega_syn/ws
    #    print omega_syn
        v_s = 50e3              # synchrotron frequency
        omega_s = v_s*2*pi
        sigma_tau = sigma_E*abs(eta_mc)/(omega_s)
        sigma_s = sigma_tau*beta*c
        Q_s = omega_s*T_rev/2/pi

    # plot beta, disp and lattice
    fig0 = plotoptic(UC, diagnostics, s, xtwiss, ytwiss, xdisp)
    if closed:
        fig1 = plotopticpars_closed(xtwiss, xdisp, ytwiss, gamma, Qx, Xx, Jx, emiteqx,
                  tau_x, Qy, Xy, Jy, E, emiteqy, tau_y, alpha_mc, eta_mc,
                  gamma_tr, Q_s, Js, sigma_E, sigma_tau, sigma_s, tau_s, U_rad,
                  P_ges, E_c, lambda_c)
    else:
        fig1 = plotopticpars_open(xtwiss, xdisp, ytwiss, gamma, E)

    # track and plot reference particle with energy offset
    if not closed:
        fig2 = plotdisptraj(s, P_UCS, E, E0, UCS, UC, diagnostics)


    if save is True:
        figsave(fig0, latt, filetype=ft)
        figsave(fig1, ''.join((latt, 'parameters')), filetype=ft)
        figsave(fig2, ''.join((latt, 'offset')), filetype=ft)
    show()
