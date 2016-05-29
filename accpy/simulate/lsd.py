#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.lsd
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import eye, dot,  trapz, pi
from .slicing import cellslice
from .rmatrices import rmatrix, UCS2R
from .tracking import (initialtwiss, tracktwiss4)
from .radiate import dipolering, synchroints
from ..lattices import bessy2
from ..visualize.plot import (plotopticpars_closed,
                              plotopticpars_open, plotoptic)


def lsd(latt, slic=int(1e4)):

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
        # one turn R matrix of ring
        M1T = eye(6)
        for i in range(N_UC):
            M1T = dot(M, M1T)

    # get sliced unit cell for finer tracking
    s, UCS, P_UCS = cellslice(UC, P_UC, slic)

    # calculate according sliced R matrix
    R = UCS2R(P_UCS, UCS, gamma)

    # track twiss and dispersion
    xtwiss, ytwiss, xdisp, xytwiss = tracktwiss4(R, P_UCS, closed, xtwiss,
                                                 ytwiss, xdisp)

    if closed:
        # tune Q_u:=1/2pi*int(ds/beta_u(s))
        Qx = N_UC*trapz(1./xtwiss[0, 0, :], s)/2/pi
        Qy = N_UC*trapz(1./ytwiss[0, 0, :], s)/2/pi
        # nat chromaticity xi_u:=1/4pi*int(k_u(s)*beta_u(s)) with k_y = - k_x
        Xx = N_UC*trapz(-UCS[4, :]*xtwiss[0, 0, :], s)/4/pi
        Xy = N_UC*trapz(UCS[4, :]*ytwiss[0, 0, :], s)/4/pi
        # calculate according ring of dipoles
        sdip, disperdip, xtwissdip, ytwissdip = \
            dipolering(s, N_UC, UD, P_UCS, UCS, xdisp, xtwiss, ytwiss, slic,
                       D_UC)
        # synchrotron integrals
        (Jx, emiteqx, tau_x, Jy, E, emiteqy, tau_y, alpha_mc, eta_mc, gamma_tr,
         Q_s, Js, sigma_E, sigma_tau, sigma_s, tau_s, U_rad, P_ges, E_c,
         lambda_c) = \
            synchroints(N_UC, s, gamma, xtwissdip, disperdip, sdip, rho, E, E0,
                        I, q, m, ytwiss)

    fig_radial = plotoptic(UC, 'radial', diagnostics, s, xtwiss, ytwiss, xdisp)
    fig_axial = plotoptic(UC, 'axial', diagnostics, s, xtwiss, ytwiss, xdisp)
    fig_dispersion = plotoptic(UC, 'dispersion', diagnostics, s, xtwiss, ytwiss, xdisp)
    fig_overview = plotoptic(UC, 'overview', diagnostics, s, xtwiss, ytwiss, xdisp)

    if closed:
        fig_pars = plotopticpars_closed(xtwiss, xdisp, ytwiss, gamma, Qx, Xx, Jx,
                                    emiteqx, tau_x, Qy, Xy, Jy, E, emiteqy,
                                    tau_y, alpha_mc, eta_mc, gamma_tr, Q_s, Js,
                                    sigma_E, sigma_tau, sigma_s, tau_s, U_rad,
                                    P_ges, E_c, lambda_c)
    else:
        fig_pars = plotopticpars_open(xtwiss, xdisp, ytwiss, gamma, E)

    return fig_radial, fig_axial, fig_dispersion, fig_overview, fig_pars
