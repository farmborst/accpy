#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.lsd
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import eye, dot, trapz, pi, nanmean, array, newaxis, hstack
from numpy.random import random_sample
from .slicing import cellslice
from .rmatrices import rmatrix, UCS2R
from .tracking import (initialtwiss, tracktwiss4)
from .radiate import dipolering, synchroints
from .particles import part2mqey
from ..lattices.bessy2 import lattice
from ..visualize.plot import (plotopticpars_closed,
                              plotopticpars_open, plotoptic)


def lsd(latt, slices, mode, particles, rounds):
#    if mode == 'trackbeta':
#        ...
#    elif mode == 'trackpart':
#        ideal = array([0, 0, 0, 0, 0, 0])     # (x, x', y, y', l, delta_p/p_0) in [mm,mrad,mm,mrad,mm,promille] ideal particle
#        start = array([1, 1, 1, 1, 1, 0])     # (x, x', y, y', l, delta_p/p_0) 1 sigma particle
#        distmean = 1e-3*ideal[newaxis, :].T   # Ideales Teilchen
#        distsigma = 1e-3*start[newaxis, :].T  # Teilchen mit 1 sigma
#        # start vectors of normally distributed ensemble up to 1 sigma particles
#        X_S = (distsigma - distmean)*random_sample(6, particles-2) + distmean
#        X_S = hstack([distmean, distsigma, X_S])

    # get parameters and unit cell of lattice
    (closed, particle, E, I, UC, diagnostics,    # always
     N_UC, HF_f, HF_V,                           # closed lattice
     xtwiss, ytwiss, xdisp) = lattice(latt)      # open lattice

    m, q, E0, gamma, P_UC = part2mqey(E, UC, particle)

    if closed:
        # one turn R matrix of unit cell
        M = eye(6)
        rho = []
        LD = []
        D_UC = 0
        for i in range(P_UC):
            UC_tmp = UC[:, i]
            # R matrices of unsliced unit cell
            M = dot(rmatrix(UC_tmp, gamma), M)
            if UC_tmp[0] == 1:
                LD.append(UC_tmp[1])
                rho.append(UC_tmp[2])
                D_UC += 1
        UC_tmp = None
        LD = nanmean(LD)
        rho = nanmean(rho)
        UD = N_UC*D_UC*LD
        xtwiss, ytwiss, xdisp = initialtwiss(M)
        # one turn R matrix of ring
        M1T = eye(6)
        for i in range(N_UC):
            M1T = dot(M, M1T)

    # get sliced unit cell for finer tracking (also get average rho)
    s, UCS, P_UCS = cellslice(UC, P_UC, slices)

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
            dipolering(s, N_UC, UD, P_UCS, UCS, xdisp, xtwiss, ytwiss, slices,
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
