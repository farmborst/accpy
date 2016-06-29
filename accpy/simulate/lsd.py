#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.lsd
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (eye, dot, trapz, pi, nanmean, array, newaxis, hstack,
                   concatenate, empty, dstack, sqrt, zeros, vstack)
from numpy.random import standard_normal
from numpy.linalg import inv
from .slicing import cellslice
from .rmatrices import rmatrix, UCS2R
from .tracking import (initialtwiss, tracktwiss4, trackparts)
from .radiate import dipolering, synchroints
from .particles import part2mqey
from ..lattices.reader import latt2py
from ..visualize.plot import (plotopticpars_closed, plottrajs,
                              plotopticpars_open, plotoptic, plotphasespace)


def oneturn(UC, P_UC, N_UC, gamma):
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
    xtwiss0, ytwiss0, xdisp0 = initialtwiss(M)
    # one turn R matrix of ring
    M1T = eye(6)
    for i in range(N_UC):
        M1T = dot(M, M1T)
    return xtwiss0, ytwiss0, xdisp0, rho, D_UC, UD, LD


def gettunes(s, xtwiss, ytwiss, N_UC):
    Qx = N_UC*trapz(1./xtwiss[0, 0, :], s)/2/pi
    Qy = N_UC*trapz(1./ytwiss[0, 0, :], s)/2/pi
    return Qx, Qy

def getchromaticity(s, xtwiss, ytwiss, N_UC, UCS):
    Xx = N_UC*trapz(-UCS[4, :]*xtwiss[0, 0, 1:], s[1:])/4/pi
    Xy = N_UC*trapz(UCS[4, :]*ytwiss[0, 0, 1:], s[1:])/4/pi
    return Xx, Xy


def lsd(closed, latt, slices, mode, particles, rounds):
    if closed:
        (particle, E, I, UC, diagnostics, N_UC,
         HF_f, HF_V) = latt2py(latt, closed)
    else:
        (particle, E, I, UC, diagnostics, N_UC,
         xtwiss0, ytwiss0, xdisp0) = latt2py(latt, closed)

    m, q, E0, gamma, P_UC = part2mqey(E, UC, particle)

    if closed:
        xtwiss0, ytwiss0, xdisp0, rho, D_UC, UD, LD = oneturn(UC, P_UC, N_UC, gamma)

    # get sliced unit cell for finer tracking
    s, UCS, P_UCS = cellslice(UC, P_UC, slices)
    # calculate according sliced R matrix
    R = UCS2R(P_UCS, UCS, gamma)

    # track twiss and dispersion
    xtwiss, ytwiss, xdisp, xytwiss = tracktwiss4(R, P_UCS, closed, xtwiss0,
                                                 ytwiss0, xdisp0)
    if closed:
        # tune Q_u:=1/2pi*int(ds/beta_u(s))
        Qx, Qy = gettunes(s, xtwiss, ytwiss, N_UC)
        # nat chromaticity xi_u:=1/4pi*int(k_u(s)*beta_u(s)) with k_y = - k_x
        Xx, Xy = getchromaticity(s, xtwiss, ytwiss, N_UC, UCS)
        # calculate according ring of dipoles
        sdip, disperdip, xtwissdip, ytwissdip = \
            dipolering(s, N_UC, UD, P_UCS, UCS, xdisp, xtwiss, ytwiss, slices,
                       D_UC)
        # synchrotron integrals
        (Cq, Jx, emiteqx, tau_x, Jy, E, emiteqy, tau_y, alpha_mc, eta_mc,
         gamma_tr, Q_s, Js, sigma_E, sigma_tau, sigma_s, tau_s, U_rad, P_ges,
         E_c, lambda_c) = \
            synchroints(N_UC, s, gamma, xtwissdip, disperdip, sdip, rho, E, E0,
                        I, q, m, ytwiss)

    if mode == 'trackbeta':
        figs = []
        figs.append(plotoptic(UC, 'radial', diagnostics, s, xtwiss, ytwiss, xdisp))
        figs.append(plotoptic(UC, 'axial', diagnostics, s, xtwiss, ytwiss, xdisp))
        figs.append(plotoptic(UC, 'dispersion', diagnostics, s, xtwiss, ytwiss, xdisp))
        figs.append(plotoptic(UC, 'overview', diagnostics, s, xtwiss, ytwiss, xdisp))
        if closed:
            figs.append(plotopticpars_closed(xtwiss, xdisp, ytwiss, gamma, Qx,
                                             Xx, Jx, emiteqx, tau_x, Qy, Xy,
                                             Jy, E, emiteqy, tau_y, alpha_mc,
                                             eta_mc, gamma_tr, Q_s, Js,
                                             sigma_E, sigma_tau, sigma_s,
                                             tau_s, U_rad, P_ges, E_c,
                                             lambda_c))
        else:
            figs.append(plotopticpars_open(xtwiss, xdisp, ytwiss, gamma, E))
    elif mode == 'trackpart':
        # [x,  x',   y,  y',   l,  delta_p/p_0]
        # [mm, mrad, mm, mrad, mm, promille]
        ideal = array([0, 0, 0, 0, 0, 0])     # Ideal particle
        start = array([1, 1, 1, 1, 1, 0])     # 1 sigma particle
        distmean = 1e-3*ideal[newaxis, :].T
        distsigma = 1e-3*start[newaxis, :].T

        # emmitanz des vorgegebenen 1-sigma teilchens (Wille 3.142)
        emittx = dot(start[:2], dot(inv(xtwiss0), start[:2]))
        emitty = dot(start[2:4], dot(inv(ytwiss0), start[2:4]))

        # Envelope E(s)=sqrt(epsilon_i*beta_i(s))
        ydisp = zeros([1, P_UCS+1])
        emit_x_beta = array([emittx*xtwiss[0, 0, :], emitty*ytwiss[0, 0, :]])
        dispdelta = (vstack([xdisp[0, :], ydisp[0, :]])*1E-3*distsigma[5])**2
        envelope = sqrt(dispdelta + emit_x_beta)

        # start vectors of normally distributed ensemble
        points = P_UCS*N_UC*rounds
        X0 = (distsigma - distmean)*standard_normal([6, particles])
        X0 = dstack([X0, empty([6, particles, points])])
        X0[:, :2, 0] = hstack([distmean, distsigma])
        X_S = [X0[:, i, :] for i in range(particles)]
        X = trackparts(R, N_UC, X_S, rounds)
        s0 = s
        envelope0 = envelope
        for i in range(1, N_UC):
            s = concatenate([s, s0[1:]+s0[-1]*i])[:]
            envelope = hstack([envelope, envelope0[:, 1:]])
        figs = plottrajs(s, X, rounds, envelope)
        figs.append(plotphasespace(s, X, rounds, xtwiss, emittx, ytwiss, emitty))
    return figs
