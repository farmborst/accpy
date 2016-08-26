# -*- coding: utf-8 -*-
''' accpy.simulate.quadscan
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (array, hstack, vstack, dstack, zeros, linspace, sqrt, cos,
                   sin, cosh, sinh, eye, dot, empty, diag)
from functools import partial
from scipy.optimize import curve_fit
from . import const
from ..visualize.plot import pltsim_quadscan, pltmeas_quadscan
from .rmatrices import UC2T


def drift(L, y):
    R = eye(6)
    R[0, 1] = R[2, 3] = L
    R[4, 5] = L/(y**2)
    return R


def quadrupole(k, L, y):
    if k == 0:
        R = drift(L, y)
    else:
        wrzlk = sqrt(abs(k))
        Omega = wrzlk*L
        coshom = cosh(Omega)
        sinhom = sinh(Omega)
        cosom = cos(Omega)
        sinom = sin(Omega)
        if k > 0:
            R = array([
            [cosom,        sinom/wrzlk, 0.,           0.,           0., 0.],
            [-wrzlk*sinom, cosom,       0.,           0.,           0., 0.],
            [0.,            0.,         coshom,       sinhom/wrzlk, 0., 0.],
            [0.,            0.,         wrzlk*sinhom, coshom,       0., 0.],
            [0.,            0.,         0.,           0.,           1., L/(y**2)],
            [0.,            0.,         0.,           0.,           0., 1.]])
        if k < 0:
            R = array([
            [coshom,       sinhom/wrzlk, 0.,           0.,          0., 0],
            [wrzlk*sinhom, coshom,       0.,           0.,          0., 0],
            [0.,           0.,           cosom,        sinom/wrzlk, 0., 0],
            [0.,           0.,           -wrzlk*sinom, cosom,       0., 0],
            [0.,           0.,           0.,           0.,          1., L/(y**2)],
            [0.,           0.,           0.,           0.,          0., 1.]])
    return R


def thinlens(k, L, y):
    K = k*L
    R = array([[1.,    0.,     0.,     0.,     0.,     0.],
               [-K,     1.,    0.,     0.,     0.,     0.],
               [0.,    0.,     1,      0.,     0.,     0.],
               [0.,    0.,     K,      1.,     0.,     0.],
               [0.,    0.,     0.,     0.,     1.,     L/(y**2)],
               [0.,    0.,     0.,     0.,     0.,     1.]])
    return R


def part2mq(particle):
    if particle == 'electron':
        m = const.me
        q = const.qe
    elif particle == 'proton':
        m = const.mp
        q = const.qe
    return m, q


def part2y(E, particle):
    m, q = part2mq(particle)
    E0 = m*const.cl**2/q
    gamma = E/E0+1
    return gamma


def Rrange2sigxy(Rrange, points, twiss0, xdisp0, epsx, epsy, epss):
    twiss = empty([4, 4, points])
    disper = empty([3, points])
    for i, R in enumerate(Rrange):
        RXY = R[:4, :4]
        RD = vstack((R[:2, [0, 1, 5]], array([[0, 0, 1]])))
        twiss[:, :, i] = dot(dot(RXY, twiss0), RXY.T)
        disper[:, [i]] = dot(RD, xdisp0)
    # sig_u^2 = eps_u*beta_u + (delta_E*D)^2
    sigx = twiss[0, 0, :]*epsx + epss*(disper[0, :])**2
    sigy = twiss[2, 2, :]*epsy
    return sigx, sigy


def simulate_quadscan(ki, kf, qL, UC, points, epsx, betx, alpx, epsy,
                      bety, alpy, epss, Dx, Dpx, energy, particle, data):
    gamma = part2y(energy, particle)
    krange = linspace(ki, kf, points)
    # get transport matrix to FOM as function of k
    Rtransfer = UC2T(UC, gamma)
    Rquad = partial(quadrupole, L=qL, y=gamma)
    Rrange = [dot(Rtransfer, Rquad(k)) for k in krange]

    # thin lens approximation
    Rtransfer2 = dot(Rtransfer, drift(qL, gamma))
    Rquad2 = partial(thinlens, L=qL, y=gamma)
    Rrange2 = [dot(Rtransfer2, Rquad2(k)) for k in krange]

    # upstream beta (twiss) matrix and dispersion vector
    gamx = (1+alpx**2)/betx
    xtwi = array([[betx, -alpx], [-alpx, gamx]])
    gamy = (1+alpy**2)/bety
    ytwi = array([[bety, -alpy], [-alpy, gamy]])
    twiss0 = hstack((xtwi, zeros((2, 2))))
    twiss0 = vstack((twiss0, hstack((zeros((2, 2)), ytwi))))
    xdisp0 = array([[Dx], [Dpx], [1]])

    # transport twiss
    sigx, sigy = Rrange2sigxy(Rrange, points, twiss0, xdisp0, epsx, epsy, epss)
    sigx2, sigy2 = Rrange2sigxy(Rrange2, points, twiss0, xdisp0, epsx, epsy, epss)
    figs = pltsim_quadscan(krange, sigx, sigy, sigx2, sigy2, data)
    return figs


def fun(x, A, B, C):
    return A*x*x + 2*B*x + C


def myfit(x, y, betagamma, qL, orientation, T, yerr=None, krange=[]):
    if len(krange) != 0:
        # trim to given range
        booleanrange = [(x >= krange[0]) & (x <= krange[1])]
        x = x[booleanrange]
        y = y[booleanrange]
    if yerr is None:
        popt, pcov = curve_fit(fun, x, y, sigma=None)
    else:
        if len(krange) != 0:
            yerr = yerr[booleanrange]
        popt, pcov = curve_fit(fun, x, y, sigma=yerr, absolute_sigma=True)
    x = linspace(x[0], x[-1], 1e3)
    fit = fun(x, *popt)
    perr = sqrt(diag(pcov))
    A = array([popt[0], perr[0]])/qL**2
    B = array([popt[1], perr[1]])/qL*orientation
    C = array([popt[2], perr[2]])
    eps, bet, alp, gam = twissfromparabola(A, B, C, T)
    epsn = betagamma*eps*1e6
    eps *= 1e9
    string = (r'$\beta=%.2f \pm %.2f$''\n'
              r'$\alpha=%.2f \pm %.2f$''\n'
              r'$\gamma=%.2f \pm %.2f$''\n'
              r'$\epsilon=(%.2f \pm %.2f)\pi\, nm\, rad$''\n'
              r'$\epsilon^*=(%.2f \pm %.2f)\pi\, \mu m\, rad$'
              % (bet[0], bet[1], alp[0], alp[1], gam[0], gam[1],
                 eps[0], eps[1], epsn[0], epsn[1]))
    return x, fit, string


def twissfromparabola(A, B, C, T):
    # make lines shorter
    T11 = T[0, 0]
    T12 = T[0, 1]
    T12s = T12*T12
    # calculate sigma matrix elements
    sig11 = A[0]/T12s
    sig11er = A[1]/T12s
    sig12 = (B[0]-T11*T12*sig11)/T12s
    sig12er = sqrt((B[1])**2+(T11*T12*sig11er)**2)/T12s
    sig22 = (C[0]-T11**2*sig11-2*T11*T12*sig12)/T12s
    sig22er = sqrt(C[1]**2-(T11**2*sig11er)**2+(2*T11*T12*sig12er)**2)/T12s
    epsilon = sqrt(sig11*sig22-sig12**2)
    epser = sqrt((sig22*sig11er)**2+(sig11*sig22er)**2+(-2*sig12*sig12er)**2)/(2*epsilon)
    beta = sig11/epsilon
    betaer = sqrt((sig11er/epsilon)**2+(-sig11*epser/epsilon**2)**2)
    alpha = -sig12/epsilon
    alphaer = sqrt((-sig12er/epsilon)**2+(sig12*epser/epsilon**2)**2)
    gamma = sig22/epsilon
    gammaer = sqrt((sig22er/epsilon)**2+(-sig22*epser/epsilon**2)**2)
    epsilon = array([epsilon, epser])
    beta = array([beta, betaer])
    alpha = array([alpha, alphaer])
    gamma = array([gamma, gammaer])
    return epsilon, beta, alpha, gamma


def measure_quadscan(figs, data, kr_fit, kr_mes, points, qL, UC, epss, Dx,
                     Dpx, energy, particle):
    if data is None:
        kqua = linspace(kr_mes[0], kr_mes[1], points)
        # sigx, sigy = measure_sigxy(kqua)
        return
    else:
        kqua = data[0]
        kx = kqua
        ky = kqua
        sigx = data[1]
        sigy = data[2]
        if len(data) > 3:
            sigxe = data[3]
            sigye = data[4]
        else:
            sigxe = None
            sigye = None

    # get transport matrix to FOM as function of k
    gamma = part2y(energy, particle)
    betagamma = sqrt(gamma**2-1)
    Rtransfer = UC2T(UC, gamma)
    Rtransfer2 = dot(Rtransfer, drift(qL, gamma))
    RX = Rtransfer2[:2, :2]
    RY = Rtransfer2[2:4, 2:4]
    Rquad = partial(quadrupole, L=qL, y=gamma)
    Rrange = [dot(Rtransfer, Rquad(k)) for k in kx]

    # subtract k-dependant dispersion from sigmas
    xdisp0 = array([[Dx], [Dpx], [1]])
    points = len(kx)
    disper = empty([3, points])
    for i, R in enumerate(Rrange):
        RD = vstack((R[:2, [0, 1, 5]], array([[0, 0, 1]])))
        disper[:, [i]] = dot(RD, xdisp0)
    # sig_u^2 = eps_u*beta_u + (delta_E*D)^2
    sigx = sigx - epss*(disper[0, :])**2

    # orientation for k<0 in focussing plane
    ox = -1
    oy = 1

    kfx, fitx, stringx = myfit(kx, sigx, betagamma, qL, ox, RX, yerr=sigxe, krange=kr_fit[:2])
    kfy, fity, stringy = myfit(ky, sigy, betagamma, qL, oy, RY, yerr=sigye, krange=kr_fit[2:])

    # transport twiss
    pltmeas_quadscan(figs,
                     kx, sigx,
                     ky, sigy,
                     kfx, fitx,
                     kfy, fity,
                     [stringx, stringy],
                     xerr=sigxe, yerr=sigye)
    return
