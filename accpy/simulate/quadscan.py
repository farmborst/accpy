# -*- coding: utf-8 -*-
''' accpy.simulate.quadscan
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (array, hstack, vstack, dstack, zeros, linspace, sqrt, cos,
                   sin, cosh, sinh, eye, dot, empty)
from functools import partial
from . import const
from ..visualize.plot import pltsim_quadscan
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
    # sig_u = eps_u*beta_u + (delta_E*D)^2
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
