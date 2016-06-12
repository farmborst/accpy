#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.simulate.tracking
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (array, empty, dstack, hstack, vstack, dot, zeros,
                   sqrt, sign, trace)
from multiprocessing import cpu_count, Pool
from functools import partial


def initialtwiss(M):
    MX = M[:2, :2]
    MY = M[2:4, 2:4]
    MD = vstack((M[:2, [0, 1, 5]], [0, 0, 1]))
    # starting values of radial twiss matrix
    tmpx = sign(MX[0, 1])*sqrt(4-(trace(MX))**2)
    beta_x = 2*MX[0, 1]/tmpx
    alph_x = (MX[0, 0]-MX[1, 1])/tmpx
    gamm_x = (1+alph_x**2)/beta_x
    xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
    # starting values of axial twiss matrix
    tmpy = sign(MY[0, 1])*sqrt(4-(trace(MY))**2)
    beta_y = 2*MY[0, 1]/tmpy
    alph_y = (MY[0, 0]-MY[1, 1])/tmpy
    gamm_y = (1+alph_y**2)/beta_y
    ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
    # starting values of radial dispersion
    dipr_x = ((MD[1, 0]*MD[0, 2])+MD[1, 2]*(1-MD[0, 0]))/(2-MD[0, 0]-MD[1, 1])
    disp_x = (MD[0, 1]*dipr_x+MD[0, 2])/(1-MD[0, 0])
    xdisp = array([[disp_x], [dipr_x], [1]])
    return xtwiss, ytwiss, xdisp


# neglecting coupling and axial dispersion
def tracktwiss2(R, P_UCS, clos, xtwiss, ytwiss, disper):
    xtwiss = dstack((xtwiss, empty([2, 2, P_UCS])))
    ytwiss = dstack((ytwiss, empty([2, 2, P_UCS])))
    disper = hstack((disper, empty([3, P_UCS])))
    for i in range(P_UCS):
        RX = R[:2, :2, i]
        RY = R[2:4, 2:4, i]
        RD = vstack((hstack((RX, R[:2, 5, i, None])), array([[0, 0, 1]])))
        xtwiss[:, :, i+1] = dot(dot(RX, xtwiss[:, :, i]), RX.T)
        ytwiss[:, :, i+1] = dot(dot(RY, ytwiss[:, :, i]), RY.T)
        disper[:, i+1] = dot(RD, disper[:, i])
    return xtwiss, ytwiss, disper


# including transverse, linear coupling and neglecting axial dispersion
def tracktwiss4(R, P_UCS, clos, xtwiss, ytwiss, disper):
    twiss = hstack((xtwiss, zeros((2, 2))))
    twiss = vstack((twiss, hstack((zeros((2, 2)), ytwiss))))
    twiss = dstack((twiss, empty([4, 4, P_UCS])))
    disper = hstack((disper, empty([3, P_UCS])))
    for i in range(P_UCS):
        RXY = R[:4, :4, i]
        RD = vstack((R[:2, [0, 1, 5], i], array([[0, 0, 1]])))
        twiss[:, :, i+1] = dot(dot(RXY, twiss[:, :, i]), RXY.T)
        disper[:, i+1] = dot(RD, disper[:, i])
    xtwiss = twiss[:2, :2, :]
    ytwiss = twiss[2:4, 2:4, :]
    xytwiss = twiss[:2, 2:4, :]
    return xtwiss, ytwiss, disper, xytwiss


# including all linear coupling and dispersion (not ready...)
def tracktwiss6(R, P_UCS, clos, xtwiss, ytwiss, disperx,
                dispery=zeros((3, 1)), xytwiss=zeros((2, 2))):
    def twiss(disper):
        disper = disper.flatten()
        beta = disper[0]
        alpha = disperx[1]/2
        gamma = (1+alpha**2)/beta
        twiss = array([[beta, -alpha], [-alpha, gamma]])
        return twiss
    xztwiss = twiss(disperx)
    yztwiss = twiss(dispery)
    ztwiss = zeros((2, 2))  # clearly not!
    beta1 = hstack((xtwiss, xytwiss, xztwiss))
    beta2 = hstack((xytwiss, ytwiss, yztwiss))
    beta3 = hstack((xztwiss, yztwiss, ztwiss))
    beta = vstack((beta1, beta2, beta3))
    beta = dstack((beta, empty([6, 6, P_UCS])))
    for i in range(P_UCS):
        Ri = R[:, :, i]
        beta[:, :, i+1] = dot(dot(Ri, beta[:, :, i]), Ri.T)
    xtwiss = beta[:2, :2, :]
    ytwiss = beta[2:4, 2:4, :]
    disper = beta[4:6, :2, :]
    disper = array([disper[0, 0, :], disper[0, 1, :], zeros(P_UCS)])
    return xtwiss, ytwiss, disper


def trackpart(X, R, P_UCS, points):
    for i in range(points):
        X[:, i+1] = dot(R[:, :, i % P_UCS], X[:, i])
    return X


def trackparts(R_UCS, N_UC, X0, rounds):
    P_UCS = R_UCS.shape[2]
    points = P_UCS*N_UC*rounds
    # parallelized computation over particles
    Ncore = cpu_count()
    trackpartN = partial(trackpart, R=R_UCS, P_UCS=P_UCS, points=points)
    pool = Pool(Ncore)
    X = pool.map(trackpartN, X0)
    pool.close()
#    pool.join()
    return X
