# -*- coding: utf-8 -*-
''' accpy.lattices.elements
author:
    Felix Kramer

functions return elementvec = ev with:
ev[0] = type of element:
    0 = drift:
        ev[1] = length of drift
    1 = dipole
        ev[1] = length of dipole
        ev[2] = bending radius
    2 = edge
        ev[1] = 0 (length of edge)
        ev[2] = bending radius
        ev[3] = edge angle
        ev[4] = 0
        ev[5] = pole shoe form factor K * pole shoe gap g
    3/4 = radial focussing/defocussing quad
        ev[1] = length of quad
        ev[4] = quadrupole strength
    5 = rotator (Rskew=Rrot(-alpha*RQ*Rrot(alpha)))
        ev[1] = 0 (length of rotator)
        ev[2] = twistangle of rotator
    6 = solenoid
        ev[1] = length of solenoid
        ev[2] = twistangle of solenoid
    7 = diagnostic (identity matrix for R)
        ev[1] = 0 (length of diagnostic)
    8 = gradientdipole
        ev[1] = length of dipole = 2 * virtual quad length
        ev[2] = bending radius
        ev[3] = virtual quad strength
'''
from __future__ import division
from numpy import zeros, concatenate


def drift(length):
    ev = zeros([6, 1])
    # ev[0] = 0
    ev[1] = length
    return ev


def uniformdipole(length, radius, edgeangle=0, Kg=0):
    ev = zeros([6, 1])
    ev[0] = 1
    ev[1] = length
    ev[2] = radius
    if edgeangle != 0:
        edge = anglededge(radius, edgeangle, Kg)
        ev = concatenate((edge, ev, edge), 1)
    return ev


def anglededge(radius, edgeangle, Kg):
    '''
    dipoles gap g
    pole shoe form factor K:
        ~ 0.5 rectangular pole shoe
        ~ 0.7 Rogowski pole shoe
    '''
    ev = zeros([6, 1])
    ev[0] = 2
    # ev[1] = 0
    ev[2] = radius
    ev[3] = edgeangle
    # ev[4] = 0
    ev[5] = Kg
    return ev


def rfquad(length, strength):
    ev = zeros([6, 1])
    ev[0] = 3
    ev[1] = length
    ev[4] = strength
    return ev


def afquad(length, strength):
    ev = zeros([6, 1])
    ev[0] = 4
    ev[1] = length
    ev[4] = strength
    return ev


def diagnostic():
    ev = zeros([6, 1])
    ev[0] = 7
    return ev


def line(*elements):
    return concatenate((elements), 1)
