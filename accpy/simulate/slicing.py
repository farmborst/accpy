#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' modul of LSD for slicing of unit cell
author:
    Felix Kramer
version:
    Created         09.11.2015
    Last Update     09.11.2015
'''
from numpy import hstack, zeros, cumsum


def cellslice(UC, P_UC, slicing):
    P_UCS = 0               # points in sliced unit cell
    UCS = zeros((6, 1))     # sliced unit cell
    for i in range(P_UC):
        if UC[0, i] in (2, 5, 7):  # no slicing of edges, rotators, diagnostics
            UCS = hstack((UCS, UC[:, i].reshape(6, 1)))
            P_UCS += 1
        else:
            UCS = hstack((UCS, UC[:, i].reshape(6, 1).repeat(slicing, 1)))
            P_UCS += slicing
    UCS[1, :] = UCS[1, :]/slicing
    s = cumsum(UCS[1, :])
    return s, UCS, P_UCS
