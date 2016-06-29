#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.lattices.reader
author:     felix.kramer(at)physik.hu-berlin.de
'''
from os import listdir
from numpy import array
header = '''
from __future__ import division
from numpy import array, zeros, append, sqrt
from .elements import drift, uniformdipole, rfquad, afquad, line, diagnostic
from ..simulate import const
diagnostics = []
'''


def latt2py(name, closed):
    '''lattice(latticename)
    input:
        latticename:    bessy2booster
                        bessy2transfer
    return:
        particle...particle type
        E       ...particle energy
        UC      ...unit cell with elementvectors:
    for closed lattices:
        rho     ...bending radius
        UD      ...total orbit length of all dipoles
        D_UC    ...dipoles per unit cell
        N_UC    ...number of unit cells
    for open lattices:
        xtwiss  ...initial radial twiss matrix
        ytwiss  ...initial axial twiss matrix
        xdisp   ...initial radial dispersion vector
        ydisp   ...initial axial dispersion vector
    '''
    l = locals()
    if closed:
        filepath = './accpy/lattices/closed/'+name+'.latt'
        code = header+open(filepath).read()
        exec(code)
        particle = l['particle']
        energy = l['energy']
        I = l['I']
        UC = l['UC']
        diagnostics = l['diagnostics']
        N_UC = l['N_UC']
        HF_f = l['HF_f']
        HF_V = l['HF_V']
        return (particle, energy, I, UC, diagnostics, N_UC,
                HF_f, HF_V)
    else:
        filepath = './accpy/lattices/open/'+name+'.latt'
        code = header+open(filepath).read()
        exec(code)
        particle = l['particle']
        energy = l['energy']
        I = l['I']
        UC = l['UC']
        diagnostics = l['diagnostics']
        N_UC = l['N_UC']
        beta_x = l['beta_x']
        alph_x = l['alph_x']
        gamm_x = l['gamm_x']
        beta_y = l['beta_y']
        alph_y = l['alph_y']
        gamm_y = l['gamm_y']
        disp_x = l['disp_x']
        dipr_x = l['dipr_x']
        xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
        ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
        xdisp = array([[disp_x], [dipr_x], [1]])
        N_UC = 1
        return (particle, energy, I, UC, diagnostics, N_UC,
                xtwiss, ytwiss, xdisp)



def lattlist():
    closedlatts = [lattfile[:-5] for lattfile in listdir('./accpy/lattices/closed') if lattfile[-5:] == '.latt']
    openlatts = [lattfile[:-5] for lattfile in listdir('./accpy/lattices/open') if lattfile[-5:] == '.latt']
    return closedlatts, openlatts
