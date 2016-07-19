# -*- coding: utf-8 -*-
''' accpy.lattices.reader
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from os import listdir
from numpy import array, zeros, append, sqrt
import numpy as np
from .elements import drift, uniformdipole, quad, line, diagnostic
from ..simulate import const


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
    try: # normal install
        path = './accpy/lattices/'
        listdir(path+'closed')
    except: # required for pyinstaller
        path = './lattices/'
    l = locals()
    diagnostics = []
    if closed:
        path += 'closed/'
        filepath = path+name+'.latt'
        code = open(filepath).read()
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
        path += 'open/'
        filepath = path+name+'.latt'
        code = open(filepath).read()
        exec(code)
        particle = l['particle']
        energy = l['energy']
        I = l['I']
        UC = l['UC']
        diagnostics = l['diagnostics']
        emit_x = l['emit_x']
        emit_y = l['emit_y']
        emit_s = l['emit_s']
        beta_x = l['beta_x']
        alph_x = l['alph_x']
        beta_y = l['beta_y']
        alph_y = l['alph_y']
        disp_x = l['disp_x']
        dipr_x = l['dipr_x']
        gamm_x = (1+alph_x**2)/beta_x
        gamm_y = (1+alph_y**2)/beta_y
        xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
        ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
        xdisp = array([[disp_x], [dipr_x], [1]])
        N_UC = 1
        return (particle, energy, I, UC, diagnostics, N_UC,
                xtwiss, ytwiss, xdisp, emit_x, emit_y, emit_s)


def txt2py(code, openclosed):
    l = locals()
    diagnostics = []
    exec(code)
    if openclosed == 'closed':
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
        particle = l['particle']
        energy = l['energy']
        I = l['I']
        UC = l['UC']
        diagnostics = l['diagnostics']
        beta_x = l['beta_x']
        alph_x = l['alph_x']
        beta_y = l['beta_y']
        alph_y = l['alph_y']
        disp_x = l['disp_x']
        dipr_x = l['dipr_x']
        gamm_x = (1+alph_x**2)/beta_x
        gamm_y = (1+alph_y**2)/beta_y
        xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
        ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
        xdisp = array([[disp_x], [dipr_x], [1]])
        N_UC = 1
        return (particle, energy, I, UC, diagnostics, N_UC,
                xtwiss, ytwiss, xdisp)



def latt2txt(name, closed):
    try:
        path = './accpy/lattices/'
        listdir(path+'closed')
    except:
        path = './lattices/'

    if closed:
        path += 'closed/'
    else:
        path += 'open/'
    filepath = path+name+'.latt'
    latticefile = open(filepath, mode='r')
    loadedfile = latticefile.read()
    latticefile.close()
    return loadedfile


def txt2latt(txt, name, openclosed):
    try:
        path = './accpy/lattices/'
        listdir(path+'closed')
    except:
        path = './lattices/'
    path += openclosed+'/'
    filepath = path+name+'.latt'
    latticefile = open(filepath, mode='w')
    latticefile.write(txt)
    latticefile.close()


def lattlist():
    try:
        path = './accpy/lattices/'
        listdir(path+'closed')
    except:
        path = './lattices/'
    closedlatts = [lattfile[:-5] for lattfile in listdir(path+'closed') if lattfile[-5:] == '.latt']
    openlatts = [lattfile[:-5] for lattfile in listdir(path+'open') if lattfile[-5:] == '.latt']
    return closedlatts, openlatts
