#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' modul of accpy.lattices containig BESSY II lattices
author:
    Felix Kramer
'''
from __future__ import division
from ..simulate import const
from numpy import array, zeros, append, pi, sqrt
from .elements import drift, uniformdipole, rfquad, afquad, line, diagnostic


def lattice(latticename):
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
    if latticename == 'bessy2injectionline':
        '''===== general ====='''
        closed = False
        particle = 'electron'
        E = 50e6
        I = 5e-3
        '''===== starting twiss parameters ====='''
        beta_x = 1
        alph_x = 1
        gamm_x = (1+alph_x**2)/beta_x
        beta_y = 1
        alph_y = 1
        gamm_y = (1+alph_y**2)/beta_y
        disp_x = 0
        dipr_x = 0
        xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
        ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
        xdisp = array([[disp_x], [dipr_x], [1]])
        '''===== drifts ====='''
        D01 = drift(1)
        D02 = drift(1)
        D03 = drift(1)
        D04 = drift(1)
        D05 = drift(1)
        D06 = drift(1)
        D07 = drift(1)
        D08 = drift(1)
        D09 = drift(1)
        D10 = drift(1)
        D11 = drift(1)
        D12 = drift(1)
        D13 = drift(1)
        D14 = drift(1)
        D15 = drift(1)
        D16 = drift(1)
        D17 = drift(1)
        D18 = drift(1)
        D19 = drift(1)
        D20 = drift(1)
        D21 = drift(1)
        D22 = drift(1)
        D23 = drift(1)
        D24 = drift(1)
        D25 = drift(1)
        D26 = drift(1)
        D27 = drift(1)
        D28 = drift(1)
        D29 = drift(1)
        D30 = drift(1)
        D31 = drift(1)
        D32 = drift(1)
        D33 = drift(1)
        '''===== fluorescent screens ====='''
        F01 = F02 = F03 = F04 = F05 = F06 = diagnostic()
        diagnostics = ['FOMZ1LI', 'FOMZ2LI', 'FOMZ2I', 'FOMZ3I', 'FOMZ5I', 'FOMZ6I']
        '''===== dipoles ====='''
        B01 = uniformdipole(.143, .143/.4373, 0, .5*.03)  # B1P1LI 81.25 A   rbend
        B02 = uniformdipole(.143, .143/.43464, 0, .5*.03)  # B1P2LI 80.63 A   rbend
        B03 = uniformdipole(.1, .1/.02322, 0, .5*.03)  # B2PI 10.34 A         sbend
        B04 = uniformdipole(.1, .1/.02322, 0, .5*.03)  # HB2PI 0 A / B2PI 10.34 A ?  sbend
        '''===== quadrupoles ====='''
        q = const.qe
        E0 = const.Ee/q
        pc = sqrt(E**2-E0**2)*q
        p = pc/const.cl
        R = p/q         # beam rigidity R = Bρ = p/q = 5.73730218421
        i2k1 = lambda I: (6.87/2.55 * I)/R  # gradient@2.55A = 6.87 T/m
        i2k2 = lambda I: (6.86/2.55 * I)/R
        i2k3 = lambda I: (1.5 * I)/R
        Q01 = afquad(.119, i2k1(-0.562000))  # Q1PLI  +
        Q02 = afquad(.119, i2k1(-0.115052))  # Q2PLI  +
        Q03 = rfquad(.119, i2k1(+0.648000))  # Q3PLI  -
        Q04 = afquad(.220, i2k2(-1.118000))  # Q4PLI  +
        Q05 = rfquad(.119, i2k1(+0.454800))  # Q5PLI  -
        Q06 = afquad(.119, i2k1(-1.479440))  # Q6PLI  +
        Q07 = rfquad(.119, i2k1(+1.398218))  # Q7PLI  -
        Q08 = rfquad(.100, i2k3(+3.661584))  # Q2P2I  -
        Q09 = rfquad(.100, i2k3(+2.500517))  # Q1P2I  -
        Q10 = rfquad(.100, i2k3(+0.000000))  # Q4PI   -
        Q11 = afquad(.100, i2k3(-0.470000))  # Q5PI   +
        Q12 = rfquad(.100, i2k3(+0.000000))  # Q4PI   -
        Q13 = rfquad(.100, i2k3(+0.600000))  # Q6PI   -
        Q14 = afquad(.100, i2k3(-0.899900))  # Q7PI   +
        Q15 = afquad(.100, i2k3(-1.170000))  # Q8PI   +
        Q16 = afquad(.100, i2k3(-0.000000))  # Q9PI   +
        Q17 = afquad(.100, i2k3(-0.676000))  # Q10PI  +
        Q18 = afquad(.100, i2k3(-0.000000))  # Q11PI  +
        Q19 = afquad(.100, i2k3(-0.000000))  # Q12PI  +
        Q20 = afquad(.100, i2k3(-0.000000))  # Q13PI  -
        Q21 = afquad(.100, i2k3(-0.000000))  # Q14PI  +
        Q22 = afquad(.100, i2k3(-0.000000))  # Q15PI  +
        Q23 = afquad(.100, i2k3(-0.000000))  # Q16PI  +
        '''===== unit cell ====='''
        UC = line(D01, Q01, D02, Q02, D03, Q03, D04, B01, D05, Q04,
                  D06, F01, D07, B02, D08, Q05, D09, Q06, D10, Q07,
                  D10, F02, D11, Q08, D12, Q09, D13, F03, D14, Q10,
                  D15, Q11, D16, Q12, D17, F04, D18, Q13, D19, Q14,
                  D20, Q15, D21, Q16, D22, Q17, D23, B03, D24, Q18,
                  D25, Q19, D26, Q20, D27, B04, D28, Q21, D29, F05,
                  D30, Q22, D31, Q23, D32, F06, D33)
    elif latticename == 'bessy2booster':
        '''===== general ====='''
        closed = True
        particle = 'electron'
        E = 1.72e9
        HF_f = 499.667e6
        HF_V = 750e3
        I = 5e-3
        N_UC = 8        # number of unit cells
        '''===== drifts ====='''
        D1 = drift(1.160)
        D2 = drift(.23)
        D3 = drift(1.1607)
        '''===== diagnostics ====='''
        diagnostics = []
        '''===== dipoles ====='''
        LD = 2.6193         # orbit length of dipole
        UD = N_UC*LD*2      # total orbit length of all dipoles
        rho = UD/2/pi       # bending radius
        phi = LD/2/rho      # edge angle of dipole
        D_UC = 2            # dipoles per unit cell
        g = 36e-3           # dipoles gap
        K = 0.7             # ~0.7 for Rogowski pole
        B = uniformdipole(LD, rho, phi, K*g)
        '''===== quadrupoles ====='''
        QF = rfquad(.3, 2.1082)
        QD = afquad(.3, -1.4658)
        '''===== unit cell ====='''
        UC = line(D3, QD, D2, B, D2, QF, D1, D1, QD, D2, B, D2, QF, D3)
    elif latticename == 'bessy2transfer':
        '''===== general ====='''
        closed = False
        particle = 'electron'
        E = 1.72e9
        I = 5e-3
        '''===== starting twiss parameters ====='''
        beta_x = 10.101
        alph_x = 3.717
        gamm_x = (1+alph_x**2)/beta_x
        beta_y = 2.717
        alph_y = -1.115
        gamm_y = (1+alph_y**2)/beta_y
        disp_x = .847
        dipr_x = -.208
        xtwiss = array([[beta_x, -alph_x], [-alph_x, gamm_x]])
        ytwiss = array([[beta_y, -alph_y], [-alph_y, gamm_y]])
        xdisp = array([[disp_x], [dipr_x], [1]])
        '''===== drifts ====='''
        D01 = drift(.075)
        D02 = drift(.075)
        D03 = drift(.257672)
        D04 = drift(.065 + .955675 + .12 + .093 + .258)
        D05 = drift(.3)
        D06 = drift(.3)
        D07 = drift(.205625)
        D08 = drift(.09 + .111875 + .1225 + .275)
        D09 = drift(.288478)
        D10 = drift(.07 + .211522 - .07)
        D11 = drift(.08 + .12 + .1)
        D12 = drift(3.1415138 + .123 + 1.676 + .123 + .201)
        D13 = drift(.07 + .141)
        D14 = drift(.3)
        D15 = drift(.3)
        D16 = drift(.637 + .123 + .125 + .125)
        D17 = drift(.06 + .23)
        D18 = drift(.25)
        D19 = drift(.08 + .12)
        D20 = drift(.085 + .115)
        D21 = drift(.095 + .10185)
        D22 = drift(.12)
        D23 = drift(.275 + .325 + .22004)
        D24 = drift(.15082)
        D25 = drift(.09)
        '''===== diagnostics ====='''
        diagnos = append([[7]], zeros([5, 1]), 0)
        F1T = F2T = F3T = F4T = F5T = F6T = F7T = F8T = diagnos
        diagnostics = ['FOMZ1T', 'FOMZ2T', 'FOMZ3T = Halo Colimator',
                       'FOMZ4T = Collimator 1', 'FOMZ5T = Collimator 2',
                       'FOMZ6T', 'FOMZ7T', 'FOMZ8T']
        '''===== dipoles ====='''
        # pole shoe form factor K (close to rectangular)
        # dipoles gap g
        LD = .7; rho = -LD/2/pi*360/6.13; phi = LD/2/rho; g = 15e-3; K = 0.5
        DP1 = uniformdipole(LD, rho, phi, K*g)
        LD = .7; rho = -LD/2/pi*360/4.38; phi = LD/2/rho; g = 15e-3; K = 0.5
        DP2 = uniformdipole(LD, rho, phi, K*g)
        LD = 1.777792; rho = LD/2/pi*360/22; phi = LD/2/rho; g = 30e-3; K = 0.5
        DP3 = uniformdipole(LD, rho, phi, K*g)
        DP4 = uniformdipole(LD, rho, phi, K*g)
        LD = 1.020001; rho = LD/2/pi*360/7.66; phi = LD/2/rho; g = 15e-3; K = 0.5
        DP5 = uniformdipole(LD, rho, phi, K*g)
        LD = .555; rho = LD/2/pi*360/3.8; phi = LD/2/rho; g = 15e-3; K = 0.5
        DP6 = uniformdipole(LD, rho, phi, K*g)
        DP7 = uniformdipole(LD, rho, phi, K*g)
        '''===== quadrupoles ====='''
        q = const.qe
        E0 = const.Ee/q
        pc = sqrt(E**2-E0**2)*q
        p = pc/const.cl
        R = p/q         # beam rigidity R = Bρ = p/q = 5.73730218421
        i2kl = lambda i: (.265410*i-.765828e-6*i**3-.239385)/R
        i2ks = lambda i: (.266015*i-.829333e-6*i**3-.233316)/R
        Q01 = afquad(.25, -i2kl(41.3949))
        Q02 = afquad(.25, -0)
        Q03 = rfquad(.25, i2kl(37.5109))
        Q04 = rfquad(.25, i2kl(23.4929))
        Q05 = afquad(.25, -i2kl(30.4454))
        Q07 = afquad(.25, -0)
        Q08 = rfquad(.25, i2kl(68.0279))
        Q09 = afquad(.25, -i2kl(55.5348))
        Q10 = afquad(.2, -i2ks(43.1419))
        Q11 = afquad(.2, -0)
        Q12 = rfquad(.2, i2ks(53.3176))
        '''===== unit cell ====='''
                # 00   01   02   03   04   05   06   07   08   09
        UC = line(DP1, D01, F1T, D02, DP2, D03, F2T, D04, Q01, D05,
                # 10   11   12   13   14   15   16   17   18   19
                  Q02, D06, Q03, D07, F3T, D08, DP3, D09, F4T, D10,
                # 20   21   22   23   24   25   26   27   28   29
                  Q04, D11, Q05, D12, F5T, D13, Q07, D14, Q08, D15,
                # 30   31   32   33   34   35   36   37   38   39
                  Q09, D16, F6T, D17, DP4, D18, Q10, D19, Q11, D20,
                # 40   41   42   43   44   45   46   47   48   49
                  Q12, D21, F7T, D22, DP5, D23, F8T, D24, DP6, D25,
                # 50
                  DP7)
    elif latticename == 'bessy2ring':   # from BESSY II parameter list (1995)
        '''===== general ====='''
        closed = True
        particle = 'electron'
        E = 1.72e9
        HF_f = 499.667e6
        HF_V = 750e3
        I = 300e-3
        N_UC = 8        # number of unit cells
        '''===== drifts ====='''
        DD = drift(2.806)
        DT = drift(2.453)
        D1 = drift(.153)
        D2 = drift(.42)
        D3 = drift(.307)
        D4 = drift(.288)
        D5 = drift(.160)
        '''===== diagnostics ====='''
        diagnostics = []
        '''===== dipoles ====='''
        LD = 0.855          # orbit length of dipole
        D_UC = 4            # dipoles per unit cell
        UD = N_UC*D_UC*LD   # total orbit length of all dipoles
        rho = UD/2/pi       # bending radius
        phi = LD/2/rho      # edge angle of dipole
        g = 30e-3           # dipoles gap
        K = 0.5             # ~0.7 for Rogowski pole
        BB = uniformdipole(LD, rho, phi, K*g)
        '''===== quadrupoles ====='''
        Q1 = rfquad(.5, 1.40503)
        Q2 = afquad(.25, -2.01497)
        Q3 = afquad(.2, -1.89757)
        Q4 = rfquad(.25, 2.4519)
        Q1T = rfquad(.5, 2.62081)
        Q2T = afquad(.25, -2.46319)
        Q5T = afquad(.2, -2.6)
        '''===== sextupoles (as drifts) ====='''
        S1 = drift(0.16)
        S2 = drift(0.16)
        S3 = drift(0.16)
        S4 = drift(0.105)
        S5 = drift(0.16)
        S6 = drift(0.16)
        '''===== unit cell ====='''
        ACHR = line(BB, D2, Q3, D3, S3, D4, Q4, D5, S4)
        LDOU = line(DD, S1, D1, Q1, D1, S2, D1, Q2, D2)
        LTRI = line(D2, Q2T, D1, S6, D1, Q1T, D1, S5, D1, Q5T, DT)
        HALF = line(LDOU, ACHR, ACHR[:, ::-1], LTRI)
        UC = line(HALF, HALF[:, ::-1])
    if closed:
        xtwiss = ytwiss = xdisp = None
    else:
        HF_f = HF_V = None
        N_UC = 1
    return (closed, particle, E, I, UC, diagnostics,    # always
            N_UC, HF_f, HF_V,                           # closed lattice
            xtwiss, ytwiss, xdisp)                      # open lattice
