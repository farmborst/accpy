#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' modul of accpy.lattices containig BESSY II lattices
author:
    Felix Kramer
'''
from ..simulate import const
from numpy import array, zeros, append, pi, concatenate, sqrt, size


def lattice(latticename):
    '''lattice(latticename)
    input:
        latticename:    bessy2booster
                        bessy2transfer
    return:
        m       ...particle mass
        q       ...particle charge
        E0      ...particle rest energy
        E       ...particle energy
        UC      ...unit cell with elementvectors:
            elementvec = ev with:
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
        m = const.me
        q = const.qe
        E0 = const.Ee/const.qe
        E = 50e6
        I = 5e-3
        gamma = E/E0+1
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
        D01, D02, D03, D04, D05, D06, D07, D08, D09, D10, D11, D12, D13, D14, \
            D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25 = \
            [zeros([6, 1]) for i in range(25)]
        D01[1] = .075
        D02[1] = .075
        D03[1] = .257672
        D04[1] = .065 + .955675 + .12 + .093 + .258
        D05[1] = .3
        D06[1] = .3
        D07[1] = .205625
        D08[1] = .09 + .111875 + .1225 + .275
        D09[1] = .288478
        D10[1] = .07 + .211522 - .07
        D11[1] = .08 + .12 + .1
        D12[1] = 3.1415138 + .123 + 1.676 + .123 + .201
        D13[1] = .07 + .141
        D14[1] = .3
        D15[1] = .3
        D16[1] = .637 + .123 + .125 + .125
        D17[1] = .06 + .23
        D18[1] = .25
        D19[1] = .08 + .12
        D20[1] = .085 + .115
        D21[1] = .095 + .10185
        D22[1] = .12
        D23[1] = .275 + .325 + .22004
        D24[1] = .15082
        D25[1] = .09
        '''===== diagnostics ====='''
        diagnos = append([[7]], zeros([5, 1]), 0)
        F1T = F2T = F3T = F4T = F5T = F6T = F7T = F8T = diagnos
        diagnostics = ['FOMZ1T', 'FOMZ2T', 'FOMZ3T = Halo Colimator',
                       'FOMZ4T = Collimator 1', 'FOMZ5T = Collimator 2',
                       'FOMZ6T', 'FOMZ7T', 'FOMZ8T']
        '''===== dipoles ====='''
        # pole shoe form factor K (close to rectangular)
        # dipoles gap g
        dipole = append([[1]], zeros([5, 1]), 0)
        edge = append([[2]], zeros([5, 1]), 0)
        LD = .7; rho = -LD/2/pi*360/6.13; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP1 = concatenate((edge, dipole, edge), 1)
        LD = .7; rho = -LD/2/pi*360/4.38; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP2 = concatenate((edge, dipole, edge), 1)
        LD = 1.777792; rho = LD/2/pi*360/22; phi = LD/2/rho; g = 30e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP3 = concatenate((edge, dipole, edge), 1)
        DP4 = DP3
        LD = 1.020001; rho = LD/2/pi*360/7.66; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP5 = concatenate((edge, dipole, edge), 1)
        LD = .555; rho = LD/2/pi*360/3.8; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP6 = concatenate((edge, dipole, edge), 1)
        DP7 = DP6
        '''===== quadrupoles ====='''
        pc = sqrt(E**2-E0**2)*q
        p = pc/const.cl
        R = p/q         # beam rigidity R = Bρ = p/q = 5.73730218421
        i2kl = lambda i: (.265410*i-.765828e-6*i**3-.239385)/R
        i2ks = lambda i: (.266015*i-.829333e-6*i**3-.233316)/R
        Q01, Q02, Q03, Q04, Q05, Q07, Q08, Q09, Q10, Q11, Q12 = \
            [zeros([6, 1]) for i in range(11)]
        Q01[[0, 1, 4]] = array([[4], [.25], [-i2kl(41.3949)]])
        Q02[[0, 1, 4]] = array([[4], [.25], [-0]])
        Q03[[0, 1, 4]] = array([[3], [.25], [i2kl(37.5109)]])
        Q04[[0, 1, 4]] = array([[3], [.25], [i2kl(23.4929)]])
        Q05[[0, 1, 4]] = array([[4], [.25], [-i2kl(30.4454)]])
        Q07[[0, 1, 4]] = array([[4], [.25], [-0]])
        Q08[[0, 1, 4]] = array([[3], [.25], [i2kl(68.0279)]])
        Q09[[0, 1, 4]] = array([[4], [.25], [-i2kl(55.5348)]])
        Q10[[0, 1, 4]] = array([[4], [.2], [-i2ks(43.1419)]])
        Q11[[0, 1, 4]] = array([[4], [.2], [-0]])
        Q12[[0, 1, 4]] = array([[3], [.2], [i2ks(53.3176)]])
        '''===== unit cell ====='''
        UC = concatenate((DP1, D01, F1T, D02, DP2, D03, F2T, D04, Q01, D05,
                             Q02, D06, Q03, D07, F3T, D08, DP3, D09, F4T, D10,
                             Q04, D11, Q05, D12, F5T, D13, Q07, D14, Q08, D15,
                             Q09, D16, F6T, D17, DP4, D18, Q10, D19, Q11, D20,
                             Q12, D21, F7T, D22, DP5, D23, F8T, D24, DP6, D25,
                             DP7), 1)
        P_UC = size(UC, 1)        # nr of elements in unit cell
    elif latticename == 'bessy2booster':
        '''===== general ====='''
        closed = True
        m = const.me
        q = const.qe
        E0 = const.Ee/const.qe
        E = 1.72e9
        HF_f = 499.667e6
        HF_V = 750e3
        I = 5e-3
        gamma = E/E0+1
        N_UC = 8        # number of unit cells
        '''===== drifts ====='''
        D1, D2, D3 = [zeros([6, 1]) for i in range(3)]
        D1[1] = 1.160
        D2[1] = .23
        D3[1] = 1.1607
        '''===== diagnostics ====='''
        diagnostics = []
        '''===== dipoles ====='''
        LD = 2.6193         # orbit length of dipole
        UD = N_UC*LD*2      # total orbit length of all dipoles
        rho = UD/2/pi    # bending radius
        phi = LD/2/rho      # edge angle of dipole
        D_UC = 2            # dipoles per unit cell
        g = 36e-3           # dipoles gap
        K = 0.7             # ~0.7 for Rogowski pole
        dipole, edge = [zeros([6, 1]) for i in range(2)]
        dipole[[0, 1, 2]] = array([[1], [LD], [rho]])
        edge[[0, 2, 3, 5]] = array([[2], [rho], [phi], [g*K]])
        B = concatenate((edge, dipole, edge), 1)
        '''===== quadrupoles ====='''
        QF, QD = [zeros([6, 1]) for i in range(2)]
        QF[[0, 1, 4]] = array([[3], [.3], [2.1082]])
        QD[[0, 1, 4]] = array([[4], [.3], [-1.57808028]])    # -1.4658
        '''===== unit cell ====='''
        UC = concatenate((
                            D3, QD, D2, B, D2, QF, D1,
                            D1, QD, D2, B, D2, QF, D3
                            ), 1)
        P_UC = size(UC, 1)        # nr of elements in unit cell
    elif latticename == 'bessy2transfer':
        '''===== general ====='''
        closed = False
        m = const.me
        q = const.qe
        E0 = const.Ee/const.qe
        E = 1.72e9
        I = 5e-3
        gamma = E/E0+1
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
        D01, D02, D03, D04, D05, D06, D07, D08, D09, D10, D11, D12, D13, D14, \
            D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25 = \
            [zeros([6, 1]) for i in range(25)]
        D01[1] = .075
        D02[1] = .075
        D03[1] = .257672
        D04[1] = .065 + .955675 + .12 + .093 + .258
        D05[1] = .3
        D06[1] = .3
        D07[1] = .205625
        D08[1] = .09 + .111875 + .1225 + .275
        D09[1] = .288478
        D10[1] = .07 + .211522 - .07
        D11[1] = .08 + .12 + .1
        D12[1] = 3.1415138 + .123 + 1.676 + .123 + .201
        D13[1] = .07 + .141
        D14[1] = .3
        D15[1] = .3
        D16[1] = .637 + .123 + .125 + .125
        D17[1] = .06 + .23
        D18[1] = .25
        D19[1] = .08 + .12
        D20[1] = .085 + .115
        D21[1] = .095 + .10185
        D22[1] = .12
        D23[1] = .275 + .325 + .22004
        D24[1] = .15082
        D25[1] = .09
        '''===== diagnostics ====='''
        diagnos = append([[7]], zeros([5, 1]), 0)
        F1T = F2T = F3T = F4T = F5T = F6T = F7T = F8T = diagnos
        diagnostics = ['FOMZ1T', 'FOMZ2T', 'FOMZ3T = Halo Colimator',
                       'FOMZ4T = Collimator 1', 'FOMZ5T = Collimator 2',
                       'FOMZ6T', 'FOMZ7T', 'FOMZ8T']
        '''===== dipoles ====='''
        # pole shoe form factor K (close to rectangular)
        # dipoles gap g
        dipole = append([[1]], zeros([5, 1]), 0)
        edge = append([[2]], zeros([5, 1]), 0)
        LD = .7; rho = -LD/2/pi*360/6.13; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP1 = concatenate((edge, dipole, edge), 1)
        LD = .7; rho = -LD/2/pi*360/4.38; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP2 = concatenate((edge, dipole, edge), 1)
        LD = 1.777792; rho = LD/2/pi*360/22; phi = LD/2/rho; g = 30e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP3 = concatenate((edge, dipole, edge), 1)
        DP4 = DP3
        LD = 1.020001; rho = LD/2/pi*360/7.66; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP5 = concatenate((edge, dipole, edge), 1)
        LD = .555; rho = LD/2/pi*360/3.8; phi = LD/2/rho; g = 15e-3; K = 0.5
        dipole[[1, 2]] = array([[LD], [rho]])
        edge[[2, 3, 5]] = array([[rho], [phi], [g*K]])
        DP6 = concatenate((edge, dipole, edge), 1)
        DP7 = DP6
        '''===== quadrupoles ====='''
        pc = sqrt(E**2-E0**2)*q
        p = pc/const.cl
        R = p/q         # beam rigidity R = Bρ = p/q = 5.73730218421
        i2kl = lambda i: (.265410*i-.765828e-6*i**3-.239385)/R
        i2ks = lambda i: (.266015*i-.829333e-6*i**3-.233316)/R
        Q01, Q02, Q03, Q04, Q05, Q07, Q08, Q09, Q10, Q11, Q12 = \
            [zeros([6, 1]) for i in range(11)]
        Q01[[0, 1, 4]] = array([[4], [.25], [-i2kl(41.3949)]])
        Q02[[0, 1, 4]] = array([[4], [.25], [-0]])
        Q03[[0, 1, 4]] = array([[3], [.25], [i2kl(37.5109)]])
        Q04[[0, 1, 4]] = array([[3], [.25], [i2kl(23.4929)]])
        Q05[[0, 1, 4]] = array([[4], [.25], [-i2kl(30.4454)]])
        Q07[[0, 1, 4]] = array([[4], [.25], [-0]])
        Q08[[0, 1, 4]] = array([[3], [.25], [i2kl(68.0279)]])
        Q09[[0, 1, 4]] = array([[4], [.25], [-i2kl(55.5348)]])
        Q10[[0, 1, 4]] = array([[4], [.2], [-i2ks(43.1419)]])
        Q11[[0, 1, 4]] = array([[4], [.2], [-0]])
        Q12[[0, 1, 4]] = array([[3], [.2], [i2ks(53.3176)]])
        '''===== unit cell ====='''   #00   01   02   03   04   05   06   07   08   09
        UC = concatenate((DP1, D01, F1T, D02, DP2, D03, F2T, D04, Q01, D05,
                            #10   11   12   13   14   15   16   17   18   19
                             Q02, D06, Q03, D07, F3T, D08, DP3, D09, F4T, D10,
                            #20   21   22   23   24   25   26   27   28   29
                             Q04, D11, Q05, D12, F5T, D13, Q07, D14, Q08, D15,
                            #30   31   32   33   34   35   36   37   38   39
                             Q09, D16, F6T, D17, DP4, D18, Q10, D19, Q11, D20,
                            #40   41   42   43   44   45   46   47   48   49
                             Q12, D21, F7T, D22, DP5, D23, F8T, D24, DP6, D25,
                            #50
                             DP7), 1)
        P_UC = size(UC, 1)        # nr of elements in unit cell
    elif latticename == 'bessy2ring':   # from BESSY II parameter list (1995)
        '''===== general ====='''
        closed = True
        m = const.me
        q = const.qe
        E0 = const.Ee/const.qe
        E = 1.72e9
        HF_f = 499.667e6
        HF_V = 750e3
        I = 300e-3
        gamma = E/E0+1
        N_UC = 8        # number of unit cells
        '''===== drifts ====='''
        DD, DT, D1, D2, D3, D4, D5 = [zeros([6, 1]) for i in range(7)]
        DD[1] = 2.806
        DT[1] = 2.453
        D1[1] = .153
        D2[1] = .42
        D3[1] = .307
        D4[1] = .288
        D5[1] = .160
        '''===== diagnostics ====='''
        diagnostics = []
        '''===== dipoles ====='''
        LD = 0.855          # orbit length of dipole
        D_UC = 4            # dipoles per unit cell
        UD = N_UC*D_UC*LD   # total orbit length of all dipoles
        rho = UD/2/pi    # bending radius
        phi = LD/2/rho      # edge angle of dipole
        g = 30e-3           # dipoles gap
        K = 0.5             # ~0.7 for Rogowski pole
        dipole, edge = [zeros([6, 1]) for i in range(2)]
        dipole[[0, 1, 2]] = array([[1], [LD], [rho]])
        edge[[0, 2, 3, 5]] = array([[2], [rho], [phi], [g*K]])
        BB = concatenate((edge, dipole, edge), 1)
        '''===== quadrupoles ====='''
        Q1, Q2, Q3, Q4, Q1T, Q2T, Q5T = [zeros([6, 1]) for i in range(7)]
        Q1[[0, 1, 4]] = array([[3], [.5], [1.40503]])
        Q2[[0, 1, 4]] = array([[4], [.25], [-2.01497]])
        Q3[[0, 1, 4]] = array([[4], [.2], [-1.89757]])
        Q4[[0, 1, 4]] = array([[3], [.25], [2.4519]])
        Q1T[[0, 1, 4]] = array([[3], [.5], [2.62081]])
        Q2T[[0, 1, 4]] = array([[4], [.25], [-2.46319]])
        Q5T[[0, 1, 4]] = array([[4], [.2], [-2.6]])
        '''===== sextupoles (as drifts) ====='''
        S1, S2, S3, S4, S5, S6 = [zeros([6, 1]) for i in range(6)]
        S1[1] = 0.16
        S2[1] = 0.16
        S3[1] = 0.16
        S4[1] = 0.105
        S5[1] = 0.16
        S6[1] = 0.16
        '''===== unit cell ====='''
        ACHR = concatenate((BB, D2, Q3, D3, S3, D4, Q4, D5, S4), 1)
        LDOU = concatenate((DD, S1, D1, Q1, D1, S2, D1, Q2, D2), 1)
        LTRI = concatenate((D2, Q2T, D1, S6, D1, Q1T, D1, S5, D1, Q5T, DT), 1)
        HALF = concatenate((LDOU, ACHR, ACHR[:, ::-1], LTRI), 1)
        UC = concatenate((HALF, HALF[:, ::-1]), 1)
        P_UC = size(UC, 1)        # nr of elements in unit cell
    if closed:
        xtwiss = ytwiss = xdisp = None
    else:
        rho = UD = D_UC = N_UC = HF_f = HF_V = None
    return m, q, E0, E, I, gamma, UC, P_UC, diagnostics, closed, \
               rho, UD, D_UC, N_UC, \
               xtwiss, ytwiss, xdisp, HF_f, HF_V
