# -*- coding: utf-8 -*-
''' accpy.simulate.rmatrices
author:     felix.kramer(at)physik.hu-berlin.de
'''
from numpy import sin, cos, tan, sinh, cosh, sqrt, eye, array, dot, empty


def rmatrix(elementvec, gamma):
    '''rmatrix(elementvec, gamma) returns the 6x6 R-matrix
    input:
        elementvec = ev of shape 5x1 with:
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
    return:
        6x6 R-matrix for tracking beam matrix or 6d particle vectors
    notice:
        RX = R[:2,:2]
        RD = [[R[:2,:2],R[:2,5]],[0,0,1]]
    '''

    def diagnostic():
        R = eye(6)
        return R

    def drift(ev, y):
        L = ev[1]    # length of element
        R = eye(6)
        R[0, 1] = R[2, 3] = L
        R[4, 5] = L/(y**2)
        return R

    def dipole(ev, y):
        L = ev[1]       # length
        rho = ev[2]     # bending radius
        alpha = L/rho
        cosalpha = cos(alpha)
        sinalpha = sin(alpha)
        R = array([
        [cosalpha,      rho*sinalpha,       0., 0., 0., rho*(1-cosalpha)],
        [-sinalpha/rho, cosalpha,           0., 0., 0., sinalpha],
        [0.,            0.,                 1., L,  0., 0],
        [0.,            0.,                 0., 1., 0., 0],
        [-sinalpha,     -rho*(1.-cosalpha), 0., 0., 1., rho*(alpha/(gamma**2)-alpha+sinalpha)],
        [0.,            0.,                 0., 0., 0., 1]
        ])
        return R

    def edge(ev):
        rho = ev[2]     # bending radius
        phi = ev[3]     # edge angle
        # distance between pole shoes g * pole shoe form faktor K
        gK = ev[5]      # K is ~0.45 for rectangular  and ~0.7 for Rogowski
        R = eye(6)
        tanphi = tan(phi)
        R[1, 0] = tanphi/rho
        if gK != 0:
            cosphi = cos(phi)
            sinphi = sin(phi)
            # Hinterberger 4.79 (exakt)
            R[3, 2] = -(tanphi-gK/rho*(1+(sinphi)**2/(cosphi**3)))/rho
            # Madx and Chao:
            # R[3,2] = -(tan(phi-gK/rho*(1+sinphi**2)/cosphi))/rho
        else:
            R[3, 2] = -tanphi/rho
        return R

    def quadrf(ev, y):
        L = ev[1]      # length
        k = ev[4]      # quadrupole strength
        if k == 0:
            R = drift(ev, y)
        else:
            wrzlk = sqrt(abs(k))
            Omega = wrzlk*L
            coshom = cosh(Omega)
            sinhom = sinh(Omega)
            cosom = cos(Omega)
            sinom = sin(Omega)
            R = array([
            [cosom,        sinom/wrzlk, 0.,           0.,           0., 0.],
            [-wrzlk*sinom, cosom,       0.,           0.,           0., 0.],
            [0.,            0.,         coshom,       sinhom/wrzlk, 0., 0.],
            [0.,            0.,         wrzlk*sinhom, coshom,       0., 0.],
            [0.,            0.,         0.,           0.,           1., L/(gamma**2)],
            [0.,            0.,         0.,           0.,           0., 1.]
            ])
        return R

    def quadaf(ev, y):
        L = ev[1]      # length
        k = ev[4]      # quadrupole strength
        if k == 0:
            R = drift(ev, y)
        else:
            wrzlk = sqrt(abs(k))
            Omega = wrzlk*L
            coshom = cosh(Omega)
            sinhom = sinh(Omega)
            cosom = cos(Omega)
            sinom = sin(Omega)
            R = array([
            [coshom,       sinhom/wrzlk, 0.,           0.,          0., 0],
            [wrzlk*sinhom, coshom,       0.,           0.,          0., 0],
            [0.,           0.,           cosom,        sinom/wrzlk, 0., 0],
            [0.,           0.,           -wrzlk*sinom, cosom,       0., 0],
            [0.,           0.,           0.,           0.,          1., L/(gamma**2)],
            [0.,           0.,           0.,           0.,          0., 1.]
            ])
        return R

    def rotator(ev):
        alpha = ev[1]      # twistangle of rotator
        sinalpha = sin(alpha)
        cosalpha = cos(alpha)
        R = array([
        [cosalpha,  0.,        sinalpha, 0.,       0., 0.],
        [0.,        cosalpha,  0.,       sinalpha, 0., 0.],
        [-sinalpha, 0.,        cosalpha, 0.,       0., 0.],
        [0.,        -sinalpha, 0.,       cosalpha, 0., 0.],
        [0.,        0.,        0.,       0.,       1., 0.],
        [0.,        0.,        0.,       0.,       0., 1.]
        ])
        return R

    def solenoid(ev, y):
        L = ev[1]           # length
        alpha = ev[2]       # twistangle of solenoid
        K = (-alpha)/(2*L)
        C = cos(K*L)
        S = sin(K*L)
        SC = C*S
        C2 = C**2
        S2 = S**2
        R = array([
        [C2,    SC/K,  SC,    S2/K, 0., 0.],
        [-K*SC, C2,    -K*S2, SC,   0., 0.],
        [-SC,   -S2/K, C2,    SC/K, 0., 0.],
        [K*S2,  -SC,   -K*SC, C2,   0., 0.],
        [0.,    0.,    0.,    0.,   1., L/(gamma**2)],
        [0.,    0.,    0.,    0.,   0., 1.]
        ])
        return R

    def gradientdipole(ev, y):
        ev[1] /= 2
        D = dipole(ev, y)
        k = ev[4]
        Q = eye(6)
        if k > 0:
            Q = quadrf(ev, y)
        else:
            Q = quadaf(ev, y)
        R = dot(dot(D, Q), D)
        return R

    q = elementvec[4]
    typ = elementvec[0]     # type of element (dipole, quadrupole, ...)
    if typ == 0:
        rmatrix = drift(elementvec, gamma)
        return rmatrix
    elif typ == 1:
        rmatrix = dipole(elementvec, gamma)
        return rmatrix
    elif typ == 2:
        rmatrix = edge(elementvec)
        return rmatrix
    elif (q == 0 and typ == 3) or (q == 0 and typ == 4):
        rmatrix = drift(elementvec, gamma)
        return rmatrix
    elif typ == 3:
        rmatrix = quadrf(elementvec, gamma)
        return rmatrix
    elif typ == 4:
        rmatrix = quadaf(elementvec, gamma)
        return rmatrix
    elif typ == 5:
        rmatrix = rotator(elementvec)
        return rmatrix
    elif typ == 6:
        rmatrix = solenoid(elementvec, gamma)
        return rmatrix
    elif typ == 7:
        rmatrix = diagnostic()
        return rmatrix
    elif typ == 8:
        rmatrix = gradientdipole(elementvec, gamma)
        return rmatrix


def rxdmatrix(elementvec):
    '''rxdmatrix(elementvec, gamma)
    input:
        elementvec = ev with:
            ev[0] = type of element
            ev[1] = length of element
            ev[2] = bending radius of dipole
            ev[3] = edge angle
            ev[4] = quadrupole strength
            ev[5] = twistangle of rotator
            ev[6] = twistangle of solenoid
    return:
        2x2 RX-matrix for tracking horizontal twiss functions
        3x3 RD-matrix for tracking dispersion
    notice:
        RX = R[0:2,0:2]
        RD = [[R[0:2,0:2],R[0:2,5]],[0,0,1]]
    '''
    def drift(ev):
        L = ev[1]    # length of element
        R = eye(3)
        R[0, 1] = L
        return R

    def dipole(ev):
        L = ev[1]       # length
        rho = ev[2]     # bending radius
        alpha = L/rho
        cosalpha = cos(alpha)
        sinalpha = sin(alpha)
        R = array([
        [cosalpha,      rho*sinalpha, rho*(1-cosalpha)],
        [-sinalpha/rho, cosalpha,     sinalpha],
        [0.,            0.,           1.]
        ])
        return R

    def edge(ev):
        rho = ev[2]      # bending radius
        phi = ev[3]      # edge angle
        R = eye(3)
        R[1, 0] = tan(phi)/rho
        return R

    def quadrf(ev):
        L = ev[1]      # length
        k = ev[4]      # quadrupole strength
        wrzlk = sqrt(abs(k))
        Omega = wrzlk*L
        cosom = cos(Omega)
        sinom = sin(Omega)
        R = array([
        [cosom,        sinom/wrzlk, 0.],
        [-wrzlk*sinom, cosom,       0.],
        [0.,            0.,         1]
        ])
        return R

    def quadaf(ev):
        L = ev[1]      # length
        k = ev[4]      # quadrupole strength
        wrzlk = sqrt(abs(k))
        Omega = wrzlk*L
        coshom = cosh(Omega)
        sinhom = sinh(Omega)
        R = array([
        [coshom,       sinhom/wrzlk, 0.],
        [wrzlk*sinhom, coshom,       0.],
        [0.,           0.,           1]
        ])
        return R

    def rotator(ev):
        alpha = ev[5]      # twistangle of rotator
        cosalpha = cos(alpha)
        R = array([
        [cosalpha,  0.,        0.],
        [0.,        cosalpha,  0.],
        [0.,        0.,        1.]
        ])
        return R

    def solenoid(ev):
        L = ev[1]           # length
        alpha = ev[5]       # twistangle of solenoid
        K = (-alpha)/(2*L)
        C = cos(K*L)
        S = sin(K*L)
        SC = C*S
        C2 = C**2
        R = array([
        [C2,    SC/K, 0.],
        [-K*SC, C2,   0.],
        [0.,    0.,   1.]
        ])
        return R

    q = elementvec[4]
    typ = elementvec[0]     # type of element (dipole, quadrupole, ...)
    if typ == 0:
        rmatrix = drift(elementvec)
        return rmatrix
    elif typ == 1:
        rmatrix = dipole(elementvec)
        return rmatrix
    elif typ == 2:
        rmatrix = edge(elementvec)
        return rmatrix
    elif (q == 0 and typ == 3) or (q == 0 and typ == 4):
        rmatrix = drift(elementvec)
        return rmatrix
    elif typ == 3:
        rmatrix = quadrf(elementvec)
        return rmatrix
    elif typ == 4:
        rmatrix = quadaf(elementvec)
        return rmatrix
    elif typ == 5:
        rmatrix = rotator(elementvec)
        return rmatrix
    elif typ == 6:
        rmatrix = solenoid(elementvec)
        return rmatrix


def UCS2R(P_UCS, UCS, gamma):
    R = empty((6, 6, P_UCS))
    for i in range(P_UCS):
            R[:, :, i] = rmatrix(UCS[:, i], gamma)
    return R
