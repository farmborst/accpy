#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.math.ode
author:     felix.kramer(at)physik.hu-berlin.de
'''
from scipy.integrate import ode


def odeint(fun, t, y0, integrator=0, method=0, rtol=0, atol=1e-12):
    s = ode(fun)
    integrators = ['vode', 'zvode', 'lsoda', 'dopri5', 'dop853']
    methods = ['adams', 'bdf']
    s.set_integrator(integrators[0],
                     method=methods[0],
                     order=10,
                     rtol=rtol,
                     atol=atol,
                     with_jacobian=False)
    t0 = t[0]
    dt = t[1]-t0
    y = [y0]
    s.set_initial_value(y0, t0)
    while s.successful() and s.t < t[-1]:
        s.integrate(s.t+dt)
        y.append(s.y)
    return y
