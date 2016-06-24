#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' accpy.math.ode
author:     felix.kramer(at)physik.hu-berlin.de
'''
from scipy.integrate import ode
from scipy.optimize import fsolve
from numpy import array


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
    y = array(y)
    return y


def euler_explicit(fun_yprime, t, y0, h):
    y = [y0]
    [y.append(y[i]+h*fun_yprime(y[i], t[i])) for i in range(len(t)-1)]
    return y


def euler_implicit(fun_yprime, t, y0, h):
    y = [y0]
    for i in range(len(t)-1):
        to_be_solved = lambda yp: y[i]+h*fun_yprime(yp, t[i+1])-yp
        y.append(fsolve(to_be_solved, y[i]))
    return y


def rungekutta2(fun_yprime, t, y0, h):
    hh = h/2.
    y = [y0]
    for i in range(len(t)-1):
        K1 = h*fun_yprime(y[i], t[i])
        K2 = h*fun_yprime(y[i] + K1/2., t[i] + hh)
        y.append(y[i]+K2)
    return y


def rungekutta4(fun_yprime, t, y0, h):
    hh = h/2.
    y = [y0]
    for i in range(len(t)-1):
        K1 = h*fun_yprime(y[i], t[i])
        K2 = h*fun_yprime(y[i] + K1/2., t[i] + hh)
        K3 = h*fun_yprime(y[i] + K2/2., t[i] + hh)
        K4 = h*fun_yprime(y[i] + K3, t[i] + h)
        y.append(y[i]+(K1+2*K2+2*K3+K4)/6.)
    return y
