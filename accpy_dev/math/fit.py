# -*- coding: utf-8 -*-
''' accpy.math.gauss
author:     felix.kramer(at)physik.hu-berlin.de
'''
from scipy.optimize import leastsq
from numpy import sqrt, diag


def CurveFit(Fun, x, ydata, guess, DFun=None):
    def Residuals(pars, args):
        x, ydata = args
        return Fun(pars, x) - ydata
    if DFun:
        popt, pcov, infodict, errmsg, success = leastsq(Residuals, guess, args=([x, ydata]), Dfun=DFun, full_output=1,
                                                       col_deriv=1)
    else:
        popt, pcov, infodict, errmsg, success = leastsq(Residuals, guess, args=([x, ydata]), full_output=1)
    rchi2 = sum((Residuals(popt, (x, ydata))**2))/(len(ydata) - len(guess))
    perr  = sqrt(diag(pcov * rchi2))
    return popt, perr


def SurfaceFit(Fun, x, ydata, guess, DFun=None):
    def Residuals(pars, args):
        x, ydata = args
        return Fun(pars, x) - ydata
    if DFun:
        popt, pcov, infodict, errmsg, success = leastsq(Residuals, guess, args=([x, ydata]), Dfun=DFun, full_output=1,
                                                       col_deriv=1)
    else:
        popt, pcov, infodict, errmsg, success = leastsq(Residuals, guess, args=([x, ydata]), full_output=1)
    rchi2 = sum((Residuals(popt, (x, ydata))**2))/(len(ydata) - len(guess))
    perr  = sqrt(diag(pcov * rchi2))
    return popt, perr