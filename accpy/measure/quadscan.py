# -*- coding: utf-8 -*-
''' accpy.measure.quadscan
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
import numpy as np
from scipy import linalg
from scipy.optimize import curve_fit
from mypy import const, UC2R


def fun(x, a, b, c):
    return a*x*x + b*x + c


def dofit(xdata, ydata, betagamma, T, dE, D, QL, orientation, yerr=None,
          krange=[]):
    # make sure that k is positive k = abs(g)/Brho_0
    xdata = abs(xdata)
    if len(krange) != 0:
        # trim to given range
        xdata = xdata[krange]
        ydata = ydata[krange]
        yerr = yerr[krange]
    if yerr is None:
        popt, pcov = curve_fit(fun, xdata, ydata, sigma=None)
    else:
        popt, pcov = curve_fit(fun, xdata, ydata, sigma=yerr, absolute_sigma=True)
    fit = fun(xdata, *popt)
    perr = np.sqrt(np.diag(pcov))
    A = np.array([popt[0], perr[0]])
    B = np.array([popt[1], perr[1]])
    C = np.array([popt[2], perr[2]])
    eps, bet, alp, gam = twissfromparabola(A, B, C, T, dE, D, QL, orientation)
    epsn = betagamma*eps*1e6
    eps *= 1e9
    A *= 1e9
    B *= 1e9
    C *= 1e9
    string = (r'$function:\;y=ax^2+bx+c$''\n'
              r'$a=(%.2f \pm %.2f)\cdot10^{-9}$''\n'
              r'$b=(%.2f \pm %.2f)\cdot10^{-9}$''\n'
              r'$c=(%.2f \pm %.2f)\cdot10^{-9}$''\n'
              r'$\beta=%.2f \pm %.2f$''\n'
              r'$\alpha=%.2f \pm %.2f$''\n'
              r'$\gamma=%.2f \pm %.2f$''\n'
              r'$\epsilon=(%.2f \pm %.2f)\pi\, nm\, rad$''\n'
              r'$\epsilon^*=(%.2f \pm %.2f)\pi\, \mu m\, rad$'
              % (A[0], A[1], B[0], B[1], C[0], C[1],
                 bet[0], bet[1], alp[0], alp[1], gam[0], gam[1],
                 eps[0], eps[1], epsn[0], epsn[1]))
    return string, fit, xdata


def fomplot(fom, quad, fig, kquad, sig2x, sig2y, ks, sxs, sxse, sys, syse,
            lorby, T, dE, D, QL, orientation, rangex, rangey, upper = True):
    if orientation == 'radialfocus':
        o = [1, 2]
    elif orientation == 'axialfocus':
        o = [2, 1]
    strx1, fitx1, xfx1 = dofit(kquad, sig2x, lorby, T, dE, D[0], QL, o[0])
    stry1, fity1, xfy1 = dofit(kquad, sig2y, lorby, T, dE, D[1], QL, o[1])
    strx2, fitx2, xfx2 = dofit(ks, sxs, lorby, T, dE, D[0], QL, o[0], yerr=sxse, krange=rangex)
    stry2, fity2, xfy2 = dofit(ks, sys, lorby, T, dE, D[1], QL, o[1], yerr=syse, krange=rangey)
    if upper:
        ax = fig.add_subplot(2, 2, 1)
        ax.set_title(r'$\sigma_x^2$ on FOMZ%iT(Q%iPT)' % (fom, quad))
        ax.set_xlabel(r'quadrupole strength k / ($1/m^2$)')
        ax.set_ylabel(r'$\sigma_x^2$ / $(m^2)$')
    #    ax.text(0.1, 0.9, strx1, horizontalalignment='left',
    #            verticalalignment='top', transform=ax.transAxes)
        ax.plot(kquad, sig2x, 'xb', label='data')
    #    ax.plot(xfx1, fitx1, '-r', label='fit')
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
        ax = fig.add_subplot(2, 2, 2)
        ax.set_title(r'$\sigma_y^2$ on FOMZ%iT(Q%iPT)' % (fom, quad))
        ax.set_xlabel(r'quadrupole strength k / ($1/m^2$)')
        ax.set_ylabel(r'$\sigma_y^2$ / $(m^2)$')
    #    ax.text(0.1, 0.9, stry1, horizontalalignment='left',
    #            verticalalignment='top', transform=ax.transAxes)
        ax.plot(kquad, sig2y, 'xb', label='data')
    #    ax.plot(xfy1, fity1, '-r', label='fit')
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
    else:
        ax = fig.add_subplot(1, 2, 1)
    if upper:
        ax = fig.add_subplot(2, 2, 3)
    else:
        ax = fig.add_subplot(1, 2, 1)
    ax.set_title(r'$\sigma_x^2$ on FOMZ%iT(Q%iPT)' % (fom, quad))
    ax.set_xlabel(r'quadrupole strength k / ($1/m^2$)')
    ax.set_ylabel(r'$\sigma_x^2$ / $(m^2)$')
    ax.text(0.3, 0.9, strx2, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.errorbar(ks, sxs, yerr=sxse, marker='.', mfc='b', ecolor='b', mec='b',
                ls='None', label='data')
    ax.plot(xfx2, fitx2, '-r', label='fit')
    ax.set_xlim([0, round(max(xfx2))+1])
    leg = ax.legend(fancybox=True, loc='lower right')
    leg.get_frame().set_alpha(0.5)
    if upper:
        ax = fig.add_subplot(2, 2, 4)
    else:
        ax = fig.add_subplot(1, 2, 2)
    ax.set_title(r'$\sigma_y^2$ on FOMZ%iT(Q%iPT)' % (fom, quad))
    ax.set_xlabel(r'quadrupole strength k / ($1/m^2$)')
    ax.set_ylabel(r'$\sigma_y^2$ / $(m^2)$')
    ax.text(0.1, 0.9, stry2, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.errorbar(ks, sys, yerr=syse, marker='.', mfc='b', ecolor='b', mec='b',
                ls='None', label='data')
    ax.plot(xfy2, fity2, '-r', label='fit')
    ax.set_xlim([0, round(max(xfy2))+1])
    leg = ax.legend(fancybox=True, loc='lower right')
    leg.get_frame().set_alpha(0.5)


def lowess(x, y, f=1. / 3., iter=5):
    """lowess(x, y, f=2./3., iter=3) -> yest
    Lowess smoother: Robust locally weighted regression.
    The lowess function fits a nonparametric regression curve to a scatterplot.
    The arrays x and y contain an equal number of elements; each pair
    (x[i], y[i]) defines a data point in the scatterplot. The function returns
    the estimated (smooth) values of y.
    The smoothing span is given by f. A larger value for f will result in a
    smoother curve. The number of robustifying iterations is given by iter. The
    function will run faster with a smaller number of iterations.
    """
    n = len(x)
    r = int(np.ceil(f * n))
    h = [np.sort(np.abs(x - x[i]))[r] for i in range(n)]
    w = np.clip(np.abs((x[:, None] - x[None, :]) / h), 0.0, 1.0)
    w = (1 - w ** 3) ** 3
    yest = np.zeros(n)
    delta = np.ones(n)
    for iteration in range(iter):
        for i in range(n):
            weights = delta * w[:, i]
            b = np.array([np.sum(weights * y), np.sum(weights * y * x)])
            A = np.array([[np.sum(weights), np.sum(weights * x)],
                          [np.sum(weights * x), np.sum(weights * x * x)]])
            beta = linalg.solve(A, b)
            yest[i] = beta[0] + beta[1] * x[i]

        residuals = y - yest
        s = np.median(np.abs(residuals))
        delta = np.clip(residuals / (6.0 * s), -1, 1)
        delta = (1 - delta ** 2) ** 2
    return yest


def parabolafit(xdata, ydata, init, sigma):
    popt, pcov = curve_fit(fun, xdata, ydata, p0=init, sigma=sigma)
    perr = np.sqrt(np.diag(pcov))
    x = np.linspace(min(xdata), max(xdata), 1e3)
    fit = fun(x, *popt)
    A = np.array([popt[0], perr[0]])
    B = np.array([popt[1], perr[1]])
    C = np.array([popt[2], perr[2]])
    return x, fit, A, B, C


def twissfromparabola(A, B, C, T, dE, D, QL, orientation):
    # K is negative in the focussing plane!
    # 1. Case: Radial Focusing X and Axial Focusing Y
    if orientation == 1:
        B *= -2
    # 2. Case: Radial Focusing Y and Axial Focusing X
    elif orientation == 2:
        B *= 2
    # make lines shorter
    T11 = T[0, 0]
    T12 = T[0, 1]
    T12s = T12*T12
    # calculate sigma matrix elements
    sig11 = A[0]/QL/QL/T12s
    sig11er = A[1]/QL/QL/T12s
    sig12 = (B[0]/2/QL-T11*T12*sig11)/T12s
    sig12er = np.sqrt((B[1]/2/QL)**2-(T11*T12*sig11er)**2)/T12s
    sig22 = (C[0]-T11**2*sig11-2*T11*T12*sig12-(dE*D)**2)/T12s
    sig22er = np.sqrt(C[1]**2-(T11**2*sig11er)**2+(2*T11*T12*sig12er)**2)/T12s
    epsilon = np.sqrt(sig11*sig22-sig12**2)
    epser = np.sqrt((sig22*sig11er)**2+(sig11*sig22er)**2+(-2*sig12*sig12er)**2)/(2*epsilon)
    beta = sig11/epsilon
    betaer = np.sqrt((sig11er/epsilon)**2+(-sig11*epser/epsilon**2)**2)
    alpha = -sig12/epsilon
    alphaer = np.sqrt((-sig12er/epsilon)**2+(sig12*epser/epsilon**2)**2)
    gamma = sig22/epsilon
    gammaer = np.sqrt((sig22er/epsilon)**2+(-sig22*epser/epsilon**2)**2)
    epsilon = np.array([epsilon, epser])
    beta = np.array([beta, betaer])
    alpha = np.array([alpha, alphaer])
    gamma = np.array([gamma, gammaer])
    return epsilon, beta, alpha, gamma


def transferlineR(start, end, I1=41.3949, I2=0, I3=37.5109, I4=23.4929,
                  I5=30.4454, I7=0, I8=68.0279, I9=55.5348, I10=43.1419,
                  I11=0, I12=53.3176):
    '''===== general ====='''
    q = const.qe
    E0 = const.Ee/const.qe
    E = 1.72e9
    gamma = E/E0+1
    pc = np.sqrt(E**2-E0**2)*q
    p = pc/const.cl
    R = p/q         # beam rigidity R = Bρ = p/q = 5.73730218421
    i2kl = lambda i: (.265410*i-.765828e-6*i**3-.239385)/R
    i2ks = lambda i: (.266015*i-.829333e-6*i**3-.233316)/R
    '''===== drifts ====='''
    D01, D02, D03, D04, D05, D06, D07, D08, D09, D10, D11, D12, D13, D14, \
        D15, D16, D17, D18, D19, D20, D21, D22, D23, D24, D25 = \
        [np.zeros([6, 1]) for i in range(25)]
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
    diagnos = np.append([[7]], np.zeros([5, 1]), 0)
    F1T = F2T = F3T = F4T = F5T = F6T = F7T = F8T = diagnos
    '''===== dipoles ====='''
    # pole shoe form factor K (close to rectangular)
    # dipoles gap g
    dipole = np.append([[1]], np.zeros([5, 1]), 0)
    edge = np.append([[2]], np.zeros([5, 1]), 0)
    LD = .7; rho = -LD/2/np.pi*360/6.13; phi = LD/2/rho; g = 15e-3; K = 0.5
    dipole[[1, 2]] = np.array([[LD], [rho]])
    edge[[2, 3, 5]] = np.array([[rho], [phi], [g*K]])
    DP1 = np.concatenate((edge, dipole, edge), 1)
    LD = .7; rho = -LD/2/np.pi*360/4.38; phi = LD/2/rho; g = 15e-3; K = 0.5
    dipole[[1, 2]] = np.array([[LD], [rho]])
    edge[[2, 3, 5]] = np.array([[rho], [phi], [g*K]])
    DP2 = np.concatenate((edge, dipole, edge), 1)
    LD = 1.777792; rho = LD/2/np.pi*360/22; phi = LD/2/rho; g = 30e-3; K = 0.5
    dipole[[1, 2]] = np.array([[LD], [rho]])
    edge[[2, 3, 5]] = np.array([[rho], [phi], [g*K]])
    DP3 = np.concatenate((edge, dipole, edge), 1)
    DP4 = DP3
    LD = 1.020001; rho = LD/2/np.pi*360/7.66; phi = LD/2/rho; g = 15e-3; K = 0.5
    dipole[[1, 2]] = np.array([[LD], [rho]])
    edge[[2, 3, 5]] = np.array([[rho], [phi], [g*K]])
    DP5 = np.concatenate((edge, dipole, edge), 1)
    LD = .555; rho = LD/2/np.pi*360/3.8; phi = LD/2/rho; g = 15e-3; K = 0.5
    dipole[[1, 2]] = np.array([[LD], [rho]])
    edge[[2, 3, 5]] = np.array([[rho], [phi], [g*K]])
    DP6 = np.concatenate((edge, dipole, edge), 1)
    DP7 = DP6
    ''' quadrupoles '''
    Q01, Q02, Q03, Q04, Q05, Q07, Q08, Q09, Q10, Q11, Q12 = \
        [np.zeros([6, 1]) for i in range(11)]
    Q01[[0, 1, 4]] = np.array([[4], [.25], [-i2kl(I1)]])
    Q02[[0, 1, 4]] = np.array([[4], [.25], [-i2kl(I2)]])
    Q03[[0, 1, 4]] = np.array([[3], [.25], [i2kl(I3)]])
    Q04[[0, 1, 4]] = np.array([[3], [.25], [i2kl(I4)]])
    Q05[[0, 1, 4]] = np.array([[4], [.25], [-i2kl(I5)]])
    Q07[[0, 1, 4]] = np.array([[4], [.25], [-i2kl(I7)]])
    Q08[[0, 1, 4]] = np.array([[3], [.25], [i2kl(I8)]])
    Q09[[0, 1, 4]] = np.array([[4], [.25], [-i2kl(I9)]])
    Q10[[0, 1, 4]] = np.array([[4], [.2], [-i2ks(I10)]])
    Q11[[0, 1, 4]] = np.array([[4], [.2], [-i2ks(I11)]])
    Q12[[0, 1, 4]] = np.array([[3], [.2], [i2ks(I12)]])
    '''===== cells ====='''
                       # 00.. 03   04   05   06.. 09   10   11   12   13
    UC = np.concatenate((DP1, D01, F1T, D02, DP2, D03, F2T, D04, Q01, D05,
                       # 14   15   16   17   18   19   20.. 23   24   25
                         Q02, D06, Q03, D07, F3T, D08, DP3, D09, F4T, D10,
                       # 26   27   28   29   30   31   32   33   34   35
                         Q04, D11, Q05, D12, F5T, D13, Q07, D14, Q08, D15,
                       # 36   37   38   39   40.. 43   44   45   46   47
                         Q09, D16, F6T, D17, DP4, D18, Q10, D19, Q11, D20,
                       # 48   49   50   51   52.. 55   56   57   58.. 61
                         Q12, D21, F7T, D22, DP5, D23, F8T, D24, DP6, D25,
                       # 62..
                         DP7), 1)
    UC = UC[:, start:end]
    print UC
    R = UC2R(UC, gamma)
    return R


def trimplot(fom, quad, fig, subplots, subplot, data, labels):
    ax = fig.add_subplot(subplots[0], subplots[1], subplot)
    ax.plot(data, '.b', label='trimmed data')
    ax.text(0.5, 0.95, r'$FOMZ%iT(Q%iPT)$' % (fom, quad),
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.locator_params(axis='x', nbins=5)
    if subplot == 1 or subplot == subplots[1]+1:
        ax.set_ylabel(labels[1])
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
    if subplot > subplots[1]:
        ax.set_xlabel(labels[0])


def trimkplot(fom, quad, fig, subplots, subplot, data, trimmin, trimmax,
              labels):
    edge = np.array([np.min(data), np.max(data)])
    ax = fig.add_subplot(subplots[0], subplots[1], subplot)
    ax.plot(data, '.b', label='data')
    ax.plot([trimmin, trimmin], [edge[0], edge[1]], '-r')
    ax.plot([trimmax, trimmax], [edge[0], edge[1]], '-r', label='trimline')
    ax.text(0.5, 0.95, r'$FOMZ%iT(Q%iPT)$' % (fom, quad),
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.locator_params(axis='x', nbins=5)
    if subplot == 1 or subplot == subplots[1]+1:
        ax.set_ylabel(labels[1])
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
    if subplot > subplots[1]:
        ax.set_xlabel(labels[0])


def trimcurplot(fom, quad, fig, subplots, subplot, data, trimmax, trimmin,
                labels):
    l = np.size(data)
    ax = fig.add_subplot(subplots[0], subplots[1], subplot)
    ax.plot(data, '.b', label='data')
    ax.plot([0, l], [trimmax, trimmax], '-r')
    ax.plot([0, l], [trimmin, trimmin], '-r', label='trimline')
    ax.text(0.5, 0.95, r'$FOMZ%iT(Q%iPT)$' % (fom, quad),
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.locator_params(axis='x', nbins=5)
    if subplot == 1 or subplot == subplots[1]+1:
        ax.set_ylabel(labels[1])
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
    if subplot > subplots[1]:
        ax.set_xlabel(labels[0])


def trimroiplot(fom, quad, fig, subplots, subplot, x, y, spline, labels):
    ax = fig.add_subplot(subplots[0], subplots[1], subplot)
    ax.plot(x, y, '.b', label='data')
    ax.text(0.5, 0.95, r'$FOMZ%iT(Q%iPT)$' % (fom, quad),
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes)
    ax.locator_params(axis='x', nbins=5)
    if subplot == 1 or subplot == subplots[1]+1:
        ax.set_ylabel(labels[1])
        leg = ax.legend(fancybox=True, loc='lower right')
        leg.get_frame().set_alpha(0.5)
    if subplot > subplots[1]:
        ax.set_xlabel(labels[0])
    if spline == 1:
        spline = lowess(x, y)
        sm = spline - np.mean(y)/3
        sp = spline + np.mean(y)/3
        ax.plot(x, spline, '-g', label='spline')
        ax.plot(x, sp, '-r')
        ax.plot(x, sm, '-r', label='trimline')
        return sm, sp
