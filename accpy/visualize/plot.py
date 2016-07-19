# -*- coding: utf-8 -*-
''' accpy.visualize.plot
author:     felix.kramer(at)physik.hu-berlin.de

'best'         : 0, (only implemented for axes legends)
'upper right'  : 1,
'upper left'   : 2,
'lower left'   : 3,
'lower right'  : 4,
'right'        : 5,
'center left'  : 6,
'center right' : 7,
'lower center' : 8,
'upper center' : 9,
'center'       : 10,
'''
from __future__ import division
from numpy import (nanmax, nanmin, concatenate, empty, linspace, array, arange,
                   abs as npabs, sqrt, sin, cos, arccos as acos, nanmean)
from itertools import product
from matplotlib.figure import Figure
from matplotlib.pyplot import cm
from matplotlib.gridspec import GridSpec
from .stringformat import SI, SId
from .lattice import drawlattice
from ..simulate import const
from ..simulate.rmatrices import UCS2R
from ..simulate.tracking import trackpart


def plot(ax, x, y, ls, xlabel, xunit, ylabel, yunit, label, col=False,
         setlim=True, rescaleX=True, rescaleY=True, xprefix=None, mx=None,
         yprefix=None, my=None, markersize=12.0, linewidth=2.0):
    if xprefix is None:
        xprefix, mx = SId(nanmean(abs(x)))
    if yprefix is None:
        yprefix, my = SId(nanmean(abs(y)))
    if rescaleX:
        x = x/mx  # carefull! numpy.ndarrays are mutable!!!
        if xunit != '':
            xunit = ' / ('+xprefix+xunit+')'
    elif xunit != '':
            xunit = ' / ('+xunit+')'
    if rescaleY:
        y = y/my
        if yunit != '':
            yunit = ' / ('+yprefix+yunit+')'
    elif yunit != '':
            yunit = ' / ('+yunit+')'
    if col is False:
        ax.plot(x, y, ls, label=label, markersize=markersize, linewidth=linewidth)
    else:
        ax.plot(x, y, ls, color=col, label=label, markersize=markersize, linewidth=linewidth)
    if xlabel != '':
        ax.set_xlabel(xlabel+xunit)
    if ylabel != '':
        ax.set_ylabel(ylabel+yunit)
    if setlim:
        epsy = (max(y)-min(y))*0.15
        ax.set_xlim([min(x), max(x)])
        ax.set_ylim([min(y)-epsy, max(y)+epsy])
    return x, y, yprefix, my


def Mplot(ax, x, ys, lss, xlabel, xunit, ylabel, yunit, labels, rescaleX=True, rescaleY=True):
    colors = getcolors(len(ys))
    xprefix, mx = SId(nanmean(x))
    yprefix, my = SId(nanmean(ys))
    if rescaleX:
        if xunit != '':
            xunit = ' / ('+xprefix+xunit+')'
    elif xunit != '':
            xunit = ' / ('+xunit+')'
    if rescaleY:
        if yunit != '':
            yunit = ' / ('+yprefix+yunit+')'
    elif yunit != '':
            yunit = ' / ('+yunit+')'
    if labels == '':
        labels = ['' for i in range(len(ys))]
    if type(x) != type([]):
        xs = [x for i in range(len(ys))]
    else:
        xs = x
    for x, y, ls, lab, col in zip(xs, ys, lss, labels, colors):
        if rescaleY:
            y = y/my
        if rescaleX:
            x = x/mx  # carefull! numpy.ndarrays are mutable!!!
        ax.plot(x, y, ls, color=col, label=lab)
    if xlabel != '':
        ax.set_xlabel(xlabel+xunit)
    if ylabel != '':
        ax.set_ylabel(ylabel+yunit)
    ax.set_xlim([min(x), max(x)])
    return mx, my


def legplot(ax, linestyles, labels, loc=0):
    colors = getcolors(len(labels))
    for ls, lab, col in zip(linestyles, labels, colors):
        ax.plot([], [], ls, label=lab, color=col)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.legend(fancybox=True, loc=loc)
    return


def getcolors(number):
    return list(cm.rainbow(linspace(0, 1, number)))


def plotoptic(UC, diagnostics, s, xtwiss, ytwiss, xdisp):
    figs = [Figure() for i in range(4)]
    ax = [figs[i].add_subplot(1, 1, 1) for i in range(4)]
    drawlattice(ax[0], UC, diagnostics, [xtwiss[0, 0, :], xtwiss[0, 1, :]], 0)
    ax[0].plot(s, xtwiss[0, 0, :], '-r', label=r'$\beta_x$')
    ax[0].plot(s, xtwiss[0, 1, :], '-c', label=r'$\alpha_x$')
    ax[0].plot([], [], '-m', label=r'$\gamma_x$')
    ax[0].set_ylabel(r'betatron function $\beta_x$ / (m)')
    ax2 = ax[0].twinx()
    ax2.plot(s, xtwiss[1, 1, :], '-m')
    ax2.set_ylabel(r'gamma function $\gamma_x$ / (m)', color='m')
    ax2.tick_params(axis='y', colors='m')
    drawlattice(ax[1], UC, diagnostics, [ytwiss[0, 0, :], ytwiss[0, 1, :]], 0)
    ax[1].plot(s, ytwiss[0, 0, :], '-b', label=r'$\beta_y$')
    ax[1].plot(s, ytwiss[0, 1, :], '-c', label=r'$\alpha_y$')
    ax[1].plot([], [], '-m', label=r'$\gamma_y$')
    ax[1].set_ylabel(r'betatron function $\beta_y$ / (m)')
    ax2 = ax[1].twinx()
    ax2.plot(s, ytwiss[1, 1, :], '-m', label=r'$\gamma_y$')
    ax2.set_ylabel(r'gamma function $\gamma_y$ / (m)', color='m')
    ax2.tick_params(axis='y', colors='m')
    drawlattice(ax[2], UC, diagnostics, [xdisp[0, :]], 0)
    ax[2].plot(s, xdisp[0, :], '-g', label=r'$D_x$')
    ax[2].plot([], [], '-m', label=r'$D_x^\prime$')
    ax[2].set_ylabel(r'dispersion function $D_x$ / (m)')
    ax2 = ax[2].twinx()
    ax2.plot(s, xdisp[1, :], '-m', label=r'$D_x^\prime$')
    ax2.set_ylabel(r'derived dispersion function $D_x^\prime$ / (m)', color='m')
    ax2.tick_params(axis='y', colors='m')
    drawlattice(ax[3], UC, diagnostics, [xtwiss[0, 0, :], ytwiss[0, 0, :]], 0)
    ax[3].plot(s, xtwiss[0, 0, :], '-r', label=r'$\beta_x$')
    ax[3].plot(s, ytwiss[0, 0, :], '-b', label=r'$\beta_y$')
    ax[3].plot([], [], '-g', label=r'$D_x$')
    ax[3].set_ylabel(r'betatron function $\beta_{x,y}$ / (m)')
    ax2 = ax[3].twinx()
    ax2.plot(s, xdisp[0, :], '-g', label=r'$D_x$')
    ax2.set_ylabel(r'dispersion function $D_x$ / (m)', color='g')
    ax2.tick_params(axis='y', colors='g')
    [ax[i].set_xlabel(r'orbit position s / (m)') for i in range(4)]
    [ax[i].set_xlim([0, nanmax(s)]) for i in range(4)]
    legs = [ax[i].legend(fancybox=True, loc='upper left') for i in range(4)]
    [legs[i].get_frame().set_alpha(0.5) for i in range(4)]
    return figs


def plotbeamsigma(UC, diagnostics, s, sigx, sigy):
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    rel = abs(nanmean(sigy)/nanmean(sigx))
    if rel > 100 or rel < 1e-2:
        drawlattice(ax, UC, diagnostics, [sigx], 0)
        ax.plot(s, sigx, '-r', label=r'$\sigma_x$')
        ax.plot([], [], '-b', label=r'$\sigma_y$')
        ax.set_ylabel(r'Beam extent $\sigma_x$ / (m)')
        ax2 = ax.twinx()
        ax2.plot(s, sigy, '-b')
        ax2.tick_params(axis='y', colors='b')
        ax2.set_ylabel(r'Beam extent $\sigma_y$ / (m)', color='b')
    else:
        drawlattice(ax, UC, diagnostics, [sigx, sigy], 0)
        ax.plot(s, sigx, '-r', label=r'$\sigma_x$')
        ax.plot(s, sigy, '-b', label=r'$\sigma_y$')
        ax.set_ylabel(r'Beam extent $\sigma_u$ / (m)')
    ax.set_xlim([min(s), max(s)])
    leg = ax.legend(fancybox=True, loc=2)
    leg.get_frame().set_alpha(0.5)
    return fig


def plotopticpars_closed(xtwiss, xdisp, ytwiss, gamma, Qx, Xx, Jx, emiteqx,
                         tau_x, Qy, Xy, Jy, E, emiteqy, tau_y, alpha_mc,
                         eta_mc, gamma_tr, Q_s, Js, sigma_E, sigma_tau,
                         sigma_s, tau_s, U_rad, P_ges, E_c, lambda_c):
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    radpars = ''.join(
        [r'------------------------------------''\n',
         r'Radial parameters''\n',
         r'------------------------------------''\n',
         r'$\beta_{x,max} = %g m$''\n' % nanmax(xtwiss[0, 0, :]),
         r'$\beta_{x,min} = %g m$''\n' % nanmin(xtwiss[0, 0, :]),
         r'$\alpha_{x,max} = %g$''\n' % nanmax(xtwiss[0, 1, :]),
         r'$\alpha_{x,min} = %g$''\n' % nanmin(xtwiss[0, 1, :]),
         r'$\gamma_{x,max} = %g$''\n' % nanmax(xtwiss[1, 1, :]),
         r'$\gamma_{x,min} = %g$''\n' % nanmin(xtwiss[1, 1, :]),
         r'$D_{x,max} = %g m$''\n' % nanmax(xdisp[0, :]),
         r'$D_{x,min} = %g m$''\n' % nanmin(xdisp[0, :]),
         r'$D_{x,max}^\prime = %g$''\n' % nanmax(xdisp[1, :]),
         r'$D_{x,min}^\prime = %g$''\n' % nanmin(xdisp[1, :]),
         r'$Q_x = %g$''\n' % Qx,
         r'$\xi_{x,nat} = %g$''\n' % Xx,
         r'$J_x = %g$''\n' % Jx,
         r'$\epsilon_x = %g \pi rad m$''\n' % emiteqx,
         r'$\tau_x = %e s$' % tau_x])
    axipars = ''.join(
        [r'------------------------------------''\n',
         r'Axial parameters''\n',
         r'------------------------------------''\n',
         r'$\beta_{y,max} = %g m$''\n' % nanmax(ytwiss[0, 0, :]),
         r'$\beta_{y,min} = %g m$''\n' % nanmin(ytwiss[0, 0, :]),
         r'$\alpha_{y,max} = %g$''\n' % nanmax(ytwiss[0, 1, :]),
         r'$\alpha_{y,min} = %g$''\n' % nanmin(ytwiss[0, 1, :]),
         r'$\gamma_{x,max} = %g$''\n' % nanmax(ytwiss[1, 1, :]),
         r'$\gamma_{x,min} = %g$''\n' % nanmin(ytwiss[1, 1, :]),
         r'$Q_y = %g$''\n' % Qy,
         r'$\xi_{y,nat} = %g$''\n' % Xy,
         r'$J_y = %g$''\n' % Jy,
         r'$\epsilon_y = %g \pi rad m$''\n' % emiteqy,
         r'$\tau_y = %e s$' % tau_y])
    lonpars = ''.join(
        [r'------------------------------------''\n',
         r'Longitudinal parameters''\n',
         r'------------------------------------''\n',
         r'$E = %g eV$''\n' % E,
         r'$\gamma_{lorentz} = %g$''\n' % gamma,
         r'$\alpha_{p} = %g $''\n' % alpha_mc,
         r'$\eta_{slip} = %g $''\n' % eta_mc,
         r'$\gamma_{tr} = %g $''\n' % gamma_tr,
         r'$Q_s = %g$''\n' % Q_s,
         r'$J_s = %g$''\n' % Js,
         r'$\sigma_{E} = %e \%% $''\n' % sigma_E,
         r'$\sigma_{\tau} = %g s$''\n' % sigma_tau,
         r'$\sigma_{s} = %g m$''\n' % sigma_s,
         r'$\tau_{s} = %e s$''\n' % tau_s,
         r'$E_{loss} = %g eV$''\n' % U_rad,
         r'$P_{rad} = %g W$''\n' % P_ges,
         r'$E_{crit} = %g eV$''\n' % E_c,
         r'$\lambda_{crit} = %g m$' % lambda_c])
    ax.text(0.05, 0.9, radpars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.text(0.35, 0.9, axipars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.text(0.65, 0.9, lonpars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    return fig


def plotopticpars_open(xtwiss, xdisp, ytwiss, gamma, E):
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    radpars = ''.join(
        [r'------------------------------------''\n',
         r'Radial parameters''\n',
         r'------------------------------------''\n',
         r'$\beta_{x,max} = %g m$''\n' % nanmax(xtwiss[0, 0, :]),
         r'$\beta_{x,min} = %g m$''\n' % nanmin(xtwiss[0, 0, :]),
         r'$\alpha_{x,max} = %g$''\n' % nanmax(xtwiss[0, 1, :]),
         r'$\alpha_{x,min} = %g$''\n' % nanmin(xtwiss[0, 1, :]),
         r'$\gamma_{x,max} = %g$''\n' % nanmax(xtwiss[1, 1, :]),
         r'$\gamma_{x,min} = %g$''\n' % nanmin(xtwiss[1, 1, :]),
         r'$D_{x,max} = %g m$''\n' % nanmax(xdisp[0, :]),
         r'$D_{x,min} = %g m$''\n' % nanmin(xdisp[0, :]),
         r'$D_{x,max}^\prime = %g$''\n' % nanmax(xdisp[1, :]),
         r'$D_{x,min}^\prime = %g$' % nanmin(xdisp[1, :])])
    axipars = ''.join(
        [r'------------------------------------''\n',
         r'Axial parameters''\n',
         r'------------------------------------''\n',
         r'$\beta_{y,max} = %g m$''\n' % nanmax(ytwiss[0, 0, :]),
         r'$\beta_{y,min} = %g m$''\n' % nanmin(ytwiss[0, 0, :]),
         r'$\alpha_{y,max} = %g$''\n' % nanmax(ytwiss[0, 1, :]),
         r'$\alpha_{y,min} = %g$''\n' % nanmin(ytwiss[0, 1, :]),
         r'$\gamma_{x,max} = %g$''\n' % nanmax(ytwiss[1, 1, :]),
         r'$\gamma_{x,min} = %g$' % nanmin(ytwiss[1, 1, :])])
    lonpars = ''.join(
        [r'------------------------------------''\n',
         r'Longitudinal parameters''\n',
         r'------------------------------------''\n',
         r'$E = %g eV$''\n' % E,
         r'$\gamma_{lorentz} = %g$' % gamma])
    ax.text(0.05, 0.9, radpars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.text(0.35, 0.9, axipars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    ax.text(0.65, 0.9, lonpars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes)
    return fig


def plotdisptraj(s, P_UCS, E, E0, UCS, UC, diagnostics):
    # measured energy dependant offset at FOMS normalized to 0 for EbE0=1
    xf1t = lambda EbE0: -.078269*EbE0 + .078269     # + .059449
    xf2t = lambda EbE0: -.241473*EbE0 + .241473     # + .229314
    xf6t = lambda EbE0: 1.174523*EbE0 - 1.174523    # - 1.196090
    xf7t = lambda EbE0: .998679*EbE0 - .998679      # - 1.018895
    xf8t = lambda EbE0: .769875*EbE0 - .769875      # - .787049
    steps = 6
    X = empty([6, P_UCS+1, steps])
    dEbE = linspace(-0.005, 0.005, steps)
    for deltaE, i in zip(dEbE, range(steps)):
        # R calculated for every energy (not necessary)
        gamma = (E+deltaE*E)/E0+1
        R = UCS2R(P_UCS, UCS, gamma)
        X0 = array([[0], [0], [0], [0], [0], [deltaE]])
        X[:, :, i] = trackpart(R, P_UCS, X0)
    fig = Figure()
    ax = fig.add_subplot(1, 1, 1)
    drawlattice(ax, UC, diagnostics, [X.flatten()], 0)
    ax.set_xlabel(r'orbit position s / (m)')
    ax.set_ylabel(r'radial displacement / (m)')
    x = [s[UCS[0, :] == 7][i] for i in [0, 1, 5, 6, 7]]
    color = iter(cm.rainbow(linspace(0, 1, steps)))
    for i in range(steps):
        c = next(color)
        EE0 = 1 + dEbE[i]
        y = array([xf1t(EE0), xf2t(EE0), xf6t(EE0), xf7t(EE0), xf8t(EE0)])
        ax.plot(x, y, 'o', c=c)
        ax.plot(s, X[0, :, i], c=c, label=r'$\Delta E/E_0=%.3f$' % dEbE[i])
    ax.plot([], [], 'ok', label=r'measured')
    leg = ax.legend(fancybox=True, loc='upper left')
    leg.get_frame().set_alpha(0.5)
    ax.set_xlim([0, nanmax(s)])
    return fig


def plottrajs(s, X, rounds, envelope):
    figs = [Figure() for i in range(7)]
    ax1 = [figs[i].add_subplot(1, 1, 1) for i in range(6)]
    GS = GridSpec(5, 3)
    ax2 = [figs[6].add_subplot(GS[:2, i]) for i in range(3)]
    ax2 += [figs[6].add_subplot(GS[3:, i]) for i in range(3)]
    ax3 = figs[6].add_subplot(GS[2, :])
    ylabs = [r'$x$ radial displacement / (mm)',
             r'$x^\prime$ radial direction deviation / (mrad)',
             r'$y$ axial displacement / (mm)',
             r'$y^\prime$ axial direction deviation / (mrad)',
             r'$l$ longitudinal displacement / (mm)',
             r'$\frac{\Delta P}{P_0}$ longitudinal momentum deviation']
    y2labs = [r'$x$ / (mm)',
              r'$x^\prime$ / (mrad)',
              r'$y$ / (mm)',
              r'$y^\prime$ / (mrad)',
              r'$l$ / (mm)',
              r'$\frac{\Delta P}{P_0}$ / \textperthousand']
    for i in range(6):
        ax1[i].set_xlabel(r'orbit position s / (m)')
        ax1[i].set_ylabel(ylabs[i])
        ax2[i].set_xlabel(r'orbit position s / (m)')
        ax2[i].set_ylabel(y2labs[i])

    color = iter(cm.rainbow(linspace(0, 1, 6)))
    labs = [r'Ideal particle',
            r'1 sigma particle',
            r'Envelope $E_{x,y}(s)=\sqrt{\epsilon_{x,y}\beta_{x,y}(s)+(\delta_ED_{x,y}(s))^2)}$',
            r'Envelope $E_{x}(s)=\sqrt{\epsilon_{x}\beta_{x}+(\delta_ED_{x}(s))^2}$',
            r'Envelope $E_{y}(s)=\sqrt{\epsilon_{y}\beta_{y}}$',
            r'Ensemble']
    for i, j in zip(range(6), [0, 3, 1, 4, 2, 5]):
        c = next(color)
        y = []
        y_ideal = []
        y_sigma = []
        for partN, traj in enumerate(X):
            for k in range(rounds):
                index = arange(len(s))+(len(s)-1)*k
                if partN > 1:
                    y.append(traj[i, index]*1e3)
                elif partN == 1:
                    y_sigma.append(traj[i, index]*1e3)
                elif partN == 0:
                    y_ideal.append(traj[i, index]*1e3)
        # ensemble trajectories
        [ax1[i].plot(s, y[l], '-', c=c) for l in range(len(y))]
        [ax2[j].plot(s, y[l], '-', c=c) for l in range(len(y))]
        ax1[i].plot([], [], '-', c=c, label=labs[5])
        ax2[j].plot([], [], '-', c=c, label=labs[5])
        # 1-sigma particle trajectories
        ax1[i].plot(s, y_sigma[0], '-b', label=labs[1])
        ax2[j].plot(s, y_sigma[0], '-b')
        [ax1[i].plot(s, y_sigma[l], '-b') for l in range(1, len(y_sigma))]
        [ax2[j].plot(s, y_sigma[l], '-b') for l in range(1, len(y_sigma))]
        # ideal particle trajectories
        ax1[i].plot(s, y_ideal[0], '-k', label=labs[0])
        if i == 0:
            ax1[i].plot([], [], '-r', label=labs[3])
        elif i == 2:
            ax1[i].plot([], [], '-r', label=labs[4])
        ax2[j].plot(s, y_ideal[0], '-k')
        leg = ax1[i].legend(fancybox=True, loc='upper right')
        leg.get_frame().set_alpha(0.5)
    ax1[0].plot(s, envelope[0, :], '-r', s, -envelope[0, :], '-r')
    ax1[2].plot(s, envelope[1, :], '-r', s, -envelope[1, :], '-r')
    ax2[0].plot(s, envelope[0, :], '-r', s, -envelope[0, :], '-r')
    ax2[2].plot(s, envelope[1, :], '-r', s, -envelope[1, :], '-r')
    legplot(ax3, ['-k', '-b', '-r'], labs[:3], loc=10)
    return figs


def plotphasespace(s, X, rounds, xtwiss, emittx, ytwiss, emitty):
    fig = Figure()
    xlabels = [r'$x$ / (mm)',
               r'$y$ / (mm)']
    ylabels = [r'$x^\prime$ / (mrad)',
               r'$y^\prime$ / (mrad)']
    titles = [r'Radial phasespace',
              r'Axial phasespace']
    axmax = []

    def roundplot(traj, ax, linestyle, label=''):
        for j in range(rounds):
            index = (len(s)-1)*j
            x = traj[i*2, index]*1e3
            y = traj[i*2+1, index]*1e3
            if j == 0:
                ax.plot(x, y, linestyle, label=label)
            else:
                ax.plot(x, y, linestyle)
            axmax.append(max(npabs([x, y])))

    GS = GridSpec(1, 5)
    ax = []
    ax.append(fig.add_subplot(GS[0, :2]))
    ax.append(fig.add_subplot(GS[0, -2:]))
    ax2 = fig.add_subplot(GS[0, 2])

    for i in range(2):
        ax[i].set_xlabel(xlabels[i])
        ax[i].set_ylabel(ylabels[i])
        ax[i].set_title(titles[i])
        for k in range(len(X)):
            traj = X[k]
            if k == 1:
                roundplot(traj, ax[i], 'xr')
            elif k == 2:
                roundplot(traj, ax[i], '.b')
            else:
                roundplot(traj, ax[i], '.b')
    axmax = max(axmax)
    for i in range(2):
        ax[i].set_xlim([-axmax, axmax])
        ax[i].set_ylim([-axmax, axmax])
    x, xp, y, yp = twissellipse(xtwiss[:, :, 0], emittx,
                                ytwiss[:, :, 0], emitty)
    ax[0].plot(x, xp, '-g')
    ax[1].plot(y, yp, '-g')
    ax[1].yaxis.tick_right()
    ax[1].yaxis.set_label_position('right')
    ax2.plot([], [], '.b', label='Ensemble')
    ax2.plot([], [], '-g', label='Twiss ellipsis')
    ax2.plot([], [], 'xr', label='1 sigma particle')
    ax2.axis('off')
    ax2.get_xaxis().set_visible(False)
    ax2.get_yaxis().set_visible(False)
    ax2.legend(loc='center')
    return fig


def twissellipse(xtwiss, emittx, ytwiss, emitty):
    def ellipse(emittance, beta, alpha, gamma):
        phi = linspace(0, 2*const.pi, 1e3)
        a = sqrt(emittance/2*(beta+gamma+sqrt((beta+gamma)**2-4)))
        b = sqrt(emittance/2*(beta+gamma-sqrt((beta+gamma)**2-4)))
        if alpha > 0:
            PHI = acos(+sqrt((beta-b/a)*emittance/(a**2-b**2)))
        else:
            PHI = acos(-sqrt((beta-b/a)*emittance/(a**2-b**2)))
        pos = a*cos(phi)*cos(PHI)+b*sin(phi)*sin(PHI)
        mom = -a*cos(phi)*sin(PHI)+b*sin(phi)*cos(PHI)
        return pos, mom
    x, xp = ellipse(emittx, xtwiss[0, 0], -xtwiss[0, 1], xtwiss[1, 1])
    y, yp = ellipse(emitty, ytwiss[0, 0], -ytwiss[0, 1], ytwiss[1, 1])
    return x, xp, y, yp


def plotramp(T, t, tt, tt2, tEgZ, tAI, tVgZ, E, EE, EEgZ, EAI, EVgZ, B, BB, loss, LL, volt,
             VV, phases, freqs, Xemitequi, Yemitequi, Semitequi, bdurequis,
             blenequis, V_HFs, Xemits, Yemits, Semits, NXemitequi, NYemitequi,
             NXemits, NYemits, bdurs, blens, t3, FF):
    def annotate(ax, xs, ys, ss, epsxs, epsys, hs, vs):
        for x, y, s, h, v, epsx, epsy in zip(xs, ys, ss, hs, vs, epsxs, epsys):
            ax.text(x+epsx, y+epsy, s, horizontalalignment=h, verticalalignment=v)

    Nfigs = 10
    legs = []
    figs = [Figure() for i in range(Nfigs)]
    ax = [figs[i].add_subplot(1, 1, 1) for i in range(6)]
    [ax[i].grid() for i in range(len(ax))]
    GS = GridSpec(2, 3)
    GS = [GS[:2, 0], GS[0, 1:], GS[1, 1:]]
    ax.append([figs[6].add_subplot(gs) for gs in GS])
    [ax.append([figs[j].add_subplot(2, 2, i) for i in range(1, 5)]) for j in range(7, 9)]
    ax.append([figs[9].add_subplot(1, 2, i) for i in range(1, 3)])
    [[ax[i][j].grid() for j in range(len(ax[i]))] for i in range(6, 10)]
    s = [r'Injection''\n(',
         r'Extraction''\n(',
         r'Maximum energy''\n(',
         r'Alternative extraction''\n(',
         r'Zero Voltage''\n(',
         r'Maximum Voltage''\n(']
    s1 = [s[0]+SI(tt[0])+'s, '+SI(EE[0]) + 'eV)',
          s[1]+SI(tt[1])+'s, '+SI(EE[1]) + 'eV)',
          s[2]+SI(tt[2])+'s, '+SI(EE[2]) + 'eV)',
          s[3]+SI(tt[3])+'s, '+SI(EE[3]) + 'eV)']
    s2 = [s[0]+SI(tt[0])+'s, '+SI(BB[0]) + 'T)',
          s[1]+SI(tt[1])+'s, '+SI(BB[1]) + 'T)',
          s[2]+SI(tt[2])+'s, '+SI(BB[2]) + 'T)',
          s[3]+SI(tt[3])+'s, '+SI(BB[3]) + 'T)']
    s3 = [s[0]+SI(tt[0])+'s, '+SI(LL[0]) + 'eV)',
          s[1]+SI(tt[1])+'s, '+SI(LL[1]) + 'eV)',
          s[2]+SI(tt[2])+'s, '+SI(LL[2]) + 'eV)',
          s[3]+SI(tt[3])+'s, '+SI(LL[3]) + 'eV)']
    s4 = [s[0]+SI(tt2[0])+'s, '+SI(VV[0]) + 'V)',
          s[1]+SI(tt2[1])+'s, '+SI(VV[1]) + 'V)',
          s[2]+SI(tt2[2])+'s, '+SI(VV[2]) + 'V)',
          s[3]+SI(tt2[3])+'s, '+SI(VV[3]) + 'V)',
          s[4]+SI(tt2[4])+'s, '+SI(VV[4]) + 'V)',
          s[5]+SI(tt2[5])+'s, '+SI(VV[5]) + 'V)']

    xlab, xunit = 'Time', 's'
    xlab2, xunit2 = 'Energy', 'eV'

    # Energy
    _, mx = SId(nanmean(abs(t)))
    epsT = array([1, -1, 0, 1])*.01*T/mx
    ha = ['left', 'right', 'center', 'left']
    va = ['top', 'bottom', 'bottom', 'bottom']
    _, _, yprefix, my = plot(ax[0], t, E, '-r', 'Time', 's', 'Energy', 'eV', 'calculated curve')
    ttn, EEn, _, _ = plot(ax[0], tt, EE, '+k', '', '', '', '', 'known points', yprefix=yprefix, my=my, setlim=False, markersize=24.0)
    legs.append(ax[0].legend(fancybox=True, loc='lower center'))
    epsY = array([-1, 0, 1, 0])*.03*max(EEn)
    annotate(ax[0], ttn, EEn, s1, epsT, epsY, ha, va)

    # Magnetic flux
    _, _, yprefix, my = plot(ax[1], t, B, '-r', 'Time', 's', 'Magnetic flux density', 'T', 'calculated curve')
    ttn, BBn, _, _ = plot(ax[1], tt, BB, '+k', '', '', '', '', 'known points', yprefix=yprefix, my=my, setlim=False, markersize=24.0)
    legs.append(ax[1].legend(fancybox=True, loc='lower center'))
    epsY = array([-1, 0, 1, 0])*.03*max(BBn)
    annotate(ax[1], ttn, BBn, s2, epsT, epsY, ha, va)

    # Energy loss
    _, _, yprefix, my = plot(ax[2], t, loss, '-r', 'Time', 's', 'Energyloss per turn', 'eV', '')
    ttn, LLn, _, _ = plot(ax[2], tt, LL, '+k', '', '', '', '', 'known points', yprefix=yprefix, my=my, setlim=False, markersize=24.0)
    epsY = array([-1, 0, 1, 0])*.03*max(LLn)
    annotate(ax[2], ttn, LLn, s3, epsT, epsY, ha, va)

    # Acceleration voltage
    epsT = array([1, -1, 1, 1, 1, -1])*.01*T/mx
    ha = ['left', 'right', 'left', 'left', 'left', 'right']
    va = ['top', 'bottom', 'bottom', 'bottom', 'bottom', 'bottom']
    _, _, yprefix, my = plot(ax[3], t, volt, '-r', 'Time', 's', 'Required acceleration voltage', 'V', '')
    ttn, VVn, _, _ = plot(ax[3], tt2, VV, '+k', '', '', '', '', 'known points', yprefix=yprefix, my=my, setlim=False, markersize=24.0)
    epsY = array([-1, 0, 1, 0, -1, 1])*.02*max(VVn)
    annotate(ax[3], ttn, VVn, s4, epsT, epsY, ha, va)

    # Synchronous phase
    labs = ['Cavity @ {0:g} kV'.format(V_HF/1e3) for V_HF in V_HFs]
    color = iter(cm.rainbow(linspace(0, 1, len(phases))))
    [plot(ax[4], t, y, '-', 'Time', 's', 'Cavity Phase', r'2$\pi$',
          l, col=next(color), setlim=False) for l, y in zip(labs, phases)]
    legs.append(ax[4].legend(fancybox=True, loc='center right'))

    # Synchrotron frequency
    colors = getcolors(len(freqs))
    [plot(ax[5], tVgZ, y/1e3, '-', 'Time', 's', 'Synchrotron frequency', r'kHz',
          l, col=c, setlim=False) for l, y, c in zip(labs, freqs, colors)]
    [ax[5].plot(t3*1e3, x/1e3, '+k', markersize=24.0) for x in FF]
    legs.append(ax[5].legend(fancybox=True, loc=1))

    # Bunchlength and duration
    labs = [r'$\delta_{{E,equilibrium}}$, $V_{{max}} = {0:g}$ kV'.format(V_HF/1e3) for V_HF in V_HFs]
    lss = ['-.' for x in range(len(bdurequis))]
    labs += [r'$\delta_{{E,0}}={0:.1f}$ \textperthousand, $V_{{max}}={1:g}$ kV'.format(y[0], V_HF/1e3) for V_HF, y in product(V_HFs, Semits)]
    lss += ['-' for x in range(len(bdurs))]
    legplot(ax[6][0], lss, labs, loc=6)
    Mplot(ax[6][1], tVgZ, bdurequis+bdurs, lss, xlab, xunit, 'Bunch duration', 's', '')
    Mplot(ax[6][2], tVgZ, blenequis+blens, lss, xlab, xunit, 'Bunch length', 'm', '')

    # Emittance
    yunit = r'm $\pi$ rad'

    # Radial Emittance
    labs = ['Equilibrium']+[r'$\epsilon_0=$ {0:g} nm rad'.format(y[0]*1e9) for y in Xemits]
    lss = ['-.']+['-' for x in range(len(Xemits))]
    ylab, ylab2 = r'$\epsilon_x$', r'$\epsilon_x^*$'
    Mplot(ax[7][0], tAI, [Xemitequi]+Xemits, lss, '', '', ylab, yunit, labs)
    Mplot(ax[7][1], EAI, [Xemitequi]+Xemits, lss, '', '', '', '', '')
    Mplot(ax[7][2], tAI, [NXemitequi]+NXemits, lss, xlab, xunit, ylab2, yunit, '')
    Mplot(ax[7][3], EAI, [NXemitequi]+NXemits, lss, xlab2, xunit2, '', '', '')
    legs.append(ax[7][0].legend(fancybox=True, loc=2))

    # Axial Emittance
    labs = ['Limit']+[r'$\epsilon_0=$ {0:g} nm rad'.format(y[0]*1e9) for y in Yemits]
    lss = ['-.']+['-' for x in range(len(Yemits))]
    ylab, ylab2 = r'$\epsilon_y$', r'$\epsilon_y^*$'
    Mplot(ax[8][0], tAI, [Yemitequi]+Yemits, lss, '', '', ylab, yunit, labs)
    Mplot(ax[8][1], EAI, [Yemitequi]+Yemits, lss, '', '', '', '', '')
    Mplot(ax[8][2], tAI, [NYemitequi]+NYemits, lss, xlab, xunit, ylab2, yunit, '')
    Mplot(ax[8][3], EAI, [NYemitequi]+NYemits, lss, xlab2, xunit2, '', '', '')
    legs.append(ax[8][0].legend(fancybox=True, loc=1))

    # Longitudinal Emittance
    labs = ['Equilibrium']+[r'$\epsilon_0=$ {} \textperthousand'.format(y[0]) for y in Semits]
    lss = ['-.']+['-' for x in range(len(Semits))]
    ylab, yunit = r'$\delta_E=\frac{\sigma_E}{E_0}$', r' \textperthousand '
    Mplot(ax[9][0], tAI, [Semitequi]+Semits, lss, xlab, xunit, ylab, yunit, labs, rescaleY=False)
    Mplot(ax[9][1], EAI, [Semitequi]+Semits, lss, xlab2, xunit2, '', '', '', rescaleY=False)
    legs.append(ax[9][0].legend(fancybox=True, loc=1))

    [leg.get_frame().set_alpha(0.5) for leg in legs]
    return figs

def pltsim_quadscan(k, sigx, sigy):
    xlabel, xunit = 'Quadrupole strength', 'm'
    ylabel1, yunit1 = r'$\sigma_x^2$', r'$mm^2$'
    ylabel2, yunit2 = r'$\sigma_y^2$', r'$mm^2$'

    figs = [Figure()]
    ax = [figs[0].add_subplot(1, 2, i) for i in range(1, 3)]
    plot(ax[0], k, sigx*1e6, '-b', xlabel, xunit, ylabel1, yunit1, '', rescaleY=False)
    plot(ax[1], k, sigy*1e6, '-b', xlabel, xunit, ylabel2, yunit2, '', rescaleY=False)
    return figs
