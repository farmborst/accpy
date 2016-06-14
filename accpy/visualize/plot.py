#!/usr/bin/python3
# -*- coding: utf-8 -*-
''' evaluation of measured data to get dispersion in transferline
author:
    Felix Kramer
version:
    Created         06.11.2015
    Last Update     06.11.2015
'''
from __future__ import division
from numpy import (size, cumsum, nanmax, nanmin, concatenate, empty, linspace,
                   array, arange, abs as npabs, sqrt, sin, cos, arccos as acos)
from matplotlib.figure import Figure as figure
from matplotlib.pyplot import cm
from matplotlib.gridspec import GridSpec
from ..simulate import const
from ..simulate.rmatrices import UCS2R
from ..simulate.tracking import trackpart


def drawlattice(ax, optic, diagnostics, ymin, ymax, height):
    ''' function for drawing given optic into figure
    inputs: ax      handle for figure axes
            optic   lattice to be drawn
            ymin    min of other functions (e.g. twiss) in same figure
            ymax    max of other functions (e.g. twiss) in same figure
            height  0 to 1 where to place lattice (1 at top)
    '''
    l = size(optic, 1)
    s = cumsum(optic[1, :])
    d = ymax-ymin
    m = ymin + d*height
    h = m + .05*d
    l = m - .05*d
    h2 = m + .025*d
    l2 = m - .025*d
    diagnostic = 0
    for i in range(size(optic, 1)):
        element = optic[0, i]
        beg = s[i] - optic[1, i]
        end = s[i]
        if element == 0:        # drift as black line
            x = [beg, end]
            y = [m, m]
            ax.plot(x, y, '-k')
        elif element == 1:      # dipole as blue box
            x = [end, beg, beg, end, end]
            y = [h, h, l, l, h]
            ax.plot(x, y, '-k')
        elif element == 2:      # edge nothing
            x = []
            y = []
        elif element == 3:      # radial focussing quad as box above
            if optic[4, i] == 0:
                x = [end, end, beg, beg, end]
                y = [l2, h2, h2, l2, l2]
                ax.plot(x, y, '--k')
            else:
                x = [end, end, beg, beg, end]
                y = [m, h, h, m, m]
                ax.plot(x, y, '-k')
                ax.text(beg+(end-beg)/2, h+0.01*d, '%.5f' % optic[4, i], rotation=90,
                        horizontalalignment='center',
                        verticalalignment='bottom')
        elif element == 4:      # axial focussing quad as box below
            if optic[4, i] == 0:
                x = [end, end, beg, beg, end]
                y = [l2, h2, h2, l2, l2]
                ax.plot(x, y, '--k')
            else:
                x = [end, end, beg, beg, end]
                y = [m, l, l, m, m]
                ax.plot(x, y, '-k')
                ax.text(beg+(end-beg)/2, m+0.01*d, '%.5f' % optic[4, i], rotation=90,
                        horizontalalignment='center',
                        verticalalignment='bottom')
        elif element == 5:      # rotator nothing
            x = []
            y = []
        elif element == 6:      # solenoid nothing
            x = []
            y = []
        elif element == 7:      # diagnostic as red vertical line
            x = [beg, end]
            y = [l, h]
            ax.plot(x, y, '-r')
            ax.text(beg, h, diagnostics[diagnostic], rotation=90,
                    horizontalalignment='center', verticalalignment='bottom')
            diagnostic += 1
    return


def plotoptic(UC, optic, diagnostics, s, xtwiss, ytwiss, xdisp):
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
    data = concatenate((xtwiss[0, 0, :], ytwiss[0, 0, :], xdisp.flatten()))
    ymin = nanmin(data)
    ymax = nanmax(data)
    ymin2 = nanmin(xdisp.flatten())
    ymax2 = nanmax(xdisp.flatten())
    if optic == 'radial':
        drawlattice(ax, UC, diagnostics, ymin, ymax, 0)
        ax.plot(s, xtwiss[0, 0, :], '-r', label=r'$\beta_x$')
        ax.plot(s, xtwiss[0, 1, :], '-c', label=r'$\alpha_x$')
        ax.plot([], [], '-m', label=r'$\gamma_x$')
        ax.set_ylabel(r'betatron function $\beta_x$ / (m)')
        ax2 = ax.twinx()
        ax2.plot(s, xtwiss[1, 1, :], '-m')
        ax2.set_ylabel(r'gamma function $\gamma_x$ / (m)', color='m')
        ax2.tick_params(axis='y', colors='m')
    if optic == 'axial':
        drawlattice(ax, UC, diagnostics, ymin, ymax, 0)
        ax.plot(s, ytwiss[0, 0, :], '-b', label=r'$\beta_y$')
        ax.plot(s, ytwiss[0, 1, :], '-c', label=r'$\alpha_y$')
        ax.plot([], [], '-m', label=r'$\gamma_y$')
        ax.set_ylabel(r'betatron function $\beta_y$ / (m)')
        ax2 = ax.twinx()
        ax2.plot(s, ytwiss[1, 1, :], '-m', label=r'$\gamma_y$')
        ax2.set_ylabel(r'gamma function $\gamma_y$ / (m)', color='m')
        ax2.tick_params(axis='y', colors='m')
    if optic == 'dispersion':
        drawlattice(ax, UC, diagnostics, ymin2, ymax2, 0)
        ax.plot(s, xdisp[0, :], '-g', label=r'$D_x$')
        ax.plot([], [], '-m', label=r'$D_x^\prime$')
        ax.set_ylabel(r'dispersion function $D_x$ / (m)')
        ax2 = ax.twinx()
        ax2.plot(s, xdisp[1, :], '-m', label=r'$D_x^\prime$')
        ax2.set_ylabel(r'derived dispersion function $D_x^\prime$ / (m)', color='m')
        ax2.tick_params(axis='y', colors='m')
    if optic == 'overview':
        drawlattice(ax, UC, diagnostics, ymin, ymax, 0)
        ax.plot(s, xtwiss[0, 0, :], '-r', label=r'$\beta_x$')
        ax.plot(s, ytwiss[0, 0, :], '-b', label=r'$\beta_y$')
        ax.plot([], [], '-g', label=r'$D_x$')
        ax.set_ylabel(r'betatron function $\beta_{x,y}$ / (m)')
        ax2 = ax.twinx()
        ax2.plot(s, xdisp[0, :], '-g', label=r'$D_x$')
        ax2.set_ylabel(r'dispersion function $D_x$ / (m)', color='g')
        ax2.tick_params(axis='y', colors='g')
    ax.set_xlabel(r'orbit position s / (m)')
    ax.set_xlim([0, nanmax(s)])
    leg = ax.legend(fancybox=True, loc='upper left')
    leg.get_frame().set_alpha(0.5)
    return fig


def plotopticpars_closed(xtwiss, xdisp, ytwiss, gamma, Qx, Xx, Jx, emiteqx,
                         tau_x, Qy, Xy, Jy, E, emiteqy, tau_y, alpha_mc,
                         eta_mc, gamma_tr, Q_s, Js, sigma_E, sigma_tau,
                         sigma_s, tau_s, U_rad, P_ges, E_c, lambda_c):
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
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
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
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
            verticalalignment='top', transform=ax.transAxes, fontsize=12)
    ax.text(0.35, 0.9, axipars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes, fontsize=12)
    ax.text(0.65, 0.9, lonpars, horizontalalignment='left',
            verticalalignment='top', transform=ax.transAxes, fontsize=12)
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
    ymin = nanmin(X.flatten())
    ymax = nanmax(X.flatten())
    fig = figure()
    ax = fig.add_subplot(1, 1, 1)
    drawlattice(ax, UC, diagnostics, ymin, ymax, 0)
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


def plottrajs(s, X, N_UC, rounds, envelope):
    figs = [figure() for i in range(7)]
    ax1 = [figs[i].add_subplot(1, 1, 1) for i in range(6)]
    ax2 = [figs[6].add_subplot(3, 3, i) for i in [1, 7, 2, 8, 3, 9]]
    ax3 = figs[6].add_subplot(3, 3, 5)
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
    labs = ['Ideal particle',
            '1 sigma particle',
            r'Envelope $E_{x,y}(s)=\sqrt{\epsilon_{x,y}\beta_{x,y}(s)+(\delta_ED_{x,y}(s))^2)}$',
            r'Envelope $E_{x}(s)=\sqrt{\epsilon_{x}\beta_{x}+(\delta_ED_{x}(s))^2}$',
            r'Envelope $E_{y}(s)=\sqrt{\epsilon_{y}\beta_{y}}$',
            'Ensemble']
    for i in range(6):
        c = next(color)
        y = []
        y_ideal = []
        y_sigma = []
        for j, traj in enumerate(X):
            for k in range(rounds):
                index = arange(len(s))+(len(s)-1)*k
                if j > 1:
                    y.append(traj[i, index]*1e3)
                elif j == 1:
                    y_sigma.append(traj[i, index]*1e3)
                else:
                    y_ideal.append(traj[i, index]*1e3)
        # ensemble trajectories
        [ax1[i].plot(s, y[l], '-', c=c) for l in range(len(y))]
        [ax2[i].plot(s, y[l], '-', c=c) for l in range(len(y))]
        ax1[i].plot([], [], '-', c=c, label=labs[5])
        ax2[i].plot([], [], '-', c=c, label=labs[5])
        # 1-sigma particle trajectories
        ax1[i].plot(s, y_sigma[0], '-b', label=labs[1])
        ax2[i].plot(s, y_sigma[0], '-b')
        [ax1[i].plot(s, y_sigma[l], '-b') for l in range(1, len(y_sigma))]
        [ax2[i].plot(s, y_sigma[l], '-b') for l in range(1, len(y_sigma))]
        # ideal particle trajectories
        ax1[i].plot(s, y_ideal[0], '-k', label=labs[0])
        if i == 0:
            ax1[i].plot([], [], '-r', label=labs[3])
        elif i == 2:
            ax1[i].plot([], [], '-r', label=labs[4])
        ax2[i].plot(s, y_ideal[0], '-k')
        leg = ax1[i].legend(fancybox=True, loc='upper right')
        leg.get_frame().set_alpha(0.5)
    ax1[0].plot(s, envelope[0, :], '-r', s, -envelope[0, :], '-r')
    ax1[2].plot(s, envelope[1, :], '-r', s, -envelope[1, :], '-r')
    ax2[0].plot(s, envelope[0, :], '-r', s, -envelope[0, :], '-r')
    ax2[2].plot(s, envelope[1, :], '-r', s, -envelope[1, :], '-r')
    ax3.plot([], [], '-k', label=labs[0])
    ax3.plot([], [], '-b', label=labs[1])
    ax3.plot([], [], '-r', label=labs[2])
    ax3.get_xaxis().set_visible(False)
    ax3.get_yaxis().set_visible(False)
    ax3.axis('off')
    leg = ax3.legend(fancybox=True, loc='center')
    return figs


def plotphasespace(s, X, rounds, xtwiss, emittx, ytwiss, emitty):
    fig = figure()
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
    x, xp, y, yp = twissellipse(xtwiss[:, :, 0], emittx, ytwiss[:, :, 0], emitty)
    ax[0].plot(x, xp, '-g')
    ax[1].plot(y, yp, '-g')
    ax2.plot([], [], '.b', label='Ensemble')
    ax2.plot([], [], '-g', label='Twiss ellipsis')
    ax2.plot([], [], 'xr', label='1 sigma particle')
    ax2.axis('off')
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
