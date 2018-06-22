# -*- coding: utf-8 -*-
''' accpy.jupyter.mechanics
author:
    Felix Kramer
'''
import numpy as np
from matplotlib.pyplot import subplots, cm, figure
from mpl_toolkits.mplot3d import axes3d, art3d


def initfigs(points, ms, bw, n=2, m=2, fs=[1920/120, 1080/120]):
    colors = cm.rainbow(np.linspace(0, 1, points))
    lines = []
    x = ['q', 't', 't', 'q']
    y = ['p', 'q', 'p', r'$E_{pot}$ & $E_{tot}$']
    fig1, axes = subplots(n, m, figsize=fs)
    axes  = [ax for subaxes in axes for ax in subaxes]
    for i, ax in enumerate(axes):
        ax.grid(True)
        ax.set_xlabel(x[i])
        ax.set_ylabel(y[i])
        lines.append([])
        for col in colors:
            lines[i].append(ax.plot([], [], '.', c=col, ms=ms)[0])
        if i in [1, 3]:
            ax.yaxis.tick_right()
            ax.yaxis.set_label_position('right')
    lines.append(axes[-1].plot([], [], '-', c=bw)[0])
    fig1.tight_layout(rect=[0, 0, 1, .9])
    fig2 = figure(figsize=fs)
    ax3D = fig2.add_subplot(111, projection='3d')
    fig2.tight_layout()
    return fig1, fig2, axes, lines, ax3D

def init(dt, T, q0, p0, points, bw, ms=1, fs=[12, 10]):
    t = np.arange(0, T, dt)
    iterations = len(t)
    q, p, U = [np.empty([iterations, points]) for i in range(3)]
    q[0, :], p[0, :] = q0, p0
    fig1, fig2, axes, lines, ax3D =  initfigs(points, ms, bw, n=2, m=2, fs=fs)
    return fig1, fig2, axes, lines, ax3D, t, iterations, q, p, U

def lims(mplotlims):
        scale = 1.021
        offset = (mplotlims[1] - mplotlims[0])*scale
        return mplotlims[1] - offset, mplotlims[0] + offset

def run(Tkin, Upot, eom, axes, lines, ax3D, t, iterations, q, p, U, dt, bw):
    for it in range(1, iterations):
        q0, p0 = q[it-1, :], p[it-1, :]
        U[it-1, :] = Upot(q0) + Tkin(p0)
        q[it, :], p[it, :] = eom(q0, p0)
    q0, p0 = q[it, :], p[it, :]
    U[it, :] = Upot(q0) + Tkin(p0)
    for i, [line0, line1, line2, line3] in enumerate(zip(lines[0], lines[1], lines[2], lines[3])):
        qdat, pdat = q[:, i], p[:, i]
        line0.set_data(qdat, pdat)
        line1.set_data(t, qdat)
        line2.set_data(t, pdat)
        line3.set_data(qdat, U[:, i])
    qmin, qmax, pmin, pmax = np.nanmin(q), np.nanmax(q), np.nanmin(p), np.nanmax(p)
    qrange = np.linspace(qmin, qmax, iterations)
    lines[4].set_data(qrange, Upot(qrange))
    for ax in axes:
        ax.relim()
        ax.autoscale_view(True,True,True)
    res = int(2e2)
    qmesh, pmesh = np.meshgrid(np.linspace(qmin, qmax, res), np.linspace(pmin, pmax, res))
    Hmesh = Tkin(pmesh) + Upot(qmesh)
    ax3D.clear()
    ax3D.plot_surface(qmesh, pmesh, Hmesh, rstride=8, cstride=8, alpha=0.3, cmap='hot')
    xlims, ylims, zlims = lims(ax3D.get_xlim()), lims(ax3D.get_ylim()), lims(ax3D.get_zlim())
    ax3D.contour(qmesh, pmesh, Hmesh, zdir='x', cmap=cm.coolwarm, offset=xlims[0])
    ax3D.contour(qmesh, pmesh, Hmesh, zdir='y', cmap=cm.coolwarm, offset=ylims[1])
    ax3D.contour(qmesh, pmesh, Hmesh, zdir='z', cmap=cm.coolwarm, offset=zlims[0])
    i = np.array([xlims[0], ylims[0], zlims[0]])
    f = np.array([xlims[0], ylims[0], zlims[1]])
    p = art3d.Poly3DCollection(np.array([[i, f]]))
    p.set_color(bw)
    ax3D.xaxis.pane.set_edgecolor(bw)
    ax3D.yaxis.pane.set_edgecolor(bw)
    ax3D.zaxis.pane.set_edgecolor(bw)
    ax3D.add_collection3d(p)
    ax3D.set_xlabel('q')
    ax3D.set_ylabel('p')
    ax3D.set_zlabel('E')
    ax3D.xaxis.pane.set_alpha(1)
    ax3D.yaxis.pane.set_alpha(1)
    ax3D.zaxis.pane.set_alpha(1)
    ax3D.xaxis.pane.fill = False
    ax3D.yaxis.pane.fill = False
    ax3D.zaxis.pane.fill = False
#     ax3D.xaxis.set_major_locator(MultipleLocator(5))
#     ax3D.yaxis.set_major_locator(MultipleLocator(30))
#     ax3D.zaxis.set_major_locator(MultipleLocator(1000))
    return