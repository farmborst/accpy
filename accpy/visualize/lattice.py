# -*- coding: utf-8 -*-
''' accpy.measure.tunes
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (size, cumsum, nanmin, nanmax)
from matplotlib.pyplot import Figure


def latticeplot(optic, diagnostics):
    ymin, ymax = 0, 100
    fig = Figure(frameon=False)
    ax = fig.add_subplot(111)
    drawlattice(ax, optic, diagnostics, [ymin, ymax], .3)
    s = cumsum(optic[1, :])
    ax.set_xlim(0, s[-1])
    ax.set_ylim(ymin, ymax)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return fig


def drawlattice(ax, optic, diagnostics, data, height):
    ''' function for drawing given optic into figure
    inputs: ax      handle for figure axes
            optic   lattice to be drawn
            ymin    min of other functions (e.g. twiss) in same figure
            ymax    max of other functions (e.g. twiss) in same figure
            height  0 to 1 where to place lattice (1 at top)
    '''
    ymin = nanmin([nanmin(x) for x in data])
    ymax = nanmax([nanmax(x) for x in data])
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


def drawlattice2d(ax, lattice):
    ''' function for drawing 2d model of given lattice into figure
    inputs: ax      handle for figure axes
            optic   lattice to be drawn
    '''
    xs = []
    ys = []

    for element in lattice:
        t = element[0]
        L = element[L]
        if t == 0:      # drift
            xs.append([0, L])
            ys.append([0, 0])
        elif t == 1:    # dipole
            xs.append(0)
            ys.append(L)
        elif t == 2:    # edge
            xs.append(0)
            ys.append(L)
        elif t == 3:    # radial focussing quad
            xs.append(0)
            ys.append(L)
        elif t == 4:    # axial focussing quad
            xs.append(0)
            ys.append(L)
        elif t == 5:    # rotator (Rskew=Rrot(-alpha*RQ*Rrot(alpha)))
            xs.append(0)
            ys.append(L)
        elif t == 6:    # solenoid
            xs.append(0)
            ys.append(L)
        elif t == 7:    # diagnostic (identity matrix for R)
            xs.append(0)
            ys.append(L)
        elif t == 8:    # gradientdipole
            xs.append(0)
            ys.append(L)
    [ax.plot(x, y) for x, y in zip(xs, ys)]
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return


def drawlattice3d(ax, lattice):
    ''' function for drawing 2d model of given lattice into figure
    inputs: ax      handle for figure axes
            optic   lattice to be drawn
    '''
    xs = []
    ys = []

    for element in lattice:
        t = element[0]
        L = element[L]
        if t == 0:      # drift
            xs.append([0, L])
            ys.append([0, 0])
        elif t == 1:    # dipole
            xs.append(0)
            ys.append(L)
        elif t == 2:    # edge
            xs.append(0)
            ys.append(L)
        elif t == 3:    # radial focussing quad
            xs.append(0)
            ys.append(L)
        elif t == 4:    # axial focussing quad
            xs.append(0)
            ys.append(L)
        elif t == 5:    # rotator (Rskew=Rrot(-alpha*RQ*Rrot(alpha)))
            xs.append(0)
            ys.append(L)
        elif t == 6:    # solenoid
            xs.append(0)
            ys.append(L)
        elif t == 7:    # diagnostic (identity matrix for R)
            xs.append(0)
            ys.append(L)
        elif t == 8:    # gradientdipole
            xs.append(0)
            ys.append(L)
    [ax.plot(x, y) for x, y in zip(xs, ys)]
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return
