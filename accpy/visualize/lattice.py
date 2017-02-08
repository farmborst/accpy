# -*- coding: utf-8 -*-
''' accpy.measure.tunes
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from numpy import (size, cumsum, nanmin, nanmax)
from matplotlib.pyplot import Figure
from ..dataio.hdf5 import confload


def latticeplot(optic, diagnostics, size=None):
    ymin, ymax = 0, 100
    fig = Figure(frameon=False)
    if size is not None:
        fig.set_figwidth(size[0])
        fig.set_figheight(size[1])
    ax = fig.add_subplot(111)
    drawlattice(ax, optic, diagnostics, [ymin, ymax], .3, checkconf=False)
    s = cumsum(optic[1, :])
    ax.set_xlim(0, s[-1])
    ax.set_ylim(ymin, ymax)
    ax.axis('off')
    ax.get_xaxis().set_visible(False)
    ax.get_yaxis().set_visible(False)
    return fig


def annotate(ax, x, y, text):
        ax.text(x, y, text, rotation=90, horizontalalignment='center',
                verticalalignment='bottom')
        return


def drawlattice(ax, optic, diagnostics, data, height, checkconf=True):
    ''' function for drawing given optic into figure
    inputs: ax      handle for figure axes
            optic   lattice to be drawn
            ymin    min of other functions (e.g. twiss) in same figure
            ymax    max of other functions (e.g. twiss) in same figure
            height  0 to 1 where to place lattice (1 at top)
    '''
    if checkconf:
        varlist, vallist = confload('./settings.conf')
        showlattice = vallist[varlist.index('showlattice')]
        showquadstrength = vallist[varlist.index('showquadstrength')]
        showdiagnostic = vallist[varlist.index('showdiagnostic')]
        showquadnr = vallist[varlist.index('showquadnr')]
    else:
        showlattice = 1
        showdiagnostic = 1
        showquadstrength = 0
        showquadnr = 1
    if showlattice == 1:
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
        quadN = 0
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
                quadN += 1
                if optic[4, i] == 0:
                    x = [end, end, beg, beg, end]
                    y = [l2, h2, h2, l2, l2]
                    ax.plot(x, y, '-k', color = '0.75')
                else:
                    x = [end, end, beg, beg, end]
                    y = [m, h, h, m, m]
                    ax.plot(x, y, '-k')
                if showquadnr and showquadstrength:
                    qtext = 'Q{0:g} {1:.5f}'.format(quadN ,optic[4, i])
                elif showquadnr:
                    qtext = 'Q{0:g}'.format(quadN)
                elif showquadstrength:
                    qtext = '{0:.5f}'.format(optic[4, i])
                else:
                    qtext = ''
                annotate(ax, beg+(end-beg)/2, h+0.01*d, qtext)
            elif element == 4:      # axial focussing quad as box below
                quadN += 1
                if optic[4, i] == 0:
                    x = [end, end, beg, beg, end]
                    y = [l2, h2, h2, l2, l2]
                    ax.plot(x, y, '-k', color = '0.75')
                else:
                    x = [end, end, beg, beg, end]
                    y = [m, l, l, m, m]
                    ax.plot(x, y, '-k')
                if showquadnr and showquadstrength:
                    qtext = 'Q{0:g}: {1:.5f}'.format(quadN ,optic[4, i])
                elif showquadnr:
                    qtext = 'Q{0:g}'.format(quadN)
                elif showquadstrength:
                    qtext = '{0:.5f}'.format(optic[4, i])
                else:
                    qtext = ''
                annotate(ax, beg+(end-beg)/2, m+0.01*d, qtext)
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
                if showdiagnostic:
                    annotate(ax, beg, h, diagnostics[diagnostic])
                diagnostic += 1
    return