# -*- coding: utf-8 -*-
"""accpy.simulate.islands
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division, print_function
from pyfftw import empty_aligned
from pyfftw.pyfftw import FFTW
from numpy import (abs as npabs, dot, roll, shape, zeros, empty, array, mean,
                   where, sort, diff, argmax, linspace, concatenate, isnan,
                   logical_or, delete, nanmin, nanmax, add, nan)
from matplotlib.mlab import dist, find
from matplotlib.cm import rainbow, ScalarMappable
from matplotlib.pyplot import tight_layout
from matplotlib.colors import Normalize



###############################################################################
# POST PROCESSING
###############################################################################

def PolyArea(x, y):
    # shoelace formula from
    # https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    return npabs(dot(x, roll(y,1)) - dot(y, roll(x, 1)))/2


def islandsloc(data, resonance, minsep=5e-3):
    data['turns'], data['particles'] = shape(data['x'])
    data['allIDs'] = range(data['particles'])
    data['A'], C, T = zeros(data['particles']), empty([data['particles'], 2]), zeros(data['particles'])
    for i in data['allIDs']:  # [0, 9, 21, 22, 92]
        x = data['x'][::3, i]
        xp = data['xp'][::3, i]
        C[i, :] = array([mean(x), mean(xp)])
        data['A'][i] = PolyArea(x, xp)
        if dist(C[i], array([0, 0])) > minsep:  # minsep – distance island to center
            T[i] = 1  # 0=normal, 1=island
    data['islandIDs'], noislandIDs = where(T == 1)[0], where(T == 0)[0]
    if len(data['islandIDs']) == 0:
        data['centerIDs'], data['enclosingIDs'] = noislandIDs, []
    else:
        Asort = sort(data['A'][noislandIDs])
        Astep = diff(Asort)
        gtisteps = Astep > max(data['A'][data['islandIDs']])
        if not gtisteps.any():
            data['centerIDs'], data['enclosingIDs'] = noislandIDs, []
        else:
            imaxstep = argmax(gtisteps)
            Amaxcent = Asort[imaxstep]
            T[noislandIDs[data['A'][noislandIDs] > Amaxcent]] = 2
            data['centerIDs'], data['enclosingIDs'] = where(T == 0)[0], where(T == 2)[0]
    return


def getmyfft(turns, frev):
    fd = linspace(0, frev/2/1e3, int(turns/2))
    dQ = linspace(0, 1/2, int(turns/2))
    fdn = concatenate((fd, -fd[::-1]))
    a = empty_aligned(turns, dtype='complex128')
    b = empty_aligned(turns, dtype='complex128')
    myfft = FFTW(a, b, threads=2)
    return dQ, fd, fdn, myfft

def getfreq(data, myfft, clip):
    return npabs(myfft(data)[1:clip])

def tunes(data, frev):
    dQ, fd, fdn, myfft = getmyfft(data['turns'], frev)
    clip = int(data['turns']/2)
    data['Q'] = array([dQ[argmax(getfreq(data['x'][:, i], myfft, clip))] for i in data['allIDs']])
#    for res in [3]:
#        N = int(data['turns']/res)
#        dQ, fd, fdn, myfft = getmyfft(N, frev)
#        for island in range(res):
#            Qstr = 'Q{}_{}'.format(res, island + 1)
#            data[Qstr] = array([dQ[argmax(getfreq(data['x'][:N, i][island::res], myfft, N))] for i in data['allIDs']])
    return

def evaltrackdat(data, resonance, frev):
    islandsloc(data, resonance)
    tunes(data, frev)
    return


###############################################################################
# VISUALIZATION
###############################################################################


def getquadsext(ax, latticename):
    quads = ['Q1_F', 'Q2_D', 'Q3D_D', 'Q4D_F', 'Q3T_D', 'Q4T_F', 'Q5T_D']
    sexts = ['S1_F', 'S2_D', 'S3D_D', 'S4D_F', 'S3T_D', 'S4T_F']

    quadstr, sextstr = '', ''
    with open(latticename, 'r') as fh:
        for line in fh:
            for quad in quads:
                if line.startswith(quad) or line.startswith('"' + quad + '"'):
                    try:
                        K1 = float(line.split('K1=')[1].split(',')[0])
                    except:
                        K1 = 0.0
                    quadstr += '{:<5}  K1={:<+.4}\n'.format(quad, K1)
            for sext in sexts:
                if line.startswith(sext) or line.startswith('"' + sext + '"'):
                    try:
                        K2 = float(line.split('K2=')[1].split(',')[0])
                    except:
                        K2 = 0.0
                    sextstr += '{:<5}  K2={:<+.4}\n'.format(sext, K2)
    props=dict(boxstyle='round', alpha=0.5)
    ax.text(1.02, 1, quadstr.rstrip('\n'), va='top', ha='left', fontproperties='monospace', bbox=props, color='r', transform=ax.transAxes)
    ax.text(1.02, .76, sextstr.rstrip('\n'), va='top', ha='left', fontproperties='monospace', bbox=props, color='g', transform=ax.transAxes)
    return

def getrdts(ax, twissdat):
    rdts = ['h11001', 'h00111', 'h20001', 'h00201', 'h10002', 'h21000', 'h30000', 'h10110', 'h10020', 'h10200','h22000',
            'h11110', 'h00220', 'h31000', 'h40000', 'h20110', 'h11200', 'h20020', 'h20200', 'h00310', 'h00400']
    string = ''
    for rdt in rdts:
        string += rdt + ' = {:<+9.4}'.format(twissdat['Re' + rdt][0])
        imp = twissdat['Im' + rdt][0]
        if npabs(imp) < 1e-10:
            string += '\n'
            continue
        if imp < 0:
            string += ' - i * {:.4}\n'.format(npabs(imp))
        else:
            string += ' + i * {:.4}\n'.format(npabs(imp))
    props=dict(boxstyle='round', alpha=0.5)
    ax.text(1.02, .55, string.rstrip('\n'), va='top', ha='left', fontproperties='monospace', bbox=props, color='y', transform=ax.transAxes)
    return

def trackplot(ax, datadict, turns=False, xy=False, fs=[16, 9], showlost=False,
              everyxturn=[0, 1]):
    x, y = xy
    colors = rainbow(linspace(0, 1, datadict['Particles'][0]))
    lost = []
    for part, col in enumerate(colors):
        i, f = everyxturn
        xdat, ydat = datadict[x][i::f, part], datadict[y][i::f, part]
        nansat = isnan(xdat)
        if nansat.any():
            lostat = where(nansat)[0][0]
            lost.append(array([part, lostat]))
            if showlost:
                ax.plot(xdat*1e3, ydat*1e3, '.', color=col)
        else:
            ax.plot(xdat*1e3, ydat*1e3, '.', color=col)
    ax.set_xlabel('x / (mm)')
    ax.set_ylabel('x\' / (mrad)')
    datadict['lost'] = array(lost)
    return

def islandsplot(ax, data, showlost=False):
    datx, daty = data['x'].copy(), data['xp'].copy()
    if showlost:
        for ID in data['islandIDs']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.r')
        for ID in data['centerIDs']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.b')
        for ID in data['enclosingIDs']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.g')
    else:
        lost = data['lost'][:, 0]
        for ID in data['islandIDs']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.r')
        for ID in data['centerIDs']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.b')
        for ID in data['enclosingIDs']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.g')
    ax.set_xlabel('x / (mm)')
    ax.set_ylabel('x\' / (mrad)')
    return

def tuneplot(ax1, ax2, data, particleIDs='allIDs', integer=1, addsub=add,
             clipint=True, showlost=False):
    particleIDs = data[particleIDs]
    Q = addsub(integer, data['Q'][particleIDs])
    if clipint:
        zeroQ = find(logical_or(logical_or(Q == 0.0, Q == 1.0), Q == 0.5))
        if len(zeroQ) > 0:  # trim reference particle with zero tune
            Q = delete(Q, zeroQ)
            particleIDs = delete(particleIDs, zeroQ)
    Qmin, Qmax = nanmin(Q), nanmax(Q)
    Qdif = Qmax - Qmin
#    fmt = '%.3f'
    if Qdif == 0.0:
        Qmin -= Qmin/1e4
        Qmax += Qmax/1e4
        Qdif = Qmax - Qmin
#        fmt = '%.5f'
    colors = rainbow((Q - Qmin) / Qdif)
    if showlost:
        for i, ID in enumerate(particleIDs):
            ax1.plot(data['x'][:, ID]*1e3, data['xp'][:, ID]*1e3, '.', c=colors[i])
    else:
        lost = data['lost']
        if len(lost) > 0:
            lost = lost[:, 0]
        for i, ID in enumerate(particleIDs):
            if ID not in lost:
                ax1.plot(data['x'][:, ID]*1e3, data['xp'][:, ID]*1e3, '.', c=colors[i])
    sm = ScalarMappable(cmap=rainbow, norm=Normalize(vmin=Qmin, vmax=Qmax))
    sm._A = []
#    cb = colorbar(sm, ax=ax2, format=fmt, ticks=[])
#    cb = colorbar(sm, ax=ax2, ticks=[])
#    cb.set_label('dQ')
    ax1.set_xlabel('x / (mm)')
    ax1.set_ylabel('x\' / (mrad)')
    for i, ID in enumerate(particleIDs):
        initialamp = data['x'][0, ID]
        ax2.plot(initialamp*1e3, Q[i], 'o', c=colors[i])
#    ax2.yaxis.set_ticklabels([])
    ax2.set_ylim([Qmin, Qmax])
    ax2.yaxis.tick_right()
    ax2.set_ylabel('dQ / (a.u.)')
    ax2.yaxis.set_label_position("right")
    ax2.set_xlabel('Initial amplitude x / (mm)')
    tight_layout()
    return
