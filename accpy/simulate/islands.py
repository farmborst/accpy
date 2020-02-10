# -*- coding: utf-8 -*-
"""accpy.simulate.islands
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import division, print_function
from pyfftw import empty_aligned
from pyfftw.pyfftw import FFTW
from numpy import (abs as npabs, dot, roll, shape, zeros, empty, array, mean,
                   where, sort, diff, argmax, linspace, concatenate, isnan, pi,
                   logical_or, delete, nanmin, nanmax, add, nan, int32, arange,
                   arctan2, argsort, sqrt, diag, sum as npsum, log2)
from scipy.optimize import curve_fit
from matplotlib.cm import rainbow, ScalarMappable, cool
from matplotlib.pyplot import tight_layout
from matplotlib.colors import Normalize



###############################################################################
# POST PROCESSING
###############################################################################

def PolyArea(x, y):
    # find "quick and dirty" center to enable angle calculation
    # cx, cy = np.mean(x), np.mean(y)  # Not so good an approximation, since points not equally spaced
    cx = nanmax(x) - npabs(nanmax(x) - nanmin(x)) / 2
    cy = nanmax(y) - npabs(nanmax(y) - nanmin(y)) / 2
    x, y = x - cx, y - cy
    # sort points given by x and y arrays by azimut anglr
    indices = argsort(arctan2(x, y))
    x, y = x[indices], y[indices]
    # calculate area using shoelace formula
    # https://stackoverflow.com/questions/24467972/calculate-area-of-polygon-given-x-y-coordinates
    # total_area = np.abs(np.dot(x, np.roll(y,1)) - np.dot(y, np.roll(x, 1)))/2
    correction = x[-1] * y[0] - y[-1] * x[0]
    main_area = dot(x[:-1], y[1:]) - dot(y[:-1], x[1:])
    total_area = 0.5*npabs(main_area + correction)
    return array([cx, cy]), total_area


def dist(a, b):
    return sqrt(npsum((a - b)**2))


def islandsloc(data, PPdata):
    PPdata['Nturns'], PPdata['Nparticles'] = shape(data['x'])
    PPdata['IDs_all'] = arange(PPdata['Nparticles'], dtype=int32)
    
    
    for sub in [''] + ['_rev' + str(i) for i in range(PPdata['resonance'])]:
        PPdata['C' + sub] = empty([PPdata['Nparticles'], 2])
        
        for s in ['A', 'fx_meas', 'fx', 'Qx', 'QxTRIBs', 'fy_meas', 'fy', 'Qy', 'QyTRIBs']:
                PPdata[s + sub] = empty(PPdata['Nparticles'])
                PPdata[s + sub].fill(nan)
    
    T = zeros(PPdata['Nparticles'])
    for i in PPdata['IDs_all']:
        PPdata['C'][i, :], PPdata['A'][i] = PolyArea(data['x'][::PPdata['resonance'], i], data['xp'][::PPdata['resonance'], i])
        if dist(PPdata['C'][i, :], array([0, 0])) > PPdata['minsep']:  # minsep – distance island to center
            T[i] = 1  # 0=normal, 1=island
    
    PPdata['IDs_isla'], noislandIDs = where(T == 1)[0], where(T == 0)[0]
    if len(PPdata['IDs_isla']) == 0:
        PPdata['IDs_core'], PPdata['IDs_encl'] = noislandIDs, []
    else:
        Asort = sort(PPdata['A'][noislandIDs])
        Astep = diff(Asort)
        gtisteps = Astep > max(PPdata['resonance'] * PPdata['A'][PPdata['IDs_isla']])
        if not gtisteps.any():
            PPdata['IDs_core'], PPdata['IDs_encl'] = noislandIDs, []
        else:
            imaxstep = argmax(gtisteps)
            Amaxcent = Asort[imaxstep]
            T[noislandIDs[PPdata['A'][noislandIDs] > Amaxcent]] = 2
            PPdata['IDs_core'], PPdata['IDs_encl'] = where(T == 0)[0], where(T == 2)[0]
                      
    PPdata['IDs_core'] = array(PPdata['IDs_core'], dtype=int32)
    PPdata['IDs_isla'] = array(PPdata['IDs_isla'], dtype=int32)
    PPdata['IDs_encl'] = array(PPdata['IDs_encl'], dtype=int32)
    
    # sort according to action
    PPdata['IDs_core'] = PPdata['IDs_core'][argsort(PPdata['A'][PPdata['IDs_core']])]
    PPdata['IDs_isla'] = PPdata['IDs_isla'][argsort(PPdata['A'][PPdata['IDs_isla']])]
    PPdata['IDs_encl'] = PPdata['IDs_encl'][argsort(PPdata['A'][PPdata['IDs_encl']])]
    
    PPdata['JxCoreMax'] = nanmax(PPdata['A'][PPdata['IDs_core']])
#     PPdata['3Jx+JxCore'] = PPdata['resonance'] * PPdata['A'] + PPdata['JxCoreMax']
    PPdata['3Jx+JxCore'] = PPdata['JxCoreMax']
    
    # add data for each island bucket
    for j in range(PPdata['resonance']):
        rev = '_rev' + str(j)
        for i in PPdata['IDs_isla']:
            PPdata['C' + rev][i, :], PPdata['A' + rev][i] = PolyArea(data['x'][j::PPdata['resonance'], i], data['xp'][j::PPdata['resonance'], i])
#         PPdata['3Jx+JxCore' + rev] = PPdata['resonance'] * PPdata['A' + str(j)] + PPdata['JxCoreMax']
        PPdata['3Jx+JxCore'] += PPdata['A' + rev]
        
    
    return


def getmyfft(turns, frev):
    fd = linspace(0, frev/2, int(turns/2))
    dQ = linspace(0, 1/2, int(turns/2))
    # get frequency vector with negative frequencies
    fdn = concatenate((fd, -fd[::-1]))
    a = empty_aligned(turns, dtype='complex128')
    b = empty_aligned(turns, dtype='complex128')
    myfft = FFTW(a, b, threads=2)
    return dQ, fd, fdn, myfft

def getfreq(data, myfft, clip):
    return npabs(myfft(data)[1:clip])

def tunes(data, PPdata):
#     dQ, fd, fdn, myfft = getmyfft(PPdata['Nturns'], PPdata['frev'])
#     clip = int(PPdata['Nturns']/2)
#     PPdata['Qx'] = array([dQ[argmax(getfreq(data['x'][:, i], myfft, clip))] for i in PPdata['IDs_all']])
#     PPdata['Qy'] = array([dQ[argmax(getfreq(data['y'][:, i], myfft, clip))] for i in PPdata['IDs_all']])
#    for res in [3]:
#        N = int(data['Nturns']/res)
#        dQ, fd, fdn, myfft = getmyfft(N, frev)
#        for island in range(res):
#            Qstr = 'Q{}_{}'.format(res, island + 1)
#            data[Qstr] = array([dQ[argmax(getfreq(data['x'][:N, i][island::res], myfft, N))] for i in PPdata['IDs_all']])
    
#     PPdata['pow2max'] = int(log2(PPdata['Nturns']))
#     PPdata['Nturns_pow2max'] = int(2**PPdata['pow2max'])
#     PPdata['Nturns_pow2max_half'] = int(PPdata['Nturns_pow2max']/2)

    Nturns = int(PPdata['Nturns'])
    if (Nturns % 2) == 0:
        # Nturns is even
        pass
    else:
        Nturns = Nturns - 1
    Nturns2 = int(Nturns / 2)
    
    fsamp = PPdata['frev']
    dQ, fd, fdn, myfft = getmyfft(Nturns, fsamp)
    for u in ['x', 'y']:
        fcalc = 'f' + u
        fmeas = fcalc + '_meas'
        Qcalc = 'Q' + u
        for p in PPdata['IDs_all']:
            # calculate FFT
            fftn = npabs(myfft(data[u][:, p][:Nturns]))

            # cancel DC part
            fftn[0] = 0

            # only positive frequencies
            fft = fftn[:Nturns2]

            # calculate peaks frequency
            ipeak = argmax(fft)
            PPdata[fmeas][p] = fd[ipeak]
            PPdata[fcalc][p] = fsamp - PPdata[fmeas][p]
            PPdata[Qcalc][p] = PPdata[fcalc][p] / fsamp
    
    if PPdata['resonance'] > 0:
        
        Nturns = int(PPdata['Nturns'] / PPdata['resonance'])
        if (Nturns % 2) == 0:
            # Nturns is even
            pass
        else:
            Nturns = Nturns - 1
        Nturns2 = int(Nturns / 2)
    
        fsamp = PPdata['frev'] / PPdata['resonance']
        dQ, fd, fdn, myfft = getmyfft(Nturns, fsamp)
        for j in range(PPdata['resonance']):
            rev = '_rev' + str(j)
            for u in ['x', 'y']:
                fcalc = 'f' + u + rev
                fmeas = 'f' + u + '_meas' + rev
                Qcalc = 'Q' + u + rev
                Qtribs = 'Q' + u + 'TRIBs' + rev
                for p in PPdata['IDs_isla']:
                    # calculate FFT
                    fftn = npabs(myfft(data['x'][j::PPdata['resonance'], p][:Nturns]))

                    # cancel DC part
                    fftn[0] = 0

                    # only positive frequencies
                    fft = fftn[:Nturns2]

                    # calculate peaks frequency
                    ipeak = argmax(fft)
                    PPdata[fmeas][p] = fd[ipeak]
                    PPdata[fcalc][p] = fsamp - fd[ipeak]
                    PPdata[Qcalc][p] = fd[ipeak] / fsamp
                    PPdata[Qtribs][p] = (2 * fsamp + fd[ipeak]) / PPdata['frev']

    return

def findlost(data, PPdata):
    lost = []
    lostIDs = []
    for part in range(data['Particles'][0]):
        nansat = isnan(data['x'][:, part])
        if nansat.any():
            lostat = where(nansat)[0][0]
            lost.append(array([part, lostat], dtype=int32))
            lostIDs.append(part)
    PPdata['lost'] = array(lost, dtype=int32)
    PPdata['IDs_lost'] = array(lostIDs, dtype=int32)
    return

def evaltrackdat(data, resonance=0, minsep=5e-3):  
    # new dict for Post Processing Results
    PPdata = {}
    
    L = data['PassLength'][0]
    PPdata['Trev'] = L/299792458
    PPdata['frev'] = 1/PPdata['Trev']/1e3   # kHz
    PPdata['resonance'] = resonance
    PPdata['minsep'] = minsep
    
    islandsloc(data, PPdata)
    tunes(data, PPdata)
    findlost(data, PPdata)
    
    data['myPP'] = PPdata
    return


###############################################################################
# VISUALIZATION
###############################################################################


def getquadsext(lines, latticename):
#    props=dict(boxstyle='round', alpha=1)
#    line1 = ax.text(1.02, 1, '', va='top', ha='left', fontproperties='monospace', bbox=props, color='r', transform=ax.transAxes)
#    line2 = ax.text(1.02, .76, '', va='top', ha='left', fontproperties='monospace', bbox=props, color='g', transform=ax.transAxes)
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
    lines[0].set_text(quadstr.rstrip('\n'))
    lines[1].set_text(sextstr.rstrip('\n'))
    return


def getrdts(line, twissdat):
#    props=dict(boxstyle='round', alpha=1)
#    line = ax.text(1.02, .55, string.rstrip('\n'), va='top', ha='left', fontproperties='monospace', bbox=props, color='y', transform=ax.transAxes)
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
    line.set_text(string.rstrip('\n'))
    return


def trackplot(ax, data, turns=False, xy=False, fs=[16, 9], showlost=False,
              everyxturn=[0, 1], ms=1):
    x, y = xy
    colors = rainbow(linspace(0, 1, data['Particles'][0]))
    IDs = data['IDs_all'].copy()
    if not showlost:
        IDs = delete(IDs, data['lostIDs'])
    for part, col in zip(IDs, colors):
        i, f = everyxturn
        xdat, ydat = data[x][i::f, part], data[y][i::f, part]
        ax.plot(xdat*1e3, ydat*1e3, '.', color=col, ms=ms)
    ax.set_xlabel(r'$x$ / (mm)')
    ax.set_ylabel(r'$x^\prime$ / (mrad)')
    return


def islandsplot(ax, data, showlost=False):
    PPdata = data['myPP']
    datx, daty = data['x'].copy(), data['xp'].copy()
    if showlost:
        for ID in PPdata['IDs_isla']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.r')
        for ID in PPdata['IDs_core']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.b')
        for ID in PPdata['IDs_encl']:
            ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.g')
    else:
        lost = data['lostIDs']
        for ID in PPdata['IDs_isla']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.r')
        for ID in PPdata['IDs_core']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.b')
        for ID in PPdata['IDs_encl']:
            if ID not in lost:
                ax.plot(datx[:, ID]*1e3, daty[:, ID]*1e3, '.g')
    ax.set_xlabel(r'$x$ / (mm)')
    ax.set_ylabel(r'$x^\prime$ / (mrad)')
    return


def tuneplot(ax1, ax2, data, particleIDs='allIDs', integer=1, addsub=add,
             clipint=True, showlost=False, QQ='Qx', ms=1, clip=[0], showfit=False):
    particleIDs = data[particleIDs]
    if not showlost:
        try:
            lost = data['lost'][:, 0]
        except:
            lost = data['lost']
        clip = concatenate([clip, lost])
    particleIDs = delete(particleIDs, clip)
    Q = addsub(integer, data[QQ][particleIDs])
    if clipint:
        zeroQ = where(logical_or(logical_or(Q == 0.0, Q == 1.0), Q == 0.5))[0]
        if len(zeroQ) > 0:  # trim reference particle with zero tune
            Q = delete(Q, zeroQ)
            particleIDs = delete(particleIDs, zeroQ)
    Qmin, Qmax = nanmin(Q), nanmax(Q)
    Qdif = Qmax - Qmin
    if Qdif == 0.0:
        Qmin -= Qmin/1e4
        Qmax += Qmax/1e4
        Qdif = Qmax - Qmin
    colors = cool((Q - Qmin) / Qdif)
    for i, ID in enumerate(particleIDs):
        ax1.plot(data['x'][:, ID]*1e3, data['xp'][:, ID]*1e3, '.', c=colors[i], ms=ms)
    if showlost:
        for ID in lost:
            ax1.plot(data['x'][:, ID]*1e3, data['xp'][:, ID]*1e3, '.', c='gray', ms=ms)
    sm = ScalarMappable(cmap=rainbow, norm=Normalize(vmin=Qmin, vmax=Qmax))
    sm._A = []
    ax1.set_xlabel(r'Position $x$ / (mm)')
    ax1.set_ylabel(r'Angle $x^\prime$ / (mrad)')
    emittance = data['A'][particleIDs]/pi
    action = emittance/2
    
    # tune shift with action
    fitfun = lambda x, a, b: a + b*x
    popt, pcov = curve_fit(fitfun, action, Q)
    perr = sqrt(diag(pcov))
    action2 = linspace(nanmin(action), nanmax(action), 1000)
    fit1 = fitfun(action2, *popt)
    print(popt[1]*1e-6*1250)

    for i, ID in enumerate(particleIDs):
        ax2.plot(action[i]*1e6, Q[i], 'o', c=colors[i], ms=ms + 1)
    if showfit:
        ax2.plot(action2*1e6, fit1, '-k', lw=1, label=r'fit with $TSWA=${:.4}$\pm${:.1} (kHz mm$^-$$^2$mrad$^-$$^2$)'.format(popt[1]*1e-6*1250, perr[1]*1e-6*1250))
#    leg = ax2.legend()
#    leg.get_frame().set_alpha(0)
    ax2.set_ylim([Qmin, Qmax])
#    ax2.yaxis.tick_right()
    ax2.set_ylabel(r'Fractional Tune d$Q$')
#    ax2.yaxis.set_label_position('right')
    ax2.set_xlabel(r'Action $J_x$ / (mm$\cdot$mrad)')
    tight_layout()
    return
