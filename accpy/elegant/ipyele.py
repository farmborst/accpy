# -*- coding: utf-8 -*-
"""accpy.elegant.ipyele
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import print_function, division, unicode_literals
from subprocess import Popen, PIPE, STDOUT
from numpy import (shape, array, max as npmax, argmax, roll, float64, core,
                   min as npmin, linspace, where, empty, nan)
from matplotlib.pylab import (plot, subplot, xlabel, ylabel, twinx, gca, xlim,
                              ylim, annotate, tight_layout)
from matplotlib.patches import Polygon, Rectangle, PathPatch
from matplotlib.path import Path
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, VPacker
from matplotlib import cm
from pylab import figure
from . import sdds
from ..simulate import const
from ..visualize.stringformat import uc
from ..dataio.hdf5 import h5save
import numpy as np


def Pelegant(filename, macro=None, verbose=False, defns=None, Ncores=4):
    if not defns:
        defns = Popen('echo $HOME', shell=True, stdout=PIPE).stdout.read().rstrip()
        defns += '/defns.rpn'
    processstring = "export RPN_DEFNS='" + defns + "' && mpiexec.hydra -n " + str(Ncores) + " Pelegant " + filename
    if macro:
        processstring += ' -macro=' + macro
    process = Popen(processstring, shell=True, stdout=PIPE, stderr=STDOUT, universal_newlines=True)
    if verbose:
        print(processstring)
        stdout = []
        for line in process.stdout:
            stdout.append(line)
            print(line)
            sys.stdout.flush()
    else:
        stdout, stderr = process.communicate()
        if not str(stderr).strip() == 'None':
            print(stderr)
        for i, line in enumerate(stdout.split('\n')):
            if 'error: ' in line:
                print('\n'.join(stdout.split('\n')[i:]))
    return stdout


def elegant(filename, macro=None, verbose=False, defns=None):
    if not defns:
        defns = Popen('echo $HOME', shell=True, stdout=PIPE).stdout.read().rstrip()
        defns += '/defns.rpn'
    processstring = "export RPN_DEFNS='" + defns + "' && elegant " + filename
    if macro:
        processstring += ' -macro=' + macro
    process = Popen(processstring, shell=True, stdout=PIPE, stderr=STDOUT, universal_newlines=True)
    if verbose:
        print(processstring)
        stdout = []
        for line in process.stdout:
            stdout.append(line)
            print(line)
            sys.stdout.flush()
    else:
        stdout, stderr = process.communicate()
        if not str(stderr).strip() == 'None':
            print(stderr)
        for i, line in enumerate(stdout.split('\n')):
            if 'error: ' in line:
                print('\n'.join(stdout.split('\n')[i:]))
    return stdout


def loadwpcols(filename):
    ''' desired output: shape(dict['x' or 'Cx']) = (turns, particles)
    with nan when particle gets lost in coo mode
    Problems:
        - when particles are lost the sddsfile does not carry the NANs around
            -> every turn the particleIDs need to be checked
        - shape of colvals opposite for tracking of
          coordinates (turns, particles) and other (1, turns)
    '''
    data = sdds.SDDS(0)
    data.load(filename)
    colkeys = data.columnName
    colvals = data.columnData
    if filename[-3:] == 'coo':
        IDsindex = colkeys.index('particleID')
        colkeys.remove('particleID')
        IDs = colvals.pop(IDsindex)
        Nturns, Nparts = len(IDs), len(IDs[0])
        tmp = empty([Nturns, Nparts])
        tmp.fill(nan)
        newcolvals = []
        for i, val in enumerate(colvals):  # loop over x,xp,y,yp,...
            newcolvals.append(tmp.copy())
            for t in range(Nturns):
                newcolvals[i][t, array(IDs[t]) - 1] = array(val[t])
    else: # filename[-3:] == 'cen' or '.twiss'
        newcolvals = [array(val).T for val in colvals]
    return data, dict(zip(colkeys, newcolvals))


def sddsload(filename, verbose=False):
    data, coldict = loadwpcols(filename)
    description = data.description

    pars = data.parameterName
    Npars = len(pars)
    parvals = data.parameterData
    pardict = dict(zip(pars, parvals))

    Ncols = len(coldict)
#    cols = data.columnName
#    Ncols = len(cols)
#    colvals = data.columnData
#    # shape of colvals opposite for tracking of coordinates (points, particles) and other (1, points)  :(
#    if shape(colvals[0])[0]==1:
#        colvals = [array(val).T for val in colvals]
#    else:
#        colvals = [array(val) for val in colvals]
#    coldict = dict(zip(cols, colvals))

    if verbose:
        # print information on loaded data
        print("SDDS description: ", )
        print()
        print(Npars, "Paramters: {0:>13}".format('shape'))
        [print('{0:<20} {1:<}'.format(key, shape(val))) for key, val in pardict.items()]
        print()
        print(Ncols, "Columns: {0:>15}".format('shape'))
        [print('{0:<20} {1:<}'.format(key, shape(val))) for key, val in coldict.items()]
        print()

    data = pardict
    data[b'description'] = description
    data.update(coldict)
    return data


def sdds2hdf5(filename, verbose=False):
    data = sdds.SDDS(0)
    data.load(filename)

    pars = data.parameterName
    parvals = data.parameterData
    datadict = dict(zip(pars, parvals))

    cols = data.columnName
    colvals = data.columnData
    # shape of colvals opposite for tracking of coordinates (points, particles) and other (1, points)  :(
    if shape(colvals[0])[0]==1:
        colvals = [array(val).T for val in colvals]
    else:
        colvals = [array(val) for val in colvals]
    datadict.update(dict(zip(cols, colvals)))
    h5save(filename, verbose, **datadict)


def eleplot(datadict, x, y, *args, **kwargs):
    try:
        selected = kwargs['sel']
        del(kwargs['sel'])
    except:
        selected = [0]
    try:
        label = kwargs['label']
        del(kwargs['label'])
    except:
        label = True
    [plot(datadict[x][:, sel], datadict[y][:, sel], *args, **kwargs) for sel in selected]
    if label:
        xlabel(x)
        ylabel(y)


def trackplot3(datadict, abscissa='Pass'):
    subplot(331)
    eleplot(datadict, abscissa, 'Cx', '.')
    subplot(334)
    eleplot(datadict, abscissa, 'Cxp', '.')
    subplot(337)
    eleplot(datadict, 'Cx', 'Cxp', '.')
    subplot(332)
    eleplot(datadict, abscissa, 'Cy', '.')
    subplot(335)
    eleplot(datadict, abscissa, 'Cyp', '.')
    subplot(338)
    eleplot(datadict, 'Cy', 'Cyp', '.')
    subplot(333)
    eleplot(datadict, abscissa, 'dCt', '.')
    subplot(336)
    eleplot(datadict, abscissa, 'Cdelta', '.')
    subplot(339)
    eleplot(datadict, 'dCt', 'Cdelta', '.')


def trackplot2(datadict, particles=8, abscissa='t'):
    subplot(321)
    eleplot(datadict, abscissa, 'x', '.b', sel=range(particles))
    subplot(323)
    eleplot(datadict, abscissa, 'xp', '.b', sel=range(particles))
    subplot(325)
    eleplot(datadict, 'x', 'xp', '.', sel=range(particles))
    subplot(322)
    eleplot(datadict, abscissa, 'y', '.', sel=range(particles))
    subplot(324)
    eleplot(datadict, abscissa, 'yp', '.', sel=range(particles))
    subplot(326)
    eleplot(datadict, 'y', 'yp', '.', sel=range(particles))


def drawlatt(ax, data, lattice, size=0.1):
    
    # get dict of nagnets and strengths
    magdict = {}
    for line in lattice.split('\n'):
        name = line.split(':')[0].replace('"', '')
        if 'kquad' in line.lower():
            try:
                K1 = float(line.split('K1=')[1].split(',')[0])
            except:
                K1 = 0.0
            magdict[name] = K1
        elif 'ksext' in line.lower():
            try:
                K2 = float(line.split('K2=')[1].split(',')[0])
            except:
                K2 = 0.0
            magdict[name] = K2

    s = data['s']
    ax.set_xlim(npmin(s), npmax(s))

    ElementTypes = data['ElementType'].astype('U13')
    ElementNames = data['ElementName'].astype('U13')

    yi, yf = ax.get_ylim()[0], ax.get_ylim()[1]
    dy = abs(yf - yi)*size
    ax.set_ylim(yi, yf + dy)
    dy = dy*.8
    yl = ax.get_ylim()[1] - dy/2  # lower edge of element

    i = 0
    while i < len(s):
        et = ElementTypes[i]

        try:
            li = where(ElementTypes[i:] != et)[0][0]
        except:
            li = len(ElementTypes[i:])

        if i == 0:
            si = s[i]
        else:
            si = s[i - 1]
        sf = s[i + li - 1]
        dx = sf - si
        i += li

        if et == 'CSBEND':
            mypatch(ax, 'yellow', si, yl, dx, dy)
        elif et == 'KQUAD':
            if magdict[ElementNames[i - li, 0]] > 0:
                mypatch(ax, 'red', si, yl, dx, dy, typ='focus')
            else:
                mypatch(ax, 'red', si, yl, dx, dy, typ='defocus')
        elif et == 'KSEXT':
            if magdict[ElementNames[i - li, 0]] > 0:
                mypatch(ax, 'green', si, yl, dx, dy, typ='sextfocus')
            else:
                mypatch(ax, 'green', si, yl, dx, dy, typ='sextdefocus')
    return


def mypatch(ax, col, si, yl, dx, dy, typ='rectangle'):
    opts = {'fc': col,
            'ec': None,
            'clip_on': False,
            'zorder': 111}
    if typ == 'rectangle':
        ax.add_patch(Rectangle(xy=(si, yl), width=dx, height=dy, **opts))
    elif typ == 'focus':
        xy = array([[si, yl + dy/2],
                    [si + dx/2, yl + dy],
                    [si + dx, yl + dy/2],
                    [si + dx/2, yl]])
        ax.add_patch(Polygon(xy, **opts))
    elif typ == 'defocus':
        xy = array([[si, yl + dy],
                    [si + dx, yl + dy],
                    [si, yl],
                    [si + dx, yl]])
        ax.add_patch(Polygon(xy, **opts))
    elif typ == 'sextfocus':
        xy = array([[si, yl],
                    [si + dx/4, yl],
                    [si + dx/2, yl + dy*0.45],
                    [si + 3*dx/4, yl],
                    [si + dx, yl],
                    [si, yl]])
        path = Path(xy, [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO])
        ax.add_patch(PathPatch(path, **opts))
        path = Path(xy + array([0, dy/2 + dy/20]), [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO])
        ax.add_patch(PathPatch(path, **opts))
#        xy = array([[si, yl + dy*0.2],
#                    [si + dx/2, yl + dy*0.4],
#                    [si + dx, yl + dy*0.2],
#                    [si + dx/2, yl]])
#        ax.add_patch(Polygon(xy + array([0, dy/2 + dy/10]), **opts))
#        xy = array([[si, yl + dy*0.4],
#                    [si + dx, yl + dy*0.4],
#                    [si, yl],
#                    [si + dx, yl]])
#        ax.add_patch(Polygon(xy, **opts))
    elif typ == 'sextdefocus':
        xy = array([[si, yl + dy*0.45],
                    [si + dx/4, yl + dy*0.45],
                    [si + dx/2, yl],
                    [si + 3*dx/4, yl + dy*0.45],
                    [si + dx, yl + dy*0.45],
                    [si, yl + dy*0.45],])
        path = Path(xy, [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO])
        ax.add_patch(PathPatch(path, **opts))
        path = Path(xy + array([0, dy/2 + dy/20]), [Path.MOVETO, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.CURVE3, Path.LINETO])
        ax.add_patch(PathPatch(path, **opts))
#        xy = array([[si, yl + dy*0.4],
#                    [si + dx, yl + dy*0.4],
#                    [si, yl],
#                    [si + dx, yl]])
#        ax.add_patch(Polygon(xy + array([0, dy/2 + dy/10]), **opts))
#        xy = array([[si, yl + dy*0.2],
#                    [si + dx/2, yl + dy*0.4],
#                    [si + dx, yl + dy*0.2],
#                    [si + dx/2, yl]])
#        ax.add_patch(Polygon(xy, **opts))
    return



def multicolorylab(ax, lablist, collist):
    boxs = [TextArea(l, textprops=dict(color=c, rotation=90, ha='right', va='bottom')) for l, c in zip(lablist, collist)]
    ybox = VPacker(children=boxs[::-1], pad=0, sep=5)
    anchored_ybox = AnchoredOffsetbox(loc=6, child=ybox, pad=-2.4, frameon=False,
                                      bbox_to_anchor=(.0, .5),
                                      bbox_transform=ax.transAxes, borderpad=.4)
    ax.add_artist(anchored_ybox)


def drawtwiss(data, lattice, ax1):
    ax2 = ax1.twinx()
    ax1.plot(data['s'], data['betax'], '-g')
    ax1.plot(data['s'], data['betay'], '-b')
    ax1.set_xlabel('s')
    multicolorylab(ax1, ['$\\beta_x$', ', ','$\\beta_y$', ' / (m)'], ['g', None, 'b', None])
    ax2.plot(data['s'], data['etax'], '-r')
    ax2.set_ylabel(r'$\eta_x$ / (m)', color='r')
    drawlatt(ax1, data, lattice)
    ax2.grid(None)
    return


def twissdata(data):
    print('General:')
    print('    L = {:}'.format(npmax(data['s'])))
    print('    E = {:}'.format(data['pCentral'][0]))
    print('    Eloss = {:}'.format(data['U0'][0]))
    print('\nTunes:')
    print('    Qx = {:.6}'.format(data['nux'][0]))
    print('    Qy = {:.6}'.format(data['nuy'][0]))
    print('\nChromaticities:')
    print('    ' + uc.greek.xi + 'x = {:.6}'.format(data['dnux|dp'][0]))
    print('    ' + uc.greek.xi + 'y = {:.6}'.format(data['dnuy|dp'][0]))
    print('\nMomentum Compaction Factor:')
    print('    ' + uc.greek.alpha + 'p = {:.5e}'.format(data['alphac'][0]))
    print('\nRadiation Damping times:')
    print('    ' + uc.greek.tau + 'x = {:.6} ms'.format(1e3*data['taux'][0]))
    print('    ' + uc.greek.tau + 'y = {:.6} ms'.format(1e3*data['tauy'][0]))
    print('    ' + uc.greek.tau + uc.greek.delta + ' = {:.6} ms'.format(1e3*data['taudelta'][0]))
    print('\nDamping partition factors:')
    print('    Jx = {:.6}'.format(data['Jx'][0]))
    print('    Jy = {:.6}'.format(data['Jy'][0]))
    print('    J' + uc.greek.delta + ' = {:.6}'.format(data['Jdelta'][0]))
    print('\nHorizontal equilibrium geometric and normalized emittances:')
    print('    ' + uc.greek.epsilon + 'x = {:.6}'.format(data['ex0'][0]))
    print('    ' + uc.greek.epsilon + 'x* = {:.6}'.format(data['enx0'][0]))
    print('\nEquilibrium fractional rms energy spread:')
    print('    ' + uc.greek.epsilon + uc.greek.delta + ' = {:.6e}'.format(data['Sdelta0'][0]))
    print('\nHigher Order::')
    print('    ' + uc.greek.xi + 'x2 = {:}'.format(data['dnux|dp2'][0]))
    print('    ' + uc.greek.xi + 'y2 = {:}'.format(data['dnuy|dp2'][0]))
    print('    ' + uc.greek.xi + 'x3 = {:}'.format(data['dnux|dp3'][0]))
    print('    ' + uc.greek.xi + 'y3 = {:}'.format(data['dnuy|dp3'][0]))
    print('    ' + uc.greek.alpha + 'p2 = {:e}'.format(data['alphac2'][0]))
    return


def tunechrom(line, data):
    # props=dict(boxstyle='round', alpha=0.5)
    # line = ax.text(9/16*0.02, 0.98, '', bbox=bboxprops, **lineprops)
    string = 'Qx = {:.4f}'.format(data['nux'][0])
    string += '\nQy = {:.4f}'.format(data['nuy'][0])
    string += '\n' + uc.greek.xi + 'x = {:.4f}'.format(data['dnux|dp'][0])
    string += '\n' + uc.greek.xi + 'y = {:.4f}'.format(data['dnuy|dp'][0])
    line.set_text(string)
    return


def twissplot(ax, data, lattice, zoom=False):
    if zoom:
        starti, endi = argmax(data['s'] >= zoom[0]), argmax(data['s'] >= zoom[1])
        clipdata = {}
        for clip in ['s', 'betax', 'betay', 'etax', 'ElementType', 'ElementName']:
            clipdata[clip] = data[clip][starti:endi]
        clipdata['description'] = data['description']
        drawtwiss(clipdata, lattice, ax)
        return
    drawtwiss(data, lattice, ax)
    return

#Dnames = ['Injection','U125','UE56','U49','UE52','UE56 + U139 (slicing)','UE112','UE49']
#Tnames = ['Landau + BAM WLS7','MPW','U41','UE49','UE46','CPMU17 + UE48 (EMIL)','PSF WLS7','Cavities']
#names = {'D': Dnames, 'T' : Tnames, 'S' :  core.defchararray.add(Dnames,core.defchararray.add(' + ',Tnames))}


def biizoom(sectyp, num):
    if sectyp == 'achromat':
        return list(array([0, 240/16]) + (num - 1)*240/16)
    elif sectyp == 'duplet':
        return list(array([3*240/32, 5*240/32]) + (num - 2)*240/8)
    elif sectyp == 'triplet':
        return list(array([240/32, 3*240/32]) + (num - 1)*240/8)


def trackplot(ax, datadict, xy=False, exclude=False, only=False, fs=[16, 9], showlost=False,
              turns=None, everyxturn=[0, 1], ms=1, color=False, lab=''):
    if xy:
        x, y = xy
        
        if only:
            IDs = np.array(only)
        else:
            IDs = datadict['allIDs'].copy()
            if exclude:
                IDs = np.delete(IDs, exclude)
            if not showlost:
                IDs = np.delete(IDs, datadict['lostIDs'])
        
        if not color:
            colors = cm.rainbow(linspace(0, 1, datadict['Particles'][0]))
        else:
            colors = np.repeat(color, len(IDs))
        
        for part, col in zip(IDs, colors):
            i, f = everyxturn
            xdat, ydat = datadict[x][:turns, part][i::f], datadict[y][:turns, part][i::f]
            ax.plot(xdat*1e3, ydat*1e3, '.', color=col, ms=ms)
        ax.plot([], [], '.', color=col, ms=ms, label=lab)
        return
    try:  # centroid watch point
        x = ['Pass', 'Pass', 'Pass', 'Pass', 'Pass', 'Pass', 'Cx', 'Cy', 'dCt']
        y = ['Cx', 'Cy', 'dCt', 'Cxp', 'Cyp', 'Cdelta', 'Cxp', 'Cyp', 'Cdelta']
        for i, (x, y) in enumerate(zip(x, y)):
            subplot(3, 3, i+1)
            if turns:
                plot(datadict[x][:turns, ], datadict[y][:turns, ], '.')
            else:
                plot(datadict[x][:, ], datadict[y][:, ], '.')
            xlabel(x)
            ylabel(y)
    except:  # coordinate watch point
        x = ['t', 't', 't', 't', 't', 't', 'x', 'y', 'dt']
        y = ['x', 'y', 'dt', 'xp', 'yp', 'p', 'xp', 'yp', 'p']
        colors = cm.rainbow(linspace(0, 1, datadict['Particles'][0]))
        for i, (x, y) in enumerate(zip(x, y)):
            subplot(3, 3, i+1)
            if turns:
                if x == 't':
                    [plot(datadict[y][:turns, part], '.', color=col) for part, col in enumerate(colors)]
                    xlabel('Pass')
                else:
                    [plot(datadict[x][:turns, part], datadict[y][:turns, part], '.', color=col) for part, col in enumerate(colors)]
                    xlabel(x)
            else:
                if x == 't':
                    [plot(datadict[y][:, part], '.', color=col) for part, col in enumerate(colors)]
                    xlabel('Pass')
                else:
                    [plot(datadict[x][:, part], datadict[y][:, part], '.', color=col) for part, col in enumerate(colors)]
                    xlabel(x)
            ylabel(y)
            tight_layout()
    return


def showbun(datadict):
    x = ['', '', '', '', '', '', 'x', 'y', 't']
    y = ['x', 'y', 't', 'xp', 'yp', 'p', 'xp', 'yp', 'p']
    colors = cm.rainbow(linspace(0, 1, datadict['Particles'][0]))
    for i, (x, y) in enumerate(zip(x, y)):
        subplot(3, 3, i+1)
        if x == '':
            [plot(datadict[y][part, ], '.', color=col) for part, col in enumerate(colors)]
            xlabel('Pass')
        else:
            [plot(datadict[x][part, ], datadict[y][part, ], '.', color=col) for part, col in enumerate(colors)]
            xlabel(x)
        ylabel(y)
    tight_layout()
    return


def mybunch(bunchname, ranges, E_mev):

    cl = const.cl
    qe = const.qe
    me = const.me
    E0 = me*cl**2/qe  # eV
    gamma = 1 + E_mev*1e6/E0
    #print('gamma = {}'.format(gamma))

    N = len(ranges['x'])
    bunch = sdds.SDDS(0)  # what does the index mean?
    bunch.setDescription('my predefined bunch', 'bunched-beam phase space')
    bunch.mode = 1  # 1 is binary, 2 is ascii

    parnames = ['Particles', 'pCentral', 'IDSlotsPerBunch']
    parsymbs = ['', 'p$bcen$n', '']
    parunits = ['', 'm$be$nc', '']
    pardescr = ['Number of particles before sampling', 'Reference beta*gamma', 'Number of particle ID slots reserved to a bunch']
    parforms = ['', '', '']
    parfixva = ['', '', '']
    partypes = [bunch.SDDS_LONG, bunch.SDDS_DOUBLE, bunch.SDDS_LONG]
    pardatas = [[N], [gamma], [N]]
    parstuff = zip(parnames, parsymbs, parunits, pardescr, parforms, partypes, parfixva, pardatas)

    for name, symbol, units, description, formatString, typ, fixedValue, data in parstuff:
        bunch.defineParameter(name, symbol, units, description, formatString, typ, fixedValue)
        bunch.setParameterValueList(name, data)

    colnames = ['x', 'xp', 'y', 'yp', 't', 'p', 'particleID']
    colsymbs = ['']*7
    colunits = ['m', '', 'm', '', 's', 'm$be$nc', '']
    coldescr = ['']*7
    colforms = ['']*7
    coltypes = [bunch.SDDS_DOUBLE]*6 + [bunch.SDDS_LONG]
    colfleng = [1]*7
    coldatas = [[list(ranges['x'])],
                [list(ranges['xp'])],
                [list(ranges['y'])],
                [list(ranges['yp'])],
                [list(ranges['t'])],
                [list(1 + ranges['p']*1e6/E0)],
                [range(1, 1 + N)]]
    colstuff = zip(colnames, colsymbs, colunits, coldescr, colforms, coltypes, colfleng, coldatas)

    for name, symbol, units, description, formatString, typ, fieldLength, data in colstuff:
        bunch.defineColumn(name, symbol, units, description, formatString, typ, fieldLength)
        bunch.setColumnValueLists(name, data)


    bunch.save(bunchname)
    return
