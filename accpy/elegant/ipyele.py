# -*- coding: utf-8 -*-
"""accpy.elegant.ipyele
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import print_function, division
from subprocess import Popen, PIPE, STDOUT
from numpy import shape, array, max as npmax, argmax, roll, float64, core, min
from matplotlib.pylab import (plot, subplot, xlabel, ylabel, twinx, gca, xlim,
                              ylim, annotate)
import matplotlib.patches as patches
from matplotlib.offsetbox import AnchoredOffsetbox, TextArea, HPacker, VPacker
from IPython.display import display, Javascript
from . import sdds
from ..visualize.stringformat import uc
from ..dataio.hdf5 import h5save


def elegant(filename):
    path = Popen('echo $HOME', shell=True, stdout=PIPE).stdout.read().rstrip()
    path += '/defns.rpn'
    processstring = "export RPN_DEFNS='" + path + "' && elegant " + filename
    process = Popen(processstring, shell=True, stdout=PIPE, stderr=STDOUT)
    return process.stdout.read()


def Pelegant(filename, Ncores=2):
    path = Popen('echo $HOME', shell=True, stdout=PIPE).stdout.read().rstrip()
    path += '/defns.rpn'
    processstring = "export RPN_DEFNS='" + path + "' && mpiexec.hydra -n " + str(Ncores) + " Pelegant " + filename
    process = Popen(processstring, shell=True, stdout=PIPE, stderr=STDOUT)
    return process.stdout.read()


def sddsload(filename, verbose=False):
    data = sdds.SDDS(0)
    data.load(filename)

    pars = data.parameterName
    Npars = len(pars)
    parvals = data.parameterData
    pardict = dict(zip(pars, parvals))

    cols = data.columnName
    Ncols = len(cols)
    colvals = data.columnData
    # shape of colvals opposite for tracking of coordinates (points, particles) and other (1, points)  :(
    if shape(colvals[0])[0]==1:
        colvals = [array(val).T for val in colvals]
    else:
        colvals = [array(val) for val in colvals]
    coldict = dict(zip(cols, colvals))


    if verbose:
        # print information on loaded data
        print("SDDS description: ", data.description)
        print()
        print(Npars, "Paramters:")
        [print('{0:} = {1:}'.format(key, val[0])) for key, val in pardict.items()]
        print()
        print(Ncols, "Columns:")
        [print('{0:} {1:}'.format(key, shape(val))) for key, val in coldict.items()]
        print()

    data = pardict
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
    [plot(datadict[x][:, sel], datadict[y][:, sel], *args, **kwargs) for sel in selected]
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


def trackplot2(datadict, abscissa='t'):
    subplot(321)
    eleplot(datadict, abscissa, 'x', '.', sel=range(8))
    subplot(323)
    eleplot(datadict, abscissa, 'xp', '.')
    subplot(325)
    eleplot(datadict, 'x', 'xp', '.')
    subplot(322)
    eleplot(datadict, abscissa, 'y', '.', sel=range(8))
    subplot(324)
    eleplot(datadict, abscissa, 'yp', '.')
    subplot(326)
    eleplot(datadict, 'y', 'yp', '.')


def twissplot(data):
    L = npmax(data['s'])
    print('L = {:}'.format(L))
    print('Qx = {:}'.format(data['nux'][0]))
    print('Qy = {:}'.format(data['nuy'][0]))
    print(uc.greek.alpha + 'p = {:e}'.format(data['alphac'][0]))

    eleplot(data, 's', 'betax', '-g')
    eleplot(data, 's', 'betay', '-b')
    twinx()
    eleplot(data, 's', 'etax', '-r')

    # Latteice graphics vertical position and size (axis coordinates!)
    lypos = gca().get_ylim()[1]
    tp = Twissplot(lypos = lypos, lysize = lypos*0.12)
    tp.axislabels(yscale=0.5)
    tp.paintlattice(data, 0, L, ec=False, fscale=2)
    xlim(0, L)


def b2twissplot(coldict, pardict):
    L = npmax(coldict['s'])
    print('L = {:}'.format(L))
    print('Qx = {:}'.format(pardict['nux'][0]))
    print('Qx = {:}'.format(pardict['nux'][0]))
    print('Qy = {:}'.format(pardict['nuy'][0]))
    print(uc.greek.alpha + 'p = {:e}'.format(pardict['alphac'][0]))

    eleplot(coldict, 's', 'betax', '-g')
    eleplot(coldict, 's', 'betay', '-b')
    twinx()
    eleplot(coldict, 's', 'etax', '-r')

    # Latteice graphics vertical position and size (axis coordinates!)
    lypos = gca().get_ylim()[1]
    tp = Twissplot(lypos = lypos, lysize = lypos*0.12)
    tp.axislabels(yscale=0.5)
    tp.paintlattice(coldict, 0, L, ec=False, fscale=2)
    xlim(0, L)



def autoscroll(threshhold):
    if threshhold==0:  # alway scroll !not good
        javastring = '''
        IPython.OutputArea.prototype._should_scroll = function(lines) {
            return true;
        }
        '''
    elif threshhold==-1:  # never scroll !not good
        javastring = '''
        IPython.OutputArea.prototype._should_scroll = function(lines) {
            return false;
        }
        '''
    else:
        javastring = 'IPython.OutputArea.auto_scroll_threshold = ' + str(threshhold)
    display(Javascript(javastring))


class Twissplot():
    Dnames = ['Injection','U125','UE56','U49','UE52','UE56 + U139 (slicing)','UE112','UE49']
    Tnames = ['Landau + BAM WLS7','MPW','U41','UE49','UE46','CPMU17 + UE48 (EMIL)','PSF WLS7','Cavities']
    names = {'D': Dnames, 'T' : Tnames, 'S' :  core.defchararray.add(Dnames,core.defchararray.add(' + ',Tnames))}

    def __init__(self, lypos = 25, lysize = 3):
        self.lypos = lypos
        self.lysize = lysize

#    def getrolled(self,s,y,fmt=None): # access lattice from -120 to 120m
#        x = array(s,dtype=float64)
#        y = array(y,dtype=float64)
#        x[x > 120] = x[x > 120] - 240.0
#        ishift = argmax(x < 0)
#        x = roll(x,-ishift)
#        y = roll(y,-ishift)
#        if fmt:
#            return x,y,fmt
#        else:
#            return x,y

    # Note:
    # It seem elegant twiss-ouput always prints the length, type and name
    # at the END of the element (!)
    #
    # Check element list
    # print np.unique(twi.ElementType)
    # ['CSBEND' 'DRIF' 'KQUAD' 'KSEXT' 'MALIGN' 'MARK' 'RECIRC' 'WATCH']
    #
    def paintlattice(self,d,s0,s1,ycenter=None,ysize=None,ec=True,labels=True,rolled=False, fscale=1.0):
        if ycenter is None:
            ycenter = self.lypos
        if ysize is None:
            ysize = self.lysize
        s = d['s']
        et = d['ElementType']
        en = d['ElementName']
        if rolled:
            s[s> 120] = s[s > 120] - 240.0
            ishift = argmax(s < 0)
            s = roll(s,-ishift)
            et = roll(et,-ishift)
            en = roll(en,-ishift)

        i0 = argmax(s >= s0)
        i1 = argmax(s >= s1)
        #print i0,i1
        start = s0
        for i in range(i0,i1 + 1):
            # save start if previous element was something else
            if i > i0:
                if et[i] != et[i-1]:
                    start = s[i-1]
            # skip if next element of same type
            if i < i1:
                if et[i] == et[i+1]:
                    continue
            end = min((s[i],s1))
            l = end - start
            #print i, s[i], et[i], en[i], '    length of element:', l
            col = 'none'
            ecol='k'
            if et[i] == 'CSBEND':
                col = 'yellow'
                ecol='black'
            if et[i] == 'KQUAD':
                col = 'red'
            if et[i] == 'KSEXT':
                col = 'green'
            if not ec:
                ecol='none'

            if col != 'none':
                gca().add_patch(patches.Rectangle((start, ycenter-0.5*ysize), l,ysize, ec=ecol, facecolor=col,clip_on=False, zorder = 101))
                if labels:
                    fs = 80 / (s1-s0) * fscale
                    if et[i] == 'KSEXT':
                        annotate(en[i], xy=(start,ycenter), xytext=(start+0.5*l, ycenter - .55*ysize),fontsize=fs,va='top',ha='center',clip_on=False, zorder = 102)
                    else:
                        annotate(en[i], xy=(start,ycenter), xytext=(start+0.5*l, ycenter + .5*ysize),fontsize=fs,va='bottom',ha='center',clip_on=False, zorder = 102)


    def axislabels(self,yscale=1,Dfac=10):
        xlabel('s / m')

        ybox3 = TextArea("       $\\eta_x / {0}".format(int(100/Dfac))+"\mathrm{cm}$", textprops=dict(color="r",rotation=90,ha='left',va='center'))
        ybox1 = TextArea("  $\\beta_y / \mathrm{m}$",     textprops=dict(color="b",rotation=90,ha='left',va='center'))
        ybox2 = TextArea("$\\beta_x / \mathrm{m}$", textprops=dict(color="g",rotation=90,ha='left',va='center'))
        ybox = VPacker(children=[ybox3, ybox1, ybox2],align="bottom", pad=0, sep=5)
        anchored_ybox = AnchoredOffsetbox(loc=8, child=ybox, pad=0., frameon=False, bbox_to_anchor=(-0.08*yscale, 0.15),
                                          bbox_transform=gca().transAxes, borderpad=0.)
        gca().add_artist(anchored_ybox)
        ylim(-1)


    def plotsection(self,d,stype,nr):
        s0 = (nr-1)*30.0 - 7.5
        if stype == 'T':
            s0 += 15.0
        if stype == 'S':
            s1 = s0 + 30.0
        else:
            s1 = s0 + 15.0

        if stype == 'D' and nr == 1:
            plot(*self.getrolled(d.s, d.betax,'g-'))
            plot(*self.getrolled(d.s, d.betay,'b-'))
            x, y, = self.getrolled(d.s, d.etax)
            plot(x,10*y,'r-')
            rolled=True
        else:
            plot(d.s, d.betax,'g-')
            plot(d.s, d.betay,'b-')
            plot(d.s, 10* array(d.etax,dtype=float64) ,'r-')
            rolled=False

        annotate(stype+'{0:0n}'.format(nr),xy=((s1+s0)/2.0,25-5), fontsize=20,ha='center',va='top',zorder=105)
        gca().yaxis.grid(alpha=0.3, zorder=0)

        #    print stype, names[stype][nr-1], names[stype]
        annotate(''+self.names[stype][nr-1]+'',xy=((s1+s0)/2.0,25-9), fontsize=8,ha='center',va='top',zorder=105)

        self.paintlattice(d,s0,s1, self.lypos, self.lysize, ec=True,rolled=rolled)
        self.axislabels()
        xlim(s0,s1)