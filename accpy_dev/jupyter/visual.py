# -*- coding: utf-8 -*-
"""accpy.elegant.jupyter
author:     felix.kramer(at)physik.hu-berlin.de
"""
from __future__ import print_function, division
from IPython.display import display, Javascript
from matplotlib import rcParams, rcdefaults
from matplotlib.pyplot import close
from jupyterthemes import jtplot
from time import strftime


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


def figsave(name, lambdaplotstandard, lambdafig):
    rcdefaults()
    jtplot.style('grade3')
    lambdaplotstandard()
    fig = lambdafig()
    name += '_' + strftime('%Y%m%d%H%M%S') + '.pdf'
    fig.savefig(name, format='pdf')
    print('saved fig as ' + name)
    close(fig)
    rcdefaults()
    jtplot.style('monokai')
    return name


def plotstandard(purpose, scalex=1, scaley=1, defaults=True):
    if defaults:
        rcdefaults()
    if purpose == 'beamer169':
        params = {'axes.labelsize': 10,
                  'axes.titlesize': 10,
                  'axes.formatter.limits': [-2, 3],
                  'axes.grid': True,
                  'figure.figsize': [4*scalex, 2.25*scaley],
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'figure.frameon': False,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 4,
                  'lines.linewidth': 1,
                  'savefig.dpi': 100,
                  'savefig.facecolor': 'white',
                  'savefig.edgecolor': 'white',
                  'savefig.format': 'pdf',
                  'savefig.bbox': 'tight',
                  'savefig.pad_inches': 0.05,
                  'text.usetex': True,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10}
    elif purpose == 'jacow':
        params = {'axes.labelsize': 8,
                  'axes.titlesize': 8,
                  'axes.formatter.limits': [-2, 3],
                  'axes.grid': True,
                  'figure.figsize': [3.5*scalex, 2.7*scaley],
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'figure.frameon': False,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 4,
                  'lines.linewidth': 1,
                  'savefig.dpi': 600,
                  'savefig.facecolor': 'white',
                  'savefig.edgecolor': 'white',
                  'savefig.format': 'pdf',
                  'savefig.bbox': 'tight',
                  'savefig.pad_inches': 0.01,
                  'text.usetex': True,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8}
    elif purpose == 'A0poster':
        params = {'axes.labelsize': 8,
                  'axes.titlesize': 8,
                  'axes.formatter.limits': [-2, 3],
                  'axes.grid': True,
                  'figure.figsize': [2.35*scalex, 2.35*scaley],
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'figure.frameon': False,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 4,
                  'lines.linewidth': 1,
                  'savefig.dpi': 600,
                  'savefig.facecolor': 'white',
                  'savefig.edgecolor': 'white',
                  'savefig.format': 'pdf',
                  'savefig.bbox': 'tight',
                  'savefig.pad_inches': 0.01,
                  'text.usetex': True,
                  'xtick.labelsize': 8,
                  'ytick.labelsize': 8}
    rcParams.update(params)
    return
    
def rctablet():
    jtplot.style('monokai')
    fs = [2048/256, 1536/256]
    rcParams['savefig.format'] = 'png'
    rcParams['savefig.dpi'] = 256
    rcParams['savefig.facecolor'] = '#232323'
    rcParams['savefig.edgecolor'] = '#232323'
    rcParams['savefig.transparent'] = False
    rcParams['nbagg.transparent'] = False
    rcParams['savefig.pad_inches'] = 0.0
    return fs

def rcpof():
    rcdefaults()
    dpi = 256
    fs = [7.12, 4.04]
    params = {'savefig.dpi'          : dpi,
              'savefig.format'       : 'png',
              'axes.formatter.limits': [-3, 3]}
    rcParams.update(params)
    return fs