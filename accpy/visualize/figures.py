# -*- coding: utf-8 -*-
''' accpy.visualize.figures
Notes:
    some useful functions for matplotlib plots
Author:
    Felix Kramer (felix.kramer@physik.hu-berlin.de)
'''
from __future__ import division
from matplotlib import rcdefaults, rcParams
from time import strftime


def plotstandards(confdict):
    rcdefaults()
    '''
    a4: 210 x 297
    from \the\textwidth conversion to inches
    width = 485.47523/72 = 6.742711527777778
    height = 9*width/16  = 3.792775234375
    '''

    width = confdict['fig_width']
    height = confdict['fig_height']
    gridon = confdict['grid']
    dpi = confdict['dpi']
    fontfamily = confdict['fontfamily']
    fontsize = confdict['fontsize']
    markersize = confdict['markersize']
    linewidth = confdict['linewidth']
    axformatterlimits = confdict['axformatterlimits']

    params = {'axes.labelsize': fontsize,
              'axes.titlesize': fontsize,
              'axes.formatter.limits': axformatterlimits,
              'axes.grid': gridon,
              'figure.figsize': [width, height],
              'figure.dpi': dpi,
              'figure.autolayout': True,
              'font.size': fontsize,
              'font.family': fontfamily,
              'legend.fontsize': fontsize,
              'lines.markersize': markersize,
              'lines.linewidth': linewidth,
              'savefig.dpi': dpi,
              'savefig.facecolor': 'white',
              'savefig.edgecolor': 'white',
              'savefig.format': 'svg',
              'savefig.bbox': 'tight',
              'savefig.pad_inches': 0.05,
              'text.usetex': True,
              'xtick.labelsize': fontsize,
              'ytick.labelsize': fontsize,
              'errorbar.capsize': 3}
    rcParams.update(params)
    return


def figsave(fig, filename, filetype='pdf'):
    ''' figsave(filename, fig)
    input:
        - fig handle to save
        - desired filename as string
    return:
        - none
    notice:
    '''
    timstamp = strftime('%Y%m%d%H%M%S')
    filename = ''.join([timstamp, '_', filename, '.', filetype])
    fig.savefig(filename,
                frameon=False,
                )
    print ('\img{'+filename+'}')
    return
