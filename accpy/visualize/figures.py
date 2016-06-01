# -*- coding: utf-8 -*-
''' accpy.visualize.figures
Notes:
    some useful functions for matplotlib plots
Author:
    Felix Kramer (felix.kramer@physik.hu-berlin.de)
'''
from matplotlib import rcdefaults, rcParams
from time import strftime

def plotstandards(formattype, scale, w=1920, h=1080):
    '''set standards for matplotlib for given formattype
        possible formattypes:
        presentation_1920x1080
        pcdisplay
        a4_portrait
        a4_landscape
        galaxytabs2(2048x1536)
        galaxytabs2(pclooking)
    '''
    def figsize(formattype, scale):
        if formattype == 'a4portrait':
            width = 485.47523/72     # \the\textwidth conversion to inches
            height = 9*width/16
        elif formattype == 'a4landscape':
            width = (297-20)/25.4
            height = (210-20)/25.4
        elif formattype == 'presentation_1920x1080':
            width = 24
            height = 13.5
        elif formattype == 'HZB_169':
            width = 24
            height = 16.4
        elif formattype == 'HZB_43':
            width = 24
            height = 16.4
        # rescaling and conversion to inches
        width = width*scale[0]
        height = height*scale[1]
        return [width, height]

    rcdefaults()
    presentations = ['presentation_1920x1080', 'HZB_169', 'HZB_43']
    if formattype in  presentations:
        params = {'figure.figsize': figsize(formattype, scale),
                  'figure.dpi': 80,
                  'figure.autolayout': True,
                  'font.size': 28,
                  'lines.linewidth': 2,
                  'savefig.dpi': 100,
                  'savefig.format': 'svg',
                  'savefig.bbox': 'tight',
                  'savefig.pad_inches': 0,
                  'axes.grid': True,
                  'axes.formatter.limits': (-2, 3)}
    elif formattype == 'pcdisplay':
        params = {'figure.figsize': (w/1920*24, h/1080*13.5),
                  'figure.dpi': 80,
                  'figure.autolayout': True,
                  'font.size': 28,
                  'font.family': 'serif',
                  'legend.fontsize': 28,
                  'lines.markersize': 3,
                  'lines.linewidth': 2,
                  'savefig.dpi': 100,
                  'savefig.pad_inches': 0,
                  'savefig.format': 'svg',
                  'text.usetex': True,
                  'xtick.labelsize': 28,
                  'ytick.labelsize': 28}
    elif formattype == 'a4_portrait':
        params = {'axes.labelsize': 10,  # fontsize for x and y labels
                  'axes.titlesize': 0,
                  'axes.formatter.limits': (-2, 3),
                  'axes.grid': True,
                  'backend': 'ps',
                  'figure.figsize': figsize(formattype, scale),
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 3,
                  'savefig.dpi': 100,
                  'savefig.pad_inches': 0,
                  'savefig.format': 'svg',
                  'text.usetex': True,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10}
    elif formattype == 'a4landscape':
        params = {'axes.labelsize': 10,  # fontsize for x and y labels
                  'axes.titlesize': 10,
                  'axes.formatter.limits': (-2, 3),
                  'axes.grid': True,
                  'backend': 'ps',
                  'figure.figsize': figsize(formattype, scale),
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 3,
                  'savefig.dpi': 100,
                  'savefig.pad_inches': 0,
                  'savefig.format': 'svg',
                  'text.usetex': True,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10}
    elif formattype == 'galaxytabs2(2048x1536)':
        params = {'axes.labelsize': 10,  # fontsize for x and y labels
                  'axes.titlesize': 10,
                  'axes.formatter.limits': (-2, 3),
                  'axes.grid': True,
                  'backend': 'ps',
                  'figure.figsize': [4*1.94, 3*1.94],
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'font.size': 10,
                  'font.family': 'serif',
                  'legend.fontsize': 10,
                  'lines.markersize': 3,
                  'savefig.dpi': 600,
                  'savefig.pad_inches': 0,
                  'savefig.format': 'svg',
                  'text.usetex': True,
                  'xtick.labelsize': 10,
                  'ytick.labelsize': 10}
    elif formattype == 'galaxytabs2(pclooking)':
        params = {'axes.labelsize': 14,  # fontsize for x and y labels
                  'axes.titlesize': 14,
                  'axes.formatter.limits': (-2, 3),
                  'axes.grid': True,
                  'backend': 'ps',
                  'figure.figsize': [4*1.94, 3*1.94],
                  'figure.dpi': 100,
                  'figure.autolayout': True,
                  'font.size': 14,
                  'font.family': 'serif',
                  'legend.fontsize': 14,
                  'lines.markersize': 3,
                  'savefig.dpi': 600,
                  'savefig.pad_inches': 0,
                  'savefig.format': 'svg',
                  'text.usetex': True,
                  'xtick.labelsize': 14,
                  'ytick.labelsize': 14}
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
