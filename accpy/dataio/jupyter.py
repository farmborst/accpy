# -*- coding: utf-8 -*-
''' accpy.dataio.jupyter
author:
    Felix Armborst
'''
import ipywidgets as widgets
from subprocess import check_output
from .hdf5 import h5load
from ..elegant.ipyele import sddsload


def loaddirs(paths, data, datapath):
    for i, path in enumerate(paths):
        print(i, path.lstrip(datapath))
        fcoos, fcens, ftwisss, fbuns, flogs = [[] for i in range(5)]
        data[path] = {
            'data_coo' : [],
            'data_cen' : [],
            'data_twi' : [],
            'data_log' : [],
            'data_bun' : [],
        }
        print('  Loaded Data Files:')
        try:
            fcoos = check_output('ls -tr ' + path + '*.coo', shell=True).decode().split('\n')[:-1]
            for fcoo in fcoos:
                trackdat = h5load(fcoo)
                data[path]['data_coo'].append(trackdat)
            print('    Coordinate Tracking: {}'.format(len(fcoos)))
        except:
            data[path]['data_coo'].append([])
            print('    Coordinate Tracking: 0')

        try:
            fcens = check_output('ls -tr ' + path + '*.cen', shell=True).decode().split('\n')[:-1]
            for fcen in fcens:
                trackdat = h5load(fcen)
                data[path]['data_cen'].append(trackdat)
            print('    Centroid Tracking: {}'.format(len(fcens)))
        except:
            data[path]['data_cen'].append([])
            print('    Centroid Tracking: 0')

        try:
            ftwisss = check_output('ls -tr ' + path + '*.twiss', shell=True).decode().split('\n')[:-1]
            for ftwiss in ftwisss:
                twissdat = h5load(ftwiss)
                data[path]['data_twi'].append(twissdat)
            print('    Twiss: {}'.format(len(ftwisss)))
        except:
            data[path]['data_twi'].append([])
            print('    Twiss: 0')
        
        try:
            fbuns = check_output('ls -tr ' + path + '*.bun', shell=True).decode().split('\n')[:-1]
            for fbun in fbuns:
                bun = sddsload(fbun)
                data[path]['data_bun'].append(bun)
            print('    Bunch: {}'.format(len(fbuns)))
        except:
            data[path]['data_bun'].append([])
            print('    Bunch: 0')

        try:
            flogs = check_output('ls -tr ' + path + '*.log.*', shell=True).decode().split('\n')[:-1]
            for flog in flogs:
                if flog[-4:] == 'hdf5':
                    log = h5load(flog)
                else:
                    with open(flog, 'r') as flog_h:
                        log = flog_h.read()
                data[path]['data_log'].append(log)
            print('    Log: {}'.format(len(flogs)))
        except:
            data[path]['data_log'].append([])
            print('    Log: 0')
    return


def dataloadmenu(datapath):
    dirs = check_output('ls -d ' + datapath + '*/', shell=True).decode().split('\n')[:-1][:-1]
    dirs_short = ['{:03}:  '.format(i) + d.lstrip(datapath) for i, d in enumerate(dirs)]
    
    
    data = {}
    paths = []
    

    print('Data in {}'.format(datapath))
    selector = widgets.SelectMultiple(
        options = dirs_short,
        value=[],
        rows=25,
        layout=widgets.Layout(width='1000px'),
        description='',
        disabled=False
    )
    
    def ClickFun(_, data=data, paths=paths, datapath=datapath):
        data.clear()
        paths.clear()
        # "linking function with output"
        with out:
            out.clear_output()
            paths += [datapath + p.split(':  ')[-1] for p in selector.value]
            loaddirs(paths, data, datapath)

    button = widgets.Button(description="Load Data")
    out = widgets.Output()
    button.on_click(ClickFun)
    return widgets.VBox([selector, button, out]), data, paths