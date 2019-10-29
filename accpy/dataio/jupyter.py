# -*- coding: utf-8 -*-
''' accpy.dataio.jupyter
author:
    Felix Armborst
'''
import ipywidgets as widgets
from .hdf5 import h5load
from subprocess import check_output


def loaddirs(paths, data, datapath):
    for i, path in enumerate(paths):
        print(i, path.lstrip(datapath))
        data[path] = {
            'data_coo' : [],
            'data_cen' : [],
            'data_twi' : [],
            'data_log' : [],
        }
        try:
            fcoos = check_output('ls -tr ' + path + '*.coo', shell=True).decode().split('\n')[:-1]
            for fcoo in fcoos:
                trackdat = h5load(fcoo)
                data[path]['data_coo'].append(trackdat)
            print('    Loaded Coordinate Tracking Data Files: {}'.format(len(fcoos)))
        except:
            data[path]['data_coo'].append([])
            print('    Loaded Coordinate Tracking Data Files: 0')

        try:
            fcens = check_output('ls -tr ' + path + '*.cen', shell=True).decode().split('\n')[:-1]
            for fcen in fcens:
                trackdat = h5load(fcen)
                data[path]['data_cen'].append(trackdat)
            print('    Loaded Centroid Tracking Data Files: {}'.format(len(fcens)))
        except:
            data[path]['data_cen'].append([])
            print('    Loaded Centroid Tracking Data Files: 0')

        try:
            ftwisss = check_output('ls -tr ' + path + '*.twiss', shell=True).decode().split('\n')[:-1]
            for ftwiss in ftwisss:
                twissdat = h5load(ftwiss)
                data[path]['data_twi'].append(twissdat)
            print('    Loaded Twiss Data Files: {}'.format(len(ftwisss)))
        except:
            data[path]['data_twi'].append([])
            print('    Loaded Twiss Data Files: 0')

        try:
            flog = check_output('ls -tr ' + path + '*.log', shell=True).decode().split('\n')[:-1][0]
            with open(flog, 'r') as flog_h:
                log = flog_h.read()
                data[path]['data_log'].append(log)
            print('    Added Log Data!')
        except:
            data[path]['data_log'].append([])
            print('    No Log Data found!')
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