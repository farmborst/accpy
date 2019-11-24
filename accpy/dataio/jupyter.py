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
        print('  Loaded Data Files:')
        data[path] = {}
        for datatype in ['coo', 'cen', 'twi', 'log', 'bun', 'lte']:
            data[path][datatype] = []
            
            # CHECK FOR CORRESPONDING FILES
            try:
                files = check_output('ls -tr ' + path + '*.' + datatype + '*', shell=True).decode().split('\n')[:-1]
            except:
                print('    *.{}: 0'.format(datatype))
                continue
                
            # REMOVE "DUPLICATE" HDF5 FILES (coo, cen, twi, log)
            hdf5 = [f[:-5] for f in files if f[-5:] == '.hdf5']
            files = [f for f in files if f not in hdf5]
            
            # LOAD FILES
            for f in files:
                ftype = f.split('.')[-1]
                if ftype in ['hdf5']:
                    data[path][datatype].append(h5load(f))
                elif ftype in ['bun']:
                    data[path][datatype].append(sddsload(f))
                elif ftype in ['log', 'lte']:
                    with open(f, 'r') as f_h:
                        data[path][datatype].append(f_h.read())

            print('    *.{}: {}'.format(datatype, len(files)))

    return


def dataloadmenu(datapath, loadtoram=True):
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
        paths[:] = []  # paths.clear() only for python 3.3+
        # "linking function with output"
        with out:
            out.clear_output()
            paths += [datapath + p.split(':  ')[-1] for p in selector.value]
            if loadtoram:
                loaddirs(paths, data, datapath)

    button = widgets.Button(description="Load Data")
    out = widgets.Output()
    button.on_click(ClickFun)
    return widgets.VBox([selector, button, out]), data, paths