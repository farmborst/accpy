# -*- coding: utf-8 -*-
''' accpy.dataio.jupyter
author:
    Felix Armborst
'''
import os
import ipywidgets as widgets
from subprocess import check_output
from .hdf5 import h5load
from ..elegant.ipyele import sddsload


def loaddirs(paths, data, datapath, progress=False):
    Npaths = len(paths)
    if progress:
            progress.max = Npaths
    for i, path in enumerate(paths):
        print('{:02d}: {:} '.format(i, path.split(datapath)[1]))
        print('  Loaded Data Files:')
        data[path] = {}

        # Get list of all files
        files = [f for f in os.listdir(path) if os.path.isfile(path + f)]

        # Find hdf5 files
        files_hdf5 = [f for f in files if f.split('.')[-1] in ['hdf5']]

        # Remove duplicate non-hdf5 (thus presumably sdds) files from list of files
        duptst = ['.'.join(f.split('.')[:-1]) for f in files_hdf5]
        files = [f for f in files if f not in duptst]

        for datatype in [['coo'], ['cen'], ['twi', 'twiss'], ['log'], ['bun'], ['lte'], ['rdt', 'srdts'], ['FPs']]:
            
            # CHECK FOR CORRESPONDING FILES
            matching = [path + f for f in files if f.split('.')[-1] in datatype]
            matching += [path + f for f in files_hdf5 if f.split('.')[-2] in datatype]
            
            # LOAD FILES
            loaded = 0
            if len(matching) == 0:
                print('    *.{}: 0'.format(datatype[0]))
                continue
            else:
                data[path][datatype[0]] = {}
                for f in matching:
                    ftype = f.split('.')[-1]
                    fshort = f.split('/')[-1]
                    if ftype in ['hdf5']:
                        data[path][datatype[0]][fshort] = h5load(f)
                        loaded += 1
                    elif ftype in ['bun']:
                        data[path][datatype[0]][fshort] = sddsload(f)
                        loaded += 1
                    elif ftype in ['log', 'lte']:
                        with open(f, 'r') as f_h:
                            data[path][datatype[0]][fshort] = f_h.read()
                        loaded += 1
                    elif ftype in ['twiss', 'coo', 'cen', 'srdts']:
                        data[path][datatype[0]][fshort] = sddsload(f)
                        loaded += 1
                    else:
                        print('    ---- WARNING: SKIPPED: ', f)

            print('    *.{}: {}'.format(datatype[0], loaded))
        
        if progress:
            progress.value = i + 1
            progress.description = '{:5.1f} %'.format(100 * float((i + 1)) / Npaths)

    return


def dataloadmenu(datapath, loadtoram=True, reverse=True, selected=False):
    button1 = widgets.Button(description="Load Selected", layout={'width': '9%'})
    button2 = widgets.Button(description="Open RUNDIR", layout={'width': '9%'})
    button3 = widgets.Button(description="Select All", layout={'width': '9%'})
    button4 = widgets.Button(description="cd ..", layout={'width': '9%'})
    progress = prg = widgets.FloatProgress(
                value=0,
                min=0,
                max=1,
                step=1,
                description='{:5.1f} %'.format(0),
                bar_style='info',
                orientation='horizontal',
                layout={'width': '55%'}
            )
    out = widgets.Output()
    acc = widgets.Accordion(title='Details', children=[out], layout=widgets.Layout(width='98%'), selected_index=None)
    
    label = widgets.Label(value='Datadirs in: ' + datapath)
    selector = widgets.SelectMultiple(
        options = [''],
        value=[],
        rows=25,
        layout=widgets.Layout(width='98%'),
        description='',
        disabled=False
    )
    
    data = {}
    paths = []
        
    def reload_path(selector, reverse=False, selected=False):
        datapath = label.value.split(': ')[1]
        dirs = os.listdir(datapath)
        dirs.sort(reverse = reverse)
        if 'atic' in dirs:
            dirs.remove('atic')
        selector.options = ['{:03}:  '.format(i) + d + '/' for i, d in enumerate(dirs)]
        if selected:
            selector.value = ['{:03}:  '.format(dirs.index(selected)) + selected + '/']
        return
    
    reload_path(selector, reverse=reverse, selected=selected)
    
    def ClickFun_LoadSelected(_, data=data, paths=paths):
        datapath = label.value.split(': ')[1]
        data.clear()
        paths[:] = []
        # "linking function with output"
        with out:
            out.clear_output()
            print('Loading from: ' + datapath)
            paths += [datapath + p.split(':  ')[-1] for p in selector.value]
            paths.sort()
            if loadtoram:
                loaddirs(paths, data, datapath, progress=progress)


    def ClickFun_OpenRundir(_, data=data, paths=paths, reverse=reverse):
        data.clear()
        paths.clear()
        p = selector.value[0]
        datapath = label.value.split(': ')[1]
        datapath += p.split(':  ')[-1]
        label.value = 'Datadirs in: ' + datapath
        reload_path(selector, reverse=reverse)
        
    
    def ClickFun_GoBack(_, data=data, paths=paths, reverse=reverse):
        data.clear()
        paths.clear()
        datapath = '/'.join(label.value.split(': ')[1].split('/')[:-2])
        label.value = 'Datadirs in: ' + datapath + '/'
        reload_path(selector, reverse=reverse)
        
    
    def ClickFun_SelectAll(_, selector=selector):
        data.clear()
        paths.clear()
        selector.value = selector.options


    button1.on_click(ClickFun_LoadSelected)
    button2.on_click(ClickFun_OpenRundir)
    button3.on_click(ClickFun_SelectAll)
    button4.on_click(ClickFun_GoBack)
    return widgets.VBox([label, selector, widgets.HBox([button1, button2, button3, button4, progress]), acc]), data, paths