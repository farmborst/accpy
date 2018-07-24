# -*- coding: utf-8 -*-
''' accpy.gui.main
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division, print_function
try:
    import Tkinter as tk
    import ttk
except:
    import tkinter as tk
    import tkinter.ttk as ttk
from .file import latticeeditor, settings, defaults
from .simulate import (gui_twisstrack, gui_parttrack, gui_ramp,
                       gui_quadscansim)
from .measure import (tunes, chromaticity, quadscanmeas, achroscan)
from .optimize import (emittex, twissmatch)
from .help import (documentation, about)
from ..visualize.figures import plotstandards
from ..dataio.hdf5 import h5load, h5save


class MainApp:
    def __init__(self, root, version, cwd):
        # add arguments to class attributes
        self.root = root
        self.version = version
        self.cwd = cwd

        # style the root window
        w = self.root.winfo_screenwidth()
        h = self.root.winfo_screenheight()
        icon = tk.PhotoImage(file=self.cwd + '/accpy/icons/' + 'icon.gif')
        self.root.geometry('{}x{}'.format(w, h))
        self.root.tk.call('wm', 'iconphoto', self.root._w, icon)
        self.root.wm_title("ACCPY gui {}".format(self.version))

        # create gui widgets
        self.mainwin = MainWindow(self.root)
        self.toolbar = ToolBar(self.root, self.cwd)
        self.menubar = MenuBar(self.root, self.mainwin, self.toolbar,
                               self.version, w, h)

        # add gui elements to root window
        self.root.config(menu=self.menubar.bar)
        self.toolbar.frame.pack(side=tk.TOP, fill=tk.X)
        self.mainwin.frame.pack(expand=True)

        # test load settings or set defaults
        try:
            confdict = h5load('./settings.hdf5')
        except:
            print('Loading Default Configuration...')
            confdict = defaults('./settings.hdf5')
            print('    ... done.')
        plotstandards(confdict)


class MenuBar:
    def __init__(self, parent, mainwin, toolbar, version, w, h):
        self.parent = parent
        self.mainwin = mainwin
        self.toolbar = toolbar
        self.version = version
        self.w = w
        self.h = h
        self.bar = tk.Menu(self.parent)
        self.FileMenu()
        self.SimulationMenu()
        self.MeasureMenu()
        self.OptimizeMenu()
        self.HelpMenu()

    def MenuMaker(self, menu, items, cmnds, name, seps=[]):
        for n, (i, c) in enumerate(zip(items, cmnds)):
            menu.add_command(label=i, command=c)
            if n in seps:
                menu.add_separator()
        self.bar.add_cascade(label=name, menu=menu)

    def SubApp(self, name, app):
        self.parent.wm_title('accpy gui - ' + name)
        # destroy all widgets in fram/tab
        for widget in self.mainwin.frame.winfo_children():
            widget.destroy()
        self.toolbar.status.set('Status: idle')
        w = self.parent.winfo_screenwidth()
        h = self.parent.winfo_screenheight()
        app(self.mainwin.frame, w, h, self.toolbar.status, self.toolbar.run)

    def FileMenu(self):
        items = ['Lattice Editor',
                 'Settings',
                 'Quit']
        cmnds = [lambda: self.SubApp(items[0], latticeeditor),
                 lambda: settings(),
                 lambda: (self.parent.quit(),
                          self.parent.destroy())]
        filemenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(filemenu, items, cmnds, 'File')

    def SimulationMenu(self):
        items = ['Betamatrix and dispersion tracking',
                 'Particle tracking',
                 'Ramp',
                 'Quadrupole scan']
        cmnds = [lambda: self.SubApp(items[0], gui_twisstrack),
                 lambda: self.SubApp(items[1], gui_parttrack),
                 lambda: self.SubApp(items[2], gui_ramp),
                 lambda: self.SubApp(items[3], gui_quadscansim)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Simulation')

    def MeasureMenu(self):
        items = ['Tunes',
                 'Chromaticity',
                 'Quadrupole scan',
                 'Achromat scan']
        cmnds = [lambda: self.SubApp(items[0], tunes),
                 lambda: self.SubApp(items[1], chromaticity),
                 lambda: self.SubApp(items[2], quadscanmeas),
                 lambda: self.SubApp(items[3], achroscan)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Measurement')

    def OptimizeMenu(self):
        items = ['Find Transverse emittance exchange section',
                 'Twiss matching']
        cmnds = [lambda: self.SubApp(items[0], emittex),
                 lambda: self.SubApp(items[0], twissmatch)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Optimization')

    def HelpMenu(self):
        items = ['Documentation',
                 'About ACCPY...']
        cmnds = [lambda: documentation(self.version, self.w, self.h),
                 lambda: about(self.version),]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Help')


class ToolBar:
    def __init__(self, parent, cwd):
        self.parent = parent
        self.frame = ttk.Frame(self.parent, relief=tk.RAISED)
        self.buttons(cwd)
        self.indicators()

    def buttons(self, cwd):
        path = cwd + '/accpy/icons/' + 'start.gif'
        self.icon = tk.PhotoImage(file=path).subsample(10, 10)
        self.run = tk.Button(self.frame, relief=tk.FLAT, image=self.icon)
        self.run.pack(side=tk.LEFT, padx=2, pady=2)

    def indicators(self):
        self.status = tk.StringVar()
        self.status.set('Status: idle')
        tk.Label(self.frame, textvariable=self.status).pack(side=tk.LEFT)


class MainWindow:
    def __init__(self, parent):
        self.parent = parent
        self.frame = ttk.Frame(self.parent)
