# -*- coding: utf-8 -*-
''' accpy.gui.main
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division, print_function
try:
    import Tkinter as tk
except:
    import tkinter as tk


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
        self.menubar = MenuBar(self.root)
        self.toolbar = ToolBar(self.root, self.cwd )
        self.mainwin = MainWindow(self.root)

        # add gui elements to root window
        self.root.config(menu=self.menubar.bar)
        self.toolbar.frame.pack(side=tk.TOP, fill=tk.X)
        self.mainwin.frame.pack(expand=True)


class MenuBar:
    def __init__(self, parent):
        self.parent = parent
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

    def FileMenu(self):
        items = ['Lattice editor',
                 'Settings',
                 'Quit']
        cmnds = [lambda: print(1),
                 lambda: print(2),
                 lambda: (self.parent.quit(),
                          self.parent.destroy())]
        filemenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(filemenu, items, cmnds, 'File')

    def SimulationMenu(self):
        items = ['Betamatrix and dispersion tracking',
                 'Particle tracking',
                 'Ramp',
                 'Quadrupole scan']
        cmnds = [lambda: print(1),
                 lambda: print(2),
                 lambda: print(3),
                 lambda: print(4)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Simulation')

    def MeasureMenu(self):
        items = ['Tunes',
                 'Chromaticity',
                 'Quadrupole scan',
                 'Achromat scan']
        cmnds = [lambda: print(1),
                 lambda: print(2),
                 lambda: print(3),
                 lambda: print(4)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Measurement')

    def OptimizeMenu(self):
        items = ['Find Transverse emittance exchange section',
                 'Twiss matching']
        cmnds = [lambda: print(1),
                 lambda: print(2)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Optimization')

    def HelpMenu(self):
        items = ['Documentation',
                 'About ACCPY...']
        cmnds = [lambda: print(1),
                 lambda: print(2)]
        simulationmenu = tk.Menu(self.bar, tearoff=0)
        self.MenuMaker(simulationmenu, items, cmnds, 'Help')


class ToolBar:
    def __init__(self, parent, cwd):
        self.parent = parent
        self.frame = tk.Frame(self.parent, relief=tk.RAISED)
        self.buttons(cwd)
        self.indicators()

    def buttons(self, cwd):
        path = cwd + '/accpy/icons/' + 'start.gif'
        icon = tk.PhotoImage(file=path).subsample(10, 10)
        self.run = tk.Button(self.parent, relief=tk.FLAT, image=icon)
        self.run.pack(side=tk.LEFT, padx=2, pady=2)

    def indicators(self):
        self.status = tk.StringVar()
        self.status.set('Status')
        tk.Label(self.parent, textvariable=self.status).pack(side=tk.RIGHT)


class MainWindow:
    def __init__(self, parent):
        self.parent = parent
        self.frame = tk.Frame(self.parent)
