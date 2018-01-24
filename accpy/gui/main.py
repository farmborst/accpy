# -*- coding: utf-8 -*-
''' accpy.gui.menu
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division, print_function
try:
    from Tkinter import (Tk, PhotoImage, Menu, TOP, X, Button, LEFT,
                         RIGHT, Label, StringVar, RAISED, FLAT)
    from ttk import Frame
except:
    from tkinter import (Tk, PhotoImage, Menu, TOP, X, Button, LEFT,
                         RIGHT, Label, StringVar, RAISED, FLAT)
    from tkinter.ttk import Frame
from .file import latticeeditor, settings, defaults
from .simulate import (gui_twisstrack, gui_parttrack, gui_ramp,
                       gui_quadscansim)
from .measure import (tunes, chromaticity, quadscanmeas, achroscan)
from .optimize import (emittex, twissmatch)
from .help import (documentation, about)
from ..visualize.figures import plotstandards
from ..dataio.hdf5 import confload


class MainApp:
    def __init__(self, parent, version, cwd):
        # add arguments to class attributes
        self.parent = parent
        self.version = version
        self.cwd = cwd

        # style the root window
        w = self.parent.winfo_screenwidth()
        h = self.parent.winfo_screenheight()
        icon = PhotoImage(file= self.cwd + '/accpy/icons/' + 'icon.gif')
        self.parent.geometry('{}x{}'.format(w, h))
        self.parent.tk.call('wm', 'iconphoto', self.parent._w, icon)
        self.parent.wm_title("ACCPY gui {}".format(self.version))

        # create gui widgets
        self.menubar = MenuBar(self.parent)
#        self.toolbar = ToolBar(self.parent)
#        self.mainwin = MainWindow(self.parent)

        # add gui elements to root window
        self.parent.config(menu=self.menubar.bar)
#        self.toolbar.frame.pack()
#        self.mainwin.frame.pack()


class MenuBar:
    def __init__(self, parent):
        self.parent = parent
        self.bar = Menu(self.parent)
        self.FileMenu()
        self.SimulationMenu()
        
    def MenuMaker(self, menu, items, cmnds):
        for i, c in zip(items, cmnds):
            menu.add_command(label=i, command=c)
            menu.add_separator()
        self.bar.add_cascade(label="File", menu=self.filemenu)
        
    def FileMenu(self):
        items = ['Lattice editor',
                 'Settings',
                 'Quit']
        cmnds = [lambda: print(1),
                 lambda: print(2),
                 lambda: print(3)]
        self.filemenu = Menu(self.bar, tearoff=0)
        self.MenuMaker(self.filemenu, items, cmnds)
        
    
    def SimulationMenu(self):
        items = ['Betamatrix and dispersion tracking',
                 'Particle tracking',
                 'Ramp',
                 'Quadrupole scan']
        cmnds = [lambda: print(1),
                 lambda: print(2),
                 lambda: print(3),
                 lambda: print(4)]
        self.simulationmenu = Menu(self.bar, tearoff=0)
        self.MenuMaker(self.filemenu, items, cmnds)

class ToolBar:
    def __init__(self, parent):
        self.parent = parent
        self.frame = Frame(self.parent, relief=RAISED)
        self.buttons()
        self.indicators()
    
    def buttons(self):
        self.run = Button(self.parent, relief=FLAT)
        self.run.pack(side=LEFT, padx=2, pady=2)
    
    def indicators(self):
        self.status = StringVar()
        self.status.set('Status')
        Label(self.parent, textvariable=self.status).pack(side=RIGHT)


class MainWindow:
    def __init__(self, parent):
        self.parent = parent
        self.frame = Frame(self.parent)