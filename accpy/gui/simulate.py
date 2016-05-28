# -*- coding: utf-8 -*-
''' accpy.gui.simulate
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as Tk
    import ttk
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    import tkinter as Tk
    import tkinter.ttk as ttk
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from threading import Thread
from .layout import tabbar
from ..simulate.lsd import lsd


def simmenu(root):
    class menu():
        @staticmethod
        def _start():
            # data plotting in new thread to keep gui (main thread&loop) responsive
            t_run = Thread(target=run(lattice.get(), tabs[1:]))
            # automatically let die with main thread -> no global stop required
            t_run.setDaemon(True)
            # start thread
            t_run.start()

        @staticmethod
        def _latt(*args):
            latt_ = latticestr.get()
            lattice.set(latt_)
    root.wm_title("accpy gui - simulations")
    tablabels = ['Menu', 'Optic', 'Parameters']
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    tabs = tabbar(root, tablabels, w, h)

    lattices = {'bessy2injectionline'   : '1',
                'bessy2booster' : '2',
                'bessy2transfer' : '3',
                'bessy2transfer': '4',
                'bessy2ring': '5',
                'bessy2bigsmallbooster' : '6'}

    lattice = Tk.StringVar()
    latticestr = Tk.StringVar()
    lattice.set('bessy2booster')
    latticestr.set(lattice.get())


    dropdown_lat = Tk.OptionMenu(tabs[0], latticestr, *lattices)
    dropdown_lat.grid(row=1, column=2)
    latticestr.trace('w', menu._latt)
    button_start = Tk.Button(master=tabs[0], text='Start', command=menu._start)
    button_start.grid(row=1, column=1)

    def run(latt, tabs, slicing=int(1e4)):
        figs = lsd(latt, slic=slicing, save=False, ft='pdf',
                   plotstandard='presentation_1920x1080', scale=[1, 1])
        for fig, tab in zip(figs, tabs):
            canvas = FigureCanvasTkAgg(fig, master=tab)
            canvas.get_tk_widget().pack(side=Tk.TOP, fill=Tk.BOTH, expand=1)
            fig.canvas.draw()


