# -*- coding: utf-8 -*-
''' accpy.gui.simulate
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as Tk
except:
    import tkinter as Tk
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from threading import Thread
from .layout import tabbar, cs_str, cs_int, cs_lab
from ..simulate.lsd import lsd
from ..visualize.figures import plotstandards

oops = ('Ooops!\n Sorry, but this feature is not ready yet...')


def twisstrack(frame):
    mode = 'trackbeta'

    def _start():
        # data plotting in new thread to keep gui (main thread&loop) responsive
        t_run = Thread(target=runtrack(w, h, mode, lattice.get(), tabs[1:],
                                       entry_slice.get()))
        # automatically let die with main thread -> no global stop required
        t_run.setDaemon(True)
        # start thread
        t_run.start()

    tablabels = ['Menu', 'Radial', 'Axial', 'Dispersion', 'Overview',
                 'Parameters']
    w, h = frame.winfo_screenwidth(), frame.winfo_screenheight()
    tabs = tabbar(frame, tablabels, w, h)

    lattices = {'bessy2injectionline': '1',
                'bessy2booster': '2',
                'bessy2transfer': '3',
                'bessy2transfer': '4',
                'bessy2ring': '5',
                'bessy2bigsmallbooster': '6'}
    lattice = cs_str('bessy2booster')

    cs_lab(tabs[0], 'Lattice', 1, 2)
    cs_lab(tabs[0], 'Nr. of slices', 1, 3)

    dropdown_lat = Tk.OptionMenu(tabs[0], lattice, *lattices)
    dropdown_lat.grid(row=2, column=2)

    entry_slice = Tk.Entry(tabs[0], textvariable=cs_int(1e3))
    entry_slice.grid(row=2, column=3)

    button_start = Tk.Button(master=tabs[0], text='Start', command=_start)
    button_start.grid(row=2, column=1)
    return


def parttrack(frame):
    mode = 'trackpart'

    def _start():
        # data plotting in new thread to keep gui (main thread&loop) responsive
        t_run = Thread(target=runtrack(w, h, mode, lattice.get(), tabs[1:],
                                       entry_slice.get(), entry_part.get(),
                                       entry_rounds.get()))
        # automatically let die with main thread -> no global stop required
        t_run.setDaemon(True)
        # start thread
        t_run.start()

    tablabels = ['Menu', 'X', 'X\'', 'Y', 'Y\'', 'Z', 'Z\'',
                 'Overview', 'Transverse phase space']
    w, h = frame.winfo_screenwidth(), frame.winfo_screenheight()
    tabs = tabbar(frame, tablabels, w, h)

    lattices = {'bessy2injectionline': '1',
                'bessy2booster': '2',
                'bessy2transfer': '3',
                'bessy2transfer': '4',
                'bessy2ring': '5',
                'bessy2bigsmallbooster': '6'}
    lattice = cs_str('bessy2booster')

    cs_lab(tabs[0], 'Lattice', 1, 2)
    cs_lab(tabs[0], 'Nr. of slices', 1, 3)
    cs_lab(tabs[0], 'Nr. of particles', 1, 4)
    cs_lab(tabs[0], 'Nr. of rounds', 1, 5)

    dropdown_lat = Tk.OptionMenu(tabs[0], lattice, *lattices)
    dropdown_lat.grid(row=2, column=2)

    entry_slice = Tk.Entry(tabs[0], textvariable=cs_int(1e3))
    entry_slice.grid(row=2, column=3)

    entry_part = Tk.Entry(tabs[0], textvariable=cs_int(1e2))
    entry_part.grid(row=2, column=4)

    entry_rounds = Tk.Entry(tabs[0], textvariable=cs_int(1e2))
    entry_rounds.grid(row=2, column=5)

    button_start = Tk.Button(master=tabs[0], text='Start', command=_start)
    button_start.grid(row=2, column=1)
    return


def runtrack(w, h, mode, latt, tabs, slices, particles=1, rounds=1):
    plotstandards('pcdisplay', [1, 1], w=w, h=h)
    figs = lsd(latt, int(slices), mode, particles=particles, rounds=rounds)
    for fig, tab in zip(figs, tabs):
        # destroy all widgets in fram/tab
        for widget in tab.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=tab)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack()
        toolbar.pack()
        canvas.draw()


def ebdynamics(frame):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def emitdynamics(frame):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def quadscansim(frame):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return
