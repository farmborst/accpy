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
from matplotlib.pyplot import close
from functools import partial
from threading import Thread
from multiprocessing import cpu_count
from time import time
from .layout import tabbar, cs_str, cs_int, cs_lab
from ..simulate.lsd import lsd
from ..visualize.figures import plotstandards

oops = ('Ooops!\n Sorry, but this feature is not ready yet...')


def twisstrack(frame):
    mode = 'trackbeta'

    def _start():
        latt = lattice.get()
        slic = int(entry_slice.get())
        go = partial(runtrack, *(status, w, h, mode, latt, tabs[1:], slic))
        # data plotting in new thread to keep gui (main thread&loop) responsive
        t_run = Thread(target=go)
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

    status = cs_lab(tabs[0], '', 3, 1)
    return


def parttrack(frame):
    mode = 'trackpart'

    def _start():
        latt = lattice.get()
        slic = int(entry_slice.get())
        prts = int(entry_part.get())
        rnds = int(entry_rounds.get())
        go = partial(runtrack, *(status, w, h, mode, latt, tabs[1:], slic, prts, rnds))
        # data plotting in new thread to keep gui (main thread&loop) responsive
        t_run = Thread(target=go)
        # automatically let die with main thread -> no global stop required
        t_run.setDaemon(True)
        # start thread
        t_run.start()

    tablabels = [' Menu ', ' X ', ' X\' ', ' Y ', ' Y\' ', ' Z ', ' Z\' ',
                 ' Overview ', ' Transverse phase space ']
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
    cs_lab(tabs[0], 'Nr. of particles (parallelized)', 1, 4)
    cs_lab(tabs[0], 'Nr. of rounds', 1, 5)

    dropdown_lat = Tk.OptionMenu(tabs[0], lattice, *lattices)
    dropdown_lat.grid(row=2, column=2)

    entry_slice = Tk.Entry(tabs[0], textvariable=cs_int(100))
    entry_slice.grid(row=2, column=3)

    entry_part = Tk.Entry(tabs[0], textvariable=cs_int(cpu_count()))
    entry_part.grid(row=2, column=4)

    entry_rounds = Tk.Entry(tabs[0], textvariable=cs_int(100))
    entry_rounds.grid(row=2, column=5)

    button_start = Tk.Button(master=tabs[0], text='Start', command=_start)
    button_start.grid(row=2, column=1)

    status = cs_lab(tabs[0], '', 3, 1)
    return


def runtrack(status, w, h, mode, latt, tabs, slices, particles=1, rounds=1):
    t0 = time()
    status.set('running...')
    plotstandards('pcdisplay', [1, 1], w=w, h=h)
    close('all')
    figs = lsd(latt, int(slices), mode, particles=int(particles),
               rounds=int(rounds))
    for fig, tab in zip(figs, tabs):
        # destroy all widgets in fram/tab and close all figures
        for widget in tab.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=tab)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack()
        toolbar.pack()
        canvas.draw()
    t1 = time() - t0
#    hh = t1//3600
#    tt = t1 % 3600
#    mm = tt//60
#    tt = tt % 60
#    ss = tt//1
#    ms = (tt % 1*1e3)//1
#    timestring = '%02i:%02i:%02i.%03i' % (hh, mm, ss, ms)
    timestring = time2str(t1)
    status.set('finished, elapsed time: '+timestring)


def time2str(t):
    dd = t//86400
    t %= 86400
    hh = t//3600
    t %= 3600
    mm = t//60
    t %= 60
    ss = t//1
    t %= 1
    ms = (t*1e3)//1
    if dd > 0:
        timestring = '%g days %g hours %g minutes %g.%g seconds' % (dd, hh, mm, ss, ms)
    elif hh > 0:
        timestring = '%g hours %g minutes %g.%g seconds' % (hh, mm, ss, ms)
    elif mm > 0:
        timestring = '%g minutes %g.%g seconds' % (mm, ss, ms)
    elif ss > 0:
        timestring = '%g.%g seconds' % (ss, ms)
    return timestring

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
