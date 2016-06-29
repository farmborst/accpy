# -*- coding: utf-8 -*-
''' accpy.gui.measure
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as Tk
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    import tkinter as Tk
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import close
from threading import Thread
from time import time
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
from ..visualize.stringformat import time2str, uc
from ..measure.tunes import measure_tunes


def showfigs(t0, status, figs, tabs):
    close('all')
    for fig, tab in zip(figs, tabs):
        # destroy all widgets in fram/tab and close all figures
        for widget in tab.winfo_children():
            widget.destroy()
        canvas = FigureCanvasTkAgg(fig, master=tab)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack()
        toolbar.pack()
        canvas.draw()
    timestring = time2str(time() - t0)
    status.set('finished, elapsed time: ' + timestring)


def runthread(status, tabs, f_simulate, argstuple):
    def run(*argstuple):
            t0 = time()
            status.set('running...')
            figs = f_simulate(*argstuple)
            showfigs(t0, status, figs, tabs[1:])
    # data plotting in new thread to keep gui (main thread&loop) responsive
    t_run = Thread(target=run, args=argstuple)
    # automatically let die with main thread -> no global stop required
    t_run.setDaemon(True)
    # start thread
    t_run.start()


oops = ('Ooops!\n Sorry, but this feature is not ready yet...')


def tunes(frame, w, h):
    def _start():
        mode = modemenu.get()
        f_HF = float(entry_f_HF.get())
        runthread(status, tabs, measure_tunes,
                  (mode, f_HF))
    def _load():
        global filename
        filename = askopenfilename()
        if filename[-5::] != '.hdf5':
            filestr.set('error: {} is not hdf5 file-type'.format(filename))
            showerror('ERROR', 'THIS IS NOT A HDF5 FILE')
        else:
            filestr.set(filename)
    def _mode(*args):
        mode = modemenu.get()
        if mode == 'From File':
            cs_button(tabs[0], 3, 1, 'Load', _load)

    tabs = cs_tabbar(frame, w, h, ['Menu', 'Radial', 'Axial', 'Longitudinal',
                                   'Overview'])

    # row 1
    cs_label(tabs[0], 1, 1, 'Mode')
    cs_label(tabs[0], 1, 2, 'Cavity frequency / MHz')

    # row 2
    modemenu = cs_dropd(tabs[0], 2, 1, ['Measurement',
                                        'From File'], action=_mode)
    entry_f_HF = cs_Dblentry(tabs[0], 2, 2, 499.667)

    # row 3
    # 3, 1Loadbutton

    # row 4
    filestr = cs_label(tabs[0], 4, 1, '')

    # last row, column
    cs_button(tabs[0], 10, 10, 'Start', _start)
    status = cs_label(tabs[0], 10, 11, '')
    return


def chromaticity(frame, w, h):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def quadscanmeas(frame, w, h):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def achroscan(frame, w, h):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return
