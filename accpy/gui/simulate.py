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
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
from ..simulate.lsd import lsd
from ..simulate.ramp import simulate_ramp

oops = ('Ooops!\n Sorry, but this feature is not ready yet...')


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
    else:
        timestring = '%g milliseconds' % (ms)
    return timestring


def runthread(go):
    # data plotting in new thread to keep gui (main thread&loop) responsive
    t_run = Thread(target=go)
    # automatically let die with main thread -> no global stop required
    t_run.setDaemon(True)
    # start thread
    t_run.start()


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


def runtrack(status, mode, latt, tabs, slices, particles=1, rounds=1):
    t0 = time()
    status.set('running...')
    figs = lsd(latt, slices, mode, particles=particles, rounds=rounds)
    showfigs(t0, status, figs, tabs)


def gui_twisstrack(frame, w, h):
    mode = 'trackbeta'

    def _start():
        latt = lattice.get()
        slic = int(entry_slice.get())
        go = partial(runtrack, *(status, mode, latt, tabs[1:], slic))
        runthread(go)

    tabs = cs_tabbar(frame, w, h, ['Menu', 'Radial', 'Axial', 'Dispersion',
                                   'Overview', 'Parameters'])

    cs_label(tabs[0], 1, 1, 'Lattice')
    cs_label(tabs[0], 1, 2, 'Nr. of slices')
    lattice = cs_dropd(tabs[0], 2, 1, ['bessy2injectionline',
                                       'bessy2booster',
                                       'bessy2transfer',
                                       'bessy2ring'])
    entry_slice = cs_Intentry(tabs[0], 2, 2, 1e3)

    cs_button(tabs[0], 3, 3, 'Start', _start)
    status = cs_label(tabs[0], 3, 4, '')
    return


def gui_parttrack(frame, w, h):
    mode = 'trackpart'

    def _start():
        latt = lattice.get()
        slic = int(entry_slice.get())
        prts = int(entry_parts.get())
        rnds = int(entry_round.get())
        go = partial(runtrack, *(status, mode, latt, tabs[1:], slic, prts, rnds))
        runthread(go)

    tabs = cs_tabbar(frame, w, h, [' Menu ', ' X ', ' X\' ', ' Y ', ' Y\' ',
                                   ' Z ', ' Z\' ', ' Overview ',
                                   ' Transverse phase space '])

    cs_label(tabs[0], 1, 1, 'Lattice')
    cs_label(tabs[0], 1, 2, 'Nr. of slices')
    cs_label(tabs[0], 1, 3, 'Nr. of particles (parallelized)')
    cs_label(tabs[0], 1, 4, 'Nr. of rounds')
    lattice = cs_dropd(tabs[0], 2, 1, ['bessy2injectionline',
                                       'bessy2booster',
                                       'bessy2transfer',
                                       'bessy2ring'])
    entry_slice = cs_Intentry(tabs[0], 2, 2, 100)
    entry_parts = cs_Intentry(tabs[0], 2, 3, cpu_count())
    entry_round = cs_Intentry(tabs[0], 2, 4, 100)

    cs_button(tabs[0], 3, 5, 'Start', _start)
    status = cs_label(tabs[0], 3, 6, '')
    return


def gui_ramp(frame, w, h):
    def _start():
        def run(T, t_inj, t_ext, text2, E_inj, E_ext, latt, points, f_HF,
                V_HFs, emitxs, emitys, emitss):
            t0 = time()
            status.set('running...')
            figs = simulate_ramp(T, t_inj, t_ext, text2, E_inj, E_ext, latt,
                                 points, f_HF, V_HFs, emitxs, emitys, emitss)
            showfigs(t0, status, figs, tabs[1:])
        points = int(entry_pnts.get())
        T_per = float(entry_Tper.get())
        t_inj = float(entry_tinj.get())
        t_ext = float(entry_text.get())
        text2 = float(entry_tex2.get())
        E_inj = float(entry_Einj.get())*1e6
        E_ext = float(entry_Eext.get())*1e6
        latt = lattice.get()
        f_HF = float(entry_f_HF.get())*1e6
        V_HFs = [float(x)*1e3 for x in entry_V_HF.get().split()]
        emitxs = [float(x)*1e-9 for x in entry_emitx.get().split()]
        emitys = [float(x)*1e-9 for x in entry_emity.get().split()]
        emitss = [float(x)*1e-3 for x in entry_emits.get().split()]
        go = partial(run, *(T_per, t_inj, t_ext, text2, E_inj, E_ext, latt,
                            points, f_HF, V_HFs, emitxs, emitys, emitss))
        runthread(go)

    tabs = cs_tabbar(frame, w, h, ['Menu', 'Energy', 'Magnetic Flux',
                                   'Energy loss', 'Acceleration voltage',
                                   'Synchronous phase',
                                   'Synchrotron frequency', 'Bunch length',
                                   'Radial Emittance', 'Axial Emittance',
                                   'Longitudinal Emittance'])
    cs_label(tabs[0], 1, 1, 'Lattice')
    cs_label(tabs[0], 1, 2, 'Acceleration Period / s')
    cs_label(tabs[0], 1, 3, 'Injection time / s')
    cs_label(tabs[0], 1, 4, 'Extraction time 1 / s')
    cs_label(tabs[0], 1, 5, 'Extraction time 2 / s')
    lattice = cs_dropd(tabs[0], 2, 1, ['bessy2booster',
                                        'bessy2ring'])
    entry_Tper = cs_Dblentry(tabs[0], 2, 2, 1e-1)
    entry_tinj = cs_Dblentry(tabs[0], 2, 3, 5518.944e-6)
    entry_text = cs_Dblentry(tabs[0], 2, 4, 38377.114e-6)
    entry_tex2 = cs_Dblentry(tabs[0], 2, 5, 57076.1e-6)

    cs_label(tabs[0], 3, 1, 'Calculation points')
    cs_label(tabs[0], 3, 2, 'Cavity frequency / MHz')
    cs_label(tabs[0], 3, 3, 'Injection energy / MeV')
    cs_label(tabs[0], 3, 4, 'Extraction energy / MeV')
    entry_pnts = cs_Intentry(tabs[0], 4, 1, 1e3)
    entry_f_HF = cs_Dblentry(tabs[0], 4, 2, 499.667)
    entry_Einj = cs_Dblentry(tabs[0], 4, 3, 52.3)
    entry_Eext = cs_Dblentry(tabs[0], 4, 4, 1720)

    pi = u'\u03c0'
    ppt = u'\u2030'

    cs_label(tabs[0], 5, 1, 'Cavity peak Voltages / kV')
    cs_label(tabs[0], 5, 3, 'Emittance @ injection')
    cs_label(tabs[0], 6, 2, 'Radial')
    cs_label(tabs[0], 7, 2, 'Axial')
    cs_label(tabs[0], 8, 2, 'Longitudinal')
    cs_label(tabs[0], 6, 4, 'nm '+pi+' rad')
    cs_label(tabs[0], 7, 4, 'nm '+pi+' rad')
    cs_label(tabs[0], 8, 4, ppt)
    entry_V_HF = cs_Strentry(tabs[0], 6, 1, '200 500 2000')
    entry_emitx = cs_Strentry(tabs[0], 6, 3, '100 200 300')
    entry_emity = cs_Strentry(tabs[0], 7, 3, '100 200 300')
    entry_emits = cs_Strentry(tabs[0], 8, 3, '.5 1 1.5')

    cs_button(tabs[0], 9, 6, 'Start', _start)
    status = cs_label(tabs[0], 9, 7, '')
    return


def gui_quadscansim(frame, w, h):
    txt = Tk.Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return
