# -*- coding: utf-8 -*-
''' accpy.gui.measure
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    from Tkinter import N, E, S, W, LabelFrame, _setit, Label, BOTH
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    from tkinter import N, E, S, W, LabelFrame, _setit, Label, BOTH
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import figure, close
from threading import Thread
from time import time
from numpy import zeros, where, sort, concatenate, array, loadtxt
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
from ..lattices.reader import lattlist, latt2py
from ..visualize.stringformat import time2str, uc
from ..measure.tunes import measure_tunes
from ..simulate.quadscan import measure_quadscan


def initfigs(tabs):
    close('all')
    figs = []
    for tab in tabs:
        # destroy all widgets in fram/tab and close all figures
        for widget in tab.winfo_children():
            widget.destroy()
        fig = figure()
        canvas = FigureCanvasTkAgg(fig, master=tab)
        figs.append(fig)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack()
        toolbar.pack()
    return figs


def runthread(status, f_simulate, argstuple):
    def run(*argstuple):
            t0 = time()
            status.set('running...')
            f_simulate(*argstuple)
            timestring = time2str(time() - t0)
            status.set('finished, elapsed time: ' + timestring)
    # data plotting in new thread to keep gui (main thread&loop) responsive
    t_run = Thread(target=run, args=argstuple)
    # automatically let die with main thread -> no global stop required
    t_run.setDaemon(True)
    # start thread
    t_run.start()


oops = ('Ooops!\n Sorry, but this feature is not ready yet...')


def tunes(frame, w, h, status, start, stop):
    def _start():
        figs = initfigs(tabs[1:4])
        filename = filestr.get()
        mode = modemenu.get()
        f_rf = float(entry_f_HF.get())
        h = float(entry_h.get())
        bunch = float(entry_bunch.get())
        steps = int(entry_steps.get())
        runthread(status, measure_tunes,
                  (figs, tunestr, mode, filename, f_rf, h, bunch, steps))

    def _load():
        filename = askopenfilename(initialdir='accpy/exampledata/')
        if filename[-5::] != '.hdf5':
            filestr.set('error: {} is not hdf5 file-type'.format(filename))
            showerror('ERROR', 'THIS IS NOT A HDF5 FILE')
        else:
            filestr.set(filename)

    def _mode(*args):
        mode = modemenu.get()
        if mode == 'From File':
            cs_button(tabs[0], 3, 1, 'Load', _load)

    start.configure(command=_start)
    tabs = cs_tabbar(frame, w, h, ['Menu', 'Radial', 'Axial', 'Longitudinal',
                                   'Overview'])

    # row 1
    cs_label(tabs[0], 1, 1, 'Mode')
    cs_label(tabs[0], 1, 2, 'Cavity frequency / MHz')
    cs_label(tabs[0], 3, 2, 'harmonic number')
    cs_label(tabs[0], 5, 2, 'bunch nr.')
    cs_label(tabs[0], 7, 2, 'steps')

    # row 2
    modemenu = cs_dropd(tabs[0], 2, 1, ['Measurement',
                                        'From File'], action=_mode)
    entry_f_HF = cs_Dblentry(tabs[0], 2, 2, 499.667)
    entry_h = cs_Dblentry(tabs[0], 4, 2, 160)
    entry_bunch = cs_Dblentry(tabs[0], 6, 2, 13)
    entry_steps = cs_Intentry(tabs[0], 8, 2, 50)

    # row 3
    # 3, 1Loadbutton

    # row 4
    filestr = cs_label(tabs[0], 4, 1, '')

    # add tunes to overview tab
    tunestr = [cs_label(tabs[4], i, 0, 'NaN', fg='green') for i in range(3)]
    return


def chromaticity(frame, w, h, status, start, stop):
    txt = Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def quadscanmeas(frame, w, h, status, start, stop):
    def _start():
        figs = initfigs([lf_OUT])
        mode = modemenu.get()
        if mode == 'From File':
            data = loadtxt(filestr.get())
            ki = None
            kf = None
            points = None
        if mode == 'Measurement':
            data = None
            ki = float(entry_ki.get())
            kf = float(entry_kf.get())
            points = int(entry_points.get())
        lattice = latticemenu.get()
        if lattice == 'drift':
            qL = float(entry_qL.get())
            driftlength = float(entry_dlen.get())
            UC = zeros([6, 1])
            UC[1] = driftlength
        else:
            quad = int(quadmenuval.get())
            screen = screenmenuval.get()
            _, _, _, UC, diagnostics, _, _, _, _, _, _, _ = latt2py(lattice, False)
            elements = UC[0, :]
            # index after selected quad (quad itself not required)
            i = sort(concatenate((where(elements==3), where(elements==4)), 1))[0][quad-1]+1
            # index before selected fom (next would be fom itself)
            fom = where(array(diagnostics) == screen)[0][0]
            f = where(elements==7)[0][fom]
            if i > f:
                showerror(title='ERROR', message='Please choose a quad before chosen screen')
                return
            qL = UC[1, i-1]
            UC = UC[:, i:f]
        kxi = float(entry_kxi.get())
        kxf = float(entry_kxf.get())
        kyi = float(entry_kyi.get())
        kyf = float(entry_kyf.get())
        kr_fit = [kxi, kxf, kyi, kyf]
        kr_mes = [ki, kf]
        epss = (float(entry_epss.get())/1e3)**2
        Dx = float(entry_Dx.get())
        Dpx = float(entry_Dpx.get())
        energy = float(entry_energy.get())*1e6
        particle = 'electron'
        runthread(status, measure_quadscan,
                  (figs, data, kr_fit, kr_mes, points, qL, UC, epss, Dx, Dpx,
                   energy, particle))

    def _load():
        filename = askopenfilename(initialdir='accpy/exampledata/')
        filestr.set(filename)

    def _mode(*args):
        mode = modemenu.get()
        if mode == 'From File':
            cs_button(lf_data, 1, 0, 'Load', _load)
        if mode == 'Measurement':
            showerror('Sorry', 'this feature is not ready')
            modemenu.set('From File')
            cs_button(lf_data, 1, 0, 'Load', _load)

    def _check(*args):
        lattice = latticemenu.get()
        def clearfigs1(tabs):
            close('all')
            for tab in tabs:
                # destroy all widgets in fram/tab and close all figures
                for widget in tab.winfo_children():
                    widget.destroy()

        def showfigs1(figs, tabs):
            clearfigs1(tabs)
            for fig, tab in zip(figs, tabs):
                canvas = FigureCanvasTkAgg(fig, master=tab)
                canvas.get_tk_widget().pack()
                canvas.draw()
        if lattice == 'drift':
            quadN.grid_remove()
            quadmenu.grid_remove()
            screenN.grid_remove()
            screenmenu.grid_remove()
            quadlab.grid()
            quadent.grid()
            driftlab.grid()
            driftent.grid()
        else:
            _, _, _, UC, diagnostics, _, _, _, _, _, _, _ = latt2py(lattice, False)
            elements = list(UC[0,:])
            quads = range(1, (elements.count(3)+elements.count(4)+1))
            quadlab.grid_remove()
            quadent.grid_remove()
            driftlab.grid_remove()
            driftent.grid_remove()
            screenmenu['menu'].delete(0, 'end')
            quadmenu['menu'].delete(0, 'end')
            [screenmenu['menu'].add_command(label=fom, command=_setit(screenmenuval, fom)) for fom in diagnostics]
            [quadmenu['menu'].add_command(label=quad, command=_setit(quadmenuval, quad)) for quad in quads]
            screenN.grid()
            quadN.grid()
            screenmenu.grid()
            quadmenu.grid()

    start.configure(command=_start)
    frame.pack(fill=BOTH, expand=1)
    lf_upbeam = LabelFrame(frame, text="Upstream longitudinal beam parameters", padx=5, pady=5)
    lf_upbeam.grid(row=0, column=0, sticky=W+E+N+S, padx=10, pady=10)
    lf_transfer = LabelFrame(frame, text="Transport matrix", padx=5, pady=5)
    lf_transfer.grid(row=1, column=0, sticky=W+E+N+S, padx=10, pady=10)
    lf_fit = LabelFrame(frame, text="Fit range", padx=5, pady=5)
    lf_fit.grid(row=0, column=1, sticky=W+E+N+S, padx=10, pady=10)
    lf_data = LabelFrame(frame, text="Data acquisition", padx=5, pady=5)
    lf_data.grid(row=1, column=1, sticky=W+E+N+S, padx=10, pady=10)
    lf_OUT = LabelFrame(frame, text="Results", padx=5, pady=5)
    lf_OUT.grid(row=2, column=0, columnspan=2, sticky=W+E+N+S, padx=10, pady=10)

    cs_label(lf_upbeam, 0, 0, 'Beam energy / MeV')
    cs_label(lf_upbeam, 1, 0, uc.delta+' / '+uc.ppt)
    cs_label(lf_upbeam, 2, 0, 'D / m')
    cs_label(lf_upbeam, 3, 0, 'D\' / rad')
    entry_energy = cs_Dblentry(lf_upbeam, 0, 1, 52)
    entry_epss = cs_Dblentry(lf_upbeam, 1, 1, 2.4)
    entry_Dx = cs_Dblentry(lf_upbeam, 2, 1, 0.)
    entry_Dpx = cs_Dblentry(lf_upbeam, 3, 1, 0.)

    cs_label(lf_transfer, 0, 0, 'Optic', sticky=W)
    _, openlatts = lattlist()
    latticemenu = cs_dropd(lf_transfer, 1, 0, ['drift'] + openlatts, action=_check)
    quadlab = cs_label(lf_transfer, 0, 1, 'Quadrupole length / m', retlab=True)[1]
    quadent = entry_qL = cs_Dblentry(lf_transfer, 0, 2, .119)
    driftlab = cs_label(lf_transfer, 1, 1, 'Driftlength / m', retlab=True)[1]
    driftent = entry_dlen = cs_Dblentry(lf_transfer, 1, 2, 1.8)
    quadN = cs_label(lf_transfer, 0, 1, 'Which quadrupole', retlab=True)[1]
    quadmenuval, quadmenu = cs_dropd(lf_transfer, 1, 1, [''], retbut=True)
    screenN = cs_label(lf_transfer, 2, 1, 'Which screen', retlab=True)[1]
    screenmenuval, screenmenu = cs_dropd(lf_transfer, 3, 1, openlatts, retbut=True)
    quadN.grid_remove()
    quadmenu.grid_remove()
    screenN.grid_remove()
    screenmenu.grid_remove()
    quadlab.grid_remove()
    quadent.grid_remove()
    driftlab.grid_remove()
    driftent.grid_remove()

    modemenu = cs_dropd(lf_data, 0, 0, ['Measurement',
                                        'From File'], action=_mode)
    filestr = cs_label(lf_data, 1,1, '')


    cs_label(lf_fit, 0, 0, 'k_x / (1/m'+uc.squared+')')
    cs_label(lf_fit, 1, 0, 'k_y / (1/m'+uc.squared+')')
    cs_label(lf_fit, 2, 0, 'Points')
    entry_kxi = cs_Dblentry(lf_fit, 0, 1, -8)
    entry_kxf = cs_Dblentry(lf_fit, 0, 2, 8)
    entry_kyi = cs_Dblentry(lf_fit, 1, 1, -8)
    entry_kyf = cs_Dblentry(lf_fit, 1, 2, 8)
    entry_points = cs_Intentry(lf_fit, 2, 1, 1000)
    return


def achroscan(frame, w, h, status, start, stop):
    txt = Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return
