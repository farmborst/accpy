# -*- coding: utf-8 -*-
''' accpy.gui.measure
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    from Tkinter import N, E, S, W, LabelFrame, _setit, Label
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    from tkinter import N, E, S, W, LabelFrame, _setit, Label
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
    figs, canvass = [], []
    for tab in tabs:
        # destroy all widgets in fram/tab and close all figures
        for widget in tab.winfo_children():
            widget.destroy()
        fig = figure()
        figs.append(fig)
        canvas = FigureCanvasTkAgg(fig, master=tab)
        canvass.append(canvas)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack()
        toolbar.pack()
        canvas.draw()
    return figs, canvass


def runthread(status, tabs, f_simulate, argstuple):
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


def tunes(frame, w, h):
    def _start():
        figs, canvass = initfigs(tabs[1:4])
        filename = filestr.get()
        mode = modemenu.get()
        f_rf = float(entry_f_HF.get())
        h = float(entry_h.get())
        bunch = float(entry_bunch.get())
        steps = int(entry_steps.get())
        runthread(status, tabs, measure_tunes,
                  (figs, canvass, tunestr, mode, filename, f_rf, h,bunch, steps))
    def _load():
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

    # last row, column
    cs_button(tabs[0], 10, 10, 'Start', _start)
    status = cs_label(tabs[0], 10, 11, '')
    return


def chromaticity(frame, w, h):
    txt = Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return


def quadscanmeas(frame, w, h):
    def _start():
        lattice = latticemenu.get()
        qL = float(entry_qL.get())
        if lattice == 'drift':
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
            UC = UC[:, i:f]
        ki = float(entry_ki.get())
        kf = float(entry_kf.get())
        points = int(entry_points.get())
        epsx = float(entry_epsx.get())/1e9
        betx = float(entry_betx.get())
        alpx = float(entry_alpx.get())
        epsy = float(entry_epsy.get())/1e9
        bety = float(entry_bety.get())
        alpy = float(entry_alpy.get())
        epss = (float(entry_epss.get())/1e3)**2
        Dx = float(entry_Dx.get())
        Dpx = float(entry_Dpx.get())
        energy = float(entry_energy.get())*1e6
        particle = 'electron'
        if filestr.get() != '':
            data = loadtxt(filestr.get())
        else:
                data = None
        if betx < 0 or bety < 0:
            showerror('ERROR', 'beta function must be positive')
            return
        runthread(status, tabs, measure_quadscan,
                  (ki, kf, qL, UC, points, epsx, betx, alpx,
                   epsy, bety, alpy, epss, Dx, Dpx, energy, particle, data))

    def _load():
        global filename
        filename = askopenfilename()
        if filename[-5::] != '.hdf5':
            #filestr.set('error: {} is not hdf5 file-type'.format(filename))
            #showerror('ERROR', 'THIS IS NOT A HDF5 FILE')
            filestr.set(filename)
        else:
            filestr.set(filename)

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


    tabs = cs_tabbar(frame, w, h, ['Menu', 'Beam extents'])
    lf_upbeam = LabelFrame(tabs[0], text="Upstream beam parameters", padx=5, pady=5)
    lf_upbeam.grid(row=1, column=0, sticky=W+E+N+S, padx=10, pady=10)
    lf_transfer = LabelFrame(tabs[0], text="Transport matrix", padx=5, pady=5)
    lf_transfer.grid(row=2, column=0, sticky=W+E+N+S, padx=10, pady=10)
    lf_quadrupole = LabelFrame(tabs[0], text="Quadrupole range", padx=5, pady=5)
    lf_quadrupole.grid(row=1, column=1, sticky=W+E+N+S, padx=10, pady=10)
    lf_data = LabelFrame(tabs[0], text="Data comparison", padx=5, pady=5)
    lf_data.grid(row=2, column=1, sticky=W+E+N+S, padx=10, pady=10)
    cs_button(tabs[0], 4, 3, 'Start', _start)
    status = cs_label(tabs[0], 4, 4, '')

    cs_label(lf_upbeam, 1, 2, uc.epsilon+' / nm rad')
    cs_label(lf_upbeam, 1, 3, uc.beta+' / m')
    cs_label(lf_upbeam, 1, 4, uc.alpha+'/ rad')
    cs_label(lf_upbeam, 2, 1, 'Radial')
    entry_epsx = cs_Dblentry(lf_upbeam, 2, 2, 202)
    entry_betx = cs_Dblentry(lf_upbeam, 2, 3, 6.37)
    entry_alpx = cs_Dblentry(lf_upbeam, 2, 4, -0.13)
    cs_label(lf_upbeam, 3, 1, 'Axial')
    entry_epsy = cs_Dblentry(lf_upbeam, 3, 2, 144)
    entry_bety = cs_Dblentry(lf_upbeam, 3, 3, 4.1)
    entry_alpy = cs_Dblentry(lf_upbeam, 3, 4, -0.51)
    cs_label(lf_upbeam, 4, 2, uc.delta+' / '+uc.ppt)
    cs_label(lf_upbeam, 4, 3, 'D / m')
    cs_label(lf_upbeam, 4, 4, 'D\' / rad')
    cs_label(lf_upbeam, 5, 1, 'Longitudinal')
    entry_epss = cs_Dblentry(lf_upbeam, 5, 2, 2.4)
    entry_Dx = cs_Dblentry(lf_upbeam, 5, 3, 0.)
    entry_Dpx = cs_Dblentry(lf_upbeam, 5, 4, 0.)

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

    cs_button(lf_data, 0, 0, 'Load', _load)
    filestr = cs_label(lf_data, 1, 0, '')

    cs_label(lf_quadrupole, 0, 0, 'Beam energy / MeV')
    entry_energy = cs_Dblentry(lf_quadrupole, 0, 1, 52)
    cs_label(lf_quadrupole, 1, 1, 'Quadrupole strength k / (1/m'+uc.squared+')')
    cs_label(lf_quadrupole, 2, 0, 'Initial (k<0 = axial focus)')
    entry_ki = cs_Dblentry(lf_quadrupole, 2, 1, -8)
    cs_label(lf_quadrupole, 3, 0, 'Final')
    entry_kf = cs_Dblentry(lf_quadrupole, 3, 1, 8)
    cs_label(lf_quadrupole, 4, 0, 'steps')
    entry_points = cs_Intentry(lf_quadrupole, 4, 1, 1000)
    return


def achroscan(frame, w, h):
    txt = Label(frame, text=oops, font=("Helvetica", 20))
    txt.pack()
    return
