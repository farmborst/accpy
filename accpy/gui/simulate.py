# -*- coding: utf-8 -*-
''' accpy.gui.simulate
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    from Tkinter import N, E, S, W, LabelFrame, _setit
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    from tkinter import N, E, S, W, LabelFrame, _setit
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
from time import time
from threading import Thread
from multiprocessing import cpu_count
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import close
from numpy import zeros, where, sort, concatenate, array, loadtxt
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
from ..lattices.reader import lattlist, latt2py
from ..visualize.stringformat import time2str, uc
from ..simulate.lsd import lsd
from ..simulate.ramp import simulate_ramp
from ..simulate.quadscan import simulate_quadscan


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


def gui_twisstrack(frame, w, h):
    def _start():
        latt = latticemenu.get()
        if latt in closedlatts:
            closed = True
        elif latt in openlatts:
            closed = False
        else:
            showerror(title='ERROR', message='Please choose a lattice')
            return
        slic = int(entry_slice.get())
        mode = 'trackbeta'
        particles = 1
        rounds = 1
        runthread(status, tabs, lsd,
                  (closed, latt, slic, mode, particles, rounds))

    tabs = cs_tabbar(frame, w, h, ['Menu', 'Radial', 'Axial', 'Dispersion',
                                   'Overview', 'Parameters', 'Beam extents'])

    cs_label(tabs[0], 1, 1, 'Lattice')
    cs_label(tabs[0], 1, 2, 'Nr. of slices')
    closedlatts, openlatts = lattlist()
    latticemenu = cs_dropd(tabs[0], 2, 1, closedlatts + openlatts)
    entry_slice = cs_Intentry(tabs[0], 2, 2, 1e3)

    cs_button(tabs[0], 3, 3, 'Start', _start)
    status = cs_label(tabs[0], 3, 4, '')
    return


def gui_parttrack(frame, w, h):
    def _start():
        latt = latticemenu.get()
        slic = int(entry_slice.get())
        if latt in closedlatts:
            closed = True
            rnds = int(entry_round.get())
        elif latt in openlatts:
            closed = False
            rnds = int(1)
        else:
            showerror(title='ERROR', message='Please choose a lattice')
            return
        mode = 'trackpart'
        prts = int(entry_parts.get())

        runthread(status, tabs, lsd,
                  (closed, latt, slic, mode, prts, rnds))
    def _check(*args):
        lattice = latticemenu.get()
        if lattice in openlatts:
            roundslabel.set('')
            entry_round.grid_remove()
        else:
            roundslabel.set('Nr. of rounds')
            entry_round.grid(row=2, column=4)




    tabs = cs_tabbar(frame, w, h, [' Menu ', ' X ', ' X\' ', ' Y ', ' Y\' ',
                                   ' Z ', ' Z\' ', ' Overview ',
                                   ' Transverse phase space '])

    cs_label(tabs[0], 1, 1, 'Lattice')
    cs_label(tabs[0], 1, 2, 'Nr. of slices')
    cs_label(tabs[0], 1, 3, 'Nr. of particles (parallelized)')
    roundslabel = cs_label(tabs[0], 1, 4, 'Nr. of rounds')
    closedlatts, openlatts = lattlist()
    latticemenu = cs_dropd(tabs[0], 2, 1, closedlatts+openlatts, action=_check)
    entry_slice = cs_Intentry(tabs[0], 2, 2, 100)
    entry_parts = cs_Intentry(tabs[0], 2, 3, cpu_count())
    entry_round = cs_Intentry(tabs[0], 2, 4, 100)
    cs_button(tabs[0], 3, 5, 'Start', _start)
    status = cs_label(tabs[0], 3, 6, '')
    return


def gui_ramp(frame, w, h):
    def _start():
        points = int(entry_pnts.get())
        T_per = float(entry_Tper.get())
        t_inj = float(entry_tinj.get())
        t_ext = float(entry_text.get())
        text2 = float(entry_tex2.get())
        E_inj = float(entry_Einj.get())*1e6
        E_ext = float(entry_Eext.get())*1e6
        latt = lattice.get()
        if latt not in closedlatts:
            showerror(title='ERROR', message='Please choose a lattice')
            return
        f_HF = float(entry_f_HF.get())*1e6
        V_HFs = [float(x)*1e3 for x in entry_V_HF.get().split()]
        emitxs = [float(x)*1e-9 for x in entry_emitx.get().split()]
        emitys = [float(x)*1e-9 for x in entry_emity.get().split()]
        emitss = [float(x)*1e-3 for x in entry_emits.get().split()]
        runthread(status, tabs, simulate_ramp,
                  (T_per, t_inj, t_ext, text2, E_inj, E_ext, latt, points,
                   f_HF, V_HFs, emitxs, emitys, emitss))

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
    closedlatts, _ = lattlist()
    lattice = cs_dropd(tabs[0], 2, 1, closedlatts)
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

    cs_label(tabs[0], 5, 1, 'Cavity peak Voltages / kV')
    cs_label(tabs[0], 5, 3, 'Emittance @ injection')
    cs_label(tabs[0], 6, 2, 'Radial')
    cs_label(tabs[0], 7, 2, 'Axial')
    cs_label(tabs[0], 8, 2, 'Longitudinal')
    cs_label(tabs[0], 6, 4, 'nm '+uc.pi+' rad')
    cs_label(tabs[0], 7, 4, 'nm '+uc.pi+' rad')
    cs_label(tabs[0], 8, 4, uc.ppt)
    entry_V_HF = cs_Strentry(tabs[0], 6, 1, '700 2000')
    entry_emitx = cs_Strentry(tabs[0], 6, 3, '200 300')
    entry_emity = cs_Strentry(tabs[0], 7, 3, '200 300')
    entry_emits = cs_Strentry(tabs[0], 8, 3, '1')

    cs_button(tabs[0], 9, 6, 'Start', _start)
    status = cs_label(tabs[0], 9, 7, '')
    return


def gui_quadscansim(frame, w, h):
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
        data = loadtxt(filestr.get())
        if betx < 0 or bety < 0:
            showerror('ERROR', 'beta function must be positive')
            return
        runthread(status, tabs, simulate_quadscan,
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
    entry_epsx = cs_Dblentry(lf_upbeam, 2, 2, 60)
    entry_betx = cs_Dblentry(lf_upbeam, 2, 3, 42)
    entry_alpx = cs_Dblentry(lf_upbeam, 2, 4, -38)
    cs_label(lf_upbeam, 3, 1, 'Axial')
    entry_epsy = cs_Dblentry(lf_upbeam, 3, 2, 20)
    entry_bety = cs_Dblentry(lf_upbeam, 3, 3, 8)
    entry_alpy = cs_Dblentry(lf_upbeam, 3, 4, 2)
    cs_label(lf_upbeam, 4, 2, uc.delta+' / '+uc.ppt)
    cs_label(lf_upbeam, 4, 3, 'D / m')
    cs_label(lf_upbeam, 4, 4, 'D\' / rad')
    cs_label(lf_upbeam, 5, 1, 'Longitudinal')
    entry_epss = cs_Dblentry(lf_upbeam, 5, 2, 0.547)
    entry_Dx = cs_Dblentry(lf_upbeam, 5, 3, 0.7)
    entry_Dpx = cs_Dblentry(lf_upbeam, 5, 4, -0.15)

    cs_label(lf_transfer, 0, 0, 'Optic', sticky=W)
    _, openlatts = lattlist()
    latticemenu = cs_dropd(lf_transfer, 1, 0, ['drift'] + openlatts, action=_check)
    quadlab = cs_label(lf_transfer, 0, 1, 'Quadrupole length / m', retlab=True)[1]
    quadent = entry_qL = cs_Dblentry(lf_transfer, 0, 2, .2)
    driftlab = cs_label(lf_transfer, 1, 1, 'Driftlength / m', retlab=True)[1]
    driftent = entry_dlen = cs_Dblentry(lf_transfer, 1, 2, 2.3)
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
    entry_energy = cs_Dblentry(lf_quadrupole, 0, 1, 1722)
    cs_label(lf_quadrupole, 1, 1, 'Quadrupole strength k / (1/m'+uc.squared+')')
    cs_label(lf_quadrupole, 2, 0, 'Initial (k<0 = axial focus)')
    entry_ki = cs_Dblentry(lf_quadrupole, 2, 1, 2)
    cs_label(lf_quadrupole, 3, 0, 'Final')
    entry_kf = cs_Dblentry(lf_quadrupole, 3, 1, 6.5)
    cs_label(lf_quadrupole, 4, 0, 'steps')
    entry_points = cs_Intentry(lf_quadrupole, 4, 1, 1000)
    return
