# -*- coding: utf-8 -*-
''' accpy.gui.simulate
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import close
from threading import Thread
from multiprocessing import cpu_count
from time import time
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
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
        latt = lattice.get()
        slic = int(entry_slice.get())
        mode = 'trackbeta'
        particles = 1
        rounds = 1
        runthread(status, tabs, lsd,
                  (latt, slic, mode, particles, rounds))

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
    def _start():
        latt = lattice.get()
        slic = int(entry_slice.get())
        mode = 'trackpart'
        prts = int(entry_parts.get())
        rnds = int(entry_round.get())
        runthread(status, tabs, lsd,
                  (latt, slic, mode, prts, rnds))

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

    cs_label(tabs[0], 5, 1, 'Cavity peak Voltages / kV')
    cs_label(tabs[0], 5, 3, 'Emittance @ injection')
    cs_label(tabs[0], 6, 2, 'Radial')
    cs_label(tabs[0], 7, 2, 'Axial')
    cs_label(tabs[0], 8, 2, 'Longitudinal')
    cs_label(tabs[0], 6, 4, 'nm '+uc.pi+' rad')
    cs_label(tabs[0], 7, 4, 'nm '+uc.pi+' rad')
    cs_label(tabs[0], 8, 4, uc.ppt)
    entry_V_HF = cs_Strentry(tabs[0], 6, 1, '200 500 2000')
    entry_emitx = cs_Strentry(tabs[0], 6, 3, '100 200 300')
    entry_emity = cs_Strentry(tabs[0], 7, 3, '100 200 300')
    entry_emits = cs_Strentry(tabs[0], 8, 3, '.5 1 1.5')

    cs_button(tabs[0], 9, 6, 'Start', _start)
    status = cs_label(tabs[0], 9, 7, '')
    return


def gui_quadscansim(frame, w, h):
    def _start():
        ki = float(entry_ki.get())
        kf = float(entry_kf.get())
        qL = float(entry_qL.get())
        points = int(entry_points.get())
        driftlength = float(entry_dlen.get())
        epsx = float(entry_epsx.get())/1e9
        betx = float(entry_betx.get())
        alpx = float(entry_alpx.get())
        epsy = float(entry_epsy.get())/1e9
        bety = float(entry_bety.get())
        alpy = float(entry_alpy.get())
        epss = float(entry_epss.get())/1e3
        Dx = float(entry_Dx.get())
        Dpx = float(entry_Dpx.get())
        energy = float(entry_energy.get())*1e6
        particle = 'electron'
        runthread(status, tabs, simulate_quadscan,
                  (ki, kf, qL, driftlength, points, epsx, betx, alpx,
                   epsy, bety, alpy, epss, Dx, Dpx, energy, particle))

    tabs = cs_tabbar(frame, w, h, ['Menu', 'Beamextent'])

    cs_label(tabs[0], 1, 2, uc.epsilon+' / nm rad')
    cs_label(tabs[0], 1, 3, uc.beta+' / m')
    cs_label(tabs[0], 1, 4, uc.alpha+'/ rad')
    cs_label(tabs[0], 2, 1, 'Radial')
    entry_epsx = cs_Dblentry(tabs[0], 2, 2, 104.17)
    entry_betx = cs_Dblentry(tabs[0], 2, 3, 4.77)
    entry_alpx = cs_Dblentry(tabs[0], 2, 4, -2.18)
    cs_label(tabs[0], 3, 1, 'Axial')
    entry_epsy = cs_Dblentry(tabs[0], 3, 2, 21.11)
    entry_bety = cs_Dblentry(tabs[0], 3, 3, 1.29)
    entry_alpy = cs_Dblentry(tabs[0], 3, 4, 1.27)
    cs_label(tabs[0], 4, 2, uc.delta+' / '+uc.ppt)
    cs_label(tabs[0], 4, 3, '1.8 / m')
    cs_label(tabs[0], 4, 4, 'D\' / rad')
    cs_label(tabs[0], 5, 1, 'Longitudinal')
    entry_epss = cs_Dblentry(tabs[0], 5, 2, 0.547)
    entry_Dx = cs_Dblentry(tabs[0], 5, 3, 1.8)
    entry_Dpx = cs_Dblentry(tabs[0], 5, 4, 0.87)

    cs_label(tabs[0], 6, 1, 'Driftlength / m')
    entry_dlen = cs_Dblentry(tabs[0], 7, 1, 2.3)
    cs_label(tabs[0], 8, 1, 'Energy / MeV')
    entry_energy = cs_Dblentry(tabs[0], 9, 1, 1722)

    cs_label(tabs[0], 6, 3, 'Quadrupole strength k / (1/m'+uc.squared+')')
    cs_label(tabs[0], 7, 2, 'Initial (k<0 = axial focus)')
    entry_ki = cs_Dblentry(tabs[0], 7, 3, -8)
    cs_label(tabs[0], 8, 2, 'Final')
    entry_kf = cs_Dblentry(tabs[0], 8, 3, 8)
    cs_label(tabs[0], 9, 2, 'steps')
    entry_points = cs_Intentry(tabs[0], 9, 3, 1000)

    cs_label(tabs[0], 6, 4, 'Quadrupole length / m')
    entry_qL = cs_Dblentry(tabs[0], 7, 4, .2)

    cs_button(tabs[0], 10, 10, 'Start', _start)
    status = cs_label(tabs[0], 10, 11, '')
    return
