# -*- coding: utf-8 -*-
''' accpy.lattices.reader
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    from Tkinter import (Text, LEFT, RIGHT, BOTH, END, NORMAL, DISABLED, N, E,
                         S, W, Scrollbar, LabelFrame, TOP, BOTTOM, RAISED,
                         Toplevel)
    from ttk import Frame
except:
    from tkinter import (Text, LEFT, RIGHT, BOTH, END, NORMAL, DISABLED, N, E,
                         S, W, Scrollbar, LabelFrame, TOP, BOTTOM, RAISED,
                         Toplevel)
    from tkinter.ttk import Frame
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import close
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry, cp_button, cs_checkbox)
from ..lattices.reader import lattlist, latt2py, latt2txt, txt2latt, txt2py
from ..visualize.lattice import latticeplot
from ..visualize.figures import plotstandards
from ..dataio.hdf5 import confload, confsave


def clearfigs(tabs):
    close('all')
    for tab in tabs:
        # destroy all widgets in fram/tab and close all figures
        for widget in tab.winfo_children():
            widget.destroy()


def showfigs(figs, tabs):
    clearfigs(tabs)
    for fig, tab in zip(figs, tabs):
        canvas = FigureCanvasTkAgg(fig, master=tab)
        toolbar = NavigationToolbar2TkAgg(canvas, tab)
        canvas.get_tk_widget().pack(fill=BOTH, expand=True)
        toolbar.pack()
        canvas.draw()


def latticeeditor(frame, w, h):
    def _save():
        txt = editor.get('1.0', END)
        txt2latt(txt, name.get(), closedmenuval.get())

    def _show():
        txt = editor.get('1.0', END)
        _, _, _, UC, diagnostics, N_UC, _, _ = txt2py(txt, closedmenuval.get())
        fig = latticeplot(UC, diagnostics)
        showfigs([fig], tabs[1:])

    def _check(*args):
        lattice = latticemenu.get()
        editor.config(state=NORMAL)
        if lattice in closedlatts:
            txt = latt2txt(lattice, True)
            _, _, _, UC, diagnostics, N_UC, _, _ = latt2py(lattice, True)
            fig = latticeplot(UC, diagnostics)
        elif lattice in openlatts:
            txt = latt2txt(lattice, False)
            _, _, _, UC, diagnostics, N_UC, _, _, _, _, _, _ = latt2py(lattice, False)
            fig = latticeplot(UC, diagnostics)
        if lattice == 'NEW LATTICE':
            txt = ''
            namelab.grid()
            name.grid()
            menulab.grid()
            closedmenu.grid()
            showbutton.grid()
            savebutton.grid()
            editor.delete('1.0', END)
            editor.insert('1.0', txt)
            clearfigs(tabs[1:])
        else:
            showfigs([fig], tabs[1:])
            editor.delete('1.0', END)
            editor.insert('1.0', txt)
            editor.config(state=DISABLED)
            closedmenu.grid_remove()
            showbutton.grid_remove()
            savebutton.grid_remove()


    tabs = cs_tabbar(frame, w, h, ['Editor', '1d Lattice View'])

    lf = Frame(tabs[0])
    lf.pack(side=LEFT, fill=BOTH, expand=False)
    rf = Frame(tabs[0])
    rf.pack(side=RIGHT, fill=BOTH, expand=True)

    closedlatts, openlatts = lattlist()
    cs_label(lf, 1, 1, 'Lattice', sticky=W)
    latticemenu = cs_dropd(lf, 2, 1, closedlatts + openlatts + ['NEW LATTICE'],
                           action=_check)

    scrollbar = Scrollbar(orient="vertical")
    editor = Text(rf, xscrollcommand=scrollbar.set)
    editor.pack(fill=BOTH, expand=True)

    namelab = cs_label(lf, 3, 1, 'Name', sticky=W, retlab=True)[1]
    name = cs_Strentry(lf, 4, 1, '', sticky=W)
    menulab = cs_label(lf, 5, 1, 'Type', sticky=W, retlab=True)[1]
    closedmenuval, closedmenu = cs_dropd(lf, 6, 1, ['open', 'closed'], sticky=W+E, retbut=True)
    showbutton = cs_button(lf, 7, 1, 'Show', _show, sticky=W+E)
    savebutton = cs_button(lf, 8, 1, 'Save', _save, sticky=W+E)
    namelab.grid_remove()
    menulab.grid_remove()
    name.grid_remove()
    closedmenu.grid_remove()
    showbutton.grid_remove()
    savebutton.grid_remove()
    return


def settings(w, h):
    def _cancel():
        r.destroy()

    def _save():
        vallist[0] = float(entry_width.get())
        vallist[1] = float(entry_height.get())
        vallist[2] = gridon.get()
        vallist[3] = float(entry_dpi.get())
        vallist[4] = dropd_fontfamily.get()
        vallist[5] = dropd_fontsize.get()
        vallist[6] = dropd_markersize.get()
        vallist[7] = dropd_linewidth.get()
        n = dropd_laxpowerlimits.get()
        m = dropd_uaxpowerlimits.get()
        vallist[8] = [n, m]
        confsave(confpath, varlist, vallist)
        plotstandards(varlist, vallist, w, h)
        r.destroy()

    confpath = './settings.conf'
    varlist, vallist = confload(confpath)

    r = Toplevel()
    w = int(w/2)
    h = int(h/2)
    r.geometry('{}x{}+{}+{}'.format(w, h, int(w/2), int(h/2)))
    r.wm_title('ACCPY Settings')

    uf = Frame(r, relief=RAISED)
    uf.pack(side=TOP, fill=BOTH, expand=True)
    lf = Frame(r)
    lf.pack(side=BOTTOM, fill=BOTH, expand=False)

    lf_figures = LabelFrame(uf, text="Figure Settings", padx=5, pady=5)
    lf_figures.grid(row=1, column=0, sticky=W+E+N+S, padx=10, pady=10)
    cs_label(lf_figures, 1, 0, 'Width / inches (25.4mm)', sticky=W)
    entry_width = cs_Dblentry(lf_figures, 2, 0, vallist[0], sticky=W+E)
    cs_label(lf_figures, 3, 0, 'Height / inches (25.4mm)', sticky=W)
    entry_height = cs_Dblentry(lf_figures, 4, 0, vallist[1], sticky=W+E)
    cs_label(lf_figures, 5, 0, 'DPI', sticky=W)
    entry_dpi = cs_Dblentry(lf_figures, 6, 0, vallist[3], sticky=W+E)
    gridon = cs_checkbox(lf_figures, 7, 0, 'Show grid', vallist[2], sticky=W)

    fontfamilys = ['serif', 'sans-serif', 'cursive', 'fantasy', 'monospace']
    cs_label(lf_figures, 8, 0, 'Fontfamily', sticky=W)
    dropd_fontfamily = cs_dropd(lf_figures, 9, 0, fontfamilys, sticky=W)
    dropd_fontfamily.set(vallist[4])
    fontsizes = range(1, 30)+range(30,80,2)
    cs_label(lf_figures, 10, 0, 'Main fontsize', sticky=W)
    dropd_fontsize = cs_dropd(lf_figures, 11, 0, fontsizes, sticky=W)
    dropd_fontsize.set(vallist[5])
    markersizes = range(1, 30)
    cs_label(lf_figures, 12, 0, 'Markersize', sticky=W)
    dropd_markersize = cs_dropd(lf_figures, 13, 0, markersizes, sticky=W)
    dropd_markersize.set(vallist[6])
    cs_label(lf_figures, 14, 0, 'Linewidth', sticky=W)
    dropd_linewidth = cs_dropd(lf_figures, 15, 0, markersizes, sticky=W)
    dropd_linewidth.set(vallist[7])
    cs_label(lf_figures, 16, 0, 'Scientific notation for ax labels, so:\n10^-n < data < 10^m', sticky=W)
    cs_label(lf_figures, 17, 0, 'n', sticky=W)
    dropd_laxpowerlimits = cs_dropd(lf_figures, 18, 0, range(-9, 1), sticky=W)
    dropd_laxpowerlimits.set(vallist[8][0])
    cs_label(lf_figures, 19, 0, 'm', sticky=W)
    dropd_uaxpowerlimits = cs_dropd(lf_figures, 20, 0, range(1, 9), sticky=W)
    dropd_uaxpowerlimits.set(vallist[8][1])


    cp_button(lf, 'Save and Close', _save, side=RIGHT, fill=BOTH, expand=True)
    cp_button(lf, 'Close', _cancel, side=RIGHT, fill=BOTH, expand=True)
    return


def defaults(confpath):
    varlist = ['fig_width',
               'fig_height',
               'grid',
               'dpi',
               'fontfamily',
               'fontsize',
               'markersize',
               'linewidth',
               'axformatterlimits']
    vallist = [19.2,
               10.8,
               1,
               100,
               'serif',
               28,
               3,
               2,
               [-2, 3]]
    confsave(confpath, varlist, vallist)
    return varlist, vallist
