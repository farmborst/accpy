# -*- coding: utf-8 -*-
''' accpy.lattices.reader
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as tk
    import ttk
except:
    import tkinter as tk
    import tkinter.ttk as ttk
from matplotlib import use
use('TkAgg')
from matplotlib.backends.backend_tkagg import (FigureCanvasTkAgg,
                                               NavigationToolbar2TkAgg)
from matplotlib.pyplot import close
from .layout import (cs_tabbar, cs_label, cs_Intentry, cs_Dblentry, cs_button,
                     cs_dropd, cs_Strentry)
from ..lattices.reader import lattlist, latt2py, latt2txt, txt2latt, txt2py
from ..visualize.lattice import latticeplot


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
        canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)
        toolbar.pack()
        canvas.draw()


def latticeeditor(frame, w, h):
    def _save():
        txt = editor.get('1.0', tk.END)
        txt2latt(txt, name.get(), closedmenuval.get())

    def _show():
        txt = editor.get('1.0', tk.END)
        _, _, _, UC, diagnostics, N_UC, _, _ = txt2py(txt, closedmenuval.get())
        fig = latticeplot(UC, diagnostics)
        showfigs([fig], tabs[1:])

    def _check(*args):
        lattice = latticemenu.get()
        editor.config(state=tk.NORMAL)
        if lattice in closedlatts:
            txt = latt2txt(lattice, True)
            _, _, _, UC, diagnostics, N_UC, _, _ = latt2py(lattice, True)
            fig = latticeplot(UC, diagnostics)
        elif lattice in openlatts:
            txt = latt2txt(lattice, False)
            _, _, _, UC, diagnostics, N_UC, _, _, _ = latt2py(lattice, False)
            fig = latticeplot(UC, diagnostics)
        if lattice == 'NEW LATTICE':
            txt = ''
            namelab.grid()
            name.grid()
            menulab.grid()
            closedmenu.grid()
            showbutton.grid()
            savebutton.grid()
            editor.delete('1.0', tk.END)
            editor.insert('1.0', txt)
            clearfigs(tabs[1:])
        else:
            showfigs([fig], tabs[1:])
            editor.delete('1.0', tk.END)
            editor.insert('1.0', txt)
            editor.config(state=tk.DISABLED)
            closedmenu.grid_remove()
            showbutton.grid_remove()
            savebutton.grid_remove()


    tabs = cs_tabbar(frame, w, h, ['Editor', '1d Lattice View'])

    lf = ttk.Frame(tabs[0])
    lf.pack(side=tk.LEFT, fill=tk.BOTH, expand=False)
    rf = ttk.Frame(tabs[0])
    rf.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True)

    closedlatts, openlatts = lattlist()
    cs_label(lf, 1, 1, 'Lattice', sticky=tk.W)
    latticemenu = cs_dropd(lf, 2, 1, closedlatts + openlatts + ['NEW LATTICE'],
                           action=_check)

    scrollbar = tk.Scrollbar(orient="vertical")
    editor = tk.Text(rf, xscrollcommand=scrollbar.set)
    editor.pack(fill=tk.BOTH, expand=True)

    namelab = cs_label(lf, 3, 1, 'Name', sticky=tk.W, retlab=True)[1]
    name = cs_Strentry(lf, 4, 1, '', sticky=tk.W)
    menulab = cs_label(lf, 5, 1, 'Type', sticky=tk.W, retlab=True)[1]
    closedmenuval, closedmenu = cs_dropd(lf, 6, 1, ['open', 'closed'], sticky=tk.W+tk.E, retbut=True)
    showbutton = cs_button(lf, 7, 1, 'Show', _show, sticky=tk.W+tk.E)
    savebutton = cs_button(lf, 8, 1, 'Save', _save, sticky=tk.W+tk.E)
    namelab.grid_remove()
    menulab.grid_remove()
    name.grid_remove()
    closedmenu.grid_remove()
    showbutton.grid_remove()
    savebutton.grid_remove()
    return


def settings(w, h):
    return
