# -*- coding: utf-8 -*-
''' accpy.gui.layout
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as Tk
    import ttk
    from tkFileDialog import askopenfilename
    from tkMessageBox import showerror
except:
    import tkinter as Tk
    import tkinter.ttk as ttk
    from tkinter.filedialog import askopenfilename
    from tkinter.messagebox import showerror
from matplotlib import use
use('TkAgg')
from threading import Thread


def tabbar(root, names, w, h):
    nb = ttk.Notebook(root, width=w, height=h)
    tabs = [ttk.Frame(nb) for i in range(len(names))]  # 5 tabs
    [nb.add(tabs[i], text=name) for i, name in enumerate(names)]
    nb.pack(expand=True)
    return tabs