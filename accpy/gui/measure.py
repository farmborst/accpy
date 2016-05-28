# -*- coding: utf-8 -*-
''' accpy.gui.measure
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


def mesmenu(root):
    root.wm_title("accpy gui - measurements")
    nb = ttk.Notebook(root)
