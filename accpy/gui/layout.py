# -*- coding: utf-8 -*-
''' accpy.gui.layout
author:     felix.kramer(at)physik.hu-berlin.de
'''
try:
    import Tkinter as Tk
    import ttk
except:
    import tkinter as Tk
    import tkinter.ttk as ttk


# create and name tabs in root
def tabbar(root, names, w, h):
    nb = ttk.Notebook(root, width=w, height=h)
    tabs = [ttk.Frame(nb) for i in range(len(names))]  # 5 tabs
    [nb.add(tabs[i], text=name) for i, name in enumerate(names)]
    nb.pack(expand=True)
    return nb, tabs


# create and set tkinter stringvar
def cs_str(name):
    svar = Tk.StringVar()
    svar.set(name)
    return svar


# create and set tkinter intvar
def cs_int(value):
    ivar = Tk.IntVar()
    ivar.set(int(value))
    return ivar


# create, position and set tkinter label
def cs_lab(root, name, r, c):
    label = Tk.Label(master=root, textvariable=cs_str(name))
    label.grid(row=r, column=c)