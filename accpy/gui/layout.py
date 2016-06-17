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
def cs_tabbar(root, w, h, names):
    nb = ttk.Notebook(root, width=w, height=h)
    tabs = [ttk.Frame(nb) for i in range(len(names))]  # 5 tabs
    [nb.add(tabs[i], text=name) for i, name in enumerate(names)]
    nb.pack()
    return tabs


# create and set tkinter stringvar
def cs_str(name):
    svar = Tk.StringVar()
    svar.set(name)
    return svar


# create and set tkinter IntVar
def cs_int(value):
    ivar = Tk.IntVar()
    ivar.set(int(value))
    return ivar


# create and set tkinter DoubleVar
def cs_dbl(value):
    dvar = Tk.DoubleVar()
    dvar.set(float(value))
    return dvar


# create, position and set tkinter label
def cs_label(root, r, c, name):
    labelstr = cs_str(name)
    label = Tk.Label(master=root, textvariable=labelstr)
    label.grid(row=r, column=c)
    return labelstr


# create, position and set IntVar label
def cs_Intentry(root, r, c, value):
    entryint = cs_int(value)
    entry = Tk.Entry(root, textvariable=entryint)
    entry.grid(row=r, column=c)
    return entry


# create, position and set DoubleVar label
def cs_Dblentry(root, r, c, value):
    entryint = cs_dbl(value)
    entry = Tk.Entry(root, textvariable=entryint)
    entry.grid(row=r, column=c)
    return entry


# create, position and set button
def cs_button(root, r, c, label, action):
    button = Tk.Button(master=root, text=label, command=action)
    button.grid(row=r, column=c)
    return button


def cs_dropd(root, r, c, options):
    startvalue = cs_str(options[0])
    dropdown = Tk.OptionMenu(root, startvalue, *options)
    dropdown.grid(row=r, column=c)
    return startvalue
