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


def gridwidget(widget, r, c, sticky):
    if sticky is None:
        widget.grid(row=r, column=c)
    else:
        widget.grid(row=r, column=c, sticky=sticky)

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
def cs_label(root, r, c, name, sticky=None, retlab=False):
    labelstr = cs_str(name)
    label = Tk.Label(master=root, textvariable=labelstr)
    gridwidget(label, r, c, sticky)
    if retlab:
        return labelstr, label
    else:
        return labelstr



# create, position and set IntVar label
def cs_Intentry(root, r, c, value, sticky=None):
    entryint = cs_int(value)
    entry = Tk.Entry(root, textvariable=entryint)
    gridwidget(entry, r, c, sticky)
    return entry


# create, position and set DoubleVar label
def cs_Dblentry(root, r, c, value, sticky=None):
    entryint = cs_dbl(value)
    entry = Tk.Entry(root, textvariable=entryint)
    gridwidget(entry, r, c, sticky)
    return entry


# create, position and set StringVar label
def cs_Strentry(root, r, c, value, sticky=None):
    entrystr = cs_str(value)
    entry = Tk.Entry(root, textvariable=entrystr)
    gridwidget(entry, r, c, sticky)
    return entry


# create, position and set button
def cs_button(root, r, c, label, action, sticky=None):
    button = Tk.Button(master=root, text=label, command=action)
    gridwidget(button, r, c, sticky)
    return button


def cs_dropd(root, r, c, options, action=None, sticky=None, retbut=False):
    startvalue = cs_str('')  # cs_str(options[0])
    dropdown = Tk.OptionMenu(root, startvalue, *options)
    gridwidget(dropdown, r, c, sticky)
    if action is not None:
        startvalue.trace('w', action)
    if retbut:
        return startvalue, dropdown
    return startvalue
