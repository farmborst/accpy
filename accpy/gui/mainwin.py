# -*- coding: utf-8 -*-
''' accpy.gui.menu
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
from .simulate import simmenu
from .measure import mesmenu


def mainwindow(root, version):
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry('{}x{}'.format(w, h))
    root.wm_title("accpy gui")
    bar = menubar(root, version)
    # add some buttons to top level menu
#    title = Tk.Label(master=root, text='Welcome to accpy gui')
#    title.pack()
    root.config(menu=bar)


def menubar(root, version):
    class menu():
        @staticmethod
        def _quit():
            root.quit()     # stops mainloop
            root.destroy()

        @staticmethod
        def _simu():
            simmenu(root)

        @staticmethod
        def _meas():
            mesmenu(root)

        @staticmethod
        def _docs():
            r = Tk.Tk()
            r.geometry('{}x{}'.format(100, 100))
            r.wm_title('Documentatioin off accpy')
            r.option_add('*font', 'Helvetica -16 bold')
            txt1 = Tk.Label(r, text='...'.format(version))
            txt1.pack()

        @staticmethod
        def _abou():
            r = Tk.Tk()
            r.wm_title('About accpy')
            r.option_add('*font', 'Helvetica -16 bold')
            txt1 = Tk.Label(r, text='accpy version {}'.format(version))
            r.option_add('*font', 'Helvetica -14')
            txt2 = Tk.Label(r, text='Created by\nfelix kramer\n(felix.kramer@physik.hu-berlin.de)')
            txt1.pack()
            txt2.pack()
    bar = Tk.Menu(root)

    filemenu = Tk.Menu(bar, tearoff=0)
    filemenu.add_separator()
    filemenu.add_command(label='New Simulation', command=menu._simu)
    filemenu.add_command(label='New Measurement', command=menu._meas)
    filemenu.add_separator()
    filemenu.add_command(label='Quit', command=menu._quit)
    bar.add_cascade(label="File", menu=filemenu)

    helpmenu = Tk.Menu(bar, tearoff=0)
    helpmenu.add_command(label='Documentation', command=menu._docs)
    filemenu.add_separator()
    helpmenu.add_command(label='About accpy...', command=menu._abou)
    bar.add_cascade(label="Help", menu=helpmenu)
    return bar
