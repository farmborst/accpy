# -*- coding: utf-8 -*-
''' accpy.gui.menu
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    import Tkinter as Tk
    from tkMessageBox import showinfo
except:
    import tkinter as Tk
    from tkinter.messagebox import showinfo
from .simulate import twisstrack, parttrack


def mainwindow(root, version):
    w, h = root.winfo_screenwidth(), root.winfo_screenheight()
    root.geometry('{}x{}'.format(w, h))
    root.wm_title("accpy gui {}".format(version))
    bar = menubar(root, version)
    root.config(menu=bar)


def menubar(root, version):
    def notready():
        showinfo('Ooops!', 'Sorry, but this feature is not ready yet...')

    # Entries for File menu
    def File_savefig():
        notready()

    def File_quit():
        root.quit()     # stops mainloop
        root.destroy()

    # Entries for Simulation menu
    def Simu_twisstrack():
        twisstrack(root)

    def Simu_parttrack():
        parttrack(root)

    # Entries for Measurement menu
    def Meas_tunes():
        notready()

    # Entries for Optimization menu
    def Opti_emittex():
        notready()

    # Entries for Help menu
    def Help_docs():
        notready()

    def Help_about():
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
    filemenu.add_command(label='Save Figure',
                         command=File_savefig)
    filemenu.add_command(label='Quit',
                         command=File_quit)
    bar.add_cascade(label="File", menu=filemenu)

    simumenu = Tk.Menu(bar, tearoff=0)
    simumenu.add_command(label='Betamatrix and dispersion tracking',
                         command=Simu_twisstrack)
    simumenu.add_command(label='Particle tracking',
                         command=Simu_parttrack)
    simumenu.add_command(label='Dynamic energy and magnetic Flux',
                         command=Simu_parttrack)
    simumenu.add_command(label='Dynamic emittances',
                         command=Simu_parttrack)
    simumenu.add_command(label='Quadrupole scan',
                         command=Simu_parttrack)
    bar.add_cascade(label="Simulation", menu=simumenu)

    measmenu = Tk.Menu(bar, tearoff=0)
    measmenu.add_command(label='Tunes',
                         command=Meas_tunes)
    measmenu.add_command(label='Chromaticity',
                         command=Meas_tunes)
    measmenu.add_command(label='Quadrupole scan',
                         command=Meas_tunes)
    measmenu.add_command(label='Achromat scan',
                         command=Meas_tunes)
    bar.add_cascade(label="Measurement",
                    menu=measmenu)

    measmenu = Tk.Menu(bar, tearoff=0)
    measmenu.add_command(label='Find Transverse emittance exchange section',
                         command=Opti_emittex)
    measmenu.add_command(label='Twiss matching',
                         command=Opti_emittex)
    bar.add_cascade(label="Optimization", menu=measmenu)

    helpmenu = Tk.Menu(bar, tearoff=0)
    helpmenu.add_command(label='Documentation',
                         command=Help_docs)
    helpmenu.add_command(label='About accpy...',
                         command=Help_about)
    bar.add_cascade(label="Help", menu=helpmenu)
    return bar
