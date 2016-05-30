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

    bar = Tk.Menu(root)

    # FILE MENU
    FM = Tk.Menu(bar, tearoff=0)
    FML = ['Save Figure',
           'Quit']

    def File_savefig():
        notready()
    FM.add_command(label=FML[0], command=File_savefig)

    def File_quit():
        root.quit()     # stops mainloop
        root.destroy()
    FM.add_command(label=FML[1], command=File_quit)
    bar.add_cascade(label="File", menu=FM)

    # SIMULATION MENU
    SM = Tk.Menu(bar, tearoff=0)
    SML = ['Betamatrix and dispersion tracking',
           'Particle tracking',
           'Dynamic energy and magnetic Flux',
           'Dynamic emittances',
           'Quadrupole scan']

    def Simu_twisstrack():
        twisstrack(root)
    SM.add_command(label=SML[0], command=Simu_twisstrack)

    def Simu_parttrack():
        parttrack(root)
    SM.add_command(label=SML[1], command=Simu_parttrack)

    def Simu_EBdynamics():
        notready()
    SM.add_command(label=SML[2], command=Simu_EBdynamics)

    def Simu_EmitDynamics():
        notready()
    SM.add_command(label=SML[3], command=Simu_EmitDynamics)

    def Simu_QuadScan():
        notready()
    SM.add_command(label=SML[4], command=Simu_QuadScan)
    bar.add_cascade(label="Simulation", menu=SM)

    # MEASUREMENT MENU
    MM = Tk.Menu(bar, tearoff=0)
    MML = ['Tunes',
           'Chromaticity',
           'Quadrupole scan',
           'Achromat scan']

    def Meas_tunes():
        notready()
    MM.add_command(label=MML[0], command=Meas_tunes)

    def Meas_chrom():
        notready()
    MM.add_command(label=MML[1], command=Meas_chrom)

    def Meas_QuadScan():
        notready()
    MM.add_command(label=MML[2], command=Meas_QuadScan)

    def Meas_AchroScan():
        notready()
    MM.add_command(label=MML[3], command=Meas_AchroScan)
    bar.add_cascade(label='Measurement', menu=MM)

    # Entries for Optimization menu
    OM = Tk.Menu(bar, tearoff=0)
    OML = ['Find Transverse emittance exchange section',
           'Twiss matching']

    def Opti_emittex():
        notready()
    OM.add_command(label=OML[0], command=Opti_emittex)

    def Opti_twissmatch():
        notready()
    OM.add_command(label=OML[1], command=Opti_twissmatch)
    bar.add_cascade(label="Optimization", menu=OM)

    # Entries for Help menu
    HM = Tk.Menu(bar, tearoff=0)
    HML = ['Documentation',
           'About accpy...']

    def Help_docs():
        notready()
    HM.add_command(label=HML[0], command=Help_docs)

    def Help_about():
        r = Tk.Tk()
        r.wm_title('About accpy')
        r.option_add('*font', 'Helvetica -16 bold')
        txt1 = Tk.Label(r, text='accpy version {}'.format(version))
        r.option_add('*font', 'Helvetica -14')
        txt2 = Tk.Label(r, text='Created by\nfelix kramer\n(felix.kramer@physik.hu-berlin.de)')
        txt1.pack()
        txt2.pack()
    HM.add_command(label=HML[1], command=Help_about)
    bar.add_cascade(label="Help", menu=HM)
    return bar
