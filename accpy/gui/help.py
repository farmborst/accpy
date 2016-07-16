# -*- coding: utf-8 -*-
''' accpy.lattices.reader
author:     felix.kramer(at)physik.hu-berlin.de
'''
from __future__ import division
try:
    from Tkinter import (Tk, Label, Listbox, Text, LEFT, RIGHT, BOTH, END,
                         NORMAL, DISABLED, Scrollbar)
    from ttk import Frame
except:
    from tkinter import (Tk, Label, Listbox, Text, LEFT, RIGHT, BOTH, END,
                         NORMAL, DISABLED, Scrollbar)
    from tkinter.ttk import Frame


def documentation(version, w, h):
    def _setdoc(evt):
        w = evt.widget
        index = int(w.curselection()[0])
        doc = w.get(index)
        textfield.config(state=NORMAL)
        textfield.delete('1.0', END)
        textfield.insert('1.0', docdic[doc])
        textfield.config(state=DISABLED)
    r = Tk()
    w = int(w/2)
    h = int(h/2)
    r.geometry('{}x{}+{}+{}'.format(w, h, int(w/2), int(h/2)))
    r.wm_title('Documentation accpy version {}'.format(version))

    lf = Frame(r)
    lf.pack(side=LEFT, fill=BOTH, expand=False)
    rf = Frame(r)
    rf.pack(side=RIGHT, fill=BOTH, expand=True)

    docmenuopts = ['General',
                   'Lattice editor']
    docmenu = Listbox(lf)
    for entry in docmenuopts:
        docmenu.insert(END, entry)
    docmenu.pack(fill=BOTH, expand=True)
    docmenu.bind('<<ListboxSelect>>', _setdoc)

    scrollbar = Scrollbar(orient="vertical")
    textfield = Text(rf, xscrollcommand=scrollbar.set)
    textfield.pack(fill=BOTH, expand=True)


def about(version, w, h):
    r = Tk()
    r.wm_title('About')
    r.option_add('*font', 'Helvetica -16 bold')
    txt1 = Label(r, text='accpy version {}'.format(version))
    r.option_add('*font', 'Helvetica -14')
    txt2 = Label(r, text='Coded by\nfelix kramer\n(felix.kramer@physik.hu-berlin.de)')
    txt1.pack()
    txt2.pack()


doc_general = '''
ACCPY
=========

#Python module for acclerator physics providing
    - linear optics particle and twiss tracking
    - energy ramp simulation for closed lattices with:
        - energy loss
        - required acceleration voltage
        - synchrotron phase and frequency
        - equilibrium emittances and bunch lengths
        - dynamic emittances and bunch lenghths
    - quadrupole scan simulation

#Requirements
    python
    matplotlib
    numpy
    scipy
    h5py

#Tested with
    python          2.7.9
    numpy           1.8.2
    scipy           0.14.0
    matplotlib      1.4.2
    h5py            2.2.1

#Installation
    install python 2.7 and matplotlib, numpy, scipy and h5py
    (on debian 8: >> apt-get install python python-matplotlib python-numpy python-scipy python-h5py
    run >> python startgui.py

#Info
    see accpy.__doc__

#FAQ (Frequently asked Questions)
'''

doc_lattice = '''
all lattice files (.latt) can be found under ~/accpy/accpy/lattices/x
    with x = closed or open
here you can edit them with an editor of your choice
you have all the powers of python including arbitrary imports at hand
'''

docdic = {'General': doc_general,
          'Lattice editor': doc_lattice}
