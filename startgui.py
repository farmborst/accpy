#!/usr/bin/python
# -*- coding: utf-8 -*-
''' accpy gui starter
author:     felix.kramer(at)physik.hu-berlin.de

grep checking:
    n: show line number
    h: hide file name
    r: recursively check all subfolders
    w: whole word match
    x: whole line match
    grep --include=\*.py -hrw '.' -e 'import'
'''

try:
   from Tkinter import Tk, PhotoImage
except:
   from tkinter import Tk, PhotoImage
from accpy.gui.main import MainApp
from gc import enable
import os


if __name__ == '__main__':
    root = Tk()
    version = 0.6
    cwd = os.getcwd()
    app = MainApp(root, version, cwd)
    enable()
    root.mainloop()
