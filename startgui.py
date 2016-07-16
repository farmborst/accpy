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
    from Tkinter import Tk, mainloop
except:
    from tkinter import Tk, mainloop
from accpy.gui.mainwin import mainwindow
from gc import enable


if __name__ == '__main__':
    root = Tk()  # create window
    version = 0.5
    mainwindow(root, version)
    enable()
    mainloop()   # start Tk mainloop
