# -*- coding: utf-8 -*-
''' accpy gui starter
author:     felix.kramer(at)physik.hu-berlin.de
'''
try:
    from Tkinter import Tk, mainloop
except:
    from tkinter import Tk, mainloop
from accpy.gui.mainwin import mainwindow
import gc
# to check for all used imports:
#   grep --include=\*.py -hrw '.' -e 'import'


if __name__ == '__main__':
    root = Tk()  # create window
    version = 0.5
    mainwindow(root, version)
    gc.enable()
    mainloop()   # start Tk mainloop
