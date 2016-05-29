# -*- coding: utf-8 -*-
''' accpy gui starter
author:     felix.kramer(at)physik.hu-berlin.de
'''
try:
    import Tkinter as Tk
except:
    import tkinter as Tk
from accpy.gui.mainwin import mainwindow

version = 0.2
root = Tk.Tk()  # create window
mainwindow(root, version)   # load toplevel menu
Tk.mainloop()   # start Tk mainloop
