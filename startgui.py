# -*- coding: utf-8 -*-
''' accpy gui starter
author:     felix.kramer(at)physik.hu-berlin.de
'''
try:
    import Tkinter as Tk
except:
    import tkinter as Tk
from accpy.gui.mainwin import mainwindow
#from multiprocessing import Process
import gc


if __name__ == '__main__':
    root = Tk.Tk()  # create window
    version = 0.2
    mainwindow(root, version)
    gc.enable()
    Tk.mainloop()   # start Tk mainloop
