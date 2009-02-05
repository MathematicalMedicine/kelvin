"""Miscellaneous program settings"""

import os
import os.path
import time
import wx
import images

def print_timing(func):
    # to use this function, you can use the @ decorator at the function
    # definition, for example:
    #
    #           @print_timing
    #           def f():
    #               ...
    #
    # this will print the timing every time the function is called.
    #
    # you can also just call the wrapper that is returned by this function
    # with the correct args, for example:
    #
    # f(foo, bar) -> print_timing(f)(foo, bar)

    def wrapper(*arg):
        t1 = time.time()
        res = func(*arg)
        t2 = time.time()
        print '%s took %0.3f s' % (func.func_name, (t2-t1))
        return res
    return wrapper

# the title of the main application window
mainWindowTitle = 'graphApp'

# how large the window is
mainWindowSize = (1024, 768)
#mainWindowSize = (800, 800)

# starting position, with (0, 0) being the upper left corner
mainWindowPosition = wx.DefaultPosition

# the size of the file format dialog box
fileformatDialogSize = (800, 600)

# file that has command line options to perform at startup
# used for debugging, and so i don't have to go through menus every time
#openCmdFile = '../data/cmd.txt'
