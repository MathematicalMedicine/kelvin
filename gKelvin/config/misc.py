"""Miscellaneous program settings"""

import os
import wx

# the title of the main application window
mainWindowTitle = 'gKelvin'

# how large the window is
#mainWindowSize = (1200, 768)
mainWindowSize = (1024, 768)

# starting position, with (0, 0) being the upper left corner
mainWindowPosition = wx.DefaultPosition

# file specified to open at startup
#openFile = '../sampledata/test.conf'
#openFile = '../sampledata/two.conf'

# what tab is focused when a new kelvin config file is opened
firstTabFocus = 0

# the size of buttons that are used to add and remove constraints and values
addRemoveButtonSize = (30, -1)

# enable or disable any kind of saving
savingEnabled = True

# color for panels that have some error
errorColor = (255, 120, 120)

# the path to the folder containing the templates used for
# the different types of new kelvin config files
# it is assumed to be in a folder called template that is
# one level above where this file lies
newConfigFilePath = os.path.dirname(__file__)
newConfigFilePath = os.path.join(newConfigFilePath, os.path.pardir)
newConfigFilePath = os.path.join(newConfigFilePath, 'template')
