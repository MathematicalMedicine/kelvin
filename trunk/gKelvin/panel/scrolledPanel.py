"""A scrolled panel with a bug fix"""

import wx
import wx.lib.scrolledpanel

from utility import util

class ScrolledPanel(wx.lib.scrolledpanel.ScrolledPanel):
    """A scrolled panel that won't jump around when pressing a child widget"""

    def __init__(self, parent):
        wx.lib.scrolledpanel.ScrolledPanel.__init__(self, parent, -1)

    def OnChildFocus(self, event):
        """In the base class, this event causes the window to focus on the 
        child, so if the child of a scrolled panel is a very large panel, then 
        the focus will jump up to the top, which is annoying.  So overwriting 
        this function to not do that.
        
        """

        event.Skip()

class TabPanel(ScrolledPanel):
    """A panel meant to be used as one of the tabs in a kelvin config file"""

    # this class has some added functions in order for it to deal with
    # having to relay information about its values changing

    def hasChanged(self):
        util.callAncestorFunction(self, 'hasChanged', self)
