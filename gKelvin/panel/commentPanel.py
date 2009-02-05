"""The panel to specify whether to insert comments in file or not"""

import wx
from utility import util

class CommentPanel(wx.Panel):
    """The panel to specify whether to insert comments or not"""

    mainLbl = 'Fully annotated config file'
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all the widgets for the panel"""

        # just a checkbox
        self.checkbox = wx.CheckBox(self, -1, label=self.mainLbl)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        # only one thing in panel, don't need to arrange
        pass

    def initialize(self):
        """Initialize the checkbox based on the file contents"""

        # the initial value depends on what the filemanager says
        self.checkbox.SetValue(self.fm.insertComments)

        # bind events
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)

    def onCheckBox(self, event):
        """Event when the checkbox is checked"""

        if self.checkbox.GetValue():
            self.fm.insertComments = True
        else:
            self.fm.insertComments = False
        util.callAncestorFunction(self, 'hasChanged')
