"""The panel to specify whether to use on marker trait loci only"""

import wx
from config import KelvinInfo
from utility import util

class TM_Panel(wx.Panel):
    """The panel to whether to use on marker trait loci only"""

    mainLbl = 'Include markers in tested trait loci'
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        temp = KelvinInfo.TMEntry
        self.tag = temp[0]

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all the widgets for the panel"""

        # just a checkbox
        self.checkbox = wx.CheckBox(self, -1, label=self.mainLbl)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        # only one widget in panel, don't need to arrange
        pass

    def initialize(self):
        """Initialize the checkbox based on the file contents"""

        if self.fm.isPresent(self.tag):
            self.checkbox.SetValue(True)
        else:
            self.checkbox.SetValue(False)

        self.onAnalysisChange(None)

        # bind events
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)

        # this panel needs to know when the analysis type changes
        # go up the parent tree until the right function is found
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)

    def onCheckBox(self, event):
        """Event when the checkbox is checked"""

        if self.checkbox.GetValue():
            # want to use TM, insert it into file
            if not self.fm.isPresent(self.tag):
                self.fm.insert(self.tag)
        else:
            # do not want to use TM tag, take it out of the file
            if self.fm.isPresent(self.tag):
                lines = self.fm.getLines(self.tag)
                for eachLine in lines:
                    self.fm.remove(eachLine)

        # need this check so it won't be called when initializing the panel
        if event is not None:
            util.callAncestorFunction(self, 'hasChanged')

    def onAnalysisChange(self, event):
        """Event when the analysis type (two point or multi point) changes"""

        # this checkbox is disabled for two point
        if self.fm.isPresent('TP'):
            self.checkbox.SetValue(False)
            self.checkbox.Disable()
            self.onCheckBox(None)
        else:
            self.checkbox.Enable(True)
