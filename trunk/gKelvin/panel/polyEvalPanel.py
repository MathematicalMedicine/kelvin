"""The panel to specify whether to use the polynomial evaluation or not"""

import wx
from config import KelvinInfo
from utility import util

class PolyEvalPanel(wx.Panel):
    """The panel to specify whether to use the polynomial evaluation or not"""

    mainLbl = 'Use polynomial evaluation'
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        temp = KelvinInfo.polynomialEvaluationEntry
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

        # only one thing in panel, don't need to arrange
        pass
        
    def initialize(self):
        """Initialize the checkbox based on the file contents"""

        if self.fm.isPresent(self.tag):
            self.checkbox.SetValue(True)
        else:
            self.checkbox.SetValue(False)

        # bind events
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)

    def onCheckBox(self, event):
        """Event when the checkbox is checked"""

        if self.checkbox.GetValue():
            # want to use polynomial, insert it in file
            if not self.fm.isPresent(self.tag):
                self.fm.insert(self.tag)
        else:
            # do not want to use polynomial, take it out of the file
            if self.fm.isPresent(self.tag):
                lines = self.fm.getLines(self.tag)
                for eachLine in lines:
                    self.fm.remove(eachLine)
        util.callAncestorFunction(self, 'hasChanged')
