"""The panel to specify liability classes"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class LiabilityClassPanel(wx.Panel):
    """The panel to specify liability classes"""

    mainLbl = 'Number of Liability Classes:'
    borderAmt = 5

    # the default value for liability classes
    default = KelvinInfo.defaultLiabilityClass

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        self.tag = KelvinInfo.liabilityClassEntry[0]

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create widgets for the panel"""

        # checkbox to see if specifying liability classes
        self.checkbox = wx.CheckBox(self, -1, label=self.mainLbl)

        # text box to specifiy number of liability classes
        self.textbox = wx.TextCtrl(self, -1)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        mainSizer.Add(self.checkbox, 0, style, self.borderAmt)
        mainSizer.Add(self.textbox, 0, style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widget values based on file contents"""

        if self.fm.isPresent(self.tag):
            line = self.fm.getLines(self.tag)[0]
            self.charge = line
            self.checkbox.SetValue(True)
            self.textbox.ChangeValue(line[1].str)
            self.textbox.charge = line[1]
        else:
            self.charge = None
            self.checkbox.SetValue(False)
            self.textbox.Disable()
            self.textbox.ChangeValue('1')
            self.textbox.charge = None

        # bind events
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)
        self.textbox.Bind(wx.EVT_TEXT, self.onTextEntry, self.textbox)

    def onCheckBox(self, event):
        """Event when the checkbox is clicked"""

        if self.checkbox.GetValue():
            # using liability classes
            if self.charge is None:
                self.charge = self.fm.insert(self.tag, self.default)
                self.textbox.charge = self.charge[1]

            # add the line back to the file
            self.fm.addLine(self.charge)

            self.textbox.Enable(True)
            self.textbox.ChangeValue(self.default)
        else:
            self.fm.remove(self.charge)
            self.textbox.ChangeValue(self.default)
            self.textbox.Disable()
            self.textbox.charge.update(self.textbox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')

        # others need to know!
        event.Skip()

    def onTextEntry(self, event):
        """Event when text is entered in the textbox"""

        self.textbox.charge.update(self.textbox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')

        # others need to know!
        event.Skip()
