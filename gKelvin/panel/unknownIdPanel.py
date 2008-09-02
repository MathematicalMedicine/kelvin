"""The panel to select the id of unknown people"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class UnknownIdPanel(wx.Panel):
    """The panel to select the id of unknown people"""

    mainLbl = 'ID for unknown individual:'
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        temp = KelvinInfo.unknownPersonId
        self.tag = temp[0]
        self.default = temp[2]

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Create all the widgets for the panel"""

        # a label with a textbox
        self.mainLabel = wx.StaticText(self, -1, self.mainLbl)
        self.unknownBox = wx.TextCtrl(self, -1)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        mainSizer.Add(self.mainLabel, 0, style, self.borderAmt)
        mainSizer.Add(self.unknownBox, 0, style, self.borderAmt)
        self.SetSizer(mainSizer)

    def initialize(self):
        """Initialize the choice boxes and number based on the file contents"""

        # check if unknown id is specified
        lines = self.fm.getLines(self.tag)

        if len(lines) > 0:
            # the textbox is responsible for something in the file
            self.unknownBox.charge = lines[0][1]
            self.unknownBox.ChangeValue(lines[0][1].str)
        else:
            # textbox is not responsible for anything
            self.unknownBox.ChangeValue(str(self.default))
            self.unknownBox.charge = None

        # bind the textbox event
        self.unknownBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.unknownBox)

    def onTextEntry(self, event):
        """Event when characters are typed into the unknown textbox"""

        if not self.unknownBox.charge:
            # the line did not exist before, add it to the file
            newLine = self.fm.insert(self.tag, self.unknownBox.GetValue())
            self.unknownBox.charge = newLine[1]
        else:
            # line already existed, simply change the value
            self.unknownBox.charge.update(self.unknownBox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')
