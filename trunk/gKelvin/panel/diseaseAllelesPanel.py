"""The panel to specify the number of disease alleles"""

import wx
from config import KelvinInfo
from utility import util

class DiseaseAllelesPanel(wx.Panel):
    """The panel to specify the number of disease alleles"""

    mainLbl = 'Number of disease alleles:'
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        temp = KelvinInfo.diseaseAlleles
        self.tag = temp[0]
        self.default = temp[2]

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all the widgets for the panel"""

        # label
        self.text = wx.StaticText(self, -1, self.mainLbl)

        # text box for value
        self.textbox = wx.TextCtrl(self, -1)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        style = wx.ALIGN_CENTER_VERTICAL | wx.RIGHT
        mainSizer.Add(self.text, 0, style, self.borderAmt)
        mainSizer.Add(self.textbox, 0, style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize the choice boxes and number based on the file contents"""

        # check if number of disease alleles are specified
        lines = self.fm.getLines(self.tag)
        if len(lines) > 0:
            # tag is present
            self.textbox.charge = lines[0][1]
            self.textbox.ChangeValue(lines[0][1].str)
        else:
            self.textbox.charge = None
            self.textbox.ChangeValue(str(self.default))

        # since this feature isn't supported, disable it for now
        self.textbox.Disable()

        # bind events
        self.textbox.Bind(wx.EVT_TEXT, self.onTextEntry, self.textbox)

    def onTextEntry(self, event):
        """Event when characters are typed into the textbox"""

        # if the texbox was not responsible for anything, that means this
        # line was not originally in the file, so have to add it to the file
        if not self.textbox.charge:
            newLine = self.fm.insert(self.tag, self.textbox.GetValue())
            self.textbox.charge = newLine[1]
        else:
            self.textbox.charge.update(self.textbox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')
