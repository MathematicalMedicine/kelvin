"""The panel to specify what type of chromosome to use"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class ChromosomePanel(wx.Panel):
    """The panel to specify what type of chromosome to use"""

    mainLbl = 'Chromosome:'
    borderAmt = 5
    choiceList = ['(pseudo) Autosomal', 'X']

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.tag = KelvinInfo.chromosomeEntry[0]

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create widgets for the panel"""

        # text
        self.mainLabel = wx.StaticText(self, -1, self.mainLbl)

        # choice box to choose the type of chromosome
        self.choice = wx.Choice(self, -1, choices=self.choiceList)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)
        
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        mainSizer.Add(self.mainLabel, 0, style, self.borderAmt)
        mainSizer.Add(self.choice, 0, style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widget values based on file contents"""

        if self.fm.isPresent(self.tag):
            # use X chromosome
            self.choice.SetSelection(1)
        else:
            self.choice.SetSelection(0)

        # bind events
        self.choice.Bind(wx.EVT_CHOICE, self.onChoice, self.choice)

    def onChoice(self, event):
        """Event when the choice box is clicked"""
    
        if self.choice.GetSelection() == 0:
            # remove all occurances of the tag
            if self.fm.isPresent(self.tag):
                lines = self.fm.getLines(self.tag)
                for eachLine in lines:
                    self.fm.remove(eachLine)
        else:
            # insert line with the tag
            if not self.fm.isPresent(self.tag):
                self.fm.insert(self.tag)

        util.callAncestorFunction(self, 'hasChanged')
