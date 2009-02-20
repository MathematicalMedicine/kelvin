"""The panel showing all the constraints"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from panel.constraintsSubpanels import PenetranceConstraintsPanel
from panel.constraintsSubpanels import ThetaConstraintsPanel
from panel.constraintsSubpanels import MeanConstraintsPanel
from panel.constraintsSubpanels import StdConstraintsPanel

class ConstraintsPanel(wx.Panel):
    """The panel to show all the different kinds of constraints together"""

    # some empty vertical space for padding
    emptyVSpace = (1,20)

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        # each type of constraint has its own subpanel
        self.penConstraintsPanel = PenetranceConstraintsPanel(self, self.fm)
        self.thetaConstraintsPanel = ThetaConstraintsPanel(self, self.fm)
        self.meanConstraintsPanel = MeanConstraintsPanel(self, self.fm)
        self.stdConstraintsPanel = StdConstraintsPanel(self, self.fm)

    def arrangePanel(self):
        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # penetrance constraints
        mainSizer.Add(self.penConstraintsPanel)

        mainSizer.Add(self.emptyVSpace)

        # theta constraints
        mainSizer.Add(self.thetaConstraintsPanel)

        mainSizer.Add(self.emptyVSpace)

        # mean constraints
        mainSizer.Add(self.meanConstraintsPanel)

        mainSizer.Add(self.emptyVSpace)

        # standard deviation constraints
        mainSizer.Add(self.stdConstraintsPanel)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        pass
