"""The choice boxes and buttons for one theta constraint"""

import wx
from config import KelvinInfo
from config import misc
from fileManager import types
from panel.constraint import ConstraintPanel

class ThetaConstraintPanel(ConstraintPanel):
    """The choice boxes and buttons for one theta constraint"""

    def __init__(self, parent, fileManager, line=None):
        ConstraintPanel.__init__(self, parent, fileManager, line)

        self.createPanel()
        self.arrangePanel()
        self.initialize()

        if self.isNew:
            self.reset()

    def createPanel(self):
        """Create widgets needed for the panel.
        For this constraint, there are only two parameters, and
        the choices for them are Tm and Tf (theta male and female)

        """

        choiceList = KelvinInfo.thetaDescriptions
        self.thetaChoice1 = wx.Choice(self, -1, choices=choiceList)
        self.thetaChoice2 = wx.Choice(self, -1, choices=choiceList)
        self.choiceBoxes.append(self.thetaChoice1)
        self.choiceBoxes.append(self.thetaChoice2)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        mainSizer.Add(self.thetaChoice1, 0, self.style, self.borderAmt)
        mainSizer.Add(self.opChoice,0, self.style, self.borderAmt)
        mainSizer.Add(self.thetaChoice2, 0, self.style, self.borderAmt)
        mainSizer.Add(self.addButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.removeButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.orButton, 0, self.style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):

        # choice box members
        self.thetaChoice1.list = KelvinInfo.thetaDescriptions
        self.thetaChoice2.list = KelvinInfo.thetaDescriptions
        self.thetaChoice1.tagToDescr = KelvinInfo.thetaTagToDescription
        self.thetaChoice2.tagToDescr = KelvinInfo.thetaTagToDescription
        self.thetaChoice1.descrToTag = KelvinInfo.thetaDescriptionToTag
        self.thetaChoice2.descrToTag = KelvinInfo.thetaDescriptionToTag

        # set responsibilities
        self.thetaChoice1.charge = self.line[0]
        self.thetaChoice2.charge = self.line[2]

        # initialize values
        description = self.thetaChoice1.tagToDescr[self.thetaChoice1.charge.str]
        index = self.thetaChoice1.list.index(description)
        self.thetaChoice1.SetSelection(index)

        description = self.thetaChoice2.tagToDescr[self.thetaChoice2.charge.str]
        index = self.thetaChoice2.list.index(description)
        self.thetaChoice2.SetSelection(index)

        # bind events
        self.thetaChoice1.Bind(wx.EVT_CHOICE, self.onChoice, self.thetaChoice1)
        self.thetaChoice2.Bind(wx.EVT_CHOICE, self.onChoice, self.thetaChoice2)

    def generateDefaultLine(self):
        """Create a default line of this constraint"""

        defaultTheta = KelvinInfo.thetaTags[0]
        defaultOp = KelvinInfo.operators[0][0]

        line = self.fm.createLine(defaultTheta, defaultOp, defaultTheta)

        line[0].isConstraint = True
        line[0].type = types.ID
        line[1].type = types.OPERATOR
        line[2].type = types.ID

        return line

