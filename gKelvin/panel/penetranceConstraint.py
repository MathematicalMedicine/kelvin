"""The choice boxes and buttons for one penetrance constraint"""

import wx
from config import KelvinInfo
from config import misc
from fileManager import types
from panel.constraint import *

class PenetranceConstraintPanel_WC(ConstraintPanel):
    """The choice boxes and buttons for one penetrance constraint that
    is within classes (the same liability classes)

    """

    def __init__(self, parent, fileManager, line=None):
        ConstraintPanel.__init__(self, parent, fileManager, line)

        self.createPanel();
        self.arrangePanel();
        self.initialize()

        if self.isNew:
            self.reset()

    def createPanel(self):
        """Create widgets needed for the panel.
        For this constraint, there are only two parameters, and
        the choices for them are genotypes (DD, Dd, dd)

        """

        # choice boxes to choose the penetrance
        self.penChoice1 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.penChoice2 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.choiceBoxes.append(self.penChoice1)
        self.choiceBoxes.append(self.penChoice2)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        mainSizer.Add(self.penChoice1, 0, self.style, self.borderAmt)
        mainSizer.Add(self.opChoice, 0, self.style, self.borderAmt)
        mainSizer.Add(self.penChoice2, 0, self.style, self.borderAmt)
        mainSizer.Add(self.addButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.removeButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.orButton, 0, self.style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Set widgets to values based on file contents"""

        # choice box members
        self.penChoice1.list = KelvinInfo.genotypes[:]
        self.penChoice2.list = KelvinInfo.genotypes[:]
        self.penChoice1.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.penChoice2.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.penChoice1.descrToTag = KelvinInfo.genotypeDescriptionToTag
        self.penChoice2.descrToTag = KelvinInfo.genotypeDescriptionToTag

        # set responsibilities
        # the token representing the first penetrance is always first
        self.penChoice1.charge = self.line[0]

        # the token representing the second penetrance is always right
        # after the operator
        self.penChoice2.charge = self.line[self.getOperatorIndex()+1]

        # initialize values
        description = self.penChoice1.tagToDescr[self.penChoice1.charge.str]
        self.penChoice1.SetSelection(self.penChoice1.list.index(description))

        description = self.penChoice2.tagToDescr[self.penChoice2.charge.str]
        self.penChoice2.SetSelection(self.penChoice2.list.index(description))

        # bind events
        self.penChoice1.Bind(wx.EVT_CHOICE, self.onChoice, self.penChoice1)
        self.penChoice2.Bind(wx.EVT_CHOICE, self.onChoice, self.penChoice2)

    def generateDefaultLine(self):
        """Create a default line of this constraint"""

        defaultPen = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]

        line = self.fm.createLine(defaultPen, defaultOp, defaultPen)

        line[0].isConstraint = True
        line[0].type = types.ID
        line[1].type = types.OPERATOR
        line[2].type = types.ID

        return line

class PenetranceConstraintPanel_AC(PenetranceConstraintPanel_WC, 
                                   ConstraintPanel_AC):
    """The choice boxes and buttons for one penetrance constraint that
    goes across classes (mutiple liability classes)
    
    """

    def __init__(self, parent, fileManager, line=None):
        PenetranceConstraintPanel_WC.__init__(self, parent, fileManager, line)
        ConstraintPanel_AC.__init__(self)

        self.createPanel_AC();
        self.arrangePanel_AC();
        self.initialize_AC();

    def createPanel_AC(self):
        """Create widgets for the panel"""

        # no extra widgets needed
        pass

    def arrangePanel_AC(self):
        """Arrange widgets in the panel"""

        # insert the new choice boxes and text boxes inside the panel
        mainSizer = self.GetSizer()

        mainSizer.Insert(1, self.classChoice1, 0, self.style, self.borderAmt)
        mainSizer.Insert(1, self.classLabel1, 0, self.style, self.borderAmt)
        mainSizer.Insert(5, self.classChoice2, 0, self.style, self.borderAmt)
        mainSizer.Insert(5, self.classLabel2, 0, self.style, self.borderAmt)

        mainSizer.Layout()
        self.Fit()

    def initialize_AC(self):
        """Set widgets to values based on file contents"""

        # nothing to do
        pass

    def generateDefaultLine(self):
        """Create a default line of this constraint"""

        defaultPen = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]
        defaultClass = KelvinInfo.defaultLiabilityClass

        line =  self.fm.createLine(defaultPen, defaultClass, defaultOp, 
                                   defaultPen, defaultClass)
        line[0].isConstraint = True
        line[0].type = types.ID
        line[1].type = types.NUMBER
        line[2].type = types.OPERATOR
        line[3].type = types.ID
        line[4].type = types.NUMBER

        return line

