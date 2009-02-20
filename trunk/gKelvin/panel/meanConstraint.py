"""The choice boxes and buttons for one mean constraint"""

import wx
from config import KelvinInfo
from config import misc
from panel.constraint import *

class MeanConstraintPanel_WC(ConstraintPanel):
    """The choice boxes and buttons for one mean constraint that
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
        """Create widgets for the panel"""

        # for this constraint, there are only two parameters, and
        # the choices for them are the same as genotypes (DD, Dd, dd)
        self.meanChoice1 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.meanChoice2 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.choiceBoxes.append(self.meanChoice1)
        self.choiceBoxes.append(self.meanChoice2)

    def arrangePanel(self):
        """Place widgets in this panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        mainSizer.Add(self.meanChoice1, 0, self.style, self.borderAmt)
        mainSizer.Add(self.opChoice, 0, self.style, self.borderAmt)
        mainSizer.Add(self.meanChoice2, 0, self.style, self.borderAmt)
        mainSizer.Add(self.addButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.removeButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.orButton, 0, self.style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Set widgets to values based on file contents"""

        # choice box members
        self.meanChoice1.list = KelvinInfo.genotypes[:]
        self.meanChoice2.list = KelvinInfo.genotypes[:]
        self.meanChoice1.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.meanChoice2.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.meanChoice1.descrToTag = KelvinInfo.genotypeDescriptionToTag
        self.meanChoice2.descrToTag = KelvinInfo.genotypeDescriptionToTag

        # set responsibilities
        # the token representing the first mean is always first
        self.meanChoice1.charge = self.line[0]

        # the token representing the second mean is always right
        # after the operator
        self.meanChoice2.charge = self.line[self.getOperatorIndex()+1]

        # initialize values
        description = self.meanChoice1.tagToDescr[self.meanChoice1.charge.str]
        self.meanChoice1.SetSelection(self.meanChoice1.list.index(description))

        description = self.meanChoice2.tagToDescr[self.meanChoice2.charge.str]
        self.meanChoice2.SetSelection(self.meanChoice2.list.index(description))

        # bind events
        self.meanChoice1.Bind(wx.EVT_CHOICE, self.onChoice, self.meanChoice1)
        self.meanChoice2.Bind(wx.EVT_CHOICE, self.onChoice, self.meanChoice2)

    def generateDefaultLine(self):
        """Create a default line of this constraint"""

        defaultMean = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]

        line = self.fm.createLine(defaultMean, defaultOp, defaultMean)

        line[0].isConstraint = True
        line[0].type = types.ID
        line[1].type = types.OPERATOR
        line[2].type = types.ID

        return line

class MeanConstraintPanel_AC(MeanConstraintPanel_WC, ConstraintPanel_AC):
    """The choice boxes and buttons for one mean constraint that
    goes across classes (mutiple liability classes)
    
    """

    def __init__(self, parent, fileManager, line=None):
        MeanConstraintPanel_WC.__init__(self, parent, fileManager, line)
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

        defaultMean = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]
        defaultClass = KelvinInfo.defaultLiabilityClass

        line =  self.fm.createLine(defaultMean, defaultClass, defaultOp, 
                                   defaultMean, defaultClass)
        line[0].isConstraint = True
        line[0].type = types.ID
        line[1].type = types.NUMBER
        line[2].type = types.OPERATOR
        line[3].type = types.ID
        line[4].type = types.NUMBER

        return line

