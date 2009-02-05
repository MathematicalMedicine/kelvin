"""The choice boxes and buttons for one standard deviation constraint"""

import wx
from config import KelvinInfo
from config import misc
from panel.constraint import *

class StdConstraintPanel_WC(ConstraintPanel):
    """The choice boxes and buttons for one standard deviation constraint that
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

        # choice boxes to choose the genotype
        self.genChoice1 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.genChoice2 = wx.Choice(self, -1, choices=KelvinInfo.genotypes)
        self.choiceBoxes.append(self.genChoice1)
        self.choiceBoxes.append(self.genChoice2)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        mainSizer.Add(self.genChoice1, 0, self.style, self.borderAmt)
        mainSizer.Add(self.opChoice, 0, self.style, self.borderAmt)
        mainSizer.Add(self.genChoice2, 0, self.style, self.borderAmt)
        mainSizer.Add(self.addButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.removeButton, 0, self.style, self.borderAmt)
        mainSizer.Add(self.orButton, 0, self.style, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Set widgets to values based on file contents"""

        # choice box members
        self.genChoice1.list = KelvinInfo.genotypes[:]
        self.genChoice2.list = KelvinInfo.genotypes[:]
        self.genChoice1.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.genChoice2.tagToDescr = KelvinInfo.genotypeTagToDescription
        self.genChoice1.descrToTag = KelvinInfo.genotypeDescriptionToTag
        self.genChoice2.descrToTag = KelvinInfo.genotypeDescriptionToTag

        # set responsibilities
        # the token representing the first genotype is always second
        self.genChoice1.charge = self.line[1]

        # the token representing the second genotype is always two
        # places after the operator
        self.genChoice2.charge = self.line[self.getOperatorIndex()+2]

        # initialize values
        description = self.genChoice1.tagToDescr[self.genChoice1.charge.str]
        self.genChoice1.SetSelection(self.genChoice1.list.index(description))

        description = self.genChoice2.tagToDescr[self.genChoice2.charge.str]
        self.genChoice2.SetSelection(self.genChoice2.list.index(description))

        # bind events
        self.genChoice1.Bind(wx.EVT_CHOICE, self.onChoice, self.genChoice1)
        self.genChoice2.Bind(wx.EVT_CHOICE, self.onChoice, self.genChoice2)

    def generateDefaultLine(self):
        """Create a default line of this constraint"""

        stdParameter = KelvinInfo.stdParameter[0]
        defaultGen = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]

        line = self.fm.createLine(stdParameter, defaultGen, defaultOp, 
                                  stdParameter, defaultGen)

        line[0].type = types.ID
        line[1].type = types.ID
        line[2].type = types.OPERATOR
        line[3].type = types.ID
        line[4].type = types.ID

        return line

class StdConstraintPanel_AC(StdConstraintPanel_WC, ConstraintPanel_AC):
    """The choice boxes and buttons for one standard deviation constraint that
    goes across classes (mutiple liability classes)
    
    """

    def __init__(self, parent, fileManager, line=None):
        StdConstraintPanel_WC.__init__(self, parent, fileManager, line)
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

        stdParameter = KelvinInfo.stdParameter[0]
        defaultGen = KelvinInfo.genotypes[0]
        defaultOp = KelvinInfo.operators[0][0]
        defaultClass = KelvinInfo.defaultLiabilityClass

        line = self.fm.createLine(stdParameter, defaultGen, defaultClass, 
                                  defaultOp, stdParameter, defaultGen, 
                                  defaultClass)

        line[0].type = types.ID
        line[1].type = types.ID
        line[2].type = types.NUMBER
        line[3].type = types.OPERATOR
        line[4].type = types.ID
        line[5].type = types.ID
        line[6].type = types.NUMBER

        return line


