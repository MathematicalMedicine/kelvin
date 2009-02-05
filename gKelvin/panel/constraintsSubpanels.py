"""This module contains panels that are used for displaying certain kinds of constraints, e.g. the panel that displays all penetrance constraints, the panel that displays all the theta constraints, etc."""

import wx
from config import KelvinInfo
from fileManager import types
from panel.scrolledPanel import ScrolledPanel
from panel.meanConstraint import *
from panel.penetranceConstraint import *
from panel.stdConstraint import *
from panel.thetaConstraint import ThetaConstraintPanel
from utility import util

class ConstraintsPanel(wx.Panel):
    """A base class for all the other types of constraints panels.
    A panel that shows a heading with a number of constraint panels under it.  
    Defines common functionality in all constraint panels.
    
    Keyword arguments for __init__():
    constraintPanelType -- a class that is the type of constraint panel that 
                           this panel uses

    getConstraintLines -- a function that gets all the lines of this constraint 
                          type from the parsetree

    indent -- a boolean to indicate whether or not to indent the list of 
              constraints in relation to the heading
        
    """

    mainLbl = "The main heading -- derived classes should override this"

    borderAmt = 5
    constraintBorder = 1
    emptyHSpace = (20,1)

    def __init__(self, parent, fileManager, constraintPanelType,
                 getConstraintLines, indent=True):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.constraintPanelType = constraintPanelType

        # whether the panel is enabled or not
        self.isEnabled = True

        # get all the lines of this constraint type
        self.lines = getConstraintLines()

        # a list of all the constraintPanels in this panel
        # this is actually a list of ConstraintLinePanels
        self.constraintPanels = []

        self.createPanel_ConstraintsPanel()
        self.arrangePanel_ConstraintsPanel(indent)
        self.initialize_ConstraintsPanel()

    def createPanel_ConstraintsPanel(self):
        """Create widgets for this panel"""

        self.mainLabel = wx.StaticText(self, -1, self.mainLbl)

        for eachLine in self.lines:
            self.constraintPanels.append(self.constraintPanelType(self,self.fm,
                                                                  eachLine))

        # if there are no lines in the file, make at least one constraint
        if len(self.lines) == 0:
            self.constraintPanels.append(self.constraintPanelType(self,self.fm))

    def arrangePanel_ConstraintsPanel(self, indent):
        """Arrange widgets in this panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        mainSizer.Add(self.mainLabel, 0, wx.ALL, self.borderAmt)

        # a horizontal box sizer that will hold empty space on the left and
        # the constraints on the right
        sizer = wx.BoxSizer(wx.HORIZONTAL)


        if indent:
            sizer.Add(self.emptyHSpace, 0, wx.EXPAND)

        # sizer for the constraints
        self.constraintSizer = wx.BoxSizer(wx.VERTICAL)
        
        for eachPanel in self.constraintPanels:
            self.constraintSizer.Add(eachPanel, 0, wx.ALL,self.constraintBorder)

        sizer.Add(self.constraintSizer)
        mainSizer.Add(sizer)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize_ConstraintsPanel(self):
        """Initialize widgets according to file contents"""

        # nothing to do
        pass

    def addConstraint(self, constraintPanel=None, line=None):
        """Add another constraint panel right underneath the given 
        constraintPanel.
        
        """

        self.Freeze()

        if constraintPanel == None:
            # use the last constraint panel
            constraintPanel = self.constraintPanels[-1]

        newConstraintPanel = self.constraintPanelType(self, self.fm, line)
        index = self.constraintPanels.index(constraintPanel)
        self.constraintPanels.insert(index+1, newConstraintPanel)
        self.constraintSizer.Insert(index+1, newConstraintPanel, 0, wx.ALL, 1)

        self.reSize()
        self.Thaw()

    def removeConstraint(self, constraintPanel):
        """Remove the given constraint panel from this panel unless it's the 
        last one there, in which case just reset it.
        
        """

        # Note: This function only removes the constraintLinePanel from
        # this panel. It doesn't modify the parsetree or anything else.
        # Right now this is only called in the 
        # ConstraintLinePanel.removeConstraint() function

        if len(self.constraintPanels) != 1:
            try:
                self.constraintPanels.remove(constraintPanel)
            except ValueError, e:
                # if there was an error, the panel was not found in the list, 
                # which means that the panel is already gone, so nothing to do
                return
            else:
                self.constraintSizer.Remove(constraintPanel)
                constraintPanel.Destroy()
                self.reSize()
        else:
            if constraintPanel.error:
                # currently in error state and is the last one in the list,
                # need to disable it then
                constraintPanel.Disable()
            constraintPanel.reset()

    def reSize(self):
        """Activate an ancestor's FitInside() function in order for scroll bars
        to function properly
        
        """

        # should be called after the panel removes or adds things
        # need to correctly update the scroll bars
        util.callAncestorTypeFunction(self, ScrolledPanel, 'FitInside')

    def Enable(self, bool=True):
        """Enables or disables panel, extra code to deal with constraint list"""

        self.isEnabled = bool
        self.mainLabel.Enable(bool)

        if bool:
            self.enableConstraintPanels()
        else:
            self.disableConstraintPanels()

    def Disable(self):
        """Disables panel"""

        self.Enable(False)

    def enableConstraintPanels(self):
        """Enables the list of constraint panels"""

        for eachPanel in self.constraintPanels:
            eachPanel.Enable()

    def disableConstraintPanels(self):
        """Disables the list of constraint panels"""

        for eachPanel in self.constraintPanels[:]:
            eachPanel.Disable()
            if eachPanel.isAllChoicesBlank() and len(self.constraintPanels) > 1:
                # not the last blank one, so remove it
                self.removeConstraint(eachPanel)

    def clear(self):
        """Remove all constraints from this panel"""

        # basically simulating pressing the remove button for each constraint
        for eachConstraintLine in self.constraintPanels[:]:
            for eachConstraint in eachConstraintLine.constraintPanels[:]:
                # Note: i'm not using util.createEvent because it's too slow,
                # since changing trait type usually results in lots 
                # of constraints being removed automatically
                eachConstraint.onRemoveButton(None)

    def clean(self):
        """Remove all blank constraints from this panel"""

        for eachConstraintLine in self.constraintPanels[:]:
            for eachConstraint in eachConstraintLine.constraintPanels[:]:
                if eachConstraint.isAllChoicesBlank():
                    # Note: i'm not using util.createEvent because it's too 
                    # slow, since changing trait type usually results in lots 
                    # of constraints being cleaned 
                    eachConstraint.onRemoveButton(None)

    def hasConstraint(self, line):
        """Given a line, returns True if that line is already represented in 
        the constraints panel. Handy when adding in defaults and you don't 
        want to input the same constraint multiple times.
        
        """

        lineLength = len(line)
        for eachConstraintLine in self.constraintPanels:
            otherLine = eachConstraintLine.line
            if otherLine is not None and lineLength == len(otherLine):
                result = True
                for i in range(lineLength):
                    if line[i].str != eachConstraintLine.line[i].str:
                        result=False
                        break
                if result == True:
                    return True
        return False

class ConstraintsWithLiabilityPanel(wx.Panel):
    """A base class for all constaints panels that can have a liability class, 
    for example the penetrance constraints.  This class should never directly 
    be used. In essence this class is a normal panel with two ConstraintsPanels.

    """

    borderAmt = ConstraintsPanel.borderAmt
    constraintBorder = ConstraintsPanel.constraintBorder
    emptyHSpace = ConstraintsPanel.emptyHSpace
    withinClassLbl = 'Within Class:'
    acrossClassLbl = 'Across Classes:'

    def __init__(self, parent, fileManager,constraintPanelType,
                 constraintLiabilityPanelType, getConstraintLines, 
                 getConstraintLiabilityLines, indent=True):

        wx.Panel.__init__(self, parent, -1)

        self.fm = fileManager

        self.createPanel_ConstraintsWithLiabilityPanel(fileManager, 
                            constraintPanelType, constraintLiabilityPanelType, 
                            getConstraintLines, getConstraintLiabilityLines)
        self.arrangePanel_ConstraintsWithLiabilityPanel(indent)
        self.intialize_ConstraintsWithLiabilityPanel()

    def createPanel_ConstraintsWithLiabilityPanel(self, fileManager,
                            constraintPanelType, constraintLiabilityPanelType, 
                            getConstraintLines, getConstraintLiabilityLines):
        """Create the headings and the two subpanels"""

        self.mainLabel = wx.StaticText(self, -1, self.mainLbl)

        # this subpanel lists all the constraints without liability classes
        self.constraintsPanel = ConstraintsPanel(self, fileManager, 
                                 constraintPanelType, getConstraintLines, False)

        # a subpanel that lists all constraints with liability classes
        self.liabilityPanel = ConstraintsPanel(self, fileManager, 
                                            constraintLiabilityPanelType, 
                                            getConstraintLiabilityLines, False)

        # need to change the labels on the subpanels
        self.constraintsPanel.mainLabel.SetLabel(self.withinClassLbl)
        self.liabilityPanel.mainLabel.SetLabel(self.acrossClassLbl)

    def arrangePanel_ConstraintsWithLiabilityPanel(self, indent):
        """Arrange widgets based on file contents"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        mainSizer.Add(self.mainLabel, 0, wx.ALL, self.borderAmt)

        # a horizontal box sizer that will hold empty space on the left 
        # since there is an indention there, and 
        # the constraints subpanels on the right
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        if indent:
            sizer.Add(self.emptyHSpace, 0, wx.EXPAND)

        # sizer for the constraints subpanels
        constraintsPanelsSizer = wx.BoxSizer(wx.VERTICAL)
        constraintsPanelsSizer.Add(self.constraintsPanel, 0)
        constraintsPanelsSizer.Add((1, 5))
        constraintsPanelsSizer.Add(self.liabilityPanel, 0)

        sizer.Add(constraintsPanelsSizer)
        mainSizer.Add(sizer)

        self.SetSizer(mainSizer)
        self.Fit()

    def intialize_ConstraintsWithLiabilityPanel(self):
        """Intiailization"""

        # need to know when the liability class changes
        util.bindLiabilityClassChange(self, self.onLiabilityChange)
        self.onLiabilityChange(None)

    def onLiabilityChange(self, event):
        """Event when the number of liability class changes"""

        # do nothing if the panel is already disabled
        if self.constraintsPanel.isEnabled is False:
            return

        # check the number of liability classes specified
        tag = KelvinInfo.liabilityClassEntry[0]
        line = self.fm.getLines(tag)
        if len(line) > 0:
            # liability class is present, get the number
            try:
                n = int(line[0][1].str)
            except ValueError, e:
                # happens when the value is blank or a bozo typed
                # a letter in instead of a number, in which case do nothing
                return
        else:
            # liability class is not present, value is default
            n = int(KelvinInfo.defaultLiabilityClass)

        if n >= 2:
            self.liabilityPanel.Enable()
        else:
            self.liabilityPanel.clear()
            self.liabilityPanel.Disable()

    def Enable(self, bool=True):
        """Enables or disables panel, extra code to disable value panels"""

        # enable or disable all labels
        self.mainLabel.Enable(bool)
        self.constraintsPanel.Enable(bool)
        self.liabilityPanel.Enable(bool)

        if bool:
            # may have to redisable something based on liability classes
            self.onLiabilityChange(None)

    def Disable(self):
        """Disables panel"""

        self.Enable(False)

    def clear(self):
        """Remove all constraints in this panel"""

        self.constraintsPanel.clear()
        self.liabilityPanel.clear()

    def clean(self):
        """Remove all blank constraints in this panel"""

        self.constraintsPanel.clean()
        self.liabilityPanel.clean()

    def addConstraint(self, line):
        """Add a constraint to this panel. Determines if it is a normal 
        constraint or one with a liability class, then add the constraint 
        to the end of the correct list.
        
        """

        # first check if the line deals with liability classes or not
        # a constraint has a liability class if one of the tokens
        # in the line is a number
        isLiability = False
        for eachToken in line:
            if eachToken.type == types.NUMBER:
                isLiability = True

        if isLiability:
            self.liabilityPanel.addConstraint(line=line)
        else:
            self.constraintsPanel.addConstraint(line=line)

    def hasConstraint(self, line):
        """Given a line, returns True if that line is already represented in 
        the constraints panel. Handy when adding in defaults and you don't 
        want to input the same constraint multiple times.
        
        """

        return self.constraintsPanel.hasConstraint(line) or \
               self.liabilityPanel.hasConstraint(line)

class PenetranceConstraintsPanel(ConstraintsWithLiabilityPanel):
    """Defines the panel used to show penetrance constraints"""

    mainLbl = 'Constraints on Penetrances'

    def __init__(self, parent, fileManager):

        # define functions that will give the correct panel types in order
        # that they be passed to the parent
        def PenConstraintPanelLine_WC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       PenetranceConstraintPanel_WC, line)

        def PenConstraintPanelLine_AC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       PenetranceConstraintPanel_AC, line)

        ConstraintsWithLiabilityPanel.__init__(self, parent, fileManager, 
                PenConstraintPanelLine_WC, PenConstraintPanelLine_AC, 
                self.getPenetranceConstraintLines, 
                self.getPenetranceLiabilityConstraintLines)

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Creates widgets for this panel"""
        # done by the parent class
        pass

    def arrangePanel(self):
        """Places widgets in panel"""
        # done by parent class
        pass

    def initialize(self):
        """Intializes widgets based on file contents"""

        # mostly done by parent class

        # need to know when the trait type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)
        util.bindTraitTypeChange(self, self.onTraitChange)
        self.checkState()

    def getPenetranceConstraintLines(self):
        tags = KelvinInfo.genotypes
        def isPenetranceConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) != 3:
                return False
            for eachTag in tags:
                if line[0].str == eachTag:
                    return True
            return False

        # penetrance constraints are only under DT (dichotomous model)
        # if QT is specified, then these constraints are not penetrance ones
        if self.fm.isPresent('QT'):
            return []

        return filter(isPenetranceConstraint, self.fm.getConstraintLines())

    def getPenetranceLiabilityConstraintLines(self):
        tags = KelvinInfo.genotypes
        def isPenetranceConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) != 5:
                return False
            for eachTag in tags:
                if line[0].str == eachTag:
                    return True
            return False

        # penetrance constraints are only under DT (dichotomous model)
        if self.fm.isPresent('QT'):
            return []

        return filter(isPenetranceConstraint, self.fm.getConstraintLines())

    def onAnalysisChange(self, event):
        """Event when the analysis type changes"""
        #self.checkState(True)
        pass

    def onTraitChange(self, event):
        """Event when the trait type changes"""
        self.checkState(True)

    def checkState(self, addDefaults=False):
        """Check if these constraints are valid under current conditions,
        add default values if they are needed.
        
        """

        # these constraints are only valid in dichotomous model
        if self.fm.isPresent('DT'):
            self.Enable()
            if addDefaults:
                self.addDefaults()
                self.clean()
        else:
            self.clear()
            self.Disable()

    def addDefaults(self):
        """Adds default constraints to the panel"""

        # any penetrance constraints must start with a genotype
        lines = []
        for eachTag in KelvinInfo.genotypes:
            lines.extend(self.fm.getDefaults(eachTag, isConstraint=True))

        for eachLine in lines:
            if not self.hasConstraint(eachLine):
                self.fm.addLine(eachLine)
                self.addConstraint(eachLine)

class ThetaConstraintsPanel(ConstraintsPanel):
    """Defines the panel used to show theta constraints"""

    mainLbl = 'Recombination Fraction (Theta) Constraints'

    def __init__(self, parent, fileManager):

        # define a function that will give the correct panel types in order
        # that they be passed to the parent
        def ThetaConstraintPanelLine(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       ThetaConstraintPanel, line)

        ConstraintsPanel.__init__(self, parent, fileManager, 
                                  ThetaConstraintPanelLine, 
                                  self.getThetaConstraints, True)

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Creates widgets for this panel"""
        # done by the parent class
        pass

    def arrangePanel(self):
        """Places widgets in panel"""
        # done by parent class
        pass

    def initialize(self):
        """Intializes widgets based on file contents"""

        # mostly done by parent class

        # need to know when the analysis type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)
        util.bindTraitTypeChange(self, self.onTraitChange)
        self.checkState()

    def getThetaConstraints(self):
        def isThetaConstraint(line):
            for eachTag in KelvinInfo.thetaTags:
                if line[0].str == eachTag:
                    return True
            return False
        return filter(isThetaConstraint, self.fm.getConstraintLines())

    def onAnalysisChange(self, event):
        """Event when the analysis type changes"""
        self.checkState(True)

    def onTraitChange(self, event):
        """Event when the trait type changes"""
        self.checkState(True)

    def checkState(self, addDefaults=False):
        """Check if these constraints are valid under current conditions"""

        # this needs to be called because need to know sex specific or 
        # sex average, and if it's two point, there is no change in the config
        # file, so can't check the state of the file directly for that info
        isSexAverage = util.callAncestorFunction(self, 'isSexAverage')

        # these constraints are only valid under two point, sex specific
        if self.fm.isPresent('TP') and not isSexAverage:
            self.Enable()
            if addDefaults:
                self.addDefaults()
                self.clean()
        else:
            self.clear()
            self.Disable()

    def addDefaults(self):
        """Adds default constraints to the panel"""

        # any theta constraints would have to deal with Tf or Tm
        lines = []
        for eachEntry in KelvinInfo.thetaConstraintParameters:
            lines.extend(self.fm.getDefaults(eachEntry[0], isConstraint=True))

        for eachLine in lines:
            if not self.hasConstraint(eachLine):
                self.fm.addLine(eachLine)
                self.addConstraint(line=eachLine)

class MeanConstraintsPanel(ConstraintsWithLiabilityPanel):
    """Defines the panel used to show mean constraints"""

    mainLbl = 'Constraints on QT Mean'

    def __init__(self, parent, fileManager):
        # define functions that will give the correct panel types in order
        # that they be passed to the parent
        def MeanConstraintPanelLine_WC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       MeanConstraintPanel_WC, line)

        def MeanConstraintPanelLine_AC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       MeanConstraintPanel_AC, line)

        ConstraintsWithLiabilityPanel.__init__(self, parent, fileManager, 
                MeanConstraintPanelLine_WC, MeanConstraintPanelLine_AC, 
                self.getMeanConstraintLines, 
                self.getMeanLiabilityConstraintLines)

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Creates widgets for this panel"""
        # done by the parent class
        pass

    def arrangePanel(self):
        """Places widgets in panel"""
        # done by parent class
        pass

    def initialize(self):
        """Intializes widgets based on file contents"""

        # mostly done by parent class

        # need to know when the analysis type and trait type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)
        util.bindTraitTypeChange(self, self.onTraitChange)
        self.checkState()

    def getMeanConstraintLines(self):
        tags = KelvinInfo.genotypes
        def isMeanConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) != 3:
                return False
            for eachTag in tags:
                if line[0].str == eachTag:
                    return True
            return False

        # mean constraints are only under QT (quantitative model)
        if not self.fm.isPresent('QT'):
            return []

        return filter(isMeanConstraint, self.fm.getConstraintLines())

    def getMeanLiabilityConstraintLines(self):
        tags = KelvinInfo.genotypes
        def isMeanConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) != 5:
                return False
            for eachTag in tags:
                if line[0].str == eachTag:
                    return True
            return False

        # mean constraints are only under QT (quantitative model)
        if not self.fm.isPresent('QT'):
            return []

        return filter(isMeanConstraint, self.fm.getConstraintLines())

    def onAnalysisChange(self, event):
        """Event when the analysis type changes"""
        #self.checkState(True)
        pass

    def onTraitChange(self, event):
        """Event when the trait type changes"""
        self.checkState(True)

    def checkState(self, addDefaults=False):
        """Check if these constraints are valid under current conditions,
        and add default values if they are needed.
        
        """

        # these constraints are only valid in quantitative model
        if self.fm.isPresent('QT'):
            self.Enable()
            if addDefaults:
                self.addDefaults()
                self.clean()
        else:
            self.clear()
            self.Disable()

    def addDefaults(self):
        """Adds default constraints to the panel"""

        # any mean constraint must start with a genotype
        lines = []
        for eachTag in KelvinInfo.genotypes:
            lines.extend(self.fm.getDefaults(eachTag, isConstraint=True))

        for eachLine in lines:
            if not self.hasConstraint(eachLine):
                self.fm.addLine(eachLine)
                self.addConstraint(eachLine)

class StdConstraintsPanel(ConstraintsWithLiabilityPanel):
    """Defines the panel used to show standard deviation constraints"""

    mainLbl = 'Constraints on QT Standard Deviation'

    def __init__(self, parent, fileManager):

        # define functions that will give the correct panel types in order
        # that they be passed to the parent
        def StdConstraintPanelLine_WC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       StdConstraintPanel_WC, line)

        def StdConstraintPanelLine_AC(parent, fileManager, line=None):
            return ConstraintLinePanel(parent, fileManager, 
                                       StdConstraintPanel_AC, line)

        ConstraintsWithLiabilityPanel.__init__(self, parent, fileManager, 
                StdConstraintPanelLine_WC, StdConstraintPanelLine_AC, 
                self.getStdConstraintLines, 
                self.getStdLiabilityConstraintLines)

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Creates widgets for this panel"""
        # done by the parent class
        pass

    def arrangePanel(self):
        """Places widgets in panel"""
        # done by parent class
        pass

    def initialize(self):
        """Intializes widgets based on file contents"""

        # mostly done by parent class

        # need to know when the trait type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)
        util.bindTraitTypeChange(self, self.onTraitChange)
        self.checkState()

    def getStdConstraintLines(self):
        tag = KelvinInfo.stdParameter[0]
        def isStdConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) == 5 and line[0].str == tag:
                return True
            return False

        return filter(isStdConstraint, self.fm.getConstraintLines())

    def getStdLiabilityConstraintLines(self):
        def isStdConstraint(line):
            line = self.fm.getFirstConstraint(line)
            if len(line) == 7:
                return True

        return filter(isStdConstraint, self.fm.getConstraintLines())

    def onAnalysisChange(self, event):
        """Event when the analysis type changes"""
        #self.checkState(True)
        pass

    def onTraitChange(self, event):
        """Event when the trait type changes"""
        self.checkState(True)

    def checkState(self, addDefaults=False):
        """Check if these constraints are valid under current conditions,
        and add default values if they are needed.
        
        """

        # these constraints are only valid in quantitative model
        if self.fm.isPresent('QT'):
            self.Enable()
            if addDefaults:
                self.addDefaults()
                self.clean()
        else:
            self.clear()
            self.Disable()

    def addDefaults(self):
        """Adds default constraints to the panel"""

        # any std constraint must start with the std parameter (P1)
        tag = KelvinInfo.stdParameter[0]
        lines = self.fm.getDefaults(tag, isConstraint=True)

        for eachLine in lines:
            if not self.hasConstraint(eachLine):
                self.fm.addLine(eachLine)
                self.addConstraint(eachLine)
