"""Contains base classes for all other constraints"""

import wx
from config import KelvinInfo
from config import misc
from fileManager import types
from panel.scrolledPanel import ScrolledPanel
from utility import util 

class ConstraintPanel(wx.Panel):
    """A base class for all other constraints. This class represents exactly 
    one constraint.  Right now it takes care of the operator choice box and 
    general intialization.
    
    """

    borderAmt = 3

    # the style of the widgets and where to place the borders
    style = wx.LEFT | wx.RIGHT | wx.TOP | wx.ALIGN_CENTER_VERTICAL

    blank = ''

    def __init__(self, parent, fileManager, line=None):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        # indicates whether the line that the constraint panel is
        # responsible for is currently inside the file
        self.inFile = False

        # indicates whether the panel is in the error state or not
        # it's in the error state when the particular constraint
        # is illegal under current file conditions
        # True if in error state, False otherwise
        self.error = False

        # a list of all the choice boxes in this constraint
        self.choiceBoxes = []

        if line is not None:
            self.line = line
            self.inFile = True
            self.isNew = False
        else:
            self.isNew = True
            self.line = self.generateDefaultLine()

        self.createPanel_ConstraintPanel();
        self.initialize_ConstraintPanel();

        if self.isNew:
            self.reset()

    def createPanel_ConstraintPanel(self):
        """Create widgets that are common to all constraints"""

        # a choice box to choose what operation the constraint is
        # (greater than, not equal to, etc.)
        self.opChoice = wx.Choice(self, -1, 
                                  choices=KelvinInfo.operatorDescriptions)
        self.choiceBoxes.append(self.opChoice)

        # add and remove constraint buttons
        buttonSize = misc.addRemoveButtonSize
        self.addButton = wx.Button(self, -1, label="+", size=buttonSize)
        self.removeButton = wx.Button(self, -1, label="-", size=buttonSize)

        # an or button to be able to or constraints together
        self.orButton = wx.Button(self, -1, label='OR', size=buttonSize)

    def initialize_ConstraintPanel(self):
        """Initialize widget values based on file contents"""

        # set responsibility and intialize operator 
        self.opChoice.list = KelvinInfo.operatorDescriptions[:]
        self.opChoice.tagToDescr = KelvinInfo.operatorTagToDescription
        self.opChoice.descrToTag = KelvinInfo.operatorDescriptionToTag
        self.opChoice.charge = [x for x in self.line \
                                            if x.type == types.OPERATOR][0]
        description = self.opChoice.tagToDescr[self.opChoice.charge.str]
        index = self.opChoice.list.index(description)
        self.opChoice.SetSelection(index)

        # bind events
        self.opChoice.Bind(wx.EVT_CHOICE, self.onChoice, self.opChoice)
        self.addButton.Bind(wx.EVT_BUTTON, self.onAddButton, self.addButton)
        self.removeButton.Bind(wx.EVT_BUTTON, self.onRemoveButton,
                                                            self.removeButton)
        self.orButton.Bind(wx.EVT_BUTTON, self.onOrButton, self.orButton)
    
    def generateDefaultLine(self):
        """Make a default line for the file"""

        # needs to be implemented by classes that derive this class
        pass

    def onChoice(self, event):
        """Event when a choice box is clicked"""

        choice = event.GetEventObject()
        index = choice.GetSelection()

        # if there is a blank in the choice list and a blank was not selected,
        # take out the blank entry
        if self.isChoiceBlank(choice) and index != 0:
            self.unBlankChoice(choice)
            index = index - 1
            choice.SetSelection(index)

        tag = choice.descrToTag[choice.list[index]]
        choice.charge.update(tag)

        # insert file if needed
        if self.inFile is False and not self.isAnyChoicesBlank():
            self.addToFile()
        self.beingUsed()

    def onAddButton(self, event):
        """Event when the add button is clicked"""

        # only able to add when not in error state
        if self.error is False:
            util.callAncestorFunction(self, 'addConstraint', self)
            self.beingUsed()

    def onRemoveButton(self, event):
        """Event when the remove button is clicked"""
        
        self.beingUsed()
        self.removeFromFile()
        util.callAncestorFunction(self, 'removeConstraint', self)

    def onOrButton(self, event):
        """Event when the OR button is pressed"""

        if self.error is False:
            util.callAncestorFunction(self, 'addConstraintToLine')
            self.beingUsed()

    def getOperatorIndex(self):
        """Returns the index at which the operator is at in the line"""

        if self.line is None:
            return

        for i in range(0, len(self.line)):
            if self.line[i].type == types.OPERATOR:
                return i

    def reset(self):
        """Takes the line out of the file and blanks out all choice boxes"""

        self.removeFromFile()

        for eachChoiceBox in self.choiceBoxes:
            self.blankChoice(eachChoiceBox)

        # reset error state
        if self.error:
            self.setError(False)

    def unBlankChoice(self, choice, selection=0):
        """Takes out blank entries in a choice box"""

        # check if there really is a blank, if so, take it out of the choices
        found = False
        index = choice.FindString(self.blank)
        while index != wx.NOT_FOUND:
            found = True
            choice.Delete(index)
            while choice.list.count(self.blank) > 0:
                choice.list.remove(self.blank)
            index = choice.FindString(self.blank)

        if found:
            choice.SetSelection(selection)
            util.createEvent(wx.EVT_CHOICE, choice, selection)

    def blankChoice(self, choice):
        """Makes the given choice box show a blank entry"""

        # check if there already is a blank
        index = choice.FindString(self.blank)
        if index == wx.NOT_FOUND:
            # add a blank entry to the list and make it selected
            choice.Insert(self.blank, 0)
            choice.SetSelection(0)
            if choice.list.count(self.blank) == 0:
                choice.list.insert(0, self.blank)
        else:
            choice.SetSelection(index)
        util.createEvent(wx.EVT_CHOICE, choice, choice.GetSelection())

    def isChoiceBlank(self, choice):
        """Returns True if there is a blank in the choice box list, 
        otherwise False
        
        """

        index = choice.FindString(self.blank)
        if index == wx.NOT_FOUND:
            return False
        else:
            return True

    def isAnyChoicesBlank(self):
        """Return True if any of the choice boxes in the constraint have a
        blank entry, returns False if all choice boxes have no blank entries
        
        """

        for eachChoiceBox in self.choiceBoxes:
            if self.isChoiceBlank(eachChoiceBox):
                return True

        return False

    def isAllChoicesBlank(self):
        """Return True if all of the choice boxes in the constraint have a
        blank entry, returns False if any choice boxes have no blank entries
        
        """

        for eachChoiceBox in self.choiceBoxes:
            if self.isChoiceBlank(eachChoiceBox) is False:
                return False

        return True

    def setError(self, bool):
        """Change the color of the panel for errors"""

        if bool:
            # there is an error
            self.SetBackgroundColour(misc.errorColor)
            self.error = True
        else:
            # there is no error, revert back to default color
            self.SetBackgroundColour(wx.NullColour)
            self.error = False
        self.Refresh()

    def addToFile(self):
        """Adds this panel's line to the file"""

        if self.inFile is True:
            # nothing to do
            return

        self.inFile = True
        self.GetParent().addToFile(self)

    def removeFromFile(self):
        """Removes this panel's line from the file"""

        if self.inFile is False:
            # nothing to do
            return

        self.inFile = False
        self.GetParent().removeFromFile(self)

    def beingUsed(self):
        """Event when something is modified in the constraint"""
        util.callAncestorFunction(self, 'hasChanged')

class ConstraintPanel_AC():
    """A base class that adds widgets required for constraints that are across 
    classes (uses liability classes).  This class is an 'addon' to the
    ConstraintsPanel, and should always be inherited from a class derived
    from ConstraintPanel.
    
    """

    # the label next to the class choice boxes
    classLbl = 'Class:'

    def __init__(self):
        self.create_ConstraintPanel_AC()
        self.intialize_ConstraintPanel_AC()
        if self.isNew:
            self.reset()

    def create_ConstraintPanel_AC(self):
        """Create widgets needed to deal with liability classes"""

        # two labels indentifying the class choice boxes
        self.classLabel1 = wx.StaticText(self, -1, self.classLbl)
        self.classLabel2 = wx.StaticText(self, -1, self.classLbl)

        # liability class choice boxes
        liabilityChoices = self.getLiabilityChoices()
        self.classChoice1 = wx.Choice(self, -1, choices=liabilityChoices)
        self.classChoice2 = wx.Choice(self, -1, choices=liabilityChoices)
        self.choiceBoxes.append(self.classChoice1)
        self.choiceBoxes.append(self.classChoice2)

    def intialize_ConstraintPanel_AC(self):
        """Intialize widgets based on file contents"""

        # other members of choice box
        liabilityChoices = self.getLiabilityChoices()
        self.classChoice1.list = liabilityChoices[:]
        self.classChoice2.list = liabilityChoices[:]
        # the dictionaries need the blank entries to work
        liabilityChoices.append(self.blank)
        dictionary = dict(zip(liabilityChoices, liabilityChoices))
        self.classChoice1.tagToDescr = dictionary
        self.classChoice1.descrToTag = dictionary
        self.classChoice2.tagToDescr = dictionary
        self.classChoice2.descrToTag = dictionary

        # set responsibilities
        # the first token specifying liability class is always
        # right before the operator of the constraint
        for i in range(0, len(self.line)):
            if self.line[i].type == types.OPERATOR:
                j = i
                break
        self.classChoice1.charge = self.line[j-1]

        # the next token specifying liability class is always last
        self.classChoice2.charge = self.line[-1]

        # intialize the values
        index = int(self.classChoice1.charge.str) - 1
        self.classChoice1.SetSelection(index)
        index = int(self.classChoice2.charge.str) - 1
        self.classChoice2.SetSelection(index)

        # bind events
        self.classChoice1.Bind(wx.EVT_CHOICE, self.onChoice, self.classChoice1)
        self.classChoice2.Bind(wx.EVT_CHOICE, self.onChoice, self.classChoice2)

        # need to know when the liabilty class changes
        util.bindLiabilityClassChange(self, self.onLiabilityChange)

    def getLiabilityChoices(self):
        """Figure out how many liability classes there are"""

        # this returns a list of strings,
        # from '1' to the number of liability classes
        tag = KelvinInfo.liabilityClassEntry[0]
        if self.fm.isPresent(tag):
            line = self.fm.getLines(tag)[0]
            try:
                list =  [str(x) for x in range(1, int(line[1].str)+1)]
            except ValueError, e:
                return []
            else:
                return list
        else:
            # no number of liability classes specified, so use default
            return [KelvinInfo.defaultLiabilityClass]

    def onLiabilityChange(self, event):
        """Event when the liabity class changes"""

        # basically update the choice list with the current number
        # of liability classes
        self.updateLiabilityChoice(self.classChoice1)
        self.updateLiabilityChoice(self.classChoice2)

    def updateLiabilityChoice(self, choice):
        """Given a choice box, update it with the current number of
        liability classes
        
        """

        if self.isChoiceBlank(choice):
            isBlank = True
        else:
            isBlank = False

        try:
            # get the current value selected
            curValue = int(choice.GetStringSelection())
        except ValueError, e:
            # the only way this could happen is if the current
            # selection for the choice is a blank
            curValue = 0

        # get current liabilities based on file contents
        curLiabilities = self.getLiabilityChoices()
        if len(curLiabilities) == 0:
            return

        if curValue > int(curLiabilities[-1]):
            # the list will go up to the current selection
            curLiabilities = [str(x) for x in range(1, curValue+1)]

        if isBlank:
            curLiabilities.insert(0, self.blank)

        if len(choice.list) < len(curLiabilities):
            # append the extra values onto the choice box
            j = len(choice.list)
            for i in range(j, len(curLiabilities)):
                choice.Insert(curLiabilities[i], i)
                choice.list.append(curLiabilities[i])
        else:
            # remove extra entries from the choice box
            j = len(choice.list)
            k = len(curLiabilities)
            for i in range(j-1, k-1, -1):
                choice.Delete(i)
                choice.list.remove(str(i))

        choice.list = curLiabilities
        dictionary = dict(zip(curLiabilities, curLiabilities))
        dictionary[self.blank] = self.blank
        choice.tagToDescr = dictionary
        choice.descrToTag = dictionary

        # set selection to previous setting
        if isBlank:
            curValue = self.blank
        index = choice.list.index(str(curValue))
        choice.SetSelection(index)

class ConstraintLinePanel(wx.Panel):
    """A class for that represents all the constraints on one line.  A panel 
    that is simply a container for single constraints.

    """

    def __init__(self, parent, fileManager, constraintType, line=None):
        """
        Keyword arguments:
        constraintType -- a class that is the type of constraint panel that 
                          this panel uses

        """

        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.line = line
        self.constraintType = constraintType
        self.error = False
        self.isEnabled = True

        # figure out if the line for this panel is in the file or not
        if line is None:
            self.inFile = False
        else:
            self.inFile = True

        # a list of all the constraint panels in this panel
        self.constraintPanels = []

        self.createPanel_ConstraintLinePanel()
        self.arrangePanel_ConstraintLinePanel()

    def createPanel_ConstraintLinePanel(self):
        """Create subpanels for the panel"""

        if self.line is not None:
            # for each constraint, make a subpanel of the correct type
            line = []
            for eachToken in self.line:
                if eachToken.type == types.SEMICOLON:
                    constraintPanel = self.constraintType(self, self.fm, line)
                    self.constraintPanels.append(constraintPanel)
                    line = []
                else:
                    line.append(eachToken)
            constraintPanel = self.constraintType(self, self.fm, line)
            self.constraintPanels.append(constraintPanel)
        else:
            # create at least one constraint of the type specified
            constraintPanel = self.constraintType(self, self.fm)
            self.constraintPanels.append(constraintPanel)
            self.line = []

    def arrangePanel_ConstraintLinePanel(self):
        """Arrange the subpanels"""

        # place subpanels along a row
        self.constraintSizer = wx.BoxSizer(wx.HORIZONTAL)

        style = wx.ALL
        for eachPanel in self.constraintPanels:
            self.constraintSizer.Add(eachPanel, 0, style, 1)

        self.SetSizer(self.constraintSizer)

    def addConstraint(self, constraintPanel):
        """Adds a whole new constraint line to the parent panel"""

        util.callAncestorFunction(self, 'addConstraint', self)

    def addConstraintToLine(self):
        """Adds a new constraint to this panel"""

        newConstraintPanel = self.constraintType(self, self.fm)
        index = len(self.constraintPanels)
        self.constraintPanels.insert(index, newConstraintPanel)
        self.constraintSizer.Add(newConstraintPanel, 0, wx.ALL, 1)

        self.Fit()
        self.GetParent().Fit()
        self.GetParent().GetParent().Fit()
        self.reSize()

        if self.isEnabled is False:
            newConstraintPanel.Disable()

    def removeConstraint(self, constraintPanel):
        """Remove the given constraint panel from this panel"""

        self.constraintPanels.remove(constraintPanel)
        self.constraintSizer.Remove(constraintPanel)
        constraintPanel.Destroy()
        self.Layout()
        self.reSize()

        if len(self.constraintPanels) == 0:
            util.callAncestorFunction(self, 'removeConstraint', self)

    def reset(self):
        """Makes one new blank constraint in this line"""

        # Note: this should only be called when there are no constraints left
        # in the constraintLinePanel, which is the case right now

        # add one new one
        self.addConstraintToLine()

    def reSize(self):
        """Activate an ancestor's FitInside() function in order for scroll bars
        to function properly
        
        """

        # should be called after the panel removes or adds things
        # need to correctly update the scroll bars
        util.callAncestorTypeFunction(self, ScrolledPanel, 'FitInside')

    def setError(self, bool):
        """Set the error states of the subpanels"""

        self.error = bool
        for eachPanel in self.constraintPanels:
            eachPanel.setError(bool)

    def Enable(self, bool=True):
        """Enables or disables panel, extra code to deal with constraint list"""

        self.isEnabled = bool
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
            eachPanel.setError(False)

    def disableConstraintPanels(self):
        """Disables the list of constraint panels"""

        for eachPanel in self.constraintPanels[:]:
            if eachPanel.inFile:
                # don't disable it, just set it to an error
                eachPanel.setError(True)
            else:
                if eachPanel.isAllChoicesBlank():
                    if len(self.constraintPanels) == 1:
                        # the last blank one, don't remove it, just disable
                        eachPanel.Disable()
                    else:
                        # not the last blank one, just remove it
                        self.removeConstraint(eachPanel)
                else:
                    # not in file, but partially filled in
                    # just set it to an error
                    eachPanel.setError(True)

    def isAllChoicesBlank(self):
        """Returns true if all choice boxes in this constraint line are blank,
        returns False otherwise
        
        """

        for eachPanel in self.constraintPanels:
            if eachPanel.isAllChoicesBlank() is False:
                return False
        return True

    def Destroy(self):
        """Event when this panel is destroyed"""

        # remove all the constraints this panel has
        for eachPanel in self.constraintPanels:
            self.removeConstraint(eachPanel)
        wx.Panel.Destroy(self)

    def addToFile(self, constraint):
        """Given a constraint, add the line the constraint is responsible 
        for to the file
        
        """

        # add a semicolon if it isn't the first constraint in this line
        if len(self.line) != 0:
            self.fm.append(self.line, ';')

        for eachToken in constraint.line:
            self.line.append(eachToken)

        # add this line to the file if necessary
        if self.inFile is False:
            self.fm.addLine(self.line)
            self.inFile = True

    def removeFromFile(self, constraint):
        """Given a constraint, remove the line the constraint is responsible 
        for from the file.
        
        """

        # figure out the starting index of where the constraint's line starts
        # and how many items to remove
        index = self.line.index(constraint.line[0])
        length = len(constraint.line)

        # if it's the first constraint and there are more on the line,
        # loop through one more item to take out a semicolon
        if index == 0 and len(self.line) > len(constraint.line):
            length += 1

        # if it's not the first constraint on the line, go back one more
        # place in order to remove a semicolon, and loop through one more item
        if index != 0:
            index -= 1
            length += 1

        for i in range(index+length-1, index-1, -1):
            self.line.pop(i)

        # remove entire line from file if nothing is left
        if len(self.line) == 0:
            self.fm.remove(self.line)
            self.inFile = False
