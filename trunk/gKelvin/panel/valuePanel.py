"""The panel to show exactly one single point or ranged value"""

import wx
from config import KelvinInfo
from config import misc
from fileManager import types
from fileManager.token import Token
from utility import util

class ValuePanel(wx.Panel):
    """The panel to show one individual range or single point value(s)"""

    beginLbl = "Begin:  "
    endLbl = "End:"
    stepLbl = "Increment:"
    pointLbl = "Values:"
    boxSize = (60, -1)
    valueTypeList = ['Value', 'Range']

    def __init__(self, parent, fileManager, parameterList, line, inFile):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.line = line
        self.paramChoiceList = [x.description for x in parameterList]
        self.tagToDescr = KelvinInfo.variableNumberTagToDescription
        self.descrToTag = KelvinInfo.variableNumberDescriptionToTag

        # indicates whether the line that the value panel is
        # responsible for is currently inside the file
        self.inFile = inFile

        self.createPanel(parameterList)
        self.arrangePanel()
        self.initialize()

    def createPanel(self, parameterList):
        """Create widgets for the panel"""

        # parameter choice box
        paramList = [x.description for x in parameterList]
        self.parameterChoice = wx.Choice(self, -1, choices=paramList)

        # choice box to choose range or point value
        self.valueTypeChoice = wx.Choice(self, -1, choices=self.valueTypeList)

        # labels to identify what the boxes are
        self.pointLabel = wx.StaticText(self, -1, self.pointLbl)
        self.beginLabel = wx.StaticText(self, -1, self.beginLbl)
        self.endLabel = wx.StaticText(self, -1, self.endLbl)
        self.stepLabel = wx.StaticText(self, -1, self.stepLbl)

        # boxes to specify values
        self.pointBox = wx.TextCtrl(self, -1, size=self.boxSize)
        self.beginBox = wx.TextCtrl(self, -1, size=self.boxSize)
        self.endBox = wx.TextCtrl(self, -1, size=self.boxSize)
        self.stepBox = wx.TextCtrl(self, -1, size=self.boxSize)

        # add and remove buttons
        buttonSize = misc.addRemoveButtonSize
        self.addButton = wx.Button(self, -1, label="+", size=buttonSize)
        self.removeButton = wx.Button(self, -1, label="-", size=buttonSize)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        # parameter choice box
        mainSizer.Add(self.parameterChoice, 0, wx.ALIGN_CENTER_VERTICAL | 
                                                                    wx.ALL, 3)

        # type choice box
        mainSizer.Add(self.valueTypeChoice, 0, wx.ALIGN_CENTER_VERTICAL | 
                                                                    wx.ALL, 3)
        mainSizer.Add((10, -1))

        # labels and textboxes
        mainSizer.Add(self.pointLabel,0,wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        mainSizer.Add(self.pointBox, 1, wx.EXPAND | wx.ALIGN_CENTER_VERTICAL | 
                                        wx.RIGHT | wx.UP | wx.DOWN, 5)

        mainSizer.Add(self.beginLabel,0,wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        mainSizer.Add(self.beginBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT
                                                        | wx.UP | wx.DOWN, 5)

        mainSizer.Add(self.endLabel, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 3)
        mainSizer.Add(self.endBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT 
                                                        | wx.UP | wx.DOWN, 5)

        mainSizer.Add(self.stepLabel, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL,3)
        mainSizer.Add(self.stepBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.RIGHT 
                                                        | wx.UP | wx.DOWN, 5)

        # add and remove buttons
        mainSizer.Add(self.addButton,0,wx.ALIGN_CENTER_VERTICAL | wx.LEFT,10)
        mainSizer.Add(self.removeButton, 0, wx.ALIGN_CENTER_VERTICAL | 
                                                       wx.LEFT | wx.RIGHT, 3)

        self.SetSizer(mainSizer)
        self.Fit()

        # need to set the minimum size as the size when the panel is longest,
        # a.k.a. whenever it's set to show a range value
        self.setPanelToRange()
        size = self.GetSize()
        self.SetMinSize(size)

    def initialize(self):
        """Initialize widget values based on file contents"""

        if self.line is not None:
            # given the tag from the line, find the right description 
            # for the tag, then initialize the parameter choice
            description = self.tagToDescr[self.line[0].str]
            index = self.paramChoiceList.index(description)
            self.parameterChoice.SetSelection(index)
            self.parameterChoice.charge = self.line[0]

            # initialize the type and value textboxes based on if the
            # line describes point values or ranged values
            # assume that the line cannot be describing both
            if self.isRange():
                self.initializeRange()
            else:
                self.initializePointValue()
        else:
            # no initial line given, set it to a point value
            self.initializeNoLine()

        # need to change the labels on the genotypes in the parameter choice 
        # based on the current trait type
        self.initializeParameterChoice()

        # there may be errors with the values depending 
        # on other values of the file
        self.checkForErrors()

        # bind some events
        self.parameterChoice.Bind(wx.EVT_CHOICE,self.onParameterChoice,
                                          self.parameterChoice)
        self.valueTypeChoice.Bind(wx.EVT_CHOICE,self.onValueTypeChoice,
                                          self.valueTypeChoice)
        self.addButton.Bind(wx.EVT_BUTTON,self.onAddButton,
                                          self.addButton)
        self.removeButton.Bind(wx.EVT_BUTTON,self.onRemoveButton,
                                          self.removeButton)
        self.beginBox.Bind(wx.EVT_TEXT,self.onTextEntry, self.beginBox)
        self.endBox.Bind(wx.EVT_TEXT,self.onTextEntry, self.endBox)
        self.stepBox.Bind(wx.EVT_TEXT,self.onTextEntry, self.stepBox)
        self.pointBox.Bind(wx.EVT_TEXT,self.onTextEntry, self.pointBox)

    def initializeNoLine(self):
        """Perform initialization when no initial line is given"""

        # this happens when the user clicks the add button to add a new value
        # make it show a point value at first 
        self.inFile = False
        self.valueTypeChoice.SetSelection(0)
        self.setPanelToValue()

        # create a mock line, just in case it ever needs 
        # to be inserted into the file
        tag = self.descrToTag[self.paramChoiceList[0]]
        self.line = self.fm.createLine(tag, '')
        self.parameterChoice.charge = self.line[0]
        self.pointBox.charge = self.line[1]
        self.beginBox.charge = Token(str='')
        self.endBox.charge = Token(str='')
        self.stepBox.charge = Token(str='')

    def initializeRange(self):
        """Perform initialization required for a ranged value line"""

        # initialize values
        self.valueTypeChoice.SetSelection(1)
        self.beginBox.ChangeValue(self.line[1].str)
        self.endBox.ChangeValue(self.line[2].str)
        self.stepBox.ChangeValue(self.line[3].str)

        # set responsibilities
        self.beginBox.charge = self.line[1]
        self.endBox.charge = self.line[2]
        self.stepBox.charge = self.line[3]
        self.pointBox.charge = Token(str='')

        self.beginLabel.Show(True)
        self.endLabel.Show(True)
        self.stepLabel.Show(True)
        self.beginBox.Show(True)
        self.endBox.Show(True)
        self.stepBox.Show(True)

        self.pointLabel.Hide()
        self.pointBox.Hide()

        self.setPanelToRange()

    def initializePointValue(self):
        """Perform initialization required for a point value line"""

        # so for a point value, multiple point values may be specified,
        # with a semicolon separating each value, for example:
        # Th 0.1; 0.2; 0.3
        # to deal with this, take all the tokens after the tag, and join
        # them all into one token, so that when text is entered in the
        # textbox, only need to update one token, let the user have the
        # responsibility of putting semicolons between values
        self.valueTypeChoice.SetSelection(0)

        # concatenate the whole line except the tag, 
        # and also remove those tokens
        s = ''
        for eachToken in self.line[1:]:
            s += eachToken.str
            if eachToken.type == types.SEMICOLON:
                s += ' '
            self.line.remove(eachToken)

        # create one token out of the concatenated line
        self.pointBox.charge = self.fm.append(self.line, s)
        self.pointBox.ChangeValue(s)

        # set other responsibilities
        self.beginBox.charge = Token('')
        self.endBox.charge = Token('')
        self.stepBox.charge = Token('')

        self.pointLabel.Show(True)
        self.pointBox.Show(True)

        self.beginLabel.Hide()
        self.endLabel.Hide()
        self.stepLabel.Hide()
        self.beginBox.Hide()
        self.endBox.Hide()
        self.stepBox.Hide()

        self.setPanelToValue()

    def initializeParameterChoice(self):
        """Based on the trait type, set the labels for the parameter choice"""

        self.updateLabels()

    def onParameterChoice(self, event):
        """Event when the parameter choice changes"""

        description = self.parameterChoice.GetStringSelection()
        tag = self.descrToTag[description]
        self.parameterChoice.charge.update(tag)
        self.hasChanged()

        # since the tag might have changed, remove then add this file back
        # so the line with the new tag gets inserted in an ok place and
        # the comment in the file doesn't lie
        if self.inFile:
            self.fm.remove(self.line)
            self.fm.addLine(self.line)

    def onTextEntry(self, event):
        """Event when something is typed into a text box"""

        textbox = event.GetEventObject()
        textbox.charge.update(textbox.GetValue())
        self.beingUsed()

    def onValueTypeChoice(self, event):
        """Event when the type of value (point value or range) changes"""

        if self.valueTypeChoice.GetSelection() == \
                                          self.valueTypeList.index('Value'):
            if self.isRange():
                # point value is chosen 
                self.setPanelToValue()
                self.setLineToValue()
        else:
            if not self.isRange():
                # range value must have been chosen
                self.setPanelToRange()
                self.setLineToRange()
        self.beingUsed()

    def onAddButton(self, event):
        """Event when the add button is clicked"""

        if hasattr(self.GetParent(), 'addValue'):
            self.GetParent().addValue(self)
        self.hasChanged()

    def onRemoveButton(self, event):
        """Event when the remove button is clicked"""

        self.removeLine()
        self.hasChanged()
        if hasattr(self.GetParent(), 'removeValue'):
            self.GetParent().removeValue(self)

    def onAnalysisChange(self):
        """When the analysis type changes, the current value may be invalid"""

        self.updateLabels()
        self.checkForErrors()

    def onTraitChange(self):
        """When the trait type changes, the current value may be invalid"""

        self.updateLabels()
        self.checkForErrors()

    def setPanelToValue(self):
        """Makes the panel show only the wigets needed to show a point value"""

        self.pointLabel.Show(True)
        self.pointBox.Show(True)

        self.beginLabel.Hide()
        self.endLabel.Hide()
        self.stepLabel.Hide()
        self.beginBox.Hide()
        self.endBox.Hide()
        self.stepBox.Hide()

        self.Layout()
        self.Fit()

    def setPanelToRange(self):
        """Makes the panel show only the wigets needed to show a range value"""

        self.beginLabel.Show(True)
        self.endLabel.Show(True)
        self.stepLabel.Show(True)
        self.beginBox.Show(True)
        self.endBox.Show(True)
        self.stepBox.Show(True)

        self.pointLabel.Hide()
        self.pointBox.Hide()

        self.Layout()
        self.Fit()

    def setLineToValue(self):
        """Changes the format of the line to a value point"""

        self.line.remove(self.beginBox.charge)
        self.line.remove(self.endBox.charge)
        self.line.remove(self.stepBox.charge)
        self.line.append(self.pointBox.charge)

    def setLineToRange(self):
        """Changes the format of the line to a range value"""

        self.line.remove(self.pointBox.charge)
        self.line.append(self.beginBox.charge)
        self.line.append(self.endBox.charge)
        self.line.append(self.stepBox.charge)

    def isRange(self):
        """Check if the line that the panel is responsible for is a range or a 
        single point value.
        
        """

        # the line is a range if three numbers are encountered without hitting
        # a semicolon, skip the first token since that's the tag
        hit = 0
        for eachToken in self.line[1:]:
            if eachToken.type == types.SEMICOLON:
                hit = 0
            else:
                hit += 1
            if hit == 3:
                return True

        return False
            
    def beingUsed(self):
        """Event when something in this panel is being modified."""

        self.insertLine()
        self.hasChanged()
        self.checkForErrors()

    def hasChanged(self):
        """Event when something changes or is clicked in the panel"""

        # the difference between this function and beingUsed() is that
        # this is called when the add and remove buttons are clicked
        util.callAncestorFunction(self, 'hasChanged')
        util.callAncestorFunction(self, 'onValuesChange', None)

    def reset(self):
        """Removes the line from the file and resets values"""

        self.removeLine()

        # set parameter choice to first selection
        self.parameterChoice.SetSelection(0)
        tag = self.descrToTag[self.paramChoiceList[0]]
        self.parameterChoice.charge.update(tag)

        # clear out all the text boxes
        self.pointBox.ChangeValue('')
        self.beginBox.ChangeValue('')
        self.endBox.ChangeValue('')
        self.stepBox.ChangeValue('')

        self.pointBox.charge.update('')
        self.beginBox.charge.update('')
        self.endBox.charge.update('')
        self.stepBox.charge.update('')

    def updateLabels(self):
        """Update parameter labels based on current trait type (DT or QT)"""

        # basically add a prefix before genotypes, either 'penetrance' 
        # for DT or 'mean' for QT, or '' for none

        curSelection = self.parameterChoice.GetSelection()

        if self.fm.isPresent('DT'):
            # using DT need to add 'Penetrance' in front of genotypes
            self.prefix = KelvinInfo.genotypePrefix[0]
        elif self.fm.isPresent('QT'):
            # using QT, need to add 'Mean' in front of genotypes
            self.prefix = KelvinInfo.genotypePrefix[1]
        else:
            # a trait type was not specified
            self.prefix = ''

        # add the prefix in
        for eachGenotype in KelvinInfo.genotypes:
            for i in range(0, len(self.paramChoiceList)):
                if self.paramChoiceList[i][-len(eachGenotype):]==eachGenotype:
                    self.paramChoiceList[i] = self.prefix + eachGenotype
                    self.parameterChoice.SetString(i, self.paramChoiceList[i])

        # make sure the previous selection is still selected
        self.parameterChoice.SetSelection(curSelection)

    def updateParameterList(self, paramList):
        """Given a new list of descriptions for parameters, set the 
        parameter choice box to the new choices.
        
        """

        # get the text of the current selection.
        # if the current selection is a genotype, need to strip
        # the selection down to the exact genotype, i.e. 'mean DD' -> 'DD'
        self.paramChoiceList = paramList
        curText = self.parameterChoice.GetStringSelection()
        for eachGenotype in KelvinInfo.genotypes:
            if curText[-len(eachGenotype):] == eachGenotype:
                curText = eachGenotype
                break
        if curText not in paramList:
            # this means the current selection is invalid for the new
            # parameter lists, so delete this entire value
            self.onRemoveButton(None)
        else:
            # re-create the choice selections, and set it back 
            # to the old selection
            self.parameterChoice.Clear()
            for eachItem in paramList:
                self.parameterChoice.Append(eachItem)

            # set the choice box to what was previously selected.
            self.parameterChoice.SetSelection(paramList.index(curText))

    def checkForErrors(self):
        """Checks if the current value should be there or not"""

        # ---------------------------------------------------------------------
        # what follows is a list of functions that are called that represents 
        # one test each. a function here should return true if there is an 
        # error, and false if there is not an error.  these functions may alter
        # other things depending on what they are checking
        # ---------------------------------------------------------------------

        def genotypeSelectedCheck(self):
            # The current selected parameter should not be a genotype 
            # unless a trait type (DT or QT) is specified.
            # DT or QT doesn't have to be in the file if MM or AM is
            # specified (Marker to marker analysis).

            # if trait type is specified, having a genotype selected is ok
            if self.hasTraitType():
                return False

            # check if current description is a genotype
            curDescription = self.parameterChoice.GetStringSelection()
            for eachItem in KelvinInfo.genotypes:
                if curDescription[-len(eachItem):] == eachItem:
                    return True

            # genotype not found, no error
            return False

        def rangeValuesCheck(self):
            # If using a range, check if increment is zero or the lower
            # limit is higher than the upper limit

            curValueType = self.valueTypeChoice.GetStringSelection()
            if curValueType == self.valueTypeList[0]:
                # using a value point instead of range, nothing to check
                return False

            try:
                if float(self.stepBox.GetValue()) == 0.0:
                    return True
                begin = float(self.beginBox.GetValue())
                end = float(self.endBox.GetValue())
                if begin >= end:
                    return True
            except ValueError, e:
                # this will happen if the text cannot be converted to a number,
                # for instance if there are letters in it
                return True

        def pointValuesCheck(self):
            # check to make sure that all the values in the texbox can be 
            # validly parsed into numbers

            curValueType = self.valueTypeChoice.GetStringSelection()
            if curValueType == self.valueTypeList[1]:
                # using a range instead of point value, nothing to check
                return False

            s = self.pointBox.GetValue()
            if s == '':
                return True
            s = s.replace(' ', '')
            for eachItem in s.split(';'):
                try:
                    if eachItem != '':
                        float(eachItem)
                except ValueError, e:
                    return True
            return False

        # ---------------------------------------------------------------------
        # list of error functions ends here, this is where they are called
        # ---------------------------------------------------------------------

        # a list of all the previous error functions so that they 
        # can be easily called in a row
        allErrorFunctions = [genotypeSelectedCheck, rangeValuesCheck, 
                             pointValuesCheck]

        for eachFunction in allErrorFunctions:
            if eachFunction(self):
                self.setError(True)
                return

        self.setError(False)

    def setError(self, bool):
        """Change the color of the panel for errors"""

        if bool:
            # there is an error
            self.SetBackgroundColour(misc.errorColor)
        else:
            # there is no error, revert back
            self.SetBackgroundColour(wx.NullColour)
        self.Refresh()

    def setTag(self, tag):
        """Change the parameter to the given tag"""

        # act as if the user clicked the parameter choice box

        # set choice box to right selection
        description = self.tagToDescr[tag]
        index = self.paramChoiceList.index(description)
        self.parameterChoice.SetSelection(index)
        description = self.parameterChoice.GetStringSelection()

        # note: i'm not using util.createEvent because when i do that,
        # the choice box is always set to the first choice, which
        # doesn't make sense since i just changed it right above.
        # maybe a bug? who knows...
        self.onParameterChoice(None)

    def hasTraitType(self):
        """Checks if the file is currently specifying a trait type or not"""

        for eachTag in KelvinInfo.traitTypeTags:
            if self.fm.isPresent(eachTag):
                return True
        return False

    def insertLine(self):
        """Inserts the line this value panel is responsible for into the file"""

        if not self.inFile:
            self.inFile = True
            self.fm.addLine(self.line)

    def removeLine(self):
        """Remove the line from the file"""

        if self.inFile:
            self.inFile = False
            self.fm.remove(self.line)
