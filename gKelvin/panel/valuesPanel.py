"""The panel to show all single point and ranged values"""

import wx
from config import KelvinInfo
import copy
from panel.scrolledPanel import ScrolledPanel
from panel.valuePanel import ValuePanel
from fileManager import types
from utility import util

class ValuesPanel(wx.Panel):
    """The panel to show all single point and ranged values"""

    paramLbl = "Parameter"
    typeLbl = "Type"
    valuesLbl = 'Note: Use semicolons to separate single values'

    valuePanelBorder = 2
    borderAmt = 5

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        # a dummy genotype line that is used when MM or AM is used
        self.dummyLine = self.fm.createLine('DD', '0.0', '0.9', '0.1')
        self.hasDummy = False

        # get all the lines in the file that are values,
        # lines with mutiple definitions are split into separate lines
        self.lines = self.getVariableNumberLines()

        # a list of the individual value panels
        self.valuePanels = []

        # get the list of valid parameters of values
        self.parameterList = self.getParameterList()

        # keep track of different states depending on what the 
        # current analysis type and trait type is.
        # each state is a list of valuePanels
        self.states = []

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Create widgets for the panel"""

        # headings
        self.paramLabel = wx.StaticText(self, -1, self.paramLbl)
        self.typeLabel = wx.StaticText(self, -1, self.typeLbl)
        self.valuesLabel = wx.StaticText(self, -1, self.valuesLbl)

        # the value panels
        tags = [x.tag for x in self.parameterList]
        for eachLine in self.lines[:]:
            if eachLine[0].str in tags:
                self.valuePanels.append(ValuePanel(self, self.fm, 
                                          self.parameterList, eachLine, True))
            else:
                # the line shouldn't be there in the first place,
                # so remove it from the file
                self.fm.remove(eachLine)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        pad1, pad2 = self.getPadding()

        labelSizer = wx.BoxSizer(wx.HORIZONTAL)
        labelSizer.Add(self.paramLabel, 0)
        labelSizer.Add((pad1, 0))
        labelSizer.Add(self.typeLabel, 0)
        labelSizer.Add((pad2, 0))
        labelSizer.Add(self.valuesLabel, 0)

        mainSizer.Add(labelSizer, 0, wx.ALL, self.borderAmt)

        for eachValuePanel in self.valuePanels:
            mainSizer.Add(eachValuePanel, 1, wx.EXPAND | wx.TOP | wx.RIGHT, 
                                                self.valuePanelBorder)

        self.sizer = mainSizer
        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widget values based on file contents"""

        # check whether the dummy line needs to be inserted or not. don't
        # worry about inserting multiple dummy lines, since at this point
        # any values that shouldn't be there should have been removed.
        if self.fm.isPresent('DT') or self.fm.isPresent('QT'):
            self.hasDummy = False
        else:
            self.fm.addLine(self.dummyLine)
            self.hasDummy = True

        self.initializeStates()

        # this panel needs to know when the trait type changes
        util.bindTraitTypeChange(self, self.onTraitChange)
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)

    def initializeStates(self):
        """Create the states"""

        def getState(analysis, trait, threshold=False):
            """Given the analysis and trait type, return the default state"""

            # the 'threshold' argument is whether you want to use the
            # threshold trait type or not

            lines = self.fm.getValueDefaults(analysis=analysis, trait=trait)
            result = []
            parameterList = self.getParameterList(analysis, trait)
            for eachLine in lines:
                panel = ValuePanel(self, self.fm, parameterList, eachLine,False)
                panel.Hide()
                panel.old_inFile = True
                result.append(panel)

            # if using threshold, add the trait type (TT) line to the state
            if threshold:
                line = copy.deepcopy(self.fm.createLine(*KelvinInfo.defaultTT))
                panel = ValuePanel(self, self.fm, parameterList, line, False)
                panel.Hide()
                panel.old_inFile = True
                result.append(panel)
            return result

        def getStateNoTrait(analysis):
            """Given the analysis type, return the default state"""

            # this function is used when using AM or MM, in which there is
            # no trait type used at all

            lines = self.fm.getValueDefaults(analysis=analysis, trait='qt')
            result = []
            parameterList = self.getParameterList(analysis, 'qt')
            for eachLine in lines:
                if KelvinInfo.genotypes.count(eachLine[0].str) > 0:
                    continue
                panel = ValuePanel(self, self.fm, parameterList, eachLine,False)
                panel.Hide()
                panel.old_inFile = True
                result.append(panel)

            return result

        # basically get the lines that are values from different defaults
        self.states.append(getState(analysis='tp', trait='dt'))
        self.states.append(getState(analysis='tp', trait='qt'))
        self.states.append(getState(analysis='mp', trait='dt'))
        self.states.append(getState(analysis='mp', trait='qt'))
        self.states.append(getState(analysis='tp', trait='qt', threshold=True))
        self.states.append(getState(analysis='mp', trait='qt', threshold=True))

        # there are also two more states to consider. when MM or AM is used,
        # there is no trait type, so that's two more states. make those states
        # by using the previous states, but removing all genotype values
        self.states.append(getStateNoTrait(analysis='tp'))
        self.states.append(getStateNoTrait(analysis='mp'))

        # to make the code easier to read, these are variables that are used
        # to index self.states
        self.tp_dt = 0
        self.tp_qt = 1
        self.mp_dt = 2
        self.mp_qt = 3
        self.tp_th = 4    # threshold trait type
        self.mp_th = 5    # threshold trait type
        self.tp_am_mm = 6 # using either AM or MM tag
        self.mp_am_mm = 7 # using either AM or MM tag

        # get the list of values currently in the file as the current state
        curState = []
        for eachPanel in self.valuePanels:
            curState.append(eachPanel)

        # replace the appropriate state with the current state
        self.curStateIndex = self.getCurStateIndex()
        self.states[self.curStateIndex] = curState
        self.curState = self.states[self.curStateIndex]

    def getCurStateIndex(self):
        """Using the current information in the file, return the current
        state's correct index.
        
        """

        tp = self.fm.isPresent('TP')
        ss = self.fm.isPresent('SS')
        sa = self.fm.isPresent('SA')
        dt = self.fm.isPresent('DT')
        qt = self.fm.isPresent('QT')
        tt = self.fm.trait_threshold

        if tp and dt:
            return self.tp_dt
        elif tp and qt:
            if tt:
                return self.tp_th
            else:
                return self.tp_qt
        elif (ss or sa) and dt:
                return self.mp_dt
        elif (ss or sa) and qt:
            if tt:
                return self.mp_th
            else:
                return self.mp_qt
        elif tp and (not dt and not qt):
            return self.tp_am_mm
        elif (ss or sa) and (not dt and not qt):
            return self.mp_am_mm

    def switchStates(self, oldState=None, newState=None):
        """Given two indices of states, switch those states around"""

        if oldState == None:
            oldState = self.curStateIndex
        if newState == None:
            newState = self.getCurStateIndex()

        # if switching to the same state, nothing to do
        if oldState == newState:
            return

        util.callAncestorFunction(self, 'hasChanged')
        self.saveState(oldState)
        self.restoreState(newState)
        self.sizer.Layout()
        self.Fit()

        self.curStateIndex = newState
        self.curState = self.states[self.curStateIndex]
        self.valuePanels = self.states[self.curStateIndex]

    def saveState(self, index=None):
        """Given an index, saves the current state in that state"""

        if index == None:
            index = self.curStateIndex

        curState = []
        for eachPanel in self.valuePanels:
            eachPanel.Hide()
            eachPanel.old_inFile = eachPanel.inFile
            eachPanel.removeLine()
            self.sizer.Remove(eachPanel)
            curState.append(eachPanel)

        self.states[index] = curState

    def restoreState(self, index=None):
        """Given an index of a state, restores the state"""

        if index == None:
            index = self.getCurStateIndex()

        for eachPanel in self.states[index]:
            eachPanel.Show()
            if eachPanel.old_inFile:
                eachPanel.insertLine()
            eachPanel.inFile = eachPanel.old_inFile
            self.sizer.Add(eachPanel, 1, wx.EXPAND | wx.TOP | wx.RIGHT, 
                                                        self.valuePanelBorder)

    def getParameterList(self, analysis=None, trait=None):
        """Returns a list of the available parameters that are valid for a 
        given value of analysis type and trait type. If they are None, use 
        the current values of the file.
        
        """

        if analysis is None or trait is None:
            tp = self.fm.isPresent('TP')
            mp = self.fm.isPresent('SS') or self.fm.isPresent('SA')
            dt = self.fm.isPresent('DT')
            qt = self.fm.isPresent('QT')
        else:
            tp = analysis == 'tp'
            mp = analysis == 'mp'
            dt = trait == 'dt'
            qt = trait == 'qt'

        result = []
        for eachEntry in KelvinInfo.variableNumberInfo:
            if eachEntry.analysisType:
                if tp and eachEntry.tp:
                    result.append(eachEntry)
                    continue
                if mp and eachEntry.mp:
                    result.append(eachEntry)
                    continue
            if eachEntry.traitType:
                if dt and eachEntry.dt:
                    result.append(eachEntry)
                    continue
                if qt and eachEntry.qt:
                    result.append(eachEntry)
                    continue
            if not eachEntry.analysisType and not eachEntry.traitType:
                result.append(eachEntry)
                continue

        return result

    def addValue(self, valuePanel, line=None):
        """Given a value panel that's already in the panel, make a new value 
        panel right after the one given
        
        """

        self.Freeze()

        inFile = line != None
        newValuePanel = ValuePanel(self,self.fm,self.parameterList,line,inFile)
        index = self.valuePanels.index(valuePanel)
        self.valuePanels.insert(index+1, newValuePanel)
        self.GetSizer().Insert(index+2, newValuePanel, 1, wx.EXPAND | wx.TOP | 
                               wx.RIGHT, self.valuePanelBorder)

        self.reSize()
        self.Thaw()

        return newValuePanel

    def removeValue(self, valuePanel):
        """Remove the given value panel"""

        if len(self.valuePanels) == 1:
            # don't remove the last one, just reset it
            self.valuePanels[0].reset()
        else:
            self.valuePanels.remove(valuePanel)
            self.GetSizer().Remove(valuePanel)
            valuePanel.Destroy()
        self.reSize()

    def reSize(self):
        """Keep going up the parent tree until you hit the scrolled panel"""

        # needs to be called when a value is added or removed.
        # need to update the scroll bars
        util.callAncestorTypeFunction(self, ScrolledPanel, 'FitInside')

    def getPadding(self):
        """Get the amount of padding needed between headings"""

        # the amount depends on the layout of a value panel, so get a
        # temporary one, get its size, and use those values for the padding
        tempValuePanel = ValuePanel(self,self.fm,self.parameterList,None,False)

        # first figure out the amount of padding between the parameter
        # label and the type label
        # need the width of the parameter label
        paramWidth, height1 = self.paramLabel.GetSizeTuple()

        # next need the distance from the beginning of the panel to
        # the beginning of the type choice
        x1, y1 = tempValuePanel.parameterChoice.GetPositionTuple()
        x2, y2 = tempValuePanel.valueTypeChoice.GetPositionTuple()
        width2 = x2 - x1

        pad1 = width2 - paramWidth

        # get distance from beginning to end of current headings
        # now need the padding between the type label and the values label
        typeWidth, height1 = self.typeLabel.GetSizeTuple()
        headingWidth = paramWidth + pad1 + typeWidth

        # get distance from the beginning of the value panel to the start
        # of the value box
        x2, y2 = tempValuePanel.pointBox.GetPositionTuple()
        width2 = x2 - x1

        # the padding for this heading
        pad2 = width2 - headingWidth


        tempValuePanel.Destroy()
        return (pad1, pad2)

    def getVariableNumberLines(self):
        """From the file, get all lines that deal with a variable number.
        There may be lines mixing single point values and range values, so 
        separate those lines out.
        
        """

        # first get the lines from the file
        lines = self.fm.getVariableNumberLines()

        result = []

        # go through each line from the file
        for eachLine in lines:
            # see if this line does not need to be expanded out
            # a line does not need to be expanded if it only
            # specifies one range, or it specifies only point values
            if not self.needsExpanding(eachLine):
                result.append(eachLine)
            else:
                splitLines = self.splitVariableNumberLine(eachLine)
                result.extend(splitLines)

                #remove the original line, and insert the expanded ones
                self.fm.remove(eachLine)
                for newLine in splitLines:
                    self.fm.addLine(newLine)

        return result

    def needsExpanding(self, line):
        """Check if a certain line needs to be expanded out into multiple 
        lines or not
        
        """

        # a line need to be expanded if it specifies at least one range 
        # and one point value, or it specifies more than one range

        # how many range values that were encountered
        range = 0

        # how many point values that were encountered
        point = 0

        # step through each token and count up how many range and point
        # values there are
        hit = 0
        for eachToken in line[1:]:
            if eachToken.type == types.SEMICOLON:
                if hit == 1:
                    point += 1
                if hit == 3:
                    range += 1
                hit = 0
            else:
                hit += 1

        # the last term isn't separated by a semicolon, so check it now
        if hit == 1:
            point += 1
        if hit == 3:
            range += 1

        # see if the line needs expanding
        if range > 1:
            return True
        if range > 0 and point > 0:
            return True

        # does not need expanding!
        return False

    def splitVariableNumberLine(self, line):
        """Given a line that needs to be expanded, split it into
        multple lines
        
        """

        # all range values need to be on a separate line
        # point values can all be on the same line at once

        result = []
        newLine = []

        # keep the tag
        tag = line[0]

        pointValueLine = [tag]
        length = len(line)

        for i in range(1, length):
            if line[i].type == types.SEMICOLON:
                continue
            if i == length-1:
                if line[i-1].type == types.SEMICOLON and \
                                           line[i].type == types.NUMBER:
                    pointValueLine.append(line[i])
                    continue
            if i == 1: 
                if line[i].type == types.NUMBER and \
                                           line[i+1].type == types.SEMICOLON:
                    pointValueLine.append(line[i])
                    pointValueLine.append(line[i+1])
                    continue
            if i > 1 and i < length-1: 
                if line[i-1].type==types.SEMICOLON and \
                                           line[i].type == types.NUMBER and \
                                           line[i+1].type == types.SEMICOLON:
                    pointValueLine.append(line[i])
                    pointValueLine.append(line[i+1])
                    continue
            if i == 1: 
                if line[i].type == types.NUMBER and \
                                           line[i+1].type == types.NUMBER:
                    # range value
                    newLine = [tag, line[i], line[i+1], line[i+2]]
                    result.append(newLine)
                    continue
            if i > 1 and i < length-2: 
                if line[i].type == types.NUMBER and \
                                           line[i+1].type == types.NUMBER and\
                                           line[i+2].type == types.NUMBER:
                    # range value
                    newLine = [tag, line[i], line[i+1], line[i+2]]
                    result.append(newLine)
                    continue

        if len(pointValueLine) > 1:
            result.insert(0, pointValueLine)
            if pointValueLine[-1].type == types.SEMICOLON:
                pointValueLine.pop(-1)

        return result

    def onTraitChange(self, event):
        """When the trait type changes, need to switch states"""

        self.parameterList = self.getParameterList()
        self.switchStates()

        for eachPanel in self.valuePanels[:]:
            eachPanel.onTraitChange()

        # now at this point eachPanel.onTraitChange() has taken out any
        # tags that don't belong to that trait type. if there is no
        # trait type (either MM or AM is specified), then the file
        # requires one genotype line (DD, Dd, or dd) in order for kelvin
        # to run properly. so in that case insert a dummy genotype line.
        if self.fm.isPresent('DT') or self.fm.isPresent('QT'):
            if self.hasDummy:
                self.fm.remove(self.dummyLine)
                self.hasDummy = False
        elif not self.hasDummy:
            self.fm.addLine(self.dummyLine)
            self.hasDummy = True

    def onAnalysisChange(self, event):
        """When the analysis type changes (two-point vs multi-point),
        need to switch states.
        
        """
        
        self.parameterList = self.getParameterList()
        self.switchStates()
        for eachPanel in self.valuePanels[:]:
            eachPanel.onAnalysisChange()

    def removeTag(self, tag):
        """Given a tag, remove all values that have that tag"""

        for eachPanel in self.valuePanels[:]:
            if eachPanel.line[0].str == tag:
                eachPanel.onRemoveButton(None)

    def hasTag(self, tag):
        """Given a tag, returns True if the tag is present in a value in 
        this panel, even if the value is not in the file yet.
        
        """

        for eachPanel in self.valuePanels[:]:
            if eachPanel.line[0].str == tag:
                return True
        return False

    def addDefaults(self):
        """Get all the default values given the current file conditions,
        and add them to the values panel if they aren't present.
        
        """

        oldLines = self.fm.getValueDefaults()
        lines = []
        # see if the lines need to be expanded out into multiple lines.
        # a lines needs to be expanded when it mixes point values and 
        # range values together
        for eachLine in oldLines:
            if not self.needsExpanding(eachLine):
                lines.append(eachLine)
            else:
                splitLines = self.splitVariableNumberLine(eachLine)
                lines.extend(splitLines)

        # add the defaults to the file and the values panel if necessary
        for eachLine in lines:
            tag = eachLine[0].str
            if self.isRangeValue(eachLine):
                # it is a range value line
                if not self.fm.isRangeValuePresent(tag):
                    self.fm.addLine(eachLine)
                    valuePanel = self.addValue(self.valuePanels[-1], eachLine)
                    valuePanel.hasChanged()
                    valuePanel.checkForErrors()
            else:
                # it is a point value line
                if not self.fm.isPointValuePresent(tag):
                    self.fm.addLine(eachLine)
                    valuePanel = self.addValue(self.valuePanels[-1], eachLine)
                    valuePanel.hasChanged()
                    valuePanel.checkForErrors()

    def isRangeValue(self, line):
        """Given a line that is a value, determine if it is a range value"""

        # Note: this assumes that the line is not mixed with range and
        # point values

        length = len(line)

        if length == 2:
            return False

        # at this point the line is at least 4 tokens. check if the third 
        # token is a semicolon or not
        if line[2].str.strip() == ';':
            return False
        else:
            return True

    def clean(self):
        """Remove all value panels that are not in the file"""

        for eachValuePanel in self.valuePanels[:]:
            if eachValuePanel.inFile == False:
                self.removeValue(eachValuePanel)

    def addLD(self):
        """Add the default line of linkage disequilibrium (LD) to the file"""

        line = self.fm.insert(*KelvinInfo.defaultLD)
        valuePanel = self.addValue(self.valuePanels[-1], line)
        valuePanel.hasChanged()
        valuePanel.checkForErrors()

    def removeLD(self):
        """Remove all occurences of linkage disequilibrium (LD) from the panel
        and from the file.

        """
        self.removeTag('LD')

    def setAnalysisTypePanel(self, panel):
        """Sets the reference to the analysis type panel"""
        self.analysisTypePanel = panel
