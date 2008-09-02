"""The panel to specify the trait type, either dichotomous or quantitative"""

import wx
from config import KelvinInfo
from config import misc
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class TraitTypePanel(wx.Panel):
    """The panel to specify the trait type, dichotomous or quantitative"""

    # how to arrange the panel, either in a row or down a column
    # labels next to choice boxes
    traitLbl = "Trait Type:"
    distributionLbl = "Distribution:"

    #labels next to text boxes
    meanLbl = "Sample Mean:"
    stdLbl = "Sample Std:"
    freedomLbl = "Degree of Freedom:"
    minLbl = "Min:"
    maxLbl = "Max:"

    blank = ' '
    borderAmt = 5

    # a flag to determine if the line that this panel is responsible for
    # is currently in the file or not
    inFile = True

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        # the text shown to describe distributions
        self.distChoiceList = KelvinInfo.distributionDescriptions
        
        # dictionaries to switch from tags to descriptions and vice versa
        self.distTagToDescr = KelvinInfo.distributionTagToDescription
        self.distDescrToTag = KelvinInfo.distributionDescriptionToTag

        # the text shown to describe trait types
        self.traitChoiceList = KelvinInfo.traitTypeDescriptions
        self.traitTags = KelvinInfo.traitTypeTags

        # dictionary to go from descriptions to tag
        self.traitDescrToTag = KelvinInfo.traitTypeDescriptionToTag
        self.traitTagToDescr = KelvinInfo.traitTypeTagToDescription

        # the index that the t-distribution will have in the choice box
        self.t_selection = KelvinInfo.distributionTags.index('T')

        # tag to use marker to marker analysis
        self.MMtag = KelvinInfo.MMEntry[0]

        # indicator to see if have already shown the warning about
        # having QT and TT but not TMIN and TMAX
        self.shownWarning = False

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all widgets for the panel"""

        # labels for trait type and distribution
        self.traitLabel = wx.StaticText(self, -1, self.traitLbl)
        self.distLabel = wx.StaticText(self,-1, self.distributionLbl)

        # choice for trait type and distrubution
        self.traitChoice = wx.Choice(self, -1, choices=self.traitChoiceList)
        self.distChoice = wx.Choice(self, -1, choices=self.distChoiceList)

        # labels for parameters
        self.meanLabel = wx.StaticText(self, -1, self.meanLbl)
        self.stdLabel = wx.StaticText(self, -1, self.stdLbl)
        self.freedomLabel = wx.StaticText(self, -1, self.freedomLbl)
        self.minLabel = wx.StaticText(self, -1, self.minLbl)
        self.maxLabel = wx.StaticText(self, -1, self.maxLbl)

        # text boxes for parameters
        self.meanBox = wx.TextCtrl(self, -1)
        self.stdBox = wx.TextCtrl(self, -1)
        self.freedomBox = wx.TextCtrl(self, -1)
        self.minBox = wx.TextCtrl(self, -1)
        self.maxBox = wx.TextCtrl(self, -1)

    def arrangePanel(self):
        """Arrange all widgets for the panel"""

        mainSizer = wx.FlexGridSizer(cols=2, hgap=self.borderAmt, 
                                             vgap=self.borderAmt)

        # add the trait label and choice
        mainSizer.Add(self.traitLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(self.traitChoice, 0, wx.ALIGN_CENTER_VERTICAL)

        # then distribution label and choice
        mainSizer.Add(self.distLabel,0,wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(self.distChoice,0,wx.ALIGN_CENTER_VERTICAL)

        mainSizer.Add(self.meanLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(self.meanBox, 0, wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(self.stdLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        mainSizer.Add(self.stdBox, 0, wx.ALIGN_CENTER_VERTICAL)

        # hide the degree of freedom and min and max boxes
        #mainSizer.Add(self.freedomLabel,0,wx.ALIGN_CENTER_VERTICAL)
        #mainSizer.Add(self.freedomBox, 0, wx.ALIGN_CENTER_VERTICAL)
        #mainSizer.Add(self.minLabel,0,wx.ALIGN_CENTER_VERTICAL)
        #mainSizer.Add(self.minBox, 0, wx.ALIGN_CENTER_VERTICAL)
        #mainSizer.Add(self.maxLabel,0,wx.ALIGN_CENTER_VERTICAL)
        #mainSizer.Add(self.maxBox, 0, wx.ALIGN_CENTER_VERTICAL)
        self.freedomLabel.Show(False)
        self.freedomBox.Show(False)
        self.minLabel.Show(False)
        self.minBox.Show(False)
        self.maxLabel.Show(False)
        self.maxBox.Show(False)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widget values based on file contents"""

        # figure out if DT or QT was specified, assume at least one of them is
        if self.fm.isPresent('DT'):
            # must be dichotomous model
            self.initializeDT()
            self.oldSelection = 0
        elif self.fm.isPresent('QT'):
            if self.fm.isPresent('TT'):
                # using the threshold model
                self.initilializeQTthreshold()
                self.oldSelection = 2
            else:
                # must be quantitative model
                self.initializeQT()
                self.oldSelection = 1
        else:
            # nothing is specified, so MM must be specified, 
            # which doesn't need a trait type
            self.intializeNoTraitType()
            self.oldSelection = -1

        # bind some events
        self.traitChoice.Bind(wx.EVT_CHOICE,self.onTraitChoice,self.traitChoice)
        self.distChoice.Bind(wx.EVT_CHOICE, self.onDistChoice, self.distChoice)
        self.meanBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.meanBox)
        self.stdBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.stdBox)
        self.freedomBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.freedomBox)
        self.minBox.Bind(wx.EVT_TEXT, self.onMinBox, self.minBox)
        self.maxBox.Bind(wx.EVT_TEXT, self.onMaxBox, self.maxBox)
        util.bindValuesChange(self, self.onValuesChange)

    def initializeDT(self):
        """Set initial conditions when DT is specified at load time"""

        self.inFile = True
        self.line = self.fm.getLines('DT')[0]
        index = self.traitTags.index('DT')
        self.traitChoice.SetSelection(index)

        # set responsibilities
        self.traitChoice.charge = self.line[0]

        # even though this is DT and there is nothing responsible for
        # these things, add it in just in case it gets changed
        # to QT later, it will make the logic alot easier
        self.distChoice.charge = self.fm.append(self.line, '')
        self.meanBox.charge = self.fm.append(self.line, '')
        self.stdBox.charge = self.fm.append(self.line, '')
        self.freedomBox.charge = self.fm.append(self.line, '')

        self.minBox.charge = self.fm.createLine('', '')
        self.maxBox.charge = self.fm.createLine('', '')
        self.minBox.inFile = False
        self.maxBox.inFile = False

        # since dichotomous doesn't need anything more parameters, 
        # disable parameter labels and text boxes
        self.disableDistributions()

    def initializeQT(self):
        """Set initial conditions when QT is specified at load time"""

        self.inFile = True
        self.line = self.fm.getLines('QT')[0]
        index = self.traitTags.index('QT')
        self.traitChoice.SetSelection(index)

        # set responsibility
        self.traitChoice.charge = self.line[0]

        # set the correct distribution
        self.distChoice.charge = self.line[1]
        description = self.distTagToDescr[self.line[1].str]
        index = self.distChoiceList.index(description)
        self.distChoice.SetSelection(index)

        # make sure the correct parameters are enabled or disabled
        # depending on distribution type
        self.meanBox.charge = self.line[2]
        self.stdBox.charge = self.line[3]

        self.meanBox.ChangeValue(self.line[2].str)
        self.stdBox.ChangeValue(self.line[3].str)

        if self.distChoice.GetSelection()==self.t_selection:
            # using degree of freedom
            self.freedomBox.charge = self.line[4]
            self.freedomBox.ChangeValue(self.line[4].str)
        else:
            # not using degree of freedom, 
            # but give it a responsiblity just in case
            self.freedomBox.charge = self.fm.append(self.line, '')

        # check if min and max values are specified. the actual tag
        # depends on if TT is specified or not
        if self.fm.isPresent('TT'):
            minTag = KelvinInfo.TMINEntry[0]
            maxTag = KelvinInfo.TMAXEntry[0]
        else:
            minTag = KelvinInfo.MINEntry[0]
            maxTag = KelvinInfo.MAXEntry[0]

        if self.fm.isPresent(minTag):
            # min tag is specified in the file
            minLine = self.fm.getLines(minTag)[0]
            self.minBox.inFile = True
        else:
            # it wasn't specified, but make a line for it anyway
            minLine = self.fm.createLine('', '')
            self.minBox.inFile = False

        if self.fm.isPresent(maxTag):
            # max tag is specified in the file
            maxLine = self.fm.getLines(maxTag)[0]
            self.maxBox.inFile = True
        else:
            # it wasn't specified, but make a line for it anyway
            maxLine = self.fm.createLine('', '')
            self.maxBox.inFile = False

        self.minBox.charge = minLine
        self.minBox.ChangeValue(self.minBox.charge[1].str)
        self.maxBox.charge = maxLine
        self.maxBox.ChangeValue(self.maxBox.charge[1].str)

        # may need to disable things based on distribution
        self.onDistChoice(None)

    def initilializeQTthreshold(self):
        """Set initial conditions when QT threshold is specified at load time"""

        # the only difference between this and normal QT is 
        # the trait choice needs to be set differently
        self.initializeQT()
        self.traitChoice.SetSelection(2)

    def intializeNoTraitType(self):
        """Set intial conditions when no trait type is specified"""

        # although no line was specified, create a line anyway just
        # in case it's needed, note that the line is not inserted
        # into the file
        self.inFile = False
        self.line = self.fm.createLine('DT', '', '', '', '')

        # set responsibilities
        self.traitChoice.charge = self.line[0]
        self.distChoice.charge = self.line[1]
        self.meanBox.charge = self.line[2]
        self.stdBox.charge = self.line[3]
        self.freedomBox.charge = self.line[4]

        # since intiailizing it to dichotomous, it doesn't need
        # the distribution information
        self.disableDistributions()

        # create lines for the min and max textboxes
        self.minBox.charge = self.fm.createLine('', '')
        self.maxBox.charge = self.fm.createLine('', '')
        self.minBox.inFile = False
        self.maxBox.inFile = False

    def onTraitChoice(self, event):
        """Event when the trait type is changed"""

        # check if using threshold trait or not
        if self.traitChoice.GetSelection() == 2:
            self.fm.trait_threshold = True
        else:
            self.fm.trait_threshold = False

        if self.traitChoice.IsEnabled() is False:
            # do nothing
            event.Skip()
            return

        util.callAncestorFunction(self, 'hasChanged')
        if self.oldSelection == self.traitChoice.GetSelection():
            # nothing changed
            return
        self.oldSelection = self.traitChoice.GetSelection()

        if self.traitChoice.GetSelection() == self.traitTags.index('DT'):
            # if dichotomous was chosen, disable distribution widgets
            self.disableDistributions()
        else:
            # quantitative must have been chosen, so enable distribution widgets
            self.enableDistributions()
            self.onDistChoice(True)

        # update the responsibility
        description = self.traitChoice.GetStringSelection()
        self.traitChoice.charge.update(self.traitDescrToTag[description])

        # since the tag might have changed, remove then add this file back
        # so the line with the new tag gets inserted in an ok place and
        # the comment in the file doesn't lie
        self.fm.remove(self.line)
        self.fm.addLine(self.line)

        # others need to know!
        event.Skip()

    def onDistChoice(self, event):
        """Event when the distribution type is changed"""

        if event != None:
            util.callAncestorFunction(self, 'hasChanged')
        if self.distChoice.IsEnabled() is False:
            # do nothing
            return

        # update the responsibility
        description = self.distChoice.GetStringSelection()
        tag = self.distDescrToTag[description]
        self.distChoice.charge.update(tag)

        # may have to disable degree of freedom textbox
        if self.distChoice.GetSelection() == self.t_selection:
            # enable the degree of freedom widgets
            self.freedomLabel.Enable(True)
            self.freedomBox.Enable(True)
        else:
            # not using degree of freedom widgets
            self.freedomLabel.Enable(False)
            self.freedomBox.Enable(False)
            self.freedomBox.ChangeValue('')
            self.updateValues()

        # put in default values for the new distribution type
        # note: event is None when this function is called for initialization
        if event != None:
            self.setDefaultValues()

    def onTextEntry(self, event):
        """Event when text is entered into any text box"""

        # update the responsibility of the box
        textbox = event.GetEventObject()
        textbox.charge.update(textbox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')

    def onMinBox(self, event):
        """Event when text has been entered in the min box"""

        # update the tag and the number next to the tag
        if self.fm.isPresent('TT'):
            # use the TMIN tag
            self.minBox.charge[0].str = KelvinInfo.TMINEntry[0]
        else:
            # use the MIN tag
            self.minBox.charge[0].str = KelvinInfo.MINEntry[0]
        self.minBox.charge[1].str = self.minBox.GetValue()

        if self.minBox.GetValue().strip() == '' and self.minBox.inFile == True:
            self.fm.remove(self.minBox.charge)
            self.minBox.inFile = False
        if self.minBox.GetValue().strip() != '' and self.minBox.inFile == False:
            self.fm.addLine(self.minBox.charge)
            self.minBox.inFile = True
        if event is not None:
            util.callAncestorFunction(self, 'hasChanged')

    def onMaxBox(self, event):
        """Event when text has been entered in the max box"""

        # update the tag and the number next to the tag
        if self.fm.isPresent('TT'):
            # use the TMAX tag
            self.maxBox.charge[0].str = KelvinInfo.TMAXEntry[0]
        else:
            # use the MAX tag
            self.maxBox.charge[0].str = KelvinInfo.MAXEntry[0]
        self.maxBox.charge[1].str = self.maxBox.GetValue()

        if self.maxBox.GetValue().strip() == '' and self.maxBox.inFile == True:
            self.fm.remove(self.maxBox.charge)
            self.maxBox.inFile = False
        if self.maxBox.GetValue().strip() != '' and self.maxBox.inFile == False:
            self.fm.addLine(self.maxBox.charge)
            self.maxBox.inFile = True
        if event is not None:
            util.callAncestorFunction(self, 'hasChanged')

    def onValuesChange(self, event):
        """Event when something in the values panel (Grid Specification) 
        changes
        
        """

        # check to make sure the min/tmin and max/tmax tags are right
        self.onMinBox(None)
        self.onMaxBox(None)

        # if QT and TT are present, tmax and tmin needs to be there
        TT = self.fm.isPresent('TT')
        if not TT:
            self.shownWarning = False
        # note: there is a very similar check to this 
        # in the file KelvinConfigFileChecker.py
        if not self.shownWarning and self.fm.isPresent('QT') and TT:
            if not self.minBox.inFile or not self.maxBox.inFile:
                s1 = "The file is configured for a multi-point analysis and "
                s2 = "has specified a trait threshold, but is missing Min "
                s3 = "and/or Max values. Please input these values in the "
                s4 = "Settings panel before running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                print message
                dlg = wx.MessageDialog(None, message, "Error!", 
                                       wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
                self.shownWarning = True

    def updateValues(self):
        """updates the responsibilities of the text boxes"""

        self.meanBox.charge.update(self.meanBox.GetValue())
        self.stdBox.charge.update(self.stdBox.GetValue())
        self.freedomBox.charge.update(self.freedomBox.GetValue())

    def setDefaultValues(self, mean=True, std=True, freedom=True):
        """Set the default values of mean, std, and degree of freedom 
        based on the current trait type if the given value is True.
        
        """

        currentTag = self.distDescrToTag[self.distChoice.GetStringSelection()]
        if currentTag == 'T':
            # t distribution currently selected
            defaults = KelvinInfo.defaultDistributionT
            if mean:
                self.meanBox.ChangeValue(defaults[0])
            if std:
                self.stdBox.ChangeValue(defaults[1])
            if freedom:
                self.freedomBox.ChangeValue(defaults[2])
        if currentTag == 'normal':
            # normal distribution is selected
            defaults = KelvinInfo.defaultDistributionNormal
            if mean:
                self.meanBox.ChangeValue(defaults[0])
            if std:
                self.stdBox.ChangeValue(defaults[1])
            if freedom:
                self.freedomBox.ChangeValue('')
        elif currentTag == 'chisq':
            # chi squared distribution is selected
            defaults = KelvinInfo.defaultDistributionChiSquared
            if mean:
                self.meanBox.ChangeValue(defaults[0])
            if std:
                self.stdBox.ChangeValue(defaults[1])
            if freedom:
                self.freedomBox.ChangeValue('')

        self.updateValues()

    def enableTrait(self):
        """Enables all widgets dealing with trait types"""

        # put the line back in the file
        if self.inFile is False:
            self.fm.addLine(self.line)
            self.inFile = True

        # enable all widgets with trait type
        self.enableDistributions()
        self.traitLabel.Enable(True)
        self.enableChoice(self.traitChoice)
        self.oldSelection = -1

    def disableTrait(self):
        """Disables all widgets dealing with trait type"""

        # take the line out of the file
        if self.inFile:
            self.fm.remove(self.line)
            self.inFile = False

        # disable all widgets dealing with trait type
        self.oldSelection = -1
        self.traitLabel.Disable()
        self.disableChoice(self.traitChoice)
        self.disableDistributions()

    def enableDistributions(self):
        """Enables all widgets dealing with the distributions"""

        self.distLabel.Enable(True)
        self.distChoice.Enable(True)
        self.meanLabel.Enable(True)
        self.meanBox.Enable(True)
        self.stdLabel.Enable(True)
        self.stdBox.Enable(True)
        self.freedomLabel.Enable(True)
        self.freedomBox.Enable(True)
        self.minBox.Enable(True)
        self.maxBox.Enable(True)

        self.enableChoice(self.distChoice)
        self.setDefaultValues()

    def disableDistributions(self):
        """Disables all widgets dealing with the distributions"""

        self.distLabel.Enable(False)
        self.distChoice.Enable(False)
        self.meanLabel.Enable(False)
        self.meanBox.Enable(False)
        self.stdLabel.Enable(False)
        self.stdBox.Enable(False)
        self.freedomLabel.Enable(False)
        self.freedomBox.Enable(False)
        self.minBox.Enable(False)
        self.maxBox.Enable(False)

        self.disableChoice(self.distChoice)
        self.distChoice.charge.update('')

        # clear out any info in the parameter boxes
        self.meanBox.ChangeValue('')
        self.stdBox.ChangeValue('')
        self.freedomBox.ChangeValue('')
        self.minBox.ChangeValue('')
        self.maxBox.ChangeValue('')
        self.updateValues()

    def enableChoice(self, choice, selection=0):
        """Enables the given choice box, and takes out blank entries"""

        choice.Enable(True)

        # check if there really is a blank, if so, take it out of the choices
        index = choice.FindString(self.blank)
        if index != wx.NOT_FOUND:
            choice.Delete(index)
            choice.SetSelection(selection)
            util.createEvent(wx.EVT_CHOICE, choice, selection)

    def disableChoice(self, choice):
        """Disables the given choice box, and makes it show a blank entry"""

        choice.Enable(False)

        # check if there already is a blank
        index = choice.FindString(self.blank)
        if index == wx.NOT_FOUND:
            # add a blank entry to the list and make it selected
            choice.Append(self.blank)
            choice.SetSelection(choice.GetCount()-1)
        else:
            choice.SetSelection(index)

        util.createEvent(wx.EVT_CHOICE, choice, choice.GetSelection())

    def Enable(self, bool=True):
        """If bool is True, enables all widgets, and inserts line into file.
        If bool is False, disables and clears all widgets, and removes line 
        from file
        
        """

        wx.Panel.Enable(self, bool)
        if bool:
            # enable all widgets, and inserts line into file
            self.enableTrait()
        else:
            # disables and clears all widgets, and removes line from file
            self.disableTrait()

    def Disable(self):
        """Disables and clears all widgets, and removes line from file"""

        self.Enable(False)

