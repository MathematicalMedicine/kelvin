"""The panel to specify affection status values"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class AffectionStatusPanel(wx.Panel):
    """The panel to specify affection status values"""

    mainLbl = 'Phenotype Codes'
    defaultLbl = 'Use defaults'
    unknownLbl = 'Unknown:'
    unaffectedLbl = 'Unaffected:'
    affectedLbl = 'Affected:'

    borderAmt = 5


    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        self.tag = KelvinInfo.affectionStatusEntry[0]

        # a flag to determine if the line that this panel is responsible for
        # is currently in the file or not
        self.inFile = True

        # flag to determine if using default values or not
        self.usingDefaults = True

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create widgets for the panel"""

        # heading
        self.mainLabel = wx.StaticText(self, -1, self.mainLbl)

        # checkbox to see if using defaults or not
        self.defaultCheckbox = wx.CheckBox(self, -1, self.defaultLbl)

        # labels for values
        self.unknownLabel = wx.StaticText(self, -1, self.unknownLbl)
        self.unaffectedLabel = wx.StaticText(self, -1, self.unaffectedLbl)
        self.affectedLabel = wx.StaticText(self, -1, self.affectedLbl)

        # text boxes for those values
        self.unknownBox = wx.TextCtrl(self, -1)
        self.unaffectedBox = wx.TextCtrl(self, -1)
        self.affectedBox = wx.TextCtrl(self, -1)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # first add the mainLabel
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL #wx.LEFT | wx.RIGHT | wx.UP
        #mainSizer.Add(self.mainLabel, 0, style, self.borderAmt)
        sizer.Add(self.mainLabel, 0, style, self.borderAmt)
        sizer.Add(self.defaultCheckbox, 0, style, self.borderAmt)

        mainSizer.Add(sizer)

        # add a line
        #mainSizer.Add(wx.StaticLine(self, style=wx.LI_HORIZONTAL), 0, 
        #              wx.EXPAND | wx.LEFT | wx.RIGHT | wx.BOTTOM, 2)


        # use a flex grid for the other labels and text boxes
        sizer = wx.FlexGridSizer(3, 2, self.borderAmt, self.borderAmt)

        sizer.Add(self.unknownLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(self.unknownBox, 1)

        sizer.Add(self.unaffectedLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(self.unaffectedBox, 1)

        sizer.Add(self.affectedLabel, 0, wx.ALIGN_CENTER_VERTICAL)
        sizer.Add(self.affectedBox, 1)

        mainSizer.Add(sizer, 0, wx.ALL, self.borderAmt)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widget values based on file contents"""

        # initialization is based on whether the AS tag was in the 
        # config file or not
        if self.fm.isPresent(self.tag):
            # AS was specified
            self.initialize_AS_present()
        else:
            # AS was not specified
            self.initialize_AS_notpresent()

        # make sure the right widgets are enable/disabled
        self.updateState()

        # bind some events
        self.unknownBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.unknownBox)
        self.unaffectedBox.Bind(wx.EVT_TEXT,self.onTextEntry,self.unaffectedBox)
        self.affectedBox.Bind(wx.EVT_TEXT, self.onTextEntry, self.affectedBox)
        self.defaultCheckbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, 
                                  self.defaultCheckbox)

        # this panel needs to know when the analysis type and trait type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)
        util.bindTraitTypeChange(self, self.onTraitChange)

    def initialize_AS_present(self):
        """Initialize widgets when AS is present in the file"""

        self.inFile = True
        self.usingDefaults = False
        self.line = self.fm.getLines(self.tag)[0]
        self.unknownBox.ChangeValue(self.line[1].str)
        self.unaffectedBox.ChangeValue(self.line[2].str)
        self.affectedBox.ChangeValue(self.line[3].str)

        # set responsibilities
        self.unknownBox.charge = self.line[1]
        self.unaffectedBox.charge = self.line[2]
        self.affectedBox.charge = self.line[3]

        self.defaultCheckbox.SetValue(False)

    def initialize_AS_notpresent(self):
        """Initialize widgets when AS is not present in the file"""

        self.inFile = False
        self.usingDefaults = True
        self.disableValues()

        # set some defaults based on what trait type was present
        if self.fm.isPresent('DT'):
            defaults = KelvinInfo.defaultAffectionDT
        else:
            defaults = KelvinInfo.defaultAffectionQT

        # even though AS was not specified, still make a line to use
        # with the widgets just in case it gets needed

        # first create a line, but don't put it in the file
        self.line = self.fm.createLine(self.tag, str(defaults[0]), 
                                   str(defaults[1]), str(defaults[2]))

        # set the default values
        self.defaultCheckbox.SetValue(True)
        self.unknownBox.ChangeValue(str(defaults[0]))
        self.unaffectedBox.ChangeValue(str(defaults[1]))
        self.affectedBox.ChangeValue(str(defaults[2]))

        # set responsibilities
        self.unknownBox.charge = self.line[1]
        self.unaffectedBox.charge = self.line[2]
        self.affectedBox.charge = self.line[3]

    def onTextEntry(self, event):
        """Event when text is entered in any of the text boxes"""

        textbox = event.GetEventObject()
        textbox.charge.update(textbox.GetValue())
        util.callAncestorFunction(self, 'hasChanged')

    def onCheckBox(self, event):
        """Event when the checkbox for defaults is checked"""

        if self.defaultCheckbox.GetValue():
            # now using defaults
            self.usingDefaults = True
        else:
            # not using defaults
            self.usingDefaults = False

        # update the state of the panel
        self.updateState()

        util.callAncestorFunction(self, 'hasChanged')

    def onAnalysisChange(self, event):
        """Event when the trait type is changed"""

        # some widgets may have to be enabled/disabled
        self.updateState()


    def onTraitChange(self, event):
        """Event when the trait type is changed"""

        if self.usingDefaults:
            # using defaults and the trait type changed, 
            # so need to change default values

            # this is assuming that a trait change immediately 
            # gets reflected in the file
            self.setDefaults()

    def enableValues(self):
        """Enables the widgets that deals with affection status values"""

        # enable the labels
        self.unknownLabel.Enable(True)
        self.unaffectedLabel.Enable(True)
        self.affectedLabel.Enable(True)

        # enable the text boxes
        self.unknownBox.Enable(True)
        self.unaffectedBox.Enable(True)
        self.affectedBox.Enable(True)

    def disableValues(self):
        """Disables the widgets that deals with affection status values"""

        # disable the labels
        self.unknownLabel.Enable(False)
        self.unaffectedLabel.Enable(False)
        self.affectedLabel.Enable(False)

        # disable the text boxes
        self.unknownBox.Enable(False)
        self.unaffectedBox.Enable(False)
        self.affectedBox.Enable(False)

    def setDefaults(self):
        """Set affection status values to their default values"""

        # the defaults depend on the current trait type (DT or QT)
        if self.fm.isPresent('DT'):
            defaults = KelvinInfo.defaultAffectionDT
        else:
            defaults = KelvinInfo.defaultAffectionQT

        self.unknownBox.ChangeValue(str(defaults[0]))
        self.unaffectedBox.ChangeValue(str(defaults[1]))
        self.affectedBox.ChangeValue(str(defaults[2]))

        self.unknownBox.charge.update(str(defaults[0]))
        self.unaffectedBox.charge.update(str(defaults[1]))
        self.affectedBox.charge.update(str(defaults[2]))

    def setToMPState(self):
        """Disable what is required when using multi-point analysis"""

        # the unaffected and affected boxes need to be disabled
        self.unaffectedLabel.Disable()
        self.unaffectedBox.Disable()
        self.affectedLabel.Disable()
        self.affectedBox.Disable()

    def usingMP(self):
        """Return True if using multi-point analysis, else False"""

        if self.fm.isPresent('SS') or self.fm.isPresent('SA'):
            return True
        else:
            return False

    def updateState(self):
        """Based on current values, enable or disable certain values"""

        usingMP = self.usingMP()
        usingDefaults = self.usingDefaults

        if usingDefaults:
            # disable everything
            self.disableValues()
            self.setDefaults()
            if self.inFile:
                self.fm.remove(self.line)
                self.inFile = False
        else:
            # not using defaults
            if self.inFile is False:
                self.fm.addLine(self.line)
                self.inFile = True

            if usingMP:
                # using multi-point
                self.enableValues()
                self.setToMPState()
            else:
                # not using multi-point, enable everything
                self.enableValues()

