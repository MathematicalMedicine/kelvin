"""The panel to specify whether to use marker to marker analysis"""

import wx
from config import KelvinInfo
from panel.traitTypePanel import TraitTypePanel
from utility import util

class MM_Panel(wx.Panel):
    """The panel to whether to use marker to marker analysis"""

    mainLbl = "Marker to marker analysis (no trait)"
    MMradioLbl = "Analyze all pairwise markers"
    AMradioLbl = "Analyze only adjacent markers"

    borderAmt = 5

    # how much the radio boxes are indented relative to the checkbox
    indent = 20

    def __init__(self, parent, fileManager, traitPanel, affectionPanel):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.traitPanel = traitPanel
        self.affectionPanel = affectionPanel
        self.MMtag = KelvinInfo.MMEntry[0]
        self.AMtag = KelvinInfo.AMEntry[0]

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all the widgets for the panel"""

        # a checkbox to see if marker to marker analysis is wanted
        self.checkbox = wx.CheckBox(self, -1, label=self.mainLbl)

        # radio buttons to determine what type of marker to marker analysis
        # MM and AM are mutually exclusive
        self.MMradio = wx.RadioButton(self, -1, self.MMradioLbl)
        self.AMradio = wx.RadioButton(self, -1, self.AMradioLbl)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        mainSizer.Add(self.checkbox, 0, wx.ALL, self.borderAmt)

        style = wx.LEFT | wx.RIGHT | wx.BOTTOM
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((self.indent, 1))
        sizer.Add(self.MMradio, 0, style, self.borderAmt)
        mainSizer.Add(sizer)

        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add((self.indent, 1))
        sizer.Add(self.AMradio, 0, style, self.borderAmt)
        mainSizer.Add(sizer)

        self.SetSizer(mainSizer)

    def initialize(self):
        """Initialize the panel based on the file contents"""

        if self.fm.isPresent(self.MMtag):
            self.checkbox.SetValue(True)
            self.MMradio.SetValue(True)
            self.AMradio.SetValue(False)
            self.onRadioButton(None)
            self.traitPanel.Disable()
            self.affectionPanel.Disable()
        elif self.fm.isPresent(self.AMtag):
            self.checkbox.SetValue(True)
            self.MMradio.SetValue(False)
            self.AMradio.SetValue(True)
            self.onRadioButton(None)
            self.traitPanel.Disable()
            self.affectionPanel.Disable()
        else:
            # nothing specified
            if self.fm.isPresent('TP'):
                # only available under two point analysis
                self.checkbox.SetValue(False)
                self.disableRadioButtons()
            else:
                self.checkbox.SetValue(False)
                self.checkbox.Disable()
                self.disableRadioButtons()

        # bind events
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)
        self.MMradio.Bind(wx.EVT_RADIOBUTTON, self.onRadioButton, self.MMradio)
        self.AMradio.Bind(wx.EVT_RADIOBUTTON, self.onRadioButton, self.AMradio)

        # need to know when analysis type changes
        util.bindAnalysisTypeChange(self, self.onAnalysisChange)

    def onCheckBox(self, event):
        """Event when the checkbox is checked"""

        if self.checkbox.GetValue():
            # using either MM or AM, therefore trait information and 
            # affection status information is not needed
            self.enableRadioButtons()
            if self.traitPanel.IsEnabled():
                self.traitPanel.Disable()
            if self.affectionPanel.IsEnabled():
                self.affectionPanel.Disable()
        else:
            # not using either MM or AM, therefore trait information and 
            # affection status information is needed
            self.disableRadioButtons()
            if not self.traitPanel.IsEnabled():
                self.traitPanel.Enable()
            if not self.affectionPanel.IsEnabled():
                self.affectionPanel.Enable()
        util.callAncestorFunction(self, 'hasChanged')

    def onRadioButton(self, event):
        """Event when a radio button is selected"""

        if self.MMradio.GetValue():
            if not self.fm.isPresent(self.MMtag):
                self.fm.insert(self.MMtag)
            lines = self.fm.getLines(self.AMtag)
            for eachLine in lines:
                self.fm.remove(eachLine)

        if self.AMradio.GetValue():
            if not self.fm.isPresent(self.AMtag):
                self.fm.insert(self.AMtag)
            lines = self.fm.getLines(self.MMtag)
            for eachLine in lines:
                self.fm.remove(eachLine)

        if event is not None:
            util.callAncestorFunction(self, 'hasChanged')

    def enableRadioButtons(self):
        """Enables all radio buttons"""

        self.MMradio.Enable()
        self.AMradio.Enable()

        # set the first option on
        self.MMradio.SetValue(True)

        # create a radio button event
        util.createEvent(wx.EVT_RADIOBUTTON, self.MMradio, 
                         self.MMradio.GetValue())

    def disableRadioButtons(self):
        """Disables all radio buttons"""

        self.MMradio.Disable()
        self.AMradio.Disable()
        self.MMradio.SetValue(False)
        self.AMradio.SetValue(False)

        # take out lines from file
        lines = self.fm.getLines(self.MMtag)
        for eachLine in lines:
            self.fm.remove(eachLine)
        lines = self.fm.getLines(self.AMtag)
        for eachLine in lines:
            self.fm.remove(eachLine)

    def onAnalysisChange(self, event):
        """Event when the analysis type changes (two point, multi point)"""

        # this panel only is available for two-point, so disable it otherwise
        if self.fm.isPresent('TP'):
            if not self.checkbox.IsEnabled():
                # re-enable it
                self.checkbox.Enable()
                self.disableRadioButtons()
        else:
            if self.checkbox.IsEnabled():
                # disable it
                self.checkbox.SetValue(False)
                self.onCheckBox(None)
                self.checkbox.Disable()
