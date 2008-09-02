"""The panel to select the analysis type, currently TP, SS, or SA"""

import wx
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from utility import util

class AnalysisTypePanel(wx.Panel):
    """The panel to select the analysis type, currently TP, SS or SA"""

    NUM_COLS = 4
    borderAmt = 5
    analysisLabel = 'Type of analysis:'
    numberLabel = 'Number of markers:'
    analysisList = ['two point', 'multi-point']
    sexList = ['sex average', 'sex specific']
    linkageList = ['linkage equilibrium', 'linkage disequilibrium']

    # the range for the spinner that does the number of markers
    spinRange = (1, 10)

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.analysisTypeTags = KelvinInfo.analysisTypeTags

        self.createPanel();
        self.arrangePanel();
        self.initialize();

    def createPanel(self):
        """Create all the widgets for the panel"""

        # type label
        self.typeLabel = wx.StaticText(self, -1, self.analysisLabel)

        # analysis type choice box
        self.analysisChoice = wx.Choice(self, -1, choices=self.analysisList)

        # sex type choice box
        self.sexChoice = wx.Choice(self, -1, choices=self.sexList)

        # linkage type choice box
        self.linkageChoice = wx.Choice(self, -1, choices=self.linkageList)

        # label for number of markers 
        self.numberLabel = wx.StaticText(self, -1, self.numberLabel)

        # spin control for the number
        self.number = wx.SpinCtrl(self, -1)
        self.number.SetRange(*self.spinRange)

    def arrangePanel(self):
        """Arrange all the widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)

        # analysis type label
        rowSizer.Add(self.typeLabel, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL, 
                     self.borderAmt)

        # three choice boxes
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        rowSizer.Add(self.analysisChoice, 0, style, self.borderAmt)
        rowSizer.Add(self.sexChoice, 0, style, self.borderAmt)
        rowSizer.Add(self.linkageChoice, 0, style, self.borderAmt)

        mainSizer.Add(rowSizer)

        # number label and box on next line
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.LEFT
        rowSizer.Add(self.numberLabel , 0, style, self.borderAmt) 
        rowSizer.Add(self.number, 0, style, self.borderAmt)

        mainSizer.Add(rowSizer)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize the choice boxes and number based on the file contents"""

        # assume that there was only one line that specified the analysis type
        if self.fm.isPresent('TP'):
            # check if TP (two point) was specified in the file
            lines = self.fm.getLines('TP')
            self.line = lines[0]
            self.analysisChoice.SetSelection(0)
            self.analysisChoice.charge = self.line[0]
            self.sexChoice.SetSelection(0)
            self.number.Enable(False)
            self.number.SetValue(2)
            self.number.charge = None
        elif self.fm.isPresent('SA'):
            # SA (multi-point sex average) was specified in the file
            lines = self.fm.getLines('SA')
            self.line = lines[0]
            self.analysisChoice.SetSelection(1)
            self.analysisChoice.charge = self.line[0]
            self.sexChoice.SetSelection(0)
            self.number.Enable(True)
            self.number.SetValue(float(self.line[1].str))
            self.number.charge = self.line[1]
        elif self.fm.isPresent('SS'):
            # SS (multi-point sex specific) was specified in the file
            lines = self.fm.getLines('SS')
            self.line = lines[0]
            self.analysisChoice.SetSelection(1)
            self.analysisChoice.charge = self.line[0]
            self.sexChoice.SetSelection(1)
            self.number.Enable(True)
            self.number.SetValue(int(self.line[1].str))
            self.number.charge = self.line[1]

        # the linkage equilibrium is independent of analysis type
        # if an LD tag was present in the file, using disequilibrium,
        # else using linkage equilibrium
        if self.fm.isPresent('LD'):
            self.linkageChoice.SetSelection(1)
        else:
            self.linkageChoice.SetSelection(0)

        # bind some events
        self.analysisChoice.Bind(wx.EVT_CHOICE, self.onAnalysisChoice, 
                                 self.analysisChoice)
        self.sexChoice.Bind(wx.EVT_CHOICE, self.onAnalysisChoice, 
                            self.sexChoice)
        self.linkageChoice.Bind(wx.EVT_CHOICE, self.onLinkageChoice, 
                                self.linkageChoice)
        self.number.Bind(wx.EVT_SPINCTRL, self.onTextEntry, self.number)
        util.bindValuesChange(self, self.onValuesChange)

    def onAnalysisChoice(self, event):
        """Event when the analysis type changes"""

        # change tag based on current values
        analysisSelection = self.analysisChoice.GetCurrentSelection()
        sexSelection = self.sexChoice.GetCurrentSelection()
        number = self.number

        if analysisSelection == self.analysisTypeTags.index('TP'):
            #TODO: do something for two point when sex choice is changed?
            #      nothing in the config file actually changes when sex
            #      choice is changed...i think because linkage equilibrium
            #      is based on other tags

            # it's two point
            self.analysisChoice.charge.update('TP')

            # disable the number box, set it to two
            self.number.Enable(False)
            self.number.SetValue(1)
            if self.number.charge:
                # if the number box had a responsibility, just blank it out
                self.number.charge.update('')

            if sexSelection == 0:
                # with two point, sex average, you can choose the linkage type
                self.linkageChoice.Enable(True)
            else:
                # with two point, sex specific, 
                # you can't choose the linkage type
                self.linkageChoice.Enable(False)

        elif analysisSelection == self.analysisTypeTags.index('SS') or \
             analysisSelection == self.analysisTypeTags.index('SA'):
            # it's multi-point

            # figure out what type of multipoint it is
            if sexSelection == 0:
                self.analysisChoice.charge.update('SA')
            else:
                self.analysisChoice.charge.update('SS')

            # if the number box had no previous responsibility 
            # a.k.a. there was no previous specification for number of markers
            # then add a new item to the line
            if not self.number.charge:
                # this is multi-point, so at least three markers are needed
                self.number.SetValue(3)
                # need to create a space for the new item in the file
                self.number.charge = self.fm.append(self.line, 
                                                    self.number.GetValue())

            # enable number box 
            self.number.Enable(True)
            if self.number.GetValue() < 2:
                self.number.SetValue(2)

            self.number.charge.update(self.number.GetValue())

            # with multi-point, only linkage equilibrium is allowed
            self.linkageChoice.SetSelection(0)
            self.linkageChoice.Disable()

        # since the tag might have changed, remove then add this file back
        # so the line with the new tag gets inserted in an ok place and
        # the comment in the file doesn't lie
        self.fm.remove(self.line)
        self.fm.addLine(self.line)

        util.callAncestorFunction(self, 'hasChanged')

        # an ancestor needs to know about this event
        event.Skip()

    def onLinkageChoice(self, event):
        """Event when user selects linkage equilibrium or disequilibrium"""

        # the way kelvin decides if it's linkage equilibrium or disequilbrium 
        # is if LD (linkage disequilibrium) values are specified in the 
        # grid specification or not. so add it or remove it from the grid
        # specification as necessary
        if self.linkageChoice.GetSelection() == 0:
            # linkage equilibrium is wanted
            self.valuesPanel.removeLD()
        if self.linkageChoice.GetSelection() == 1:
            # linkage disequilibrium is wanted
            if not self.fm.isPresent('LD'):
                self.valuesPanel.addLD()

        # other things needs to know if this value has changed
        util.callAncestorFunction(self, 'hasChanged')

    def onTextEntry(self, event):
        """Event when the number of markers change"""

        self.number.charge.update(self.number.GetValue())
        util.callAncestorFunction(self, 'hasChanged')

    def onValuesChange(self, event):
        """Event when something in the values panel changes"""

        # when this happens, check if linkage disequilibrium matches with
        # the choice in the linkage equilibrium choice box. the user may
        # have removed all LD values, or added LD values in.
        if self.fm.isPresent('LD'):
            self.linkageChoice.SetSelection(1)
        else:
            self.linkageChoice.SetSelection(0)

    def setValuesPanel(self, valuesPanel):
        """Give this panel a reference to the values panel."""
        self.valuesPanel = valuesPanel

    def isSexAverage(self):
        """Returns true if sex average is currently chosen, else False"""

        if self.sexChoice.GetCurrentSelection() == 0:
            return True
        else:
            return False
