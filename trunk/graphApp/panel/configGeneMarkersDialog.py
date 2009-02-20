"""A dialog box to show and edit gene marker options"""

import copy

import matplotlib.pyparsing
import wx
import wx.lib.colourselect as colorsel

import colors
from config import misc, graphInfo
from panel.scrolledPanel import ScrolledPanel

class ConfigGeneMarkerDialog(wx.Dialog):
    """A dialog box to hold the gene marker config panel"""

    title = 'Configure Gene Markers'

    def __init__(self, mainFrame):

        style = wx.DEFAULT_DIALOG_STYLE
        size = (100,100)
        pre = wx.PreDialog()
        pre.Create(mainFrame, -1, self.title, wx.DefaultPosition, size, style)
        self.PostCreate(pre)

        self.parentFrame = mainFrame

        scrollpanel = ScrolledPanel(self)
        self.configPanel = ConfigGeneMarkerPanel(scrollpanel, mainFrame)

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        scrollSizer = wx.BoxSizer(wx.VERTICAL)
        scrollSizer.Add(self.configPanel, 1, wx.EXPAND | wx.LEFT | wx.TOP, 20)
        scrollpanel.SetSizer(scrollSizer)
        scrollpanel.SetAutoLayout(1)
        scrollpanel.SetupScrolling()

        # have to manually set the size of the scrollpanel
        size = self.configPanel.GetSize()
        x = size[0] + 60
        y = min(size[1] + 30, misc.mainWindowSize[1] - 50)
        scrollpanel.SetMinSize((x, y))
        scrollpanel.Fit()

        scrollSizer = wx.BoxSizer(wx.VERTICAL)
        scrollSizer.Add(scrollpanel, 1, wx.EXPAND)
        mainSizer.Add(scrollpanel, 0, wx.EXPAND)

        # ok and cancel buttons
        self.okButton = wx.Button(self, wx.ID_OK)
        self.cancelButton = wx.Button(self, wx.ID_CANCEL)
        self.okButton.SetDefault()

        # add the ok and cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        btnsizer.AddButton(self.okButton)
        btnsizer.AddButton(self.cancelButton)
        btnsizer.Realize()
        mainSizer.Add(btnsizer, 0,  wx.ALIGN_RIGHT | wx.ALL, 10)

        self.SetSizer(mainSizer)
        self.Fit()
        self.CenterOnParent()

        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)
        self.cancelButton.Bind(wx.EVT_BUTTON, self.onCancel, self.cancelButton)
        self.Bind(wx.EVT_CLOSE, self.onCloseWindow)

    def onOK(self, event):
        """Event when the ok button is pressed"""
        self.Destroy()

    def onCancel(self, event):
        """Event when the cancel button is pressed"""

        # need to restore all old settings
        self.configPanel.restoreBackup()
        self.Destroy()

    def onCloseWindow(self, event):
        """Event when the user presses the close window icon"""

        # treat this event as if they pressed the cancel button
        self.onCancel(None)

class ConfigGeneMarkerPanel(wx.Panel):
    """A panel to show and edit gene marker options"""

    headingsLbl = ['Name', 'Size', 'Color', 'Symbol']

    def __init__(self, parent, mainFrame):
        wx.Panel.__init__(self, parent, -1)

        self.mainFrame = mainFrame

        self.createWidgets()
        self.arrange()
        self.initialize()
        self.createBackup()

    def createWidgets(self):
        """Create all widgets needed for this panel"""

        # create all the headings
        self.headings = []
        for eachLbl in self.headingsLbl:
            self.headings.append(wx.StaticText(self, -1, eachLbl))

        # create widgets for each marker group
        self.groupWidgets = []
        for eachGroup in self.mainFrame.geneMarkerGroups:
            self.groupWidgets.append(GroupWidgets(self, self.mainFrame, 
                                                  eachGroup))

    def arrange(self):
        """Arrange widgets for this panel"""

        cols = 4
        rows = len(self.mainFrame.geneMarkerGroups)
        sizer = wx.FlexGridSizer(rows, cols, 5, hgap=15)

        # add headings
        style = wx.ALIGN_CENTER_HORIZONTAL
        for eachHeading in self.headings:
            sizer.Add(eachHeading, 0, style)

        # arrange widgets for each group
        for eachRow in self.groupWidgets:
            eachRow.arrange(sizer)

        self.SetSizer(sizer)

    def initialize(self):
        """Set initial values of the widgets"""
        pass

    def createBackup(self):
        """Backup values that may change"""

        for eachGroup in self.mainFrame.geneMarkerGroups:
            eachGroup.backup()

    def restoreBackup(self):
        """Restore backup values"""

        for eachGroup in self.mainFrame.geneMarkerGroups:
            eachGroup.restore()

        self.update()

    def update(self):
        """Updates the the graph"""
        self.mainFrame.draw()

class GroupWidgets():
    """A class that takes care of all the widgets for one gene marker group"""
    
    def __init__(self, parent, mainFrame, group):
        self.parent = parent
        self.mainFrame = mainFrame
        self.group = group
        self.createWidgets()
        self.initialize()

    def createWidgets(self):
        """Create all widgets needed for one line"""

        parent = self.parent

        # text indicating the name of the group
        self.nameBox = wx.TextCtrl(parent, -1, self.group.name)

        # a color button
        color = colors.hex2rgb(self.group.color)
        self.colorButton = colorsel.ColourSelect(parent, -1, '', color)

        # choice box for marker symbol
        choices = graphInfo.lineSymbolDescriptions
        self.markerSymbolChoice = wx.Choice(parent, -1, choices=choices)

        # spinner for symbol size
        self.symbolSizeSpin = wx.SpinCtrl(parent, -1)

    def arrange(self, sizer):
        """Arrange widgets for the markers"""

        style = wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL
        sizer.Add(self.nameBox, 0, style | wx.LEFT, 3)
        sizer.Add(self.symbolSizeSpin, 0, style)
        sizer.Add(self.colorButton, 0, style)
        sizer.Add(self.markerSymbolChoice, 0, style)

    def initialize(self):
        """Initialize widgets to values"""

        # intialize choice boxes
        index = graphInfo.lineSymbolTags.index(self.group.symbol)
        self.markerSymbolChoice.SetSelection(index)

        # initialize spinners
        self.symbolSizeSpin.SetRange(0, 100)
        self.symbolSizeSpin.SetValue(self.group.size)

        # bind events
        self.nameBox.Bind(wx.EVT_TEXT, self.onNameText, self.nameBox)
        self.colorButton.Bind(colorsel.EVT_COLOURSELECT, self.onColor, 
                                                         self.colorButton)
        self.markerSymbolChoice.Bind(wx.EVT_CHOICE, self.onMarkerSymbolChoice, 
                                                    self.markerSymbolChoice)
        self.symbolSizeSpin.Bind(wx.EVT_SPINCTRL, self.onSymbolSize, 
                                                  self.symbolSizeSpin)

    def onNameText(self, event):
        """Event when user types in the group name text box"""

        newName = self.nameBox.GetValue()
        self.group.name = newName
        self.update()

    def onColor(self, event):
        """Event when a color is selected"""

        color = colors.hexcolor(event.GetValue())
        self.group.set_color(color)
        self.update()

    def onMarkerSymbolChoice(self,event):
        """Event when a marker symbol was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineSymbolDescrToTag[event.GetString()]
        self.group.set_marker(tag)
        self.update()

    def onSymbolSize(self, event):
        """Event when the symbol size spinner changes value"""

        self.group.set_markersize(self.symbolSizeSpin.GetValue())
        self.update()

    def update(self):
        self.parent.update()
