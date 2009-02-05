"""A dialog box to show and edit line sets options"""

import wx
import wx.lib.colourselect as colorsel

import colors
from config import misc, graphInfo
from panel.scrolledPanel import ScrolledPanel

class ConfigLineSetsDialog(wx.Dialog):
    """A dialog box to hold the line sets config panel"""

    title = 'Configure Line Sets'

    def __init__(self, mainFrame):

        style = wx.DEFAULT_DIALOG_STYLE
        size = (100,100)
        pre = wx.PreDialog()
        pre.Create(mainFrame, -1, self.title, wx.DefaultPosition, size, style)
        self.PostCreate(pre)

        self.parentFrame = mainFrame

        scrollpanel = ScrolledPanel(self)
        self.configPanel = ConfigLineSetsPanel(scrollpanel, mainFrame)

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

class ConfigLineSetsPanel(wx.Panel):
    """A panel to show and edit line set options"""

    headingsLbl = ['Name', 'Color', 'Line Style', 'Thickness', 'Symbol', 
                   'Symbol Size']

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

        # create widgets for each line set
        self.groupWidgets = []
        for eachSet in self.mainFrame.lines.lineSets:
            self.groupWidgets.append(GroupWidgets(self, self.mainFrame, 
                                                  eachSet))

    def arrange(self):
        """Arrange widgets for this panel"""

        cols = 6
        rows = len(self.mainFrame.lines.lineSets)
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
        # GroupWidgets already does this
        pass

    def createBackup(self):
        """Backup values that may change"""

        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.backup()

    def restoreBackup(self):
        """Restore backup values"""

        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.restore()

        self.update()

    def update(self):
        """Updates the the graph"""

        self.mainFrame.updateLegend(False)
        self.mainFrame.draw()

class GroupWidgets():
    """A class that takes care of all the widgets for one line set"""
    
    def __init__(self, parent, mainFrame, lineSet):
        self.parent = parent
        self.mainFrame = mainFrame
        self.lineSet = lineSet
        self.createWidgets()
        self.initialize()

    def createWidgets(self):
        """Create all widgets needed for one line set"""

        parent = self.parent

        # text indicating the name of the line set
        self.nameBox = wx.TextCtrl(parent, -1, self.lineSet.name)

        # a color button
        color = colors.hex2rgb(self.lineSet.color)
        self.colorButton = colorsel.ColourSelect(parent, -1, '', color)

        # choice box for line style
        choices = graphInfo.lineStyleDescriptions
        self.lineStyleChoice = wx.Choice(parent, -1, choices=choices)

        # spinner for line thickness
        self.thickSpin = wx.SpinCtrl(parent, -1)

        # choice box for symbol
        choices = graphInfo.lineSymbolDescriptions
        self.lineSymbolChoice = wx.Choice(parent, -1, choices=choices)

        # spinner for symbol size
        self.symbolSizeSpin = wx.SpinCtrl(parent, -1)

    def arrange(self, sizer):
        """Arrange widgets"""

        style = wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL
        sizer.Add(self.nameBox, 0, style | wx.LEFT, 3)
        sizer.Add(self.colorButton, 0, style)
        sizer.Add(self.lineStyleChoice, 0, style)
        sizer.Add(self.thickSpin, 0, style)
        sizer.Add(self.lineSymbolChoice, 0, style)
        sizer.Add(self.symbolSizeSpin, 0, style)

    def initialize(self):
        """Initialize widgets to values"""

        # intialize choice boxes
        index = graphInfo.lineStyleTags.index(self.lineSet.linestyle)
        self.lineStyleChoice.SetSelection(index)
        index = graphInfo.lineSymbolTags.index(self.lineSet.marker)
        self.lineSymbolChoice.SetSelection(index)

        # initialize spinners
        self.thickSpin.SetRange(0, 100)
        self.thickSpin.SetValue(self.lineSet.thickness)
        self.symbolSizeSpin.SetRange(0, 100)
        self.symbolSizeSpin.SetValue(self.lineSet.markersize)

        # bind events
        self.nameBox.Bind(wx.EVT_TEXT, self.onNameText, self.nameBox)
        self.colorButton.Bind(colorsel.EVT_COLOURSELECT, self.onColor, 
                                                             self.colorButton)
        self.lineStyleChoice.Bind(wx.EVT_CHOICE, self.onLineStyleChoice, 
                                                         self.lineStyleChoice)
        self.thickSpin.Bind(wx.EVT_SPINCTRL, self.onThickness, self.thickSpin)
        self.lineSymbolChoice.Bind(wx.EVT_CHOICE, self.onLineSymbolChoice, 
                                                        self.lineSymbolChoice)
        self.symbolSizeSpin.Bind(wx.EVT_SPINCTRL, self.onSymbolSize, 
                                                          self.symbolSizeSpin)

    def onNameText(self, event):
        """Event when user types in the name text box"""

        newName = self.nameBox.GetValue()
        self.lineSet.name = newName
        self.update()

    def onColor(self, event):
        """Event when a color is selected"""

        color = colors.hexcolor(event.GetValue())
        self.lineSet.set_color(color)
        self.update()

    def onThickness(self, event):
        """Event when the thickness spinner changes"""

        self.lineSet.set_linewidth(self.thickSpin.GetValue())
        self.update()

    def onLineStyleChoice(self,event):
        """Event when a line style was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineStyleDescrToTag[event.GetString()]
        self.lineSet.set_linestyle(tag)
        self.update()

    def onLineSymbolChoice(self,event):
        """Event when a line symbol was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineSymbolDescrToTag[event.GetString()]
        self.lineSet.set_marker(tag)
        self.update()

    def onSymbolSize(self, event):
        """Event when the symbol size spinner changes value"""

        self.lineSet.set_markersize(self.symbolSizeSpin.GetValue())
        self.update()

    def update(self):
        self.parent.update()
