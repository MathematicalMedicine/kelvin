"""A dialog box to show and edit graph configuration options"""

import copy

import matplotlib.pyparsing
import wx
import wx.lib.colourselect as colorsel

from colors import mpl_color, hexcolor, hex2rgb
from config import misc, graphInfo
import panel.mainFrame_parts.annotation as ann
from panel.scrolledPanel import ScrolledPanel

class ConfigDialog(wx.Dialog):
    """A dialog box to hold the config panel"""

    title = 'Configure Graph'

    def __init__(self, mainFrame):

        style = wx.DEFAULT_DIALOG_STYLE
        size = (100,100)
        pre = wx.PreDialog()
        pre.Create(mainFrame, -1, self.title, wx.DefaultPosition, size, style)
        self.PostCreate(pre)

        self.parentFrame = mainFrame

        scrollpanel = ScrolledPanel(self)
        self.configPanel = ConfigPanel(scrollpanel, mainFrame)

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        scrollSizer = wx.BoxSizer(wx.VERTICAL)
        scrollSizer.Add(self.configPanel, 1, wx.EXPAND)
        scrollSizer.Add((100,100))
        scrollpanel.SetSizer(scrollSizer)
        scrollpanel.SetAutoLayout(1)
        scrollpanel.SetupScrolling()

        # have to manually set the size of the scrollpanel
        size = self.configPanel.GetSize()
        x = size[0] + 100
        y = min(size[1], misc.mainWindowSize[1] - 100)
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
        self.configPanel.restoreBackupConfig()
        self.Destroy()

    def onCloseWindow(self, event):
        """Event when the user presses the close window icon"""

        # treat this event as if they pressed the cancel button
        self.onCancel(None)

class ConfigPanel(wx.Panel):
    """A panel to show and edit graph configuration options"""

    TitleLbl = 'Title:'
    xAxisLbl = 'X Label:'
    yAxisLbl = 'Y Label:'
    titleSizeLbl = 'Title Text Size'
    axisSizeLbl = 'Axis Text Size'
    tickSizeLbl = 'Tick Mark Text Size'
    legendChkBxLbl = 'Show Legend'
    legendFrameChkBxLbl = 'Show Legend Frame'
    legendPosLbl = 'Legend Location:'
    legendShowSetsLbl = 'List Line Sets'
    legendShowIndvLinesLbl = 'List Individual Lines'
    annPPLPrecLbl = 'Annotation PPL Precision:'
    defaultLbl = 'Default'
    userSpecifiedLbl = 'Custom'
    annPosPrecLbl = 'Annotation Position Precision'
    annColorChkBxLbl = 'Keep annotations the same color as its line'
    duplicateNumChkBxLbl ='Suppress duplicate chromosome numbers in line labels'

    title_flag = 'title_flag'
    xaxis_flag = 'x_axis_flag'
    yaxis_flag = 'y_axis_flag'

    def __init__(self, parent, mainFrame):
        wx.Panel.__init__(self, parent, -1)

        self.mainFrame = mainFrame

        # only show the lines panel if there are lines in the graph
        if len(self.mainFrame.lines) > 0:
            self.showLines = True
        else:
            self.showLines = False

        self.createWidgets()
        self.arrange()
        self.initialize()
        self.createBackupConfig()

    def createWidgets(self):
        """Create all widgets needed for this panel"""

        # text labels
        self.titleLabel = wx.StaticText(self, -1, self.TitleLbl)
        self.xAxisLabel = wx.StaticText(self, -1, self.xAxisLbl)
        self.yAxisLabel = wx.StaticText(self, -1, self.yAxisLbl)

        # text boxes
        self.titleBox = wx.TextCtrl(self, -1, self.mainFrame.config.title)
        self.xAxisBox = wx.TextCtrl(self, -1, self.mainFrame.config.xlabel)
        self.yAxisBox = wx.TextCtrl(self, -1, self.mainFrame.config.ylabel)

        # button that selects background color
        axes = self.mainFrame.figure.gca()
        bgcolor = mpl_color(axes.get_axis_bgcolor(), default=(255,255,255))
        self.bgcolorButton = colorsel.ColourSelect(self, -1, 'Background Color',
                                                   bgcolor)

        # labels and spinners for text sizes
        self.titleSizeLabel = wx.StaticText(self, -1, self.titleSizeLbl)
        self.axisSizeLabel = wx.StaticText(self, -1, self.axisSizeLbl)
        self.titleSizeSpin = wx.SpinCtrl(self, -1)
        self.axisSizeSpin = wx.SpinCtrl(self, -1)

        # widgets for annotation ppl precision
        self.annPPLPrecLabel = wx.StaticText(self, -1, self.annPPLPrecLbl)
        self.annPPLPrecSpin = wx.SpinCtrl(self, -1)
        self.defaultRadio = wx.RadioButton(self, -1, self.defaultLbl, 
                                                              style=wx.RB_GROUP)
        self.userSpecifiedRadio = wx.RadioButton(self, -1,self.userSpecifiedLbl)

        # widgets for annotation position precision
        self.annPosPrecLabel = wx.StaticText(self, -1, self.annPosPrecLbl)
        self.annPosPrecSpin = wx.SpinCtrl(self, -1)

        # labels and checkboxes for legend options
        self.showLegend = wx.CheckBox(self, -1, self.legendChkBxLbl)
        self.lineSetRadio = wx.RadioButton(self, -1, self.legendShowSetsLbl, 
                                                              style=wx.RB_GROUP)
        self.indvLinesRadio=wx.RadioButton(self, -1,self.legendShowIndvLinesLbl)
        self.showLegendFrame = wx.CheckBox(self, -1, self.legendFrameChkBxLbl)
        self.legendPosLabel = wx.StaticText(self, -1, self.legendPosLbl)
        self.legendPosChoice = wx.Choice(self, -1, 
                                         choices=graphInfo.legendPositions)

        # checkmark for the option to display annotations as all the same color
        # or to make annotations the same color as its line
        self.annColorChkbx = wx.CheckBox(self, -1, self.annColorChkBxLbl)

        # checkmark for the option to not allow duplicate numbers in the
        # labels for the lines at the bottom of the graph
        self.duplicateNumChkBx = wx.CheckBox(self, -1,self.duplicateNumChkBxLbl)

        # a panel for line styles
        if self.showLines:
            self.linesPanel = LinesPanel(self, self.mainFrame)

    def arrange(self):
        """Arrange widgets of this panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # title and axis labels and text boxes
        sizer = wx.FlexGridSizer(3, 2, 5, 5)
        sizer.SetFlexibleDirection(wx.HORIZONTAL)
        sizer.AddGrowableCol(1)

        # this row has text boxes to specify title, x-label, y-label
        style = wx.ALIGN_CENTER_VERTICAL 
        style2 = wx.ALIGN_CENTER_VERTICAL | wx.EXPAND | wx.ALL
        sizer.Add(self.titleLabel, 0, style)
        sizer.Add(self.titleBox, 1, style2, 2)
        sizer.Add(self.xAxisLabel, 0, style)
        sizer.Add(self.xAxisBox, 1, style2, 2)
        sizer.Add(self.yAxisLabel, 0, style)
        sizer.Add(self.yAxisBox, 1, style2, 2)
        mainSizer.Add(sizer, 1, wx.EXPAND | wx.LEFT | wx.RIGHT | wx.TOP, 10)

        # next row is background color, title and axis text sizes
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        sizer.Add(self.bgcolorButton, 0, style, 5)
        sizer.Add(self.titleSizeLabel, 0, style, 5)
        sizer.Add(self.titleSizeSpin, 0, style, 5)
        sizer.Add(self.axisSizeLabel, 0, style, 5)
        sizer.Add(self.axisSizeSpin, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 10)

        # next two rows is legend options
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        sizer.Add(self.showLegend, 0, style, 5)
        sizer.Add(self.lineSetRadio, 0, style, 5)
        sizer.Add(self.indvLinesRadio, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        sizer.Add(self.showLegendFrame, 0, style, 5)
        sizer.Add(self.legendPosLabel, 0, style, 5)
        sizer.Add(self.legendPosChoice, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)

        # next row is annotation options
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        sizer.Add(self.annPPLPrecLabel, 0, style, 5)
        sizer.Add(self.defaultRadio, 0, style, 5)
        sizer.Add(self.userSpecifiedRadio, 0, style, 5)
        sizer.Add(self.annPPLPrecSpin, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)

        # next row is other annotation options
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        sizer.Add(self.annPosPrecLabel, 0, style, 5)
        sizer.Add(self.annPosPrecSpin, 0, style, 5)
        sizer.Add(self.annColorChkbx, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)

        # next row is duplicate numbers option
        sizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        sizer.Add(self.duplicateNumChkBx, 0, style, 5)
        mainSizer.Add(sizer, 0, wx.ALIGN_CENTER_HORIZONTAL | wx.TOP, 5)

        # add the lines panel
        if self.showLines:
            style = wx.ALIGN_CENTER_HORIZONTAL | wx.TOP
            mainSizer.Add(self.linesPanel, 0, style, 10)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widgets of this panel"""

        def handle(arg):
            """Creates text handlers for event watching"""
            return lambda event: self.onText(event, arg)

        config = self.mainFrame.config

        self.titleSizeSpin.SetRange(5, 100)
        self.titleSizeSpin.SetValue(config.titleSize)
        self.axisSizeSpin.SetRange(5, 100)
        self.axisSizeSpin.SetValue(config.axisSize)

        pos = config.legendPosition
        self.legendPosChoice.SetStringSelection(pos)
        self.showLegend.SetValue(config.showLegend)
        self.lineSetRadio.SetValue(config.showLineSetsInLegend)
        self.indvLinesRadio.SetValue(not config.showLineSetsInLegend)
        self.showLegendFrame.SetValue(config.showLegendFrame)

        self.defaultRadio.SetValue(config.defaultAnnPPLPrec)
        self.userSpecifiedRadio.SetValue(not config.defaultAnnPPLPrec)
        self.annPPLPrecSpin.SetRange(2, 100)
        self.annPPLPrecSpin.SetValue(config.annPPLPrec)
        self.annPosPrecSpin.SetRange(0, 100)
        self.annPosPrecSpin.SetValue(config.annPosPrec)
        self.annColorChkbx.SetValue(not config.isAnnSameColor)
        self.duplicateNumChkBx.SetValue(not config.allowDuplicates)

        # bind events
        self.titleSizeSpin.Bind(wx.EVT_SPINCTRL, self.onTitleSpin, 
                                self.titleSizeSpin)
        self.axisSizeSpin.Bind(wx.EVT_SPINCTRL, self.onAxisSpin, 
                                self.axisSizeSpin)

        self.titleBox.Bind(wx.EVT_TEXT, handle(self.title_flag), self.titleBox)
        self.xAxisBox.Bind(wx.EVT_TEXT, handle(self.xaxis_flag), self.xAxisBox)
        self.yAxisBox.Bind(wx.EVT_TEXT, handle(self.yaxis_flag), self.yAxisBox)

        self.bgcolorButton.Bind(colorsel.EVT_COLOURSELECT, self.onBgColor, 
                                                           self.bgcolorButton)
        self.showLegend.Bind(wx.EVT_CHECKBOX, self.onShowLegend,self.showLegend)
        self.lineSetRadio.Bind(wx.EVT_RADIOBUTTON, self.onLineSetRadio, 
                               self.lineSetRadio)
        self.indvLinesRadio.Bind(wx.EVT_RADIOBUTTON, self.onIndvLinesRadio, 
                                 self.indvLinesRadio)
        self.showLegendFrame.Bind(wx.EVT_CHECKBOX, self.onShowLegendFrame,
                                  self.showLegendFrame)
        self.legendPosChoice.Bind(wx.EVT_CHOICE, self.onLegendPosChoice, 
                                  self.legendPosChoice)
        self.defaultRadio.Bind(wx.EVT_RADIOBUTTON, self.onDefaultRadio, 
                               self.defaultRadio)
        self.userSpecifiedRadio.Bind(wx.EVT_RADIOBUTTON,
                                     self.onUserSpecifiedRadio,
                                     self.userSpecifiedRadio)
        self.annPPLPrecSpin.Bind(wx.EVT_SPINCTRL, self.onAnnPPLPrecSpin, 
                                 self.annPPLPrecSpin)
        self.annPosPrecSpin.Bind(wx.EVT_SPINCTRL, self.onAnnPosPrecSpin, 
                                 self.annPosPrecSpin)
        self.annColorChkbx.Bind(wx.EVT_CHECKBOX, self.onAnnColorChkBx, 
                                self.annColorChkbx)
        self.duplicateNumChkBx.Bind(wx.EVT_CHECKBOX, self.onDuplicateNumChkBx, 
                                self.duplicateNumChkBx)

    def onText(self, event, which):
        """Event when text is typed in a text box.  The 'which' argument
        determines which box was typed in.
        
        """

        s = event.GetString()
        if '"' not in s:
            t = r"%s" % s
        elif "'" not in s:
            t = r'%s' % s
        else:
            t = r"""%s""" % s

        if which == self.title_flag:
            self.mainFrame.config.title = t
        elif which == self.xaxis_flag:
            self.mainFrame.config.xlabel = t
        elif which == self.yaxis_flag:
            self.mainFrame.config.ylabel = t

        try:
            # putting this in a try block because the draw() command in 
            # self.update() could result in an error when someone is typing 
            # in a TeX expression, and it's not complete yet. so it reads
            # in the incomplete expression as an error
            self.update()
        except ValueError, e:
            pass
        except matplotlib.pyparsing.ParseException, e:
            pass

    def onTitleSpin(self, event):
        """Event when the title size spinner has changed value"""

        value = self.titleSizeSpin.GetValue()
        self.mainFrame.config.titleSize = value
        self.update()

    def onAxisSpin(self, event):
        """Event when the axis size spinner has changed value"""

        value = self.axisSizeSpin.GetValue()
        self.mainFrame.config.axisSize = value
        self.update()

    def onBgColor(self, event):
        """Event when a background color is selected"""

        axes = self.mainFrame.figure.gca()
        color = hexcolor(event.GetValue())
        self.mainFrame.config.bgcolor = color
        axes.set_axis_bgcolor(color)

        self.update()

    def onShowLegend(self, event):
        """Event when the show legend checkbox is clicked"""

        self.mainFrame.config.showLegend = self.showLegend.GetValue()
        self.update()

    def onLineSetRadio(self, event):
        """Event when the user clicks on the list line sets for legend display"""

        if self.lineSetRadio.GetValue():
            self.mainFrame.config.showLineSetsInLegend = True
        self.update()

    def onIndvLinesRadio(self, event):
        """Event when the user clicks on the list individual lines for 
        legend display.
        
        """

        if self.indvLinesRadio.GetValue():
            self.mainFrame.config.showLineSetsInLegend = False
        self.update()

    def onShowLegendFrame(self, event):
        """Event when the show legend frame checkbox is clicked"""

        self.mainFrame.config.showLegendFrame = self.showLegendFrame.GetValue()
        self.update()

    def onLegendPosChoice(self, event):
        """Event when the legend position is changed"""

        self.mainFrame.config.legendPosition = \
                                       self.legendPosChoice.GetStringSelection()
        self.update()

    def onDefaultRadio(self, event):
        """Event when the default annotation ppl precision radio button is
        clicked
        
        """

        if self.defaultRadio.GetValue():
            self.mainFrame.config.defaultAnnPPLPrec = True
        self.update(updateAnnotations=True)

    def onUserSpecifiedRadio(self, event):
        """Event when the user specified annotation ppl precision radio button
        is clicked
        
        """

        if self.userSpecifiedRadio.GetValue():
            self.mainFrame.config.defaultAnnPPLPrec = False
        self.update(updateAnnotations=True)

    def onAnnPPLPrecSpin(self, event):
        """Event when the annotatino position precision spinner 
        has changed value
        
        """

        value = self.annPPLPrecSpin.GetValue()
        self.mainFrame.config.annPPLPrec = value
        self.update(updateAnnotations=True)

    def onAnnPosPrecSpin(self, event):
        """Event when the annotatino position precision spinner 
        has changed value
        
        """

        value = self.annPosPrecSpin.GetValue()
        self.mainFrame.config.annPosPrec = value
        self.update(updateAnnotations=True)

    def onAnnColorChkBx(self, event):
        """Event when the change annotation color checkbox is clicked"""

        value = not self.annColorChkbx.GetValue()
        self.mainFrame.config.isAnnSameColor = value
        self.update(updateAnnotations=True)

    def onDuplicateNumChkBx(self, event):
        """Event when the allow duplicate line number checkbox is clicked"""

        value = not self.duplicateNumChkBx.GetValue()
        self.mainFrame.config.allowDuplicates = value
        for eachGroup in self.mainFrame.lines.getGroups():
            eachGroup.allowDuplicates = value
            eachGroup.redoLabel()
        self.mainFrame.drawLabels()
        self.update(updateAnnotations=True)

    def createBackupConfig(self):
        """Create a backup copy of the config object"""

        self.backupConfig = copy.copy(self.mainFrame.graph.config)

        # also backup all lines and line groups
        for eachLine in self.mainFrame.lines:
            eachLine.backup()
        for eachGroup in self.mainFrame.lines.getGroups():
            eachGroup.backup()

    def restoreBackupConfig(self):
        """Use the backup config for the settings"""

        self.mainFrame.graph.config = self.backupConfig
        self.mainFrame.config = self.backupConfig
        for eachLine in self.mainFrame.lines:
            eachLine.restore()
        for eachGroup in self.mainFrame.lines.getGroups():
            eachGroup.restore()

        # line names may have changed
        self.mainFrame.changeLineName()

        self.update(updateAnnotations=True)

    def update(self, updateAnnotations=False):
        """Updates the the graph"""

        self.mainFrame.config.update()
        self.mainFrame.updateLegend(False)
        if updateAnnotations:
            ann.updateAnnotations(self.mainFrame)
        self.mainFrame.updateInfoText()
        self.mainFrame.draw()

class LinesPanel(wx.Panel):
    """A panel that shows a list of all lines and ways to edit line styles"""

    headingsLbl = ['Label', 'Color', 'Line Style', 'Thickness', 'Symbol', 
                   'Symbol Size']

    def __init__(self, parent, mainFrame):
        wx.Panel.__init__(self, parent, -1)

        self.mainFrame = mainFrame

        self.createWidgets()
        self.arrange()
        self.initialize()

    def createWidgets(self):
        """Create widgets for this panel"""

        # create all the headings
        self.headings = []
        for eachLbl in self.headingsLbl:
            self.headings.append(wx.StaticText(self, -1, eachLbl))

        # create widgets for each line in the graph
        self.lineWidgets = []
        for eachLine in self.mainFrame.lines:
            self.lineWidgets.append(LineWidgets(self, self.mainFrame, eachLine))

    def arrange(self):
        """Arrange widgets for this panel"""

        cols = 6
        rows = len(self.mainFrame.lines)
        sizer = wx.FlexGridSizer(rows, cols, 5, hgap=15)

        # add headings
        style = wx.ALIGN_CENTER_HORIZONTAL
        for eachHeading in self.headings:
            sizer.Add(eachHeading, 0, style)

        # arrange widgets for each line in the graph
        for eachRow in self.lineWidgets:
            eachRow.arrange(sizer)

        self.SetSizer(sizer)

    def initialize(self):
        """Set initial values of the widgets"""
        pass

    def update(self, updateAnnotations=False):
        self.GetParent().update(updateAnnotations)

class LineWidgets():
    """A class that takes care of all the widgets for one line"""
    
    def __init__(self, parent, mainFrame, line):
        self.parent = parent
        self.mainFrame = mainFrame
        self.line = line
        self.group = mainFrame.lines.getGroup(line.name)
        self.createWidgets()
        self.initialize()

    def createWidgets(self):
        """Create all widgets needed for one line"""

        parent = self.parent

        # text indicating the name of the chromosome
        self.nameBox = wx.TextCtrl(parent, -1, self.line.name)

        # a color button
        color = mpl_color(self.line.color)
        self.colorButton = colorsel.ColourSelect(parent, -1, '', color)

        # choice box for line style
        choices = graphInfo.lineStyleDescriptions
        self.lineStyleChoice = wx.Choice(parent, -1, choices=choices)

        # spinner for line thickness
        self.thickSpin = wx.SpinCtrl(parent, -1)

        # choice box for line symbol
        choices = graphInfo.lineSymbolDescriptions
        self.lineSymbolChoice = wx.Choice(parent, -1, choices=choices)

        # spinner for symbol size
        self.symbolSizeSpin = wx.SpinCtrl(parent, -1)

    def arrange(self, sizer):
        """Arrange widgets for the line"""

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
        index = graphInfo.lineStyleTags.index(self.line.linestyle)
        self.lineStyleChoice.SetSelection(index)
        index = graphInfo.lineSymbolTags.index(self.line.marker)
        self.lineSymbolChoice.SetSelection(index)

        # initialize spinners
        self.thickSpin.SetRange(0, 100)
        self.thickSpin.SetValue(self.line.thickness)
        self.symbolSizeSpin.SetRange(0, 100)
        self.symbolSizeSpin.SetValue(self.line.markersize)

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
        """Event when user types in the line name text box"""

        newName = self.nameBox.GetValue()
        self.mainFrame.changeLineName(self.line, self.group, newName)
        self.update()

    def onColor(self, event):
        """Event when a color is selected"""

        color = hexcolor(event.GetValue())
        self.line.set_color(color)
        self.update(updateAnnotations=True)

    def onThickness(self, event):
        """Event when the thickness spinner changes"""

        self.line.set_linewidth(self.thickSpin.GetValue())
        self.update()

    def onLineStyleChoice(self,event):
        """Event when a line style was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineStyleDescrToTag[event.GetString()]
        self.line.set_linestyle(tag)
        self.update()

    def onLineSymbolChoice(self,event):
        """Event when a line symbol was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineSymbolDescrToTag[event.GetString()]
        self.line.set_marker(tag)
        self.update()

    def onSymbolSize(self, event):
        """Event when the symbol size spinner changes value"""

        self.line.set_markersize(self.symbolSizeSpin.GetValue())
        self.update()

    def update(self, updateAnnotations=False):
        self.parent.update(updateAnnotations)
