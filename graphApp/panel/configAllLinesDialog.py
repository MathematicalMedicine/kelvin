"""A dialog box to show and edit options of all lines at once"""

import wx
import wx.lib.colourselect as colorsel

from colors import mpl_color, hexcolor
from config import graphInfo

class ConfigAllLinesDialog(wx.Dialog):
    """A dialog box to let the user configure all lines at once"""

    title = 'Configure All Lines'

    def __init__(self, mainFrame):

        style = wx.DEFAULT_DIALOG_STYLE
        size = (100,100)
        pre = wx.PreDialog()
        pre.Create(mainFrame, -1, self.title, wx.DefaultPosition, size, style)
        self.PostCreate(pre)

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        self.configPanel = ConfigPanel(self, mainFrame)
        mainSizer.Add(self.configPanel, 0, wx.ALL, 20)

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
        self.configPanel.restoreLineValues()
        self.Destroy()

    def onCloseWindow(self, event):
        """Event when the user presses the close window icon"""

        # treat this event as if they pressed the cancel button
        self.onCancel(None)

class ConfigPanel(wx.Panel):
    """A panel to show and edit line configuration options"""

    headingsLbl = ['Color', 'Line Style', 'Thickness', 'Symbol', 'Symbol Size']

    def __init__(self, parent, mainFrame):
        wx.Panel.__init__(self, parent, -1)

        self.mainFrame = mainFrame

        self.createWidgets()
        self.arrange()
        self.initialize()
        self.backupLineValues()

    def createWidgets(self):
        """Create all widgets needed for this panel"""

        # create all the headings
        self.headings = []
        for eachLbl in self.headingsLbl:
            self.headings.append(wx.StaticText(self, -1, eachLbl))

        # a color button
        color = mpl_color('b')
        self.colorButton = colorsel.ColourSelect(self, -1, '', color)

        # choice box for line style
        choices = graphInfo.lineStyleDescriptions
        self.lineStyleChoice = wx.Choice(self, -1, choices=choices)

        # spinner for line thickness
        self.thickSpin = wx.SpinCtrl(self, -1)

        # choice box for line symbol
        choices = graphInfo.lineSymbolDescriptions
        self.lineSymbolChoice = wx.Choice(self, -1, choices=choices)

        # spinner for symbol size
        self.symbolSizeSpin = wx.SpinCtrl(self, -1)

    def arrange(self):
        """Arrange widgets of this panel"""

        cols = 5
        rows = 1
        mainSizer = wx.FlexGridSizer(rows, cols, 5, hgap=15)

        # add headings
        style = wx.ALIGN_CENTER_HORIZONTAL
        for eachHeading in self.headings:
            mainSizer.Add(eachHeading, 0, style)

        # arrange widgets about line options
        style = wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL
        mainSizer.Add(self.colorButton, 0, style)
        mainSizer.Add(self.lineStyleChoice, 0, style)
        mainSizer.Add(self.thickSpin, 0, style)
        mainSizer.Add(self.lineSymbolChoice, 0, style)
        mainSizer.Add(self.symbolSizeSpin, 0, style)

        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Initialize widgets of this panel"""

        # since the lines in the whole graph could have different options,
        # search through them, and show the most used options
        if len(self.mainFrame.lines) > 1:
            options = self.getMostUsedOptions()
        else:
            options = {}
            options['color'] = 'b'
            options['linestyle'] = '-'
            options['thickness'] = 2
            options['symbol'] = 'None'
            options['symbol_size'] = 6

        # initialize widget values
        self.colorButton.SetColour(mpl_color(options['color']))
        index = graphInfo.lineStyleTags.index(options['linestyle'])
        self.lineStyleChoice.SetSelection(index)
        index = graphInfo.lineSymbolTags.index(options['symbol'])
        self.lineSymbolChoice.SetSelection(index)
        self.thickSpin.SetRange(0, 100)
        self.thickSpin.SetValue(options['thickness'])
        self.symbolSizeSpin.SetRange(0, 100)
        self.symbolSizeSpin.SetValue(options['symbol_size'])

        # bind events
        self.colorButton.Bind(colorsel.EVT_COLOURSELECT, self.onColor, 
                                                             self.colorButton)
        self.lineStyleChoice.Bind(wx.EVT_CHOICE, self.onLineStyleChoice, 
                                                         self.lineStyleChoice)
        self.thickSpin.Bind(wx.EVT_SPINCTRL, self.onThickness, self.thickSpin)
        self.lineSymbolChoice.Bind(wx.EVT_CHOICE, self.onLineSymbolChoice, 
                                                        self.lineSymbolChoice)
        self.symbolSizeSpin.Bind(wx.EVT_SPINCTRL, self.onSymbolSize, 
                                                          self.symbolSizeSpin)

    def getMostUsedOptions(self):
        """Return a dict of the most used options"""

        def increment(d, line, attrName):
            # increment the score of a certain attribute in a dict
            attr = getattr(line, attrName)
            if d.has_key(attr):
                d[attr] += 1
            else:
                d[attr] = 1

        def getMaxKey(d):
            # given a dict, return key with the highest value
            maxKey = None
            maxValue = None
            for key, value in d.iteritems():
                if maxKey is None:
                    maxKey = key
                if maxValue is None:
                    maxValue = value
                if value > maxValue:
                    maxKey = key
                    maxValue = value

            return maxKey

        # count up how many times options show up
        colors = {}
        linestyle = {}
        thickness = {}
        symbol = {}
        symbol_size = {}
        for eachLine in self.mainFrame.lines:
            increment(colors, eachLine, 'color')
            increment(linestyle, eachLine, 'linestyle')
            increment(thickness, eachLine, 'thickness')
            increment(symbol, eachLine, 'marker')
            increment(symbol_size, eachLine, 'markersize')

        # collect the most used options
        results = {}
        results['color'] = getMaxKey(colors)
        results['linestyle'] = getMaxKey(linestyle)
        results['thickness'] = getMaxKey(thickness)
        results['symbol'] = getMaxKey(symbol)
        results['symbol_size'] = getMaxKey(symbol_size)

        return results

    def onColor(self, event):
        """Event when a color is selected"""

        color = hexcolor(event.GetValue())
        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.set_color(color)
        self.update()

    def onThickness(self, event):
        """Event when the thickness spinner changes"""

        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.set_linewidth(self.thickSpin.GetValue())
        self.update()

    def onLineStyleChoice(self,event):
        """Event when a line style was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineStyleDescrToTag[event.GetString()]
        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.set_linestyle(tag)
        self.update()

    def onLineSymbolChoice(self,event):
        """Event when a line symbol was chosen"""
        
        # change the description to the right command
        tag = graphInfo.lineSymbolDescrToTag[event.GetString()]
        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.set_marker(tag)
        self.update()

    def onSymbolSize(self, event):
        """Event when the symbol size spinner changes value"""

        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.set_markersize(self.symbolSizeSpin.GetValue())
        self.update()

    def backupLineValues(self):
        """Create a backup of the present line values""" 

        # backup line values
        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.backup()

    def restoreLineValues(self):
        """Restore the line options to before anything was changed"""

        for eachSet in self.mainFrame.lines.lineSets:
            eachSet.restore()

        self.update()

    def update(self):
        """Updates the the graph"""

        self.mainFrame.updateLegend(False)
        self.mainFrame.draw()
