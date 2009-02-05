"""Dialog box asking about format of gene marker file"""

import os
import wx
import wx.lib.colourselect as colorsel

import colors
from config import graphInfo
from config import misc
from widgets.geneMarkers import loadGeneMarkerFile

def getDefaultParameters(filename):
    """Returns a dictionary containing the default parameters needed to
    load a gene marker file.
    
    """

    info = {}

    # the filename containing the gene markers
    info['filename'] = filename

    # the column that indicated the chromosome number
    info['chr_col'] = 0

    # the column that has the gene marker position
    info['pos_col'] = 1

    # name given to this group of gene markers
    info['name'] = ''

    # the size of the markers
    info['size'] = graphInfo.markers_size

    # the user specified color, in hex format
    info['color'] = graphInfo.markers_color

    # the symbol of the markers
    info['symbol'] = graphInfo.markers_symbol

    return info

class InsertGeneMarkerDialog(wx.Dialog):
    """A dialog box that determines the format of a file containing gene marker
    positions.  The dialog box asks the user to choose what column the
    chromosome number is listed under and what column the gene positions are 
    listed under.  

    """

    TITLE = 'Please Specify Gene Marker File Format'
    previewLbl = 'Preview '
    chrLbl = 'Chromosome Column:'
    posLbl = 'Gene Position Column:'
    nameLbl = 'Name of Gene Marker'
    sizeLbl = 'Size'
    colorLbl = 'Color'
    symbolLbl = 'Symbol'

    numPreviewLines = 5

    def __init__(self, parent, ID, filename, size=misc.fileformatDialogSize, 
                 pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE):

        self.parent = parent
        self.filename = filename
        self.info = getDefaultParameters(self.filename)
        self.previewLbl += 'for ' + os.path.basename(filename) + ':' 
        self.columns = self.getNumColumns()

        pre = wx.PreDialog()
        pre.Create(parent, ID, self.TITLE, pos, size, style)
        self.PostCreate(pre)

        self.createWidgets()
        self.arrange()
        self.initialize()

    def createWidgets(self):
        """Creates widgets for the dialog box"""

        # the ok and cancel buttons
        self.okButton = wx.Button(self, wx.ID_OK)
        self.okButton.SetDefault()
        self.cancelButton = wx.Button(self, wx.ID_CANCEL)

        # some text labels
        self.previewLabel = wx.StaticText(self, -1, self.previewLbl)
        self.chrLabel = wx.StaticText(self, -1, self.chrLbl)
        self.posLabel = wx.StaticText(self, -1, self.posLbl)
        self.nameLabel = wx.StaticText(self, -1, self.nameLbl)
        self.sizeLabel = wx.StaticText(self, -1, self.sizeLbl)
        self.colorLabel = wx.StaticText(self, -1, self.colorLbl)
        self.symbolLabel = wx.StaticText(self, -1, self.symbolLbl)

        # this text box shows a preview of the file, a couple of lines
        self.previewBox = wx.StaticText(self, -1, label=self.getPreviewText())

        # choices boxes to pick the column for chr, pos, 
        columnChoices = [str(x) for x in range(1, self.columns+1, 1)]
        self.chrChoice = wx.Choice(self, -1, choices=columnChoices)
        self.posChoice = wx.Choice(self, -1, choices=columnChoices)

        # text boxe for the name of the markers
        self.nameBox = wx.TextCtrl(self, -1)

        # spinner for symbol size
        self.sizeSpin = wx.SpinCtrl(self, -1)
        
        # button to determine the color
        color = colors.hex2rgb(self.info['color'])
        self.colorButton = colorsel.ColourSelect(self, -1, '', color)

        # choice box for marker symbol
        choices = graphInfo.lineSymbolDescriptions
        self.symbolChoice = wx.Choice(self, -1, choices=choices)

    def arrange(self):
        """Arranges widgets"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # everything is added to this sizer because padding is needed
        # around the whole dialog box, so this sizer is added to the
        # main sizer with some padding
        otherSizer = wx.BoxSizer(wx.VERTICAL)

        # a preview of the file is at the top
        otherSizer.Add(self.previewLabel, 0, wx.ALL, 5)
        otherSizer.Add(self.previewBox, 1, wx.EXPAND | wx.ALL, 5)

        # horizontal line as divider
        otherSizer.Add(wx.StaticLine(self, style=wx.LI_HORIZONTAL), 0, 
                       wx.EXPAND | wx.ALL, 5)

        # this sizer contains all the choice boxes and radio boxes,
        # which need to be centered horizontally on the dialog box
        optionsSizer = wx.BoxSizer(wx.VERTICAL)

        # next is one line containing the chr and position column choice boxes
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        chr_pos_sizer = wx.BoxSizer(wx.HORIZONTAL)
        chr_pos_sizer.Add(self.chrLabel, 0, style, 5)
        chr_pos_sizer.Add(self.chrChoice, 0, style, 5)
        chr_pos_sizer.Add(self.posLabel, 0, style, 5)
        chr_pos_sizer.Add(self.posChoice, 0, style, 5)
        optionsSizer.Add(chr_pos_sizer, 0)

        # some empty space
        optionsSizer.Add((1, 10))

        # next is headings and the widgets to control options below headings
        cols = 4
        rows = 2
        gridSizer = wx.FlexGridSizer(rows, cols, vgap=7, hgap=15)

        # add headings
        style = wx.ALIGN_CENTER_HORIZONTAL
        gridSizer.Add(self.nameLabel, 0, style)
        gridSizer.Add(self.sizeLabel, 0, style)
        gridSizer.Add(self.colorLabel, 0, style)
        gridSizer.Add(self.symbolLabel, 0, style)

        # add widgets
        style = wx.ALIGN_CENTER_HORIZONTAL | wx.ALIGN_CENTER_VERTICAL
        gridSizer.Add(self.nameBox, 0, style)
        gridSizer.Add(self.sizeSpin, 0, style)
        gridSizer.Add(self.colorButton, 0, style)
        gridSizer.Add(self.symbolChoice, 0, style)
        optionsSizer.Add(gridSizer, 0, wx.LEFT | wx.TOP, 5)

        otherSizer.Add(optionsSizer, 0, wx.ALIGN_CENTER_HORIZONTAL)

        # some empty space
        otherSizer.Add((1, 20))

        # add the ok and cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        btnsizer.AddButton(self.okButton)
        btnsizer.AddButton(self.cancelButton)
        btnsizer.Realize()

        style = wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | wx.ALL
        otherSizer.Add(btnsizer, 0, style, 5)
        mainSizer.Add(otherSizer, 0, wx.ALL, 10)
        self.SetSizer(mainSizer)
        self.Fit()

    def initialize(self):
        """Give widgets a default value"""

        if self.chrChoice.GetCount() <= self.info['chr_col']:
            self.info['chr_col'] = self.chrChoice.GetCount()-1
        self.chrChoice.SetSelection(self.info['chr_col'])

        if self.posChoice.GetCount() <= self.info['pos_col']:
            self.info['pos_col'] = self.posChoice.GetCount()-1
        self.posChoice.SetSelection(self.info['pos_col'])

        self.sizeSpin.SetValue(self.info['size'])
        index = graphInfo.lineSymbolTags.index(self.info['symbol'])
        self.symbolChoice.SetSelection(index)

        # bind events
        self.chrChoice.Bind(wx.EVT_CHOICE, self.onChrChoice, self.chrChoice)
        self.posChoice.Bind(wx.EVT_CHOICE, self.onPosChoice, self.posChoice)
        self.nameBox.Bind(wx.EVT_TEXT, self.onNameBox, self.nameBox)
        self.sizeSpin.Bind(wx.EVT_SPINCTRL, self.onSize, self.sizeSpin)
        self.colorButton.Bind(colorsel.EVT_COLOURSELECT, self.onColor, 
                                                             self.colorButton)
        self.symbolChoice.Bind(wx.EVT_CHOICE, self.onSymbolChoice, 
                               self.symbolChoice)
        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)

    def onChrChoice(self, event):
        """Event when the chr column choice is selected"""
        self.info['chr_col'] = int(self.chrChoice.GetStringSelection()) - 1

    def onPosChoice(self, event):
        """Event when the position column choice is selected"""
        self.info['pos_col'] = int(self.posChoice.GetStringSelection()) - 1

    def onNameBox(self, event):
        """Event when the name text box is used"""
        self.info['name'] = self.nameBox.GetValue()

    def onSize(self, event):
        """Event when the symbol size spinner changes value"""
        self.info['size'] = self.sizeSpin.GetValue()

    def onColor(self, event):
        """Event when a color is selected"""
        self.info['color'] = colors.hexcolor(event.GetValue())

    def onSymbolChoice(self, event):
        """Event when the symbol choice is selected"""
        selection = self.symbolChoice.GetStringSelection()
        tag = graphInfo.lineSymbolDescrToTag[selection]
        self.info['symbol'] = tag

    def onOK(self, event):
        """Event when the ok button is pressed"""

        try:
            loadGeneMarkerFile(self.parent, self.info)
        except Exception, e:
            m = "Error loading gene markers: %s" % str(e)
            print m
            dlg = wx.MessageDialog(self, m,"Error!", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            event.Skip()

    def getPreviewText(self):
        """Returns the first couple of lines of the file as one string"""

        f = open(self.filename, 'r')
        s = ''
        lines = f.readlines()
        i = 0
        j = 0
        while i < self.numPreviewLines and j < len(lines):
            stripped = lines[j].strip()
            if len(stripped) > 0 and stripped[0] != '#':
                s += lines[j]
                i += 1
            j += 1

        f.close()
        return s

    def getNumColumns(self):
        """Look at one line in the file and determine the number of columns"""

        f = open(self.filename, 'r')
        numColumns = 0
        for eachLine in f.readlines():
            split = eachLine.split()
            if len(split) > 0 and split[0][0] != '#':
                numColumns = len(eachLine.split())
                break

        f.close()
        return numColumns

