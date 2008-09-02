"""Dialog box asking about the format of the data file to be opened"""

import os
import wx
import wx.lib.colourselect as colorsel

from colors import mpl_color, hexcolor
from config import graphInfo
from config import misc

class FileFormatDialog(wx.Dialog):
    """A dialog box that determines the format of opened files.  The dialog box 
    asks the user to choose what column the chromosome position is listed under
    and what column the ppl values are listed under.  It also asks if the file 
    specifies one or more chromosomes, and if more than one, what column 
    that's listed in

    """

    TITLE = 'Please Specify File Format'
    previewLbl = 'Preview '
    posLbl = 'Position Column:'
    pplLbl = 'PPL Column:'
    chrLbl = 'This file has data for:'
    chrNumLbl = '     Number of chromosome:'
    chrColLbl = 'Column that indicates chromosome number:'
    singleLbl = 'Single chromosome'
    multipleLbl = 'Multiple chromosomes'
    allChrLbl = 'Import all chromosomes in file'
    specificChrLbl = 'Import only specified chromosomes (comma separated)'
    overlapLbl = 'Chromosome to begin overlap with:'
    sameColorLbl = 'Keep all imported lines the same color'
    lineSetLbl = 'Name of this set of lines'
    prefixLbl = 'Prefix chromosome number with'
    postfixLbl = 'Postfix chromosome number with'
    metadataLbl = 'Column(s) with extra data'
    headerLbl = 'Number of header lines:'

    numPreviewLines = 5

    def __init__(self, parent, ID, info, file, filename, chrNames=None,
                 size=misc.fileformatDialogSize, pos=wx.DefaultPosition, 
                 style=wx.DEFAULT_DIALOG_STYLE):

        self.parent = parent
        self.file = file
        self.info = info
        self.previewLbl += 'for ' + filename + ':' 

        # flag to indicate if the chromosome being added will need to be
        # overlapped with existing ones in the graph
        self.overlap = info['overlap']
        self.chrNames = chrNames

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
        self.posLabel = wx.StaticText(self, -1, self.posLbl)
        self.pplLabel = wx.StaticText(self, -1, self.pplLbl)
        self.chrLabel = wx.StaticText(self, -1, self.chrLbl)
        self.chrNumLabel = wx.StaticText(self, -1, self.chrNumLbl)
        self.chrColLabel = wx.StaticText(self, -1, self.chrColLbl)
        self.lineSetLabel = wx.StaticText(self, -1, self.lineSetLbl)
        self.prefixLabel = wx.StaticText(self, -1, self.prefixLbl)
        self.postfixLabel = wx.StaticText(self, -1, self.postfixLbl)
        self.metadataLabel = wx.StaticText(self, -1, self.metadataLbl)
        self.headerLabel = wx.StaticText(self, -1, self.headerLbl)

        # this text box shows a preview of the file, a couple of lines
        self.previewBox = wx.StaticText(self, -1, label=self.getPreviewText())

        # choices boxes to pick the column for pos, ppl,
        # and chromosome number or chromosome column
        columnChoices = [str(x) for x in range(1, self.columns+1, 1)]
        self.posChoice = wx.Choice(self, -1, choices=columnChoices)
        self.pplChoice = wx.Choice(self, -1, choices=columnChoices)

        # spinner to specify how many header lines there are.
        # header lines are ignored when processng the file
        self.headerSpin = wx.SpinCtrl(self, -1)
        
        # there might not be a column to specify chromosome number
        columnChoices[0:0] = ['None']
        self.chrColChoice = wx.Choice(self, -1, choices=columnChoices)

        # a human has 23 chromosomes, for now...
        chrChoices = [str(x) for x in range(1, 24, 1)]
        self.chrNumChoice = wx.Choice(self, -1, choices=chrChoices)

        # radio buttons to choose whether the file is about a single
        # chromosome or multiple chromsomes
        self.singleRadio = wx.RadioButton(self, -1, self.singleLbl)
        self.multipleRadio = wx.RadioButton(self, -1, self.multipleLbl)

        # radio buttons to choose whether to extract all chromosomes or
        # only specific ones
        self.allChrRadio = wx.RadioButton(self, -1, self.allChrLbl, 
                                          style=wx.RB_GROUP)
        self.specificChrRadio = wx.RadioButton(self, -1, self.specificChrLbl)

        # text box to specify which chromosome to extract
        self.chrBox = wx.TextCtrl(self, -1, size=(200, -1))

        # checkbox to determine if using same color or not
        self.colorCheckBox = wx.CheckBox(self, -1, self.sameColorLbl)

        textBoxSize = (200, -1)

        # text box to input a line set name
        self.lineSetNameBox = wx.TextCtrl(self, -1, size=textBoxSize)

        # text boxes to input a prefix or postfix for chromosome label
        self.prefixBox = wx.TextCtrl(self, size=textBoxSize)
        self.postfixBox = wx.TextCtrl(self, size=textBoxSize)

        # text box to input columns for metadata
        self.metadataBox = wx.TextCtrl(self, size=textBoxSize)

        # button to determine the color of lines
        color = mpl_color(self.info['color'])
        self.colorButton = colorsel.ColourSelect(self, -1, '', color)

        # text and choice box to determine overlapping chromosomes
        if self.overlap:
            self.overlapLabel = wx.StaticText(self, -1, self.overlapLbl)
            self.overlapChoice = wx.Choice(self, -1, choices=self.chrNames)

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
        # which need to be centered on the dialog box
        optionsSizer = wx.BoxSizer(wx.VERTICAL)

        # next is one line containg the position and ppl column, and
        # spinner to specify number of header lines
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        pos_ppl_sizer = wx.BoxSizer(wx.HORIZONTAL)
        pos_ppl_sizer.Add(self.posLabel, 0, style, 5)
        pos_ppl_sizer.Add(self.posChoice, 0, style, 5)
        pos_ppl_sizer.Add(self.pplLabel, 0, style, 5)
        pos_ppl_sizer.Add(self.pplChoice, 0, style, 5)
        pos_ppl_sizer.Add(self.headerLabel, 0, style, 5)
        pos_ppl_sizer.Add(self.headerSpin, 0, style, 5)
        optionsSizer.Add(pos_ppl_sizer, 0)

        # this area asks what chromosomes are specified
        chrSizer = wx.BoxSizer(wx.HORIZONTAL)
        chrSizer.Add(self.chrLabel, 0, wx.ALL, 5)

        # sizer for the radio buttons
        choiceSizer = wx.BoxSizer(wx.VERTICAL)

        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        rowSizer.Add(self.singleRadio, 0, wx.ALIGN_CENTER_VERTICAL)
        rowSizer.Add(self.chrNumLabel, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 10)
        rowSizer.Add(self.chrNumChoice, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 10)
        choiceSizer.Add(rowSizer)

        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        rowSizer.Add(self.multipleRadio, 0, wx.ALIGN_CENTER_VERTICAL)
        rowSizer.Add(self.chrColLabel, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 10)
        rowSizer.Add(self.chrColChoice, 0, wx.ALIGN_CENTER_VERTICAL|wx.LEFT, 10)
        choiceSizer.Add(rowSizer)

        # sizer for the other two radio buttons
        colSizer = wx.BoxSizer(wx.VERTICAL)
        colSizer.Add(self.allChrRadio, 0, wx.ALL, 3)

        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        rowSizer.Add(self.specificChrRadio, 0)
        rowSizer.Add(self.chrBox, 0, wx.LEFT | wx.ALIGN_CENTER_VERTICAL, 10)
        colSizer.Add(rowSizer, 0, wx.ALL, 3)
        choiceSizer.Add(colSizer, 0, wx.LEFT, 15)

        chrSizer.Add(choiceSizer)
        optionsSizer.Add(chrSizer, 0)

        # sizer for some options that require text boxes
        textBoxSizer = wx.FlexGridSizer(rows=3, cols=2, vgap=0, hgap=0)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL

        # widgets about line set name
        textBoxSizer.Add(self.lineSetLabel, 0, style, 5)
        textBoxSizer.Add(self.lineSetNameBox, 0, style, 5)

        # text boxes asking for prefix and postfix
        textBoxSizer.Add(self.prefixLabel, 0, style, 5)
        textBoxSizer.Add(self.prefixBox, 0, style, 5)
        textBoxSizer.Add(self.postfixLabel, 0, style, 5)
        textBoxSizer.Add(self.postfixBox, 0, style, 5)

        # add in the text box about metadata
        textBoxSizer.Add(self.metadataLabel, 0, style, 5)
        textBoxSizer.Add(self.metadataBox, 0, style, 5)

        optionsSizer.Add(textBoxSizer)

        # checkbox asking about same color or not
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        rowSizer.Add(self.colorCheckBox, 0, style, 5)
        rowSizer.Add(self.colorButton, 0, style, 5)
        optionsSizer.Add(rowSizer)

        if self.overlap:
            # add a line of widgets asking where to start overlap
            overlapSizer = wx.BoxSizer(wx.HORIZONTAL)
            style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
            overlapSizer.Add(self.overlapLabel, 0, style, 5)
            overlapSizer.Add(self.overlapChoice, 0, style, 5)
            optionsSizer.Add(overlapSizer)

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
        
        # default values for info
        self.setDefaultValues()

        self.posChoice.SetSelection(self.info['pos_col']-1)
        self.pplChoice.SetSelection(self.info['ppl_col']-1)
        self.headerSpin.SetValue(self.info['header_lines'])
        self.chrColChoice.SetSelection(self.info['chr_col'])
        self.chrNumChoice.SetSelection(0)

        self.singleRadio.SetValue(False)
        self.multipleRadio.SetValue(True)

        self.allChrRadio.SetValue(True)
        self.specificChrRadio.SetValue(False)

        self.colorCheckBox.SetValue(True)

        self.chrNumChoice.Disable()
        self.chrBox.Disable()

        # bind events
        self.posChoice.Bind(wx.EVT_CHOICE, self.onPosChoice, self.posChoice)
        self.pplChoice.Bind(wx.EVT_CHOICE, self.onPPLChoice, self.pplChoice)
        self.headerSpin.Bind(wx.EVT_SPINCTRL, self.onHeaderSpin,self.headerSpin)
        self.singleRadio.Bind(wx.EVT_RADIOBUTTON, self.onSingleMultRadioButton, 
                              self.singleRadio)
        self.multipleRadio.Bind(wx.EVT_RADIOBUTTON,self.onSingleMultRadioButton,
                                self.multipleRadio)
        self.chrNumChoice.Bind(wx.EVT_CHOICE, self.onChrNumChoice, 
                               self.chrNumChoice)
        self.chrColChoice.Bind(wx.EVT_CHOICE, self.onChrColChoice, 
                               self.chrColChoice)
        self.allChrRadio.Bind(wx.EVT_RADIOBUTTON,self.onAllSpecRadioButton,
                              self.allChrRadio)
        self.specificChrRadio.Bind(wx.EVT_RADIOBUTTON,self.onAllSpecRadioButton,
                                   self.specificChrRadio)
        self.chrBox.Bind(wx.EVT_TEXT, self.onChrText, self.chrBox)
        self.colorCheckBox.Bind(wx.EVT_CHECKBOX, self.onColorCheck, 
                                self.colorCheckBox)
        self.colorButton.Bind(colorsel.EVT_COLOURSELECT, self.onColor, 
                                                             self.colorButton)
        self.lineSetNameBox.Bind(wx.EVT_TEXT, self.onSetLineSetNameText, 
                                                          self.lineSetNameBox)
        self.prefixBox.Bind(wx.EVT_TEXT, self.onPrefixText, self.prefixBox)
        self.postfixBox.Bind(wx.EVT_TEXT, self.onPostfixText, self.postfixBox)
        self.metadataBox.Bind(wx.EVT_TEXT, self.onMetadataText,self.metadataBox)
        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)

        if self.overlap:
            self.overlapChoice.SetSelection(0)
            self.overlapChoice.Bind(wx.EVT_CHOICE, self.onOverlapChoice, 
                                    self.overlapChoice)

    def setDefaultValues(self):
        """Set default values for info"""

        if self.columns > 2:
            self.info['pos_col'] = 2
            self.info['ppl_col'] = 3
            self.info['chr_col'] = 1
        elif self.columns == 2:
            self.info['pos_col'] = 1
            self.info['ppl_col'] = 2
            self.info['chr_col'] = 0
        else:
            # number of columns is 1? good luck with that.
            self.info['pos_col'] = 1
            self.info['ppl_col'] = 1
            self.info['chr_col'] = 0

        self.info['multiple_chr'] = True
        self.info['chr_num'] = 0
        self.info['import_all'] = True
        self.info['import_list'] = []
        self.info['overlap_index'] = 0
        self.info['header_lines'] = 1
        self.info['same_color'] = True
        self.info['color'] = hexcolor(mpl_color(self.info['color']))
        self.info['line_set_name'] = ''
        self.info['prefix'] = ''
        self.info['postfix'] = ''
        self.info['metadata'] = []

    def onPosChoice(self, event):
        """Event when the position column choice is selected"""
        self.info['pos_col'] = int(self.posChoice.GetStringSelection())

    def onPPLChoice(self, event):
        """Event when the ppl column choice is selected"""
        self.info['ppl_col'] = int(self.pplChoice.GetStringSelection())

    def onHeaderSpin(self, event):
        """Event when the spinner for number of header lines is clicked"""
        self.info['header_lines'] = self.headerSpin.GetValue()

    def onSingleMultRadioButton(self, event):
        """Event when the single or multiple chromosome 
        radio buttons are clicked
        
        """

        if self.singleRadio.GetValue():
            # single chromosome were chosen
            self.chrNumChoice.Enable()
            self.chrColChoice.Disable()
            self.allChrRadio.Disable()
            self.specificChrRadio.Disable()
            self.chrBox.Disable()

            self.info['multiple_chr'] = False
            self.info['chr_col'] = 0
            self.info['chr_num'] = int(self.chrNumChoice.GetStringSelection())
        else:
            # multiple chromosome were chosen
            self.chrColChoice.Enable()
            self.chrNumChoice.Disable()
            self.allChrRadio.Enable()
            self.specificChrRadio.Enable()
            if self.allChrRadio.GetValue() == False and \
                    self.specificChrRadio.GetValue() == False:
                self.allChrRadio.SetValue(True)
                self.onAllSpecRadioButton(None)

            self.info['multiple_chr'] = True
            val = self.chrColChoice.GetStringSelection()
            if val.isdigit():
                self.info['chr_col'] = int(val)
            else:
                # the no column choice must have been selected
                self.info['chr_col'] = 0

    def onAllSpecRadioButton(self, event):
        """Event when the all or specific chromosome radio buttons are hit"""

        if self.allChrRadio.GetValue():
            self.chrBox.Disable()
            self.info['import_all'] = True
        else:
            # specific chromosome radio button was chosen
            self.chrBox.Enable()
            self.info['import_all'] = False

    def onChrText(self, event):
        """Event when the text box to specify chromosome is typed in"""
        # don't need to do anything, will extract this when ok is pressed
        pass

    def onColor(self, event):
        """Event when a color is selected"""
        self.info['color'] = hexcolor(event.GetValue())

    def onMetadataText(self, event):
        """Event when the text box to metadata columns is typed in"""
        # don't need to do anything, will extract this when ok is pressed
        pass

    def onSetLineSetNameText(self, event):
        """Event when the line set name text box is typed in"""
        self.info['line_set_name'] = self.lineSetNameBox.GetValue()

    def onPrefixText(self, event):
        """Event when the prefix text box is typed in"""
        self.info['prefix'] = self.prefixBox.GetValue()

    def onPostfixText(self, event):
        """Event when the postfix text box is typed in"""
        self.info['postfix'] = self.postfixBox.GetValue()

    def onOK(self, event):
        """Event when the ok button is pressed"""

        # parse out some parameters
        if not self.parseChrBox():
            return
        if not self.parseMetadataBox():
            return

        event.Skip()

    def onChrNumChoice(self, event):
        """Event when the chromosome number choice is selected"""
        self.info['chr_num'] = int(self.chrNumChoice.GetStringSelection())

    def onChrColChoice(self, event):
        """Event when the chromosome column choice is selected"""

        val = self.chrColChoice.GetStringSelection()
        if val.isdigit():
            self.info['chr_col'] = int(self.chrColChoice.GetStringSelection())
        else:
            self.info['chr_col'] = 0

    def onOverlapChoice(self, event):
        """Event when the overlap choice is selected """
        self.info['overlap_index'] = self.overlapChoice.GetSelection()

    def onColorCheck(self, event):
        """Event when the color checkbox is clicked"""
        self.info['same_color'] = self.colorCheckBox.GetValue()

    def parseChrBox(self):
        """Parse out the text in the box to choose specific chromosomes"""

        if self.multipleRadio.GetValue() and self.specificChrRadio.GetValue():
            # this needs to be cleared out just in case it fails 
            # and the user tries again
            self.info['import_list'] = []

            # take out all spaces, then split it using commas as delmiters
            s = self.chrBox.GetValue()
            s = s.replace(' ', '')
            s = s.split(',')
            for eachValue in s:
                if not eachValue.isdigit():
                    # show some error message
                    m1 = eachValue + ' is not a valid number.'
                    m2 = 'Please use commas to separate numbers.'
                    m = m1 + '\n' + m2
                    title = 'Error in choosing specific chromosomes'
                    dlg = wx.MessageDialog(self, m, title, wx.OK|wx.ICON_ERROR)
                    dlg.ShowModal()
                    dlg.Destroy()
                    return False
                else:
                    self.info['import_list'].append(int(eachValue))

        # parsing is ok
        return True

    def parseMetadataBox(self):
        """Parse out text in box that chooses columns that have extra data"""

        def showError(message):
            title = 'Error in choosing columns with extra data'
            dlg = wx.MessageDialog(self, message, title, wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()


        if self.metadataBox.GetValue() != '':
            # this needs to be cleared out just in case it fails 
            # and the user tries again
            self.info['metadata'] = []

            # take out all spaces, then split it using commas as delimiters
            s = self.metadataBox.GetValue()
            s = s.replace(' ', '')
            s = s.split(',')
            for eachValue in s:
                if not eachValue.isdigit():
                    # show some error message
                    m1 = eachValue + ' is not a valid number.'
                    m2 = 'Please use commas to separate numbers.'
                    m = m1 + '\n' + m2
                    showError(m)
                    return False
                elif int(eachValue) > self.columns:
                    # show some other error message
                    m1 = eachValue + ' is too large a column value.'
                    m2 = str(self.columns) + ' is the maximum column value.'
                    m = m1 + '\n' + m2
                    showError(m)
                    return False
                elif int(eachValue) < 1:
                    # show some other other error message
                    m1 = eachValue + ' is too small a column value.'
                    m2 = '1 is the minimum column value.'
                    m = m1 + '\n' + m2
                    showError(m)
                    return False
                else:
                    self.info['metadata'].append(int(eachValue))

        # parsing is ok
        return True

    def getPreviewText(self):
        """Returns the first couple of lines of the file as one string"""

        s = ''
        self.file.seek(0)
        lines = self.file.readlines()
        i = 0
        j = 0
        while i < self.numPreviewLines and j < len(lines):
            stripped = lines[j].strip()
            if len(stripped) > 0 and stripped[0] != '#':
                s += lines[j]
                i += 1
            j += 1

        self.file.seek(0)
        return s

    def getNumColumns(self):
        """Look at one line in the file and determine the number of columns"""

        for eachLine in self.file.readlines():
            split = eachLine.split()
            if len(split) > 0 and split[0][0] != '#':
                self.file.seek(0)
                return len(eachLine.split())

        # if execution gets to this point then that means the file is blank
        return None
