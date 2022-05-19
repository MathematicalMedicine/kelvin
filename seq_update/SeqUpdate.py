#!/usr/bin/env python
"""
This is a front-end to the Calculate Sequentially Updated PPL program that
accompanies Kelvin, implemented as a dialog box that assembles a series of
command line arguments based on user choices.

"""

# Author: Sang-Cheol Seok
# Converted to standalone form by Jo Valentine-Cooper
#
# Copyright (C) 2011, 2015, 2022 Mathematical Medicine LLC
# This program is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the Free
# Software Foundation, either version 3 of the License, or (at your option)
# any later version.
# You should have received a copy of the GNU General Public License along
# with this program. If not, see <https://www.gnu.org/licenses/>.

import os
import wx
import subprocess


class SequpApp(wx.App):
    
    def OnInit(self):
        info = SeqUpdatePPLInfo()
        info['proceed'] = False
        dlg = SeqUpdatePPLDialog(info)
        self.SetTopWindow(dlg)
        return True

class SeqUpdatePPLDialog(wx.Frame):
    """A dialog box that options for calc_updated_ppl """
    
    #Labels
    inputFilesLbl = 'Input files (space or comma separated): '
    
    optionLbl = 'Options'
    multiPtLbl = 'Multipoint '
    twoPtLbl = 'Twopoint '
    
    fixedGridLbl = 'From Fixed grid '
    sexSpecificLbl = 'Sex Specific '
    comMarkerNamesLbl = 'Suppress comparing marker names across input files'
    priorProbLbl = 'set linkage prior probability to '
    cutoffLbl = 'set small-Theta cutoff to '
    weightLbl = 'set small-Theta weight to '
    mapfileLbl = 'use mapfile to order markers :'
    partfileLbl ='write updated Bayes Ratios to partfile'
    
    outFileLbl = 'Output file name: '
    searchLbl = 'Search'
    
    readFromDlg = False #flag to show inputfile names are read from dlg
    
    def __init__(self, info):
        
        wx.Frame.__init__(self, None, wx.ID_ANY, title='Please Specify Options')
        
        self.info = info
        
        self.createWidgets()
        self.arrange()
        self.initialize()
        self.checkForProgram()
    
    def checkForProgram(self):
        def which(program):
            def is_exe(fpath):
                return os.path.exists(fpath) and os.access(fpath, os.X_OK)
            
            def ext_candidates(fpath):
                yield fpath
                for ext in os.environ.get("PATHEXT", "").split(os.pathsep):
                    yield fpath + ext
            
            fpath, fname = os.path.split(program)
            if fpath:
                if is_exe(program):
                    return program
            else:
                for path in os.environ["PATH"].split(os.pathsep):
                    exe_file = os.path.join(path, program)
                    for candidate in ext_candidates(exe_file):
                        if is_exe(candidate):
                            return candidate
            
            return None
        
        def _program_check_report():
            if returnvalue == None:
                dlg = wx.MessageDialog(self,
                        'Kelviz uses {} as the main computational engine '
                        'to calculate sequentially updated PPL. Please make it '
                        'callable. That is, download kelvin from '
                        'kelvin.mathmed.org and Set PATH for it in '
                        'seq_update/'.format(seqBinary),
                        'Error when Calculating Sequentially updated PPL  ',
                        style=wx.OK)
                if dlg.ShowModal() == wx.ID_CANCEL:
                    # there is a patently ludicrous bug in wxPython 3 on OSX
                    # that will automatically close dialog boxes opened in
                    # frame init with ID_CANCEL, whether or not the CANCEL
                    # button even exists! so if that's the response we got,
                    # then it came from this bug and the user did not see our
                    # dialog.
                    _program_check_report()
                else:
                    self.Destroy()
            else:
                self.Show()
        
        seqBinary = 'calc_updated_ppl'
        returnvalue = which(seqBinary)
        _program_check_report()
    
    
    def createWidgets(self):
        """Creates widgets for the dialog box"""
        
        # the ok and cancel buttons
        self.okButton = wx.Button(self, wx.ID_OK)
        self.okButton.SetDefault()
        self.cancelButton = wx.Button(self, wx.ID_CANCEL)
        
        self.inputFilesButton = wx.Button(self, -1, label = self.searchLbl)
        self.outfileButton = wx.Button(self, -1, label = self.searchLbl)
        self.mapfileButton = wx.Button(self, -1, label = 'Search')
        # text box to input inputfiles
        self.inputFilesBox = wx.TextCtrl(self, -1, size=(400,-1))
        
        # some text labels
        self.inputFilesLabel = wx.StaticText(self, -1, self.inputFilesLbl)
        self.multiPtLabel = wx.StaticText(self, -1, self.multiPtLbl)
        self.twoPtLabel = wx.StaticText(self, -1, self.twoPtLbl)
        self.optionLabel = wx.StaticText(self, -1, self.optionLbl)
        
        self.garbageBox = wx.TextCtrl(self, -1, size=(400,-1))
        
        self.fixedLabel = wx.StaticText(self, -1, self.fixedGridLbl)
        self.sexSpecificLabel = wx.StaticText(self, -1, self.sexSpecificLbl)
        self.comMarkerNamesLabel = wx.StaticText(self, -1,
                self.comMarkerNamesLbl)
        self.priorProbLabel = wx.StaticText(self, -1, self.priorProbLbl)
        self.cutoffLabel = wx.StaticText(self, -1, self.cutoffLbl)
        self.weightLabel = wx.StaticText(self, -1, self.weightLbl)
        self.mapfileLabel = wx.StaticText(self, -1, self.mapfileLbl)
        self.partfileLabel = wx.StaticText(self, -1, self.partfileLbl)
        
        self.outFileLabel = wx.StaticText(self, -1, self.outFileLbl)
        
        
        # radio buttons to choose 2pt or multipt analysis
        self.twoRadio = wx.RadioButton(self, -1, self.twoPtLbl)
        self.multiRadio = wx.RadioButton(self, -1, self.multiPtLbl)
        
        # checkbox to determine if using fixed grid or not
        self.fixedGridCheckBox = wx.CheckBox(self, -1, self.fixedGridLbl)
        
        # checkbox to determine if using sex specific analysis
        self.sexSpecificCheckBox = wx.CheckBox(self, -1, self.sexSpecificLbl)
        
        # checkbox to determine if suppressing comparing marker names accross
        # input files
        self.comMarkerCheckBox = wx.CheckBox(self, -1, self.comMarkerNamesLbl)
        
        
        textBoxSize = (100, -1)
        
        # text box to input linkage prior probability
        self.priorProbBox = wx.TextCtrl(self, -1, size=textBoxSize)
        
        # text boxes to input a cutoff and a weight for small-Theta
        self.cutoffBox = wx.TextCtrl(self, -1, size=textBoxSize)
        self.weightBox = wx.TextCtrl(self, -1, size=textBoxSize)
        
        # text box to input mapfile
        self.mapfileBox = wx.TextCtrl(self, -1, size=textBoxSize)
        
        # text box to input mapfile
        self.partfileBox = wx.TextCtrl(self, -1, size=textBoxSize)
        
        textBoxSize = (200, -1)
        # text box to input mapfile
        self.outfileBox = wx.TextCtrl(self, -1, size=textBoxSize)
        self.garbage2Box = wx.TextCtrl(self, -1, size=(430,30))
                # To cover the garbages showing on the Dialog box!    
        self.garbage2Box.Disable()
    
    def arrange(self):
        """Arranges widgets"""
        
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        
        # everything is added to this sizer because padding is needed
        # around the whole dialog box, so this sizer is added to the
        # main sizer with some padding
        otherSizer = wx.BoxSizer(wx.VERTICAL)
        # some empty space
        otherSizer.Add((1, 20))
        
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        rowSizer.Add(self.inputFilesLabel, 0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        rowSizer.Add(self.inputFilesButton, 0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        otherSizer.Add(rowSizer)
        otherSizer.Add(self.inputFilesBox, 0,
                wx.ALIGN_CENTER_VERTICAL | wx.ALL, 5)
        
        # horizontal line as divider
        divider_line(otherSizer, self)
        
        otherSizer.Add(self.optionLabel, 0, wx.LEFT, 5)
        
        # this sizer contains all the choice boxes and radio boxes,
        # which need to be centered on the dialog box
        optionsSizer = wx.BoxSizer(wx.VERTICAL)
        
        # next is one line containg the position and ppl column, and
        # spinner to specify number of header lines
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        
        # sizer for the radio buttons
        choiceSizer = wx.BoxSizer(wx.VERTICAL)
        choiceSizer.Add(self.multiRadio,0)
        choiceSizer.Add(self.twoRadio,0)
        choiceSizer.Add(self.fixedGridCheckBox, 0, wx.LEFT, 15)
        choiceSizer.Add(self.sexSpecificCheckBox, 0, wx.LEFT, 15)
        choiceSizer.Add(self.comMarkerCheckBox, 0, wx.LEFT, 15)
        
        optionsSizer.Add(choiceSizer)
        
        # sizer for some options that require text boxes
        textBoxSizer = wx.FlexGridSizer(cols=2, vgap=0, hgap=0)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        
        # text boxes asking for prefix and postfix
        textBoxSizer.Add(self.priorProbLabel, 0, style, 5)
        textBoxSizer.Add(self.priorProbBox, 0, style, 5)
        textBoxSizer.Add(self.cutoffLabel, 0, style, 5)
        textBoxSizer.Add(self.cutoffBox, 0, style, 5)
        textBoxSizer.Add(self.weightLabel, 0, style, 5)
        textBoxSizer.Add(self.weightBox, 0, style, 5)
        
        textBoxSizer.Add(self.mapfileLabel, 0, style, 5)
        textBoxSizer.Add(self.mapfileBox, 0, style, 5)
        textBoxSizer.Add(self.partfileLabel, 0, style, 5)
        textBoxSizer.Add(self.partfileBox, 0, style, 5)
        
        optionsSizer.Add(textBoxSizer, 0, wx.LEFT, 15)
        otherSizer.Add(optionsSizer, 0, wx.ALIGN_CENTER_HORIZONTAL)
        
        # horizontal line as divider
        divider_line(otherSizer, self)
        
        # checkbox asking about same color or not
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        rowSizer.Add(self.outFileLabel, 0, style, 5)
        rowSizer.Add(self.outfileButton, 0, style, 5)
        otherSizer.Add(rowSizer)
        otherSizer.Add(self.outfileBox, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALL,
                5)
        
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
        
        self.twoRadio.SetValue(False)
        self.multiRadio.SetValue(True)
        
        self.outfileBox.SetValue(self.info['outfile'])
                #change it to ppl.out for twoPoint analysis
        
        #All options for two point is disabled
        self.fixedGridCheckBox.Disable()
        self.sexSpecificCheckBox.Disable()
        self.comMarkerCheckBox.Disable()
        self.priorProbBox.Disable()
        self.cutoffBox.Disable()
        self.weightBox.Disable()
        self.mapfileBox.Disable()
        self.partfileBox.Disable()
        
        
        # bind events
        self.inputFilesBox.Bind(wx.EVT_TEXT, self.onInputfiles,
                self.inputFilesBox)
        self.twoRadio.Bind(wx.EVT_RADIOBUTTON, self.onSingleMultRadioButton,
                self.twoRadio)
        self.multiRadio.Bind(wx.EVT_RADIOBUTTON,self.onSingleMultRadioButton,
                self.multiRadio)
        
        self.fixedGridCheckBox.Bind(wx.EVT_CHECKBOX, self.onFixedGridCheck,
                self.fixedGridCheckBox)
        self.sexSpecificCheckBox.Bind(wx.EVT_CHECKBOX, self.onSexSpecificCheck,
                self.fixedGridCheckBox)
        self.comMarkerCheckBox.Bind(wx.EVT_CHECKBOX, self.onComMarkerCheck,
                self.fixedGridCheckBox)
        self.priorProbBox.Bind(wx.EVT_TEXT, self.onPrior, self.priorProbBox)
        self.cutoffBox.Bind(wx.EVT_TEXT, self.onCutoff, self.cutoffBox)
        self.weightBox.Bind(wx.EVT_TEXT, self.onWeight, self.weightBox)
        self.mapfileBox.Bind(wx.EVT_TEXT, self.onMapfile, self.mapfileBox)
        self.partfileBox.Bind(wx.EVT_TEXT, self.onPartfile, self.partfileBox)
        self.outfileBox.Bind(wx.EVT_TEXT, self.onOutfile, self.outfileBox)
        
        
        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)
        self.cancelButton.Bind(wx.EVT_BUTTON, self.onCancel, self.cancelButton)
        self.inputFilesButton.Bind(wx.EVT_BUTTON, self.onOpenInFiles,
                self.inputFilesButton)
        self.outfileButton.Bind(wx.EVT_BUTTON, self.onOpenOutFile,
                self.outfileButton)
    
    def onOpenInFiles(self,event):
        
        self.inputFilesBox.SetEditable(True)
        dlg = wx.FileDialog(self, "Open File", os.getcwd(),
                style=wx.OPEN | wx.FD_MULTIPLE)
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        
        filenames = dlg.GetPaths()
        self.info['inputfiles'] = filenames
        #print  filenames
        fls = ', '.join(filenames)
        self.inputFilesBox.SetValue(fls)
        self.readFromDlg = True
        self.inputFilesBox.SetEditable(False)
        #print "on onOpenInFiles filenames are " + fls
    
    def onOpenOutFile(self,event):
        
        default_file = 'br.out'
        self.outfileBox.SetEditable(True)
        # show save dialog to get filename
        dlg = wx.FileDialog(self, "Save as...", os.getcwd(), default_file,
                style=wx.SAVE | wx.OVERWRITE_PROMPT)
        
        if dlg.ShowModal() != wx.ID_OK:
            dlg.Destroy()
            return
        
        filename = dlg.GetPath()
        if os.path.exists(filename):
            os.remove(filename)
        self.outfileBox.SetValue(filename)
        self.outfileBox.SetEditable(False)
        #print "onOpenoutFilesfilename are " + filename
    
    def onFixedGridCheck(self, event):
        """Event when the color checkbox is clicked"""
        self.info['fixed_grid'] = self.fixedGridCheckBox.GetValue()
    
    def onSexSpecificCheck(self, event):
        """Event when the color checkbox is clicked"""
        self.info['sex_specific'] = self.sexSpecificCheckBox.GetValue()
    
    def onComMarkerCheck(self, event):
        """Event when the color checkbox is clicked"""
        self.info['comMarker'] = self.comMarkerCheckBox.GetValue()
    
    def onInputfiles(self,event):
        """Event when the text box to specify input files"""
        # don't need to do anything, will extract this when ok is pressed
        pass
    
    def onPrior(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        self.info['priorProb'] =  self.priorProbBox.GetValue()
    
    def onCutoff(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        self.info['cutoff'] =  self.cutoffBox.GetValue()
    
    def onWeight(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        self.info['weight'] =  self.weightBox.GetValue()
    
    def onMapfile(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        self.info['mapfile'] = self.mapfileBox.GetValue()
    
    def onPartfile(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        self.info['partfile'] = self.partfileBox.GetValue()
    
    def onOutfile(self, event):
        """Event when the text box to specify prior probability distribution"""
        # don't need to do anything, will extract this when ok is pressed
        pass
        #self.info['outfile'] = self.outfileBox.GetValue()
    
    
    def onSingleMultRadioButton(self, event):
        """Event when the single or multiple chromosome
        radio buttons are clicked
        
        """
        if self.twoRadio.GetValue():
            # two point were chosen
            self.info['twoPoint'] = True
            self.info['outfile'] = 'ppl.out'
            self.outfileBox.SetValue(self.info['outfile'])
            self.twoRadio.SetValue(True)
            self.multiRadio.SetValue(False)
            self.fixedGridCheckBox.Enable()
            self.sexSpecificCheckBox.Enable()
            self.comMarkerCheckBox.Enable()
            self.priorProbBox.Enable()
            self.cutoffBox.Enable()
            self.weightBox.Enable()
            self.mapfileBox.Enable()
            self.partfileBox.Enable()
        
        else:
            # multi point were chosen
            self.info['outfile'] = 'br.out'
            self.outfileBox.SetValue(self.info['outfile'])
            self.info['twoPoint'] = False
            self.twoRadio.SetValue(False)
            self.multiRadio.SetValue(True)
            self.fixedGridCheckBox.Disable()
            self.sexSpecificCheckBox.Disable()
            self.comMarkerCheckBox.Disable()
            self.priorProbBox.Disable()
            self.cutoffBox.Disable()
            self.weightBox.Disable()
            self.mapfileBox.Disable()
            self.partfileBox.Disable()
    
    
    def onOK(self, event):
        """Event when the ok button is pressed"""
        
        # parse out some parameters
        if self.readFromDlg:
            #print "outfiles read already from dlg so leaving now"
            #print 'Outfile name read from dlg is ' +  self.info['outfile']
            self.info['proceed'] = True
            event.Skip()
            return
        
        if self.parseInputfilesBox():
            
            self.info['outfile'] = self.outfileBox.GetValue()
            if len(self.info['outfile']) == 0:
                s = "No output filename"
                dlg = wx.MessageDialog(self, s, "Error! in output filename", wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
            
            else:
                #print 'Outfile name is ' +  self.info['outfile']
                self.info['proceed'] = True
                event.Skip()
        
        calc_seq_updated_ppl(info)
    
    def onCancel(self, event):
        self.Destroy()
    
    
    def parseInputfilesBox(self):
        """Parse out text in box for input files
        Multiple files can be specified using spaces or commas.
        
        """
        
        def showError(message):
            title = 'Error in parsing input files'
            dlg = wx.MessageDialog(self, message, title, wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        
        #print "In parsing inputfiles"
        
        if self.inputFilesBox.GetValue() != '':
            # this needs to be cleared out just in case it fails
            # and the user tries again
            self.info['inputfiles'] = []
            
            # replace commas with spaces, and then split it
            s = self.inputFilesBox.GetValue()
            s = s.replace(',', ' ')
            s = s.split()
            for eachValue in s:
                #print eachValue
                #eachValue.replace(' ', '\ ')
                self.info['inputfiles'].append(eachValue)
        
        else:
            showError('No inputfiles')
            return False
        
        # parsing is ok
        return True

class SeqUpdatePPLInfo(dict):
    # class that stores information about the format of a text file
    
    # this is a dictionary describing the parameters needed to
    # call calc_updated_ppl.
    # --- parameters ---
    # inputfiles : files to be used fo seq. updating.
    # twoPoint : two point analysis or multipoint analysis
    # options for twoPoint
    #    fixed_grid :
    #    sex_specific:
    #    comMarker:
    #    priorProb: probability distribution
    #    cutoff:  small-Theta cutoff
    #    weight:  small-Theta weight
    #    mapfile: to order markers
    #    partfile: updated Bayes Ratios
    # outfile: file to write PPLs
    #
    # Default it multi point
    
    def __init__(self):
        self['twoPoint'] = False
        self['fixed_grid'] = False
        self['sex_specific'] = False
        self['comMarker'] = False
        self['priorProb'] = ''
        self['cutoff'] = ''
        self['weight'] = ''
        self['mapfile'] = ''
        self['partfile'] = ''
        self['outfile'] = 'br.out'
        self['inputfiles'] = []


def divider_line(sizer, parent):
    """Adds to our sizer a semistandard horizontal dividing line as used in
    several of our dialogs.
    
    """
    
    sizer.Add(wx.StaticLine(parent, style=wx.LI_HORIZONTAL), 0, 
            wx.EXPAND | wx.ALL, 5)


def calc_seq_updated_ppl(info):
    """Brings up the Calculate Sequentially Updated PPL dialog, assembles a
    series of command-line arguments with the info from same, and calls the
    calc_updated_ppl program.
    
    """
    
    if info['proceed']:
        
        #print info
        args = []
        if info['twoPoint']:
            if info['fixed_grid']:
                args.append('-o')
            if info['sex_specific']:
                args.append('-s')
            if info['comMarker']:
                args.append('-r')
            if info['priorProb']:
                args.append('-p')
                args.append(info['priorProb'])
            if info['cutoff']:
                args.append('-c')
                args.append(info['cutoff'])
            if info['weight']:
                args.append('-w')
                args.append(info['weight'])
            if len(info['mapfile']) > 0:
                args.append('-M')
                args.append(info['mapfile'])
            if len(info['partfile']) > 0:
                args.append('-O')
                args.append(info['partfile'])
        
        else:
            args.append('-m')
        
        args.append('-R')
        
        
        outFile = str(info['outfile'])
        files = [str(x) for x in info['inputfiles']]
        #print files
        
        #Check if the outputfile already exist.
        if os.path.exists(outFile):
            wx.MessageBox('Output file {} already exists. Please provide'
                    ' different name '.format(outFile),
                    'Error when Calculating Sequentially updated PPL  ',
                    style=wx.OK)
            return
        
        commandline = [seqBinary] + args + [outFile] + files
        fixedCommandline = []
        for eachTerm in commandline:
            eachTerm = str(eachTerm)
            eachTerm = eachTerm.strip()
            #eachTerm = eachTerm.replace(' ','\ ')
            fixedCommandline.append(eachTerm)
        
        
        commandlineString = ' '.join(fixedCommandline)
        #print "commandline is " + commandlineString
        
        subprocess.call(fixedCommandline)

if __name__ == "__main__":
    asdf = SequpApp()
    asdf.MainLoop()
    
