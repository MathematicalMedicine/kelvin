"""The panel to edit Kelvin config files"""

import os
import wx
import wx.lib.flatnotebook as FNoteBook
from config import misc
from fileManager.KelvinConfigFileChecker import KelvinConfigFileChecker
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from panel.affectionStatusPanel import AffectionStatusPanel
from panel.analysisTypePanel import AnalysisTypePanel
from panel.chromosomePanel import ChromosomePanel
from panel.commentPanel import CommentPanel
from panel.constraintsPanel import ConstraintsPanel
from panel.diseaseAllelesPanel import DiseaseAllelesPanel
from panel.fileNamesPanel import FileNamesPanel
from panel.liabilityClassPanel import LiabilityClassPanel
from panel.MM_Panel import MM_Panel
from panel.polyEvalPanel import PolyEvalPanel
from panel.scrolledPanel import TabPanel
from panel.simpleTextFile import SimpleTextFilePanel
from panel.TM_Panel import TM_Panel
from panel.traitTypePanel import TraitTypePanel
#from panel.unknownIdPanel import UnknownIdPanel
from panel.valuesPanel import ValuesPanel

class KelvinConfigPanel(wx.Panel):
    """The panel used to edit Kelvin config files"""

    borderAmt = 5
    settingsTabLbl = 'Settings'
    constraintsTabLbl = 'Constraints'
    valuesTabLbl = 'Grid Specification'

    def __init__(self, parent, defaultsManager, fileName='', isNew=False):
        # the fourth argument, isNew, specifies whether the filename is
        # one of the templates used when starting a new config file from
        # scratch, or the filename is actually a file on the disk

        wx.Panel.__init__(self, parent, -1) 
        if isNew:
            self.fileName = ''
        else:
            self.fileName = fileName
        self.fm = KelvinConfigFileManager(fileName, isNew, defaultsManager)

        # a list of functions that will be called when the trait type changes
        self.doWhenTraitChanges = []

        # a list of functions that will be called when the analysis type changes
        self.doWhenAnalysisTypeChanges = []

        # a list of functions that will be called when 
        # the number of liabilty classes change
        self.doWhenLiabilityClassChanges = []

        # a list of functions that will be called when 
        # something in the values panel changes
        self.doWhenValuesChange = []

        # add the main notebook to the panel
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.book = self.createNotebook()
        mainSizer.Add(self.book, 1, wx.EXPAND)
        self.SetSizer(mainSizer)

        # add tabs to notebook
        self.addSettingsTab()
        self.addConstraintsTab()
        self.addValuesTab()

        # panels may need references to other panels, so send them around
        self.setPanelReferences()

        self.book.SetSelection(misc.firstTabFocus)

        # bind events
        self.Bind(wx.EVT_CHOICE, self.onAnalysisTypeChange, 
                  self.analysisTypePanel.analysisChoice)
        self.Bind(wx.EVT_CHOICE, self.onAnalysisTypeChange, 
                  self.analysisTypePanel.sexChoice)
        self.Bind(wx.EVT_CHOICE, self.onTraitTypeChange, 
                  self.traitTypePanel.traitChoice)
        self.Bind(wx.EVT_TEXT, self.onLiabilityClassChange, 
                  self.liabilityPanel.textbox)

        # check the file for serious errors and not-so-serious warnings
        self.fileChecker = KelvinConfigFileChecker(self, self.fm, self.fileName)
        self.checkForErrors()

        # at this point i know the file hasn't been modified, but
        # the filemanager may think it's been modified because of the 
        # way it's been opened, so reset the flag
        # for example, if a line has point and range values and they
        # need to be split into separate lines, the filemanager will
        # think it's been modified
        self.fm.isModified = False
        self.resetTabNames()

    def createNotebook(self):
        book = FNoteBook.FlatNotebook(self, -1)
        style = book.GetWindowStyleFlag()

        # remove old tabs style
        mirror = ~(FNoteBook.FNB_VC71 | FNoteBook.FNB_VC8 | 
                   FNoteBook.FNB_FANCY_TABS)
        style &= mirror

        # set style
        #style |= FNoteBook.FNB_VC8
        style |= FNoteBook.FNB_VC71
        style |= FNoteBook.FNB_NO_X_BUTTON
        style |= FNoteBook.FNB_NO_NAV_BUTTONS
        style |= FNoteBook.FNB_NODRAG

        book.SetWindowStyleFlag(style)

        book.SetNonActiveTabTextColour((0.5, 0.5, 0.5))

        return book

    def createSettingsPanel(self, parent):
        panel = TabPanel(parent)

        # create the subpanels
        self.analysisTypePanel = AnalysisTypePanel(panel, self.fm)
        self.traitTypePanel = TraitTypePanel(panel, self.fm)
        self.chromosomePanel = ChromosomePanel(panel, self.fm)
        #self.unknownIdPanel = UnknownIdPanel(panel, self.fm)
        self.affectionStatusPanel = AffectionStatusPanel(panel, self.fm)
        self.liabilityPanel = LiabilityClassPanel(panel, self.fm)
        self.diseaseAllelesPanel = DiseaseAllelesPanel(panel, self.fm)
        self.polyEvalPanel = PolyEvalPanel(panel, self.fm)
        self.TM_Panel = TM_Panel(panel, self.fm)
        self.MM_Panel = MM_Panel(panel, self.fm, self.traitTypePanel, 
                                 self.affectionStatusPanel)
        self.commentPanel = CommentPanel(panel, self.fm)
        self.fileNamesPanel = FileNamesPanel(panel, self.fm)

        # arrange them
        # this is the main sizer
        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        # sizer for the left side of the panel
        sizer = wx.BoxSizer(wx.VERTICAL)

        # analysis type panel
        # TODO: change the border later, right now there 
        # is no border on the right...
        sizer.Add(self.analysisTypePanel, 0, wx.LEFT | wx.TOP | wx.BOTTOM, 
                  self.borderAmt)

        sizer.Add(self.chromosomePanel, 0, wx.LEFT, self.borderAmt)

        # marker to marker analysis panel
        sizer.Add(self.MM_Panel, 0, wx.LEFT | wx.UP, self.borderAmt)

        # trait type panel
        # TODO: change the border amount to make it more uniform like
        sizer.Add(self.traitTypePanel, 0, wx.ALL, self.borderAmt*2)

        # affection status panel
        sizer.Add(self.affectionStatusPanel, 0, wx.ALL, self.borderAmt)

        # add settings on the left
        mainSizer.Add(sizer, 0, wx.ALL, self.borderAmt)

        # add a divider
        mainSizer.Add(wx.StaticLine(panel, style=wx.LI_VERTICAL), 0, 
                      wx.EXPAND | wx.TOP | wx.BOTTOM, self.borderAmt)

        # sizer for right side
        sizer = wx.BoxSizer(wx.VERTICAL)

        # unknown person id panel
        #sizer.Add(self.unknownIdPanel, 0, wx.ALL, self.borderAmt)

        # disease alleles panel
        sizer.Add(self.diseaseAllelesPanel, 0, wx.ALL, self.borderAmt*2)

        # liability classes panel
        sizer.Add(self.liabilityPanel, 0, wx.ALL, self.borderAmt)

        # polynomial evaluation panel
        sizer.Add(self.polyEvalPanel, 0, wx.ALL, self.borderAmt*2)

        # TM directive
        sizer.Add(self.TM_Panel, 0, wx.ALL, self.borderAmt*2)

        # comment option checkbox
        sizer.Add(self.commentPanel, 0, wx.ALL, self.borderAmt*2)

        # horizontal divider
        sizer.Add(wx.StaticLine(panel, style=wx.LI_HORIZONTAL), 0, 
                      wx.EXPAND | wx.RIGHT, self.borderAmt)

        # add filenames panel
        sizer.Add(self.fileNamesPanel, 1, wx.EXPAND | wx.ALL, self.borderAmt)

        mainSizer.Add(sizer, 1)

        panel.SetSizer(mainSizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        return panel

    def addSettingsTab(self):
        self.settingsPanel = self.createSettingsPanel(self.book)
        self.book.AddPage(self.settingsPanel, self.settingsTabLbl, True)

    def addConstraintsTab(self):
        panel = TabPanel(self.book)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.constraintsPanel = ConstraintsPanel(panel, self.fm)
        style = wx.ALL
        sizer.Add(self.constraintsPanel, 0, style, self.borderAmt*3)
        panel.SetSizer(sizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()
        self.book.AddPage(panel, self.constraintsTabLbl, True)

    def addValuesTab(self):
        panel = TabPanel(self.book)
        sizer = wx.BoxSizer(wx.VERTICAL)
        self.valuesPanel = ValuesPanel(panel, self.fm)
        style = wx.ALL
        sizer.Add(self.valuesPanel, 0, style, self.borderAmt*3)
        panel.SetSizer(sizer)
        panel.SetAutoLayout(1)
        panel.SetupScrolling()

        # make it invisible to the user
        panel.Show(False)
        #self.book.AddPage(panel, self.valuesTabLbl, True)

    def setPanelReferences(self):
        """Give references of some panels to ones that need them"""
        self.analysisTypePanel.setValuesPanel(self.valuesPanel)
        self.valuesPanel.setAnalysisTypePanel(self.analysisTypePanel)

    def save(self):
        self.fm.save()
        self.postSave()

    def saveAs(self, fileName):
        self.fileName = fileName
        self.fileChecker.fileName = fileName
        self.fm.saveAs(fileName)
        self.postSave()

    def postSave(self):
        """Things to do after saving a file"""

        # check for errors and show warnings about them
        self.fileNamesPanel.checkFileExistence()
        self.checkForErrors()
        self.resetTabNames()

    def resetTabNames(self):
        """Set Tab names to normal"""

        # if there are any astericks in front of tab names because
        # something was changed in that tab, take out asterick
        numPages = self.book.GetPageCount()
        for i in range(numPages):
            text = self.book.GetPageText(i)
            if text[0] == '*':
                text = text[1:]
                self.book.SetPageText(i, text)

    def hasChanged(self, panel):
        """Given a panel where some value in the panel has changed, 
        show that the tab has changed.
        
        """

        # change the name on the tab so that the name starts with an asterick,
        # which indicates that there are unsaved changes in that tab
        i = self.book.GetPageIndex(panel)
        if i is not None and i >= 0:
            text = self.book.GetPageText(i)
            if text[0] != '*':
                self.book.SetPageText(i, '*' + text)
            self.fm.isModified = True

    def checkForErrors(self):
        """Checks the file for errors and warnings."""

        self.fileChecker.checkFileForWarnings()
        self.fileChecker.checkFileForErrors()

    def onAnalysisTypeChange(self, event):
        """Call all functions that previously registered 
        to be notified of an analysis type change.
        
        """
        
        for eachFunction in self.doWhenAnalysisTypeChanges[:]:
            try:
                eachFunction(event)
            except wx.PyDeadObjectError, e:
                # delete the function
                self.doWhenAnalysisTypeChanges.remove(eachFunction)

    def onTraitTypeChange(self, event):
        """Call all functions that previously registered to be 
        notified of a trait type change
        
        """

        for eachFunction in self.doWhenTraitChanges[:]:
            try:
                eachFunction(event)
            except wx.PyDeadObjectError, e:
                # delete the function
                self.doWhenTraitChanges.remove(eachFunction)

    def onLiabilityClassChange(self, event):
        """Call all functions that previously registered to be 
        notified of a liability class change
        
        """
        
        for eachFunction in self.doWhenLiabilityClassChanges[:]:
            try:
                eachFunction(event)
            except wx.PyDeadObjectError, e:
                # delete the function
                self.doWhenLiabilityClassChanges.remove(eachFunction)

    def onValuesChange(self, event):
        """Call all functions that previously registered to be notified
        when something in the values panel changes.
        
        """
        
        for eachFunction in self.doWhenValuesChange[:]:
            try:
                eachFunction(event)
            except wx.PyDeadObjectError, e:
                # delete the function
                self.doWhenValuesChange.remove(eachFunction)

    def bindAnalysisTypeChange(self, function):
        """Register a function to be called when the analysis type changes"""
        self.doWhenAnalysisTypeChanges.append(function)

    def bindTraitTypeChange(self, function):
        """Register a function to be called when the trait type changes"""
        self.doWhenTraitChanges.append(function)

    def bindLiabilityClassChange(self, function):
        """Register a function to be called when the liability class changes"""
        self.doWhenLiabilityClassChanges.append(function)

    def bindValuesChange(self, function):
        """Register a function to be called when the values panel changes"""
        self.doWhenValuesChange.append(function)

    def isSexAverage(self):
        """Returns true if sex average is currently chosen, else False if sex
        specific is chosen. Asks the analysis type panel about it. Used in the
        theta constraints panel.
        
        """

        return self.analysisTypePanel.isSexAverage()
