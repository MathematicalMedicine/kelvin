"""The main window of the application"""

import os
import sys
import wx
import wx.lib.flatnotebook as FNoteBook
from config import misc
from config import KelvinInfo
from utility import menuBar
from panel.KelvinConfig import KelvinConfigPanel
from panel.simpleTextFile import SimpleTextFilePanel
from fileManager.KelvinConfigDefaults import KelvinConfigDefaults
from fileManager.parseError import ParseError
import subprocess

class MainFrame(wx.Frame):
    """The frame (window) that the main application uses"""

    # used to determine what files to filter when opening and saving
    # kelvin config files
    kelvinConfigWildcard = "Kelvin config file (*.conf)|*.conf|" + \
                           "All files (*.*)|*.*"

    def __init__(self):
        wx.Frame.__init__(self, None, -1, title=misc.mainWindowTitle, 
                          size = misc.mainWindowSize,
                          pos = misc.mainWindowPosition)

        # used to keep track of default values of a kelvin config file
        self.kelvinConfigDefaults = KelvinConfigDefaults()

        # the last directory used when opening and saving files
        self.last_dir = os.getcwd()

        # create some components
        self.SetTitle(misc.mainWindowTitle)
        self.mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.menuBar = menuBar.createMenuBar(self, self.menuData())
        self.statusBar = self.CreateStatusBar()
        self.createToolBar()
        self.createNotebook()
        self.SetSizer(self.mainSizer)
        self.Refresh()

        # create an icon if the file is present
        try:
            path = os.path.join(sys.path[0], 'gKelvin.ico')
            file = open(path)
        except IOError, e:
            # file not found, do nothing
            pass
        else:
            # file found, set the icon 
            file.close()
            self.SetIcon(wx.Icon(path, wx.BITMAP_TYPE_ICO))

        # bind some main events
        self.Bind(wx.EVT_CLOSE, self.onCloseWindow)
        self.Bind(wx.EVT_KEY_DOWN, self.onKeyPress)
        self.Bind(FNoteBook.EVT_FLATNOTEBOOK_PAGE_CLOSING, self.onCloseTab)

        # check for file to open immediately
        if hasattr(misc, 'openFile') and misc.openFile:
            self.openKelvinConfigFile(misc.openFile)

        # check if files were given in the command line
        if len(sys.argv) > 1:
            # assume all arguments are kelvin config filenames
            for eachFilename in sys.argv[1:]:
                self.openKelvinConfigFile(eachFilename)
            self.book.SetSelection(0)

    def menuData(self):
        """Returns a list that defines the structure of the menu bar"""
        kelvinConfigText = "Create New Kelvin config File"
        return [("File", (
                    ("&New Kelvin Config File\tCtrl+N", kelvinConfigText, 
                                                        self.onNewKelvinConfig),
                    ("&Open Kelvin Config File\tCtrl+O", 
                     "Open Kelvin Config File", self.onOpen),
                    ("&Open Other", "Open File", self.onOpenOther),
                    ("&Save\tCtrl+S", "Save File", self.onSave),
                    ("Save &As", "Save File As", self.onSaveAs),
                    ("&Close", "Close File", self.onClose),
                    ("", "", ""),
                    ("&Quit\tCtrl+Q", "Quit", self.onCloseWindow)))]

                # since we are giving this out to other people and i don't
                # think anyone uses this, take out the run command
                #("Run", (
                #    ("Run &Kelvin\tCtrl+R", "Run Kelvin with config File", 
                #                                             self.onRunKelvin),
                #   ("Run PedCheck Levels 0 and 1", "Run PedCheck Levels 0 & 1",
                #                                           self.onPedCheck_01),
                #    ("Run PedCheck Level 2", "Run PedCheck Level 2", 
                #                                        self.onPedCheck_2)))]

    def createNotebook(self):
        self.book = FNoteBook.FlatNotebook(self, -1)
        self.book.SetNonActiveTabTextColour((0.5, 0.5, 0.5))
        style = self.book.GetWindowStyleFlag()
        style |= FNoteBook.FNB_NO_NAV_BUTTONS
        style |= FNoteBook.FNB_VC8
        style |= FNoteBook.FNB_NODRAG
        style |= FNoteBook.FNB_NO_X_BUTTON
        self.book.SetWindowStyleFlag(style)
        self.mainSizer.Add(self.book, 1, wx.EXPAND)

    def createToolBar(self):
        self.toolbar = self.CreateToolBar(wx.TB_HORIZONTAL)

        # get some bitmaps
        new_bmp = wx.ArtProvider.GetBitmap(wx.ART_NEW, wx.ART_TOOLBAR)
        open_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_OPEN, wx.ART_TOOLBAR)
        save_bmp = wx.ArtProvider.GetBitmap(wx.ART_FILE_SAVE, wx.ART_TOOLBAR)
        #run_bmp = wx.ArtProvider.GetBitmap(wx.ART_GO_FORWARD, wx.ART_TOOLBAR)

        # add buttoms
        newId = wx.NewId()
        openId = wx.NewId()
        saveId = wx.NewId()
        #runId = wx.NewId()
        self.toolbar.AddLabelTool(newId, "New", new_bmp, shortHelp="New", 
                                  longHelp="New Kelvin Config File")
        self.toolbar.AddLabelTool(openId, "Open", open_bmp, shortHelp="Open", 
                                  longHelp="Open Kelvin Config File")
        self.toolbar.AddLabelTool(saveId, "Save", save_bmp, shortHelp="Save", 
                                  longHelp="Save Current File")
        #self.toolbar.AddLabelTool(runId, "Run Kelvin", run_bmp, 
        #                          shortHelp="Run Kelvin", 
        #                          longHelp="Run Kelvin on current file")

        self.Bind(wx.EVT_TOOL, self.onNewKelvinConfig, id=newId)
        self.Bind(wx.EVT_TOOL, self.onOpen, id=openId)
        self.Bind(wx.EVT_TOOL, self.onSave, id=saveId)
        #self.Bind(wx.EVT_TOOL, self.onRunKelvin, id=runId)

        # render toolbar
        self.toolbar.Realize()

    def onNewKelvinConfig(self, event):
        """Event when a new kelvin config file is created"""

        # ask user what type of new file to use
        choices = [entry[0] for entry in KelvinInfo.newConfigFile]
        dlg = wx.SingleChoiceDialog(self, 'Please specify type of config file',
                                    'Kelvin Config File Type', choices,
                                    wx.CHOICEDLG_STYLE)
        val = dlg.ShowModal()
        if val == wx.ID_OK:
            self.Freeze()
            i = dlg.GetSelection()
            self.openKelvinConfigFile(KelvinInfo.newConfigFile[i][1], True)
            self.Thaw()

        dlg.Destroy()

    def onKeyPress(self, event):
        if event.GetKeyCode() == 27:
            # esc key was pressed
            self.onCloseWindow(event)

    def onOpen(self, event):
        dlg = wx.FileDialog(self, "Open File", self.last_dir, style=wx.OPEN, 
                            wildcard=self.kelvinConfigWildcard)
        if dlg.ShowModal() == wx.ID_OK:
            self.Freeze()
            self.openKelvinConfigFile(dlg.GetPath())
            self.Thaw()
            self.last_dir = os.path.abspath(os.path.dirname(dlg.GetPath()))
        dlg.Destroy()

    def onOpenOther(self, event):
        dlg = wx.FileDialog(self, "Open File", self.last_dir, style=wx.OPEN)
                            
        if dlg.ShowModal() == wx.ID_OK:
            self.openSimpleFile(dlg.GetPath())
            self.last_dir = os.path.abspath(os.path.dirname(dlg.GetPath()))
        dlg.Destroy()

    def onSave(self, event):
        # get the currently selected page
        currentPanel = self.book.GetCurrentPage()
        if not currentPanel: 
            # no file open at all
            return

        if not currentPanel.fileName:
            self.onSaveAs(event)
        else:
            try:
                currentPanel.save()
            except IOError, e:
                s = "Error saving file '%s' : %s" % \
                        (currentPanel.fileName, e.strerror)
                dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()

    def onSaveAs(self, event):
        dlg = wx.FileDialog(self, "Save as...", self.last_dir,
                            style=wx.SAVE | wx.OVERWRITE_PROMPT, 
                            wildcard = self.kelvinConfigWildcard)
        if dlg.ShowModal() == wx.ID_OK:
            fileName = dlg.GetPath()

            # if the new filename doesn't have the right extension, add it
            if not os.path.splitext(fileName)[1]:
                fileName = fileName + '.conf'

            # get the currently selected page and save it
            currentPanel = self.book.GetCurrentPage()
            try:
                currentPanel.saveAs(fileName)
                self.last_dir = os.path.abspath(os.path.dirname(dlg.GetPath()))
            except IOError, e:
                s = "Error saving file '%s' : %s" % \
                                    (currentPanel.fileName, e.strerror)
                dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()

            # set title of tab to new filename
            baseFileName = os.path.basename(fileName)
            currentIndex = self.book.GetPageIndex(self.book.GetCurrentPage())
            self.book.SetPageText(currentIndex, baseFileName)
        dlg.Destroy()

    def onClose(self, event):
        """Event when the close menu item is chosen"""

        # close the tab currently in focus
        self.book.DeletePage(self.book.GetSelection())

    def onCloseTab(self, event):
        """Event when a tab is closed"""

        # if the file has been modified, ask the user if want to save

        # get the page being deleted
        page = self.book.GetCurrentPage()

        if type(page) == KelvinConfigPanel:
            if page.fm.isModified:
                # the file has been modified, ask if want to save or not
                name = os.path.basename(page.fm.fileName)
                title = "The file %s has been modified" % name
                s="Do you want to save changes to this document before closing?"
                dlg = wx.MessageDialog(self, s, title, wx.YES_NO | wx.CANCEL |
                                                            wx.ICON_ERROR)
                button = dlg.ShowModal()
                dlg.Destroy()

                if button == wx.ID_CANCEL:
                    event.Veto()
                    return
                elif button == wx.ID_YES:
                    self.onSave(None)
                elif button == wx.ID_NO:
                    pass

    def onCloseWindow(self, event):
        # individually close each open tab to see if any files need saving
        n = self.book.GetPageCount()
        for i in range(n):
            self.book.DeletePage(self.book.GetSelection())

        # if any tabs are left keep the window open
        if self.book.GetPageCount() > 0:
            event.Veto()
        else:
            self.Destroy()

    def onRunKelvin(self, event):
        """Event when Kelvin runs"""

        # get the currently selected page
        page = self.book.GetCurrentPage()
        if type(page) != KelvinConfigPanel:
            s = "Kelvin can only run when a kelvin config file is selected."
            dlg = wx.MessageDialog(None, s, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return
        if page.fileName is None:
            s = "Please save the file before running Kelvin"
            dlg = wx.MessageDialog(None, s, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return

        # check for any errors
        a = page.fileChecker.checkFileForErrors()
        b = page.fileChecker.checkFileForWarnings()
        if a or b:
            s = "Errors found, Kelvin run aborted!"
            print s
            dlg = wx.MessageDialog(None, s, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return


        program = KelvinInfo.kelvinRunCommand
        args = ' ' + os.path.basename(page.fileName)
        dir = os.path.dirname(os.path.abspath(page.fileName))
        fileName = os.path.basename(page.fileName)
        p = subprocess.Popen(program + args, shell=True, cwd=dir)

        # show a progress bar
        max = 100
        dlg = wx.ProgressDialog("Running kelvin...",
                                "Running file: " + fileName,
                                maximum = max,
                                parent=page,
                                style = wx.PD_CAN_ABORT
                                | wx.PD_ELAPSED_TIME
                                | wx.PD_SMOOTH)

        keepGoing = True
        pressedCancel = False
        count = 1

        while keepGoing:
            if count % 1000 == 0:
                (keepGoing, skip) = dlg.Pulse()
            if keepGoing is False:
                pressedCancel = True
            if p.poll() is not None:
                keepGoing = False
            count += 1

        dlg.Destroy()

        if pressedCancel is True:
            # cancel was presssed, need to manually terminate process
            subprocess.call(['kill', '-9', str(p.pid)])

    def onPedCheck_01(self, event):
        """Run pedcheck levels 0 an 1 on current kelvin config file"""

        def getCommand(datafile, pedfile):
            """Return the command needed to execute pedcheck level 0 & 1"""
            return 'pedcheck -d %s -p %s' % (datafile, pedfile)

        self.runPedCheck(getCommand)

    def onPedCheck_2(self, event):
        """Run pedcheck level 2 on current kelvin config file"""

        def getCommand(datafile, pedfile):
            """Return the command needed to execute pedcheck level 2"""
            return 'pedcheck -d %s -p %s -2' % (datafile, pedfile)

        self.runPedCheck(getCommand)

    def runPedCheck(self, getCommand):
        """Run pedcheck using the given function to get the right command"""

        def getFileName(tag):
            """Gets the filename given a tag of the filename type"""
            lines = page.fm.getLines(tag)
            if len(lines) > 0:
                # the file was specified
                return lines[0][1].str
            else:
                # the file was not specified, use the default
                for eachItem in KelvinInfo.fileNameInfo:
                    if eachItem.tag == tag:
                        return eachItem.default

        page = self.book.GetCurrentPage()
        if type(page) != KelvinConfigPanel:
            # can't run pedcheck on this, show error message
            s = "Pedcheck can only run when a kelvin config file is selected."
            dlg = wx.MessageDialog(None, s, "Error", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return

        datafile = getFileName('DF')
        pedfile = getFileName('PD')
        dir = os.path.dirname(os.path.abspath(page.fileName))
        p = subprocess.Popen(getCommand(datafile, pedfile), shell=True, cwd=dir)

    def openKelvinConfigFile(self, fileName='', isNew=False):
        """Open a Kelvin config file and add a tab for that file"""

        # the third argument, isNew, specifies whether the filename is
        # one of the templates used when starting a new config file from
        # scratch, or the filename is actually a file on the disk

        try:
            panel = KelvinConfigPanel(self.book, self.kelvinConfigDefaults, 
                                      fileName, isNew)
        except IOError, e:
            str = "Error opening file '%s' : %s" % (fileName, e.strerror)
            dlg = wx.MessageDialog(self, str, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        except ParseError, e:
            dlg = wx.MessageDialog(self, e.message, "Error!", 
                                   wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            if isNew:
                # this means it was a new file
                self.book.AddPage(panel, 'Untitled', True)
            else:
                baseFileName = os.path.basename(fileName)
                self.book.AddPage(panel, baseFileName, True)

    def openSimpleFile(self, fileName):
        """Given a filename, open a textbox for it"""
        try:
            panel = SimpleTextFilePanel(self.book, fileName)
        except IOError, e:
            s = "Error opening file '%s' : %s" % (fileName, e.strerror)
            dlg = wx.MessageDialog(self, s, "Error!", 
                                   wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            baseFileName = os.path.basename(fileName)
            self.book.AddPage(panel, baseFileName, True)
