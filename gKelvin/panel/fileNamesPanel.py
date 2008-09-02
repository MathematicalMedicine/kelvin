"""Panel to display and alter the filenames specified"""

import os
import wx
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager
from config import KelvinInfo
from config import misc
from utility import util

class FileNamesPanel(wx.Panel):
    """The panel that displays and alters all filenames specified"""

    lbl1 = 'File Type'
    lbl2 = 'File Name'
    borderAmt = 5
    fileNameBorder = 0

    def __init__(self, parent, fileManager):
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Create widgets for the panel"""

        self.label1 = wx.StaticText(self, -1, self.lbl1)
        self.label2 = wx.StaticText(self, -1, self.lbl2)

        # create a filename panel for each filetype
        self.fileNamePanels = []
        for eachEntry in KelvinInfo.fileNameEntry:
            fileNamePanel = FileNamePanel(self, self.fm, eachEntry)
            self.fileNamePanels.append(fileNamePanel)

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # sizer to go along the row for the headings
        sizer = wx.BoxSizer(wx.HORIZONTAL)

        #add headings
        sizer.Add(self.label1, 0, wx.ALL, self.borderAmt)
        sizer.Add(self.label2, 0, wx.ALL, self.borderAmt)

        mainSizer.Add(sizer)

        # add filename panels
        for eachPanel in self.fileNamePanels:
            mainSizer.Add(eachPanel, 1, wx.ALL | wx.EXPAND,self.fileNameBorder)

        self.SetSizer(mainSizer)
        self.alignTextBoxes()
        self.alignHeadings()

    def alignTextBoxes(self):
        """Align the text boxes in the panel"""

        # in order for the text boxes to line up flush on the left side, 
        # find the maximum size checkbox, then set all checkboxes to that size
        maxSize = [0, 0]
        for eachPanel in self.fileNamePanels:
            size = eachPanel.checkbox.GetSizeTuple()
            maxSize[0] = max(maxSize[0], size[0])
            maxSize[1] = max(maxSize[1], size[1])

        for eachPanel in self.fileNamePanels:
            eachPanel.checkbox.SetMinSize(maxSize)

        self.Layout()

    def alignHeadings(self):
        """Make sure the headings are in the right position"""

        # first get the size of a checkbox, then make the first heading the
        # same size as the checkbox
        size = self.fileNamePanels[0].checkbox.GetSizeTuple()
        size = list(size)

        # need to account for differences in borders
        size[0] -= (self.borderAmt - self.fileNamePanels[0].borderAmt)*2
        self.label1.SetMinSize(size)
        self.Layout()

    def initialize(self):
        """Intialize widgets based on file contents"""
        # nothing to do
        pass

    def checkFileExistence(self):
        """Checks if all filenames specified actually exist"""

        for eachPanel in self.fileNamePanels:
            eachPanel.checkFileExistence()

class FileNamePanel(wx.Panel):
    """The panel that shows a single file name"""

    borderAmt = 3

    def __init__(self, parent, fileManager, tagEntry): 
        """tagEntry is a list detailing the tag, description, default value, and
        whether this file is required or not
        
        """
        
        wx.Panel.__init__(self, parent, -1)
        self.fm = fileManager
        self.tag = tagEntry[0]
        self.mainLbl = tagEntry[1]
        self.default = tagEntry[2]
        self.isRequired = tagEntry[3]

        self.createPanel()
        self.arrangePanel()
        self.initialize()

    def createPanel(self):
        """Creates the widgets for this panel"""

        # a checkbox to see if want to specify filename or use default
        self.checkbox = wx.CheckBox(self, -1, self.mainLbl)

        # a textbox to input filenames
        self.textbox = wx.TextCtrl(self, -1)

        # browse and open button
        self.browseButton = wx.Button(self, -1, label='Browse')
        self.openButton = wx.Button(self, -1, label='Open')

    def arrangePanel(self):
        """Arrange widgets in the panel"""

        mainSizer = wx.BoxSizer(wx.HORIZONTAL)

        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL 
        mainSizer.Add(self.checkbox, 0, style, self.borderAmt)
        mainSizer.Add(self.textbox, 1, style | wx.EXPAND, self.borderAmt)
        mainSizer.Add(self.browseButton, 0, style, self.borderAmt)
        mainSizer.Add(self.openButton, 0, style, self.borderAmt)

        self.SetSizer(mainSizer)

    def initialize(self):
        """Intialize widgets based on file contents"""

        # figure out responsibilities
        # see if line with tag exists
        if self.fm.isPresent(self.tag):
            self.line = self.fm.getLines(self.tag)[0]
            self.checkbox.SetValue(True)
        else:
            # create a default line just in case it's needed
            self.line = self.fm.createLine(self.tag, self.default)
            self.checkbox.SetValue(False)
            self.textbox.Disable()

        self.textbox.charge = self.line[1]
        self.textbox.ChangeValue(self.textbox.charge.str)

        # bind events
        self.browseButton.Bind(wx.EVT_BUTTON, self.onBrowse, self.browseButton)
        self.openButton.Bind(wx.EVT_BUTTON, self.onOpen, self.openButton)
        self.checkbox.Bind(wx.EVT_CHECKBOX, self.onCheckBox, self.checkbox)
        self.textbox.Bind(wx.EVT_TEXT, self.onTextEntry, self.textbox)

        # check file existance if necessary
        self.checkFileExistence()

    def onCheckBox(self, event):
        """Event when the checkbox is clicked"""

        if self.checkbox.GetValue() is True:
            # user wants to specify file
            self.fm.addLine(self.line)
            self.textbox.Enable()
        else:
            # user wants to use the default
            self.fm.remove(self.line)
            self.textbox.Disable()
            self.textbox.ChangeValue(self.default)
            self.onTextEntry(None)
        util.callAncestorFunction(self, 'hasChanged')

    def onBrowse(self, event):
        """Event when the browse button is pressed"""
        if self.textbox.IsEnabled():
            # start the dialog box in the same folder as the kelvin config file 
            startFolder = os.path.dirname(os.path.abspath(self.fm.fileName))
            dlg = wx.FileDialog(self, "Choose File", startFolder,style=wx.OPEN)
            if dlg.ShowModal() == wx.ID_OK:
                fileName = dlg.GetPath()

                # only show the relative path to the file in the textbox
                relativePath = self.relpath(self.fm.fileName, fileName)
                self.textbox.ChangeValue(relativePath)
                self.onTextEntry(None)
            dlg.Destroy()

    def onOpen(self, event):
        """Event when the open button is pressed"""

        # get the filename, and get an ancestor to take care of it
        fileName = self.textbox.GetValue()
        if not os.path.isabs(fileName):
            # fileName is a relative address, make it an absolute one by joining
            path = os.path.dirname(self.fm.fileName)
            fileName = os.path.join(path, fileName)

        # it with the path of the file being opened
        util.callAncestorFunction(self, 'openSimpleFile', fileName)

    def onTextEntry(self, event):
        """Event when text is entered in the textbox"""

        self.textbox.charge.update(self.textbox.GetValue())
        self.checkFileExistence()
        util.callAncestorFunction(self, 'hasChanged')

    def checkFileExistence(self):
        """If the file is required, check if it exists or not"""

        if not self.isRequired:
            return

        # don't check if there is no filename (it's a new file)
        if not self.fm.fileName:
            return

        # make sure the path is relative to the filename
        path = os.path.normpath(os.path.abspath(self.fm.fileName))
        path = os.path.dirname(path)
        path = os.path.join(path, self.textbox.GetValue())

        if not os.path.exists(path):
            self.setError(True)
        else:
            self.setError(False)

    def setError(self, bool):
        """Put this panel in an error state or out of an error state"""

        # change the background color if in error
        if bool:
            self.SetBackgroundColour(misc.errorColor)
        else:
            self.SetBackgroundColour(wx.NullColour)

        self.Refresh()

    def relpath(self, file1, file2):
        """Given two full file paths (with filenames included), return a path
        that is equivalent to path2 but is relative to path1
        
        """

        # first make sure both of them are absolute and normalized
        path1 = os.path.normpath(os.path.abspath(file1))
        path2 = os.path.normpath(os.path.abspath(file2))

        # get the directory path
        dirpath1 = os.path.dirname(path1) + os.sep
        dirpath2 = os.path.dirname(path2) + os.sep

        # get the common prefix between the two paths
        common = os.path.commonprefix([dirpath1, dirpath2])

        if common == os.sep:
            # the only common path is the very top level
            # totally different paths, nothing to do
            return file2

        # strip out the common prefix of the paths
        spath1 = path1.replace(common, '', 1)
        spath2 = path2.replace(common, '', 1)

        # split the paths into different folders to see 
        # the number of folders in the path
        split1 = spath1.split(os.sep)
        split2 = spath2.split(os.sep)

        for i in range(0, len(split1)-1):
            # need to go one level up of path2
            spath2 = os.path.pardir + os.path.sep + spath2

        return spath2

