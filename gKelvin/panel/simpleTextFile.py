import wx
from config import misc

class SimpleTextFilePanel(wx.Panel):
    """A panel that takes a filename and shows the contents of that file."""

    borderAmt = 20

    def __init__(self, parent, fileName):
        wx.Panel.__init__(self, parent, -1) 
        self.fileName = fileName

        # open the file
        self.file = open(self.fileName, 'r')

        if self.file == None:
            return

        sizer = wx.BoxSizer(wx.VERTICAL)

        # make a text control
        self.textbox = wx.TextCtrl(self, -1, self.file.read(), 
                                   style=wx.TE_MULTILINE)

        sizer.Add(self.textbox, 1, wx.EXPAND | wx.ALL, self.borderAmt)
        self.SetSizer(sizer)
        #self.Fit()

        self.file.close()

    def openFile(self):
            return file

    def save(self, fileName=''):
        if misc.savingEnabled:
            if fileName == '':
                fileName = self.fileName

            self.file = open(self.fileName, 'w')
            self.file.write(self.textbox.GetValue())
            self.file.close()

    def saveAs(self, fileName):
        self.fileName = fileName
        self.save()
