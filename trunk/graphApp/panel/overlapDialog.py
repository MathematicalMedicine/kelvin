"""Dialog boxes used in the overlap menu"""

import wx
from config import misc

class OverlapExistingChrDialog(wx.Dialog):
    """A dialog used when the user wants to overlap existing lines with
    another existing line in the graph
    
    """

    # Note: source chromosomes mean the chromosomes that get moved
    # and destination chromosome refers to the chromosome whose location
    # the source ones will move to

    title = 'Please Specify Chromosomes to Overlap'
    sourceLbl = 'Chromosomes that will be moved'
    destLbl = 'Chromosome to overlap with'

    size = misc.fileformatDialogSize
    checkBoxListSize = wx.DefaultSize

    def __init__(self, parent, id, info, chrNames):

        self.file = file
        self.info = info

        # all the names of the other chromosomes in the graph
        self.chrNames = chrNames

        pre = wx.PreDialog()
        style = wx.DEFAULT_DIALOG_STYLE
        pre.Create(parent, id, self.title, wx.DefaultPosition, self.size, style)
        self.PostCreate(pre)

        self.createWidgets()
        self.arrange()
        self.intialize()
        self.CenterOnParent()

    def createWidgets(self):
        """Creates widgets for the dialog box"""

        # the ok and cancel buttons
        self.okButton = wx.Button(self, wx.ID_OK)
        self.okButton.SetDefault()
        self.cancelButton = wx.Button(self, wx.ID_CANCEL)

        # some text labels
        self.sourceLabel = wx.StaticText(self, -1, self.sourceLbl)
        self.destLabel = wx.StaticText(self, -1, self.destLbl)

        # a checklistbox to determine the source chromosome
        self.sourceListBox = wx.CheckListBox(self, -1, (80,50), 
                                           self.checkBoxListSize, self.chrNames)

        # choices boxes to pick the destination chromosome
        self.destChoice = wx.Choice(self, -1, choices=self.chrNames)

    def arrange(self):
        """Arranges widgets"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # everything is added to this sizer because padding is needed
        # around the whole dialog box, so this sizer is added to the
        # main sizer with some padding
        otherSizer = wx.BoxSizer(wx.VERTICAL)

        style = wx.ALL 
        otherSizer.Add(self.sourceLabel, 0, style, 5)
        otherSizer.Add(self.sourceListBox, 0, style, 5)
        otherSizer.Add((1, 20))
        otherSizer.Add(self.destLabel, 0, style, 5)
        otherSizer.Add(self.destChoice, 0, style, 5)

        # some empty space
        otherSizer.Add((1, 20))

        # add the ok and cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        btnsizer.AddButton(self.okButton)
        btnsizer.AddButton(self.cancelButton)
        btnsizer.Realize()

        otherSizer.Add(btnsizer, 0, wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_RIGHT | 
                      wx.BOTTOM | wx.TOP, 5)
        mainSizer.Add(otherSizer, 0, wx.ALL, 10)
        self.SetSizer(mainSizer)
        self.Fit()

    def intialize(self):
        """Give widgets a default value"""

        self.info[0] = []
        self.info[1] = 0
        self.destChoice.SetSelection(0)

        # bind events
        self.destChoice.Bind(wx.EVT_CHOICE, self.onDestChoice, self.destChoice)
        self.sourceListBox.Bind(wx.EVT_CHECKLISTBOX, self.onSourceCheck, 
                                self.sourceListBox)
        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)

    def onDestChoice(self, event):
        """Event when the destination choice box is changed"""
        self.info[1] = self.destChoice.GetSelection()

    def onSourceCheck(self, event):
        """Event when the user checks or unchecks the checklistbox"""

        index = event.GetSelection()
        if self.sourceListBox.IsChecked(index):
            # add the selection to the list
            self.info[0].append(index)
        else:
            # must have unchecked, remove from list
            self.info[0].remove(index)

    def onOK(self, event):
        """Event when the ok button is pressed"""

        # make sure that one of the source chromosomes isn't the same 
        # as the destination chromosome
        if self.info[1] in self.info[0]:
            s1 = "One of the chromosomes to be moved is the same as "
            s2 = 'the chromosome to overlap with'
            message = s1 + s2
            dlg = wx.MessageDialog(self, message,"Error!",wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            # everything's ok, unshow this dialog box
            event.Skip()
