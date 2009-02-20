"""Dialog box asking user about how to annotate above or below a threshold"""

import wx

from config import misc

class AnnotateThresholdDialog(wx.Dialog):
    """A dialog box asking the user what value the threshold is and whether
    to annotate peaks above the threshold or troughs below the threshold

    """

    TITLE = 'Please Specify Threshold Parameters'
    thresholdLbl = 'Threshold Value:'
    annotatePeaksLbl = 'Annotate peaks above threshold'
    annotateTroughsLbl = 'Annotate valleys below threshold'

    def __init__(self, parent, ID, info, size=misc.fileformatDialogSize, 
                 pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE):

        self.info = info

        pre = wx.PreDialog()
        pre.Create(parent, ID, self.TITLE, pos, size, style)
        self.PostCreate(pre)

        self.createWidgets()
        self.arrange()
        self.intialize()

    def createWidgets(self):
        """Creates widgets for the dialog box"""

        # the ok and cancel buttons
        self.okButton = wx.Button(self, wx.ID_OK)
        self.okButton.SetDefault()
        self.cancelButton = wx.Button(self, wx.ID_CANCEL)

        # label for text box
        self.thresholdLabel = wx.StaticText(self, -1, self.thresholdLbl)

        # this text box is where to input the threshold value
        self.thresholdBox = wx.TextCtrl(self, -1)

        # radio buttons to choose whether the to annotate the peaks or troughs
        self.peakRadio = wx.RadioButton(self, -1, self.annotatePeaksLbl)
        self.troughRadio = wx.RadioButton(self, -1, self.annotateTroughsLbl)

    def arrange(self):
        """Arranges widgets"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        # everything is added to this sizer because padding is needed
        # around the whole dialog box, so this sizer is added to the
        # main sizer with some padding
        otherSizer = wx.BoxSizer(wx.VERTICAL)

        # the first row has the threshold text and box
        rowSizer = wx.BoxSizer(wx.HORIZONTAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        rowSizer.Add(self.thresholdLabel, 0, style, 10)
        rowSizer.Add(self.thresholdBox, 0, style, 10)
        otherSizer.Add(rowSizer)

        # some empty space
        otherSizer.Add((1, 20))

        radioSizer = wx.BoxSizer(wx.VERTICAL)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALL
        radioSizer.Add(self.peakRadio, 0, style, 5)
        radioSizer.Add(self.troughRadio, 0, style, 5)
        otherSizer.Add(radioSizer)

        # some more empty space
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

    def intialize(self):
        """Give widgets a default value"""
        
        # default values for info
        self.info['threshold'] = 0.0
        self.info['annotate_peaks'] = True

        self.thresholdBox.SetValue(str(self.info['threshold']))
        self.peakRadio.SetValue(self.info['annotate_peaks'])
        self.troughRadio.SetValue(not self.info['annotate_peaks'])

        # bind events
        self.thresholdBox.Bind(wx.EVT_TEXT, self.onThresholdText, 
                                                        self.thresholdBox)
        self.peakRadio.Bind(wx.EVT_RADIOBUTTON, self.onRadioBtn, self.peakRadio)
        self.troughRadio.Bind(wx.EVT_RADIOBUTTON, self.onRadioBtn, 
                                                        self.troughRadio)

    def onThresholdText(self, event):
        """Event when text is entered in the threshold box"""

        try:
            self.info['threshold'] = float(self.thresholdBox.GetValue())
        except ValueError, e:
            self.info['threshold'] = 0.0

    def onRadioBtn(self, event):
        """Event when a radio button is clicked"""

        if self.peakRadio.GetValue():
            self.info['annotate_peaks'] = True
        else:
            self.info['annotate_peaks'] = False
