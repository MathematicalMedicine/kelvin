"""Dialog box asking user what to set the view limits to"""

import wx

class ViewLimitsDialog(wx.Dialog):
    """A dialog box asking the user what to set the view limits to.
    y-axis limits are handled by just a number, but x-axis values
    are specified by a chromosome and the position value of that chromosome.
    
    """

    TITLE = 'Please Specify x-axis and y-axis limits'
    yminLbl = 'ymin:'
    ymaxLbl = 'ymax:'
    xminLbl = 'xmin:'
    xmaxLbl = 'xmax:'
    chrLbl = '   chr:'
    rangeLbl = '   range:'

    def __init__(self, parent, info, size=wx.DefaultSize,
                 pos=wx.DefaultPosition, style=wx.DEFAULT_DIALOG_STYLE):

        self.info = info
        self.mainFrame = parent

        pre = wx.PreDialog()
        pre.Create(parent, -1, self.TITLE, pos, size, style)
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

        # labels for the dialog box
        self.xminLabel = wx.StaticText(self, -1, self.xminLbl)
        self.xmaxLabel = wx.StaticText(self, -1, self.xmaxLbl)
        self.yminLabel = wx.StaticText(self, -1, self.yminLbl)
        self.ymaxLabel = wx.StaticText(self, -1, self.ymaxLbl)
        self.chrLabel1 = wx.StaticText(self, -1, self.chrLbl)
        self.chrLabel2 = wx.StaticText(self, -1, self.chrLbl)
        self.rangeLabel1 = wx.StaticText(self, -1, self.rangeLbl)
        self.rangeLabel2 = wx.StaticText(self, -1, self.rangeLbl)

        # text boxes for the actual values
        self.xminBox = wx.TextCtrl(self, -1)
        self.xmaxBox = wx.TextCtrl(self, -1)
        self.yminBox = wx.TextCtrl(self, -1)
        self.ymaxBox = wx.TextCtrl(self, -1)

        # choice boxes to select what chromosome to use for the x-values
        choices = self.mainFrame.lines.getNames()
        self.xminChrChoice = wx.Choice(self, -1, choices=choices)
        self.xmaxChrChoice = wx.Choice(self, -1, choices=choices)

        # two labels that show the possible range for each chromosome
        self.xminRangeLabel = wx.StaticText(self, -1)
        self.xmaxRangeLabel = wx.StaticText(self, -1)

    def arrange(self):
        """Arranges widgets"""

        mainSizer = wx.BoxSizer(wx.VERTICAL)

        sizer = wx.FlexGridSizer(4, 6, 5, 5)
        style = wx.ALIGN_CENTER_VERTICAL | wx.ALIGN_CENTER_HORIZONTAL

        sizer.Add(self.xminLabel, 0, style)
        sizer.Add(self.xminBox, 0, style)
        sizer.Add(self.chrLabel1, 0, style)
        sizer.Add(self.xminChrChoice, 0, style)
        sizer.Add(self.rangeLabel1, 0, style)
        sizer.Add(self.xminRangeLabel, 0, style)

        sizer.Add(self.xmaxLabel, 0, style)
        sizer.Add(self.xmaxBox, 0, style)
        sizer.Add(self.chrLabel2, 0, style)
        sizer.Add(self.xmaxChrChoice, 0, style)
        sizer.Add(self.rangeLabel2, 0, style)
        sizer.Add(self.xmaxRangeLabel, 0, style)

        sizer.Add(self.yminLabel, 0, style)
        sizer.Add(self.yminBox, 0, style)
        sizer.Add((1, 1))
        sizer.Add((1, 1))
        sizer.Add((1, 1))
        sizer.Add((1, 1))

        sizer.Add(self.ymaxLabel, 0, style)
        sizer.Add(self.ymaxBox, 0, style)
        sizer.Add((1, 1))
        sizer.Add((1, 1))
        sizer.Add((1, 1))
        sizer.Add((1, 1))

        mainSizer.Add(sizer, 0, wx.ALL, 15)

        # add the ok and cancel buttons
        btnsizer = wx.StdDialogButtonSizer()
        btnsizer.AddButton(self.okButton)
        btnsizer.AddButton(self.cancelButton)
        btnsizer.Realize()
        mainSizer.Add(btnsizer, 0, wx.ALIGN_RIGHT | wx.ALL, 10)

        self.SetSizer(mainSizer)
        self.Fit()

    def intialize(self):
        """Give widgets a default value"""

        # Note: the 'index' attribute that is being assigned indicates
        # the index of info for which that widget is responsible for

        # the y-values are easy to initialize
        self.yminBox.SetValue(str(round(self.info[2], 3)))
        self.yminBox.index = 2
        self.ymaxBox.SetValue(str(round(self.info[3], 3)))
        self.ymaxBox.index = 3

        # the x-values need to be converted to the right chromosome, etc.
        self.init_x(self.xminBox, self.xminChrChoice, self.xminRangeLabel, 0)
        self.init_x(self.xmaxBox, self.xmaxChrChoice, self.xmaxRangeLabel, 1)

        # the choice boxes need to know the other wigets in its row
        self.xminChrChoice.textbox = self.xminBox
        self.xminChrChoice.rangeLabel = self.xminRangeLabel
        self.xmaxChrChoice.textbox = self.xmaxBox
        self.xmaxChrChoice.rangeLabel = self.xmaxRangeLabel

        # this is because the range labels will change size
        self.Fit()

        # bind events
        self.yminBox.Bind(wx.EVT_TEXT, self.onText, self.yminBox)
        self.ymaxBox.Bind(wx.EVT_TEXT, self.onText, self.ymaxBox)
        self.xminBox.Bind(wx.EVT_TEXT, self.onText, self.xminBox)
        self.xmaxBox.Bind(wx.EVT_TEXT, self.onText, self.xmaxBox)
        self.xminChrChoice.Bind(wx.EVT_CHOICE, self.onChoice,self.xminChrChoice)
        self.xmaxChrChoice.Bind(wx.EVT_CHOICE, self.onChoice,self.xmaxChrChoice)
        self.okButton.Bind(wx.EVT_BUTTON, self.onOK, self.okButton)

    def init_x(self, textbox, choice, rangeLabel, index):
        """Initializes the widgets for the x-axis values of one row
        (either xmin or xmax).
        
        """

        # find the group that overlaps the x-value
        x = self.info[index]
        group = self.mainFrame.lines.getGroupByX(x)

        if group is None:
            # this means the point isn't in any group. either the point is
            # before the first group or after the last group. if before,
            # use the first group, if after, use the last group.
            groups = self.mainFrame.lines.getGroups()
            if x <= groups[0].xmin:
                group = groups[0]
            else:
                group = groups[-1]

        # use the first line in the group for initialization
        line = group.lines[0]

        # set widgets to correct values
        val = x - group.xmin + group.start
        val = max(val, group.start)
        val = min(val, group.end)
        textbox.SetValue(str(round(val, 4)))
        choice.SetStringSelection(line.name)
        rangeLabel.SetLabel(str(line.pos[0]) + '-' + str(line.pos[-1]))

        # keep track of some values
        textbox.index = index
        textbox.offset = group.xmin

    def onText(self, event):
        """Event when the value in a textbox changes"""

        textbox = event.GetEventObject()
        try:
            value = float(textbox.GetValue())
        except ValueError, e:
            value = None
        self.info[textbox.index] = value

    def onChoice(self, event):
        """Event when a choice box is selected"""

        choice = event.GetEventObject()
        index = choice.GetSelection()
        line = self.mainFrame.lines.getLine(index)

        # update the range label
        choice.rangeLabel.SetLabel(str(line.pos[0]) + '-' + str(line.pos[-1]))

        # update the textbox if its current value is outside the line interval
        x = float(choice.textbox.GetValue())
        if x < line.pos[0]:
            choice.textbox.SetValue(str(round(line.pos[0], 4)))
        if x > line.pos[-1]:
            choice.textbox.SetValue(str(round(line.pos[-1], 4)))

    def onOK(self, event):
        """Event when user clicks the ok button"""

        # do lots of error checking
        hasError = False

        # make sure all info in text boxes are numbers
        hasError = hasError or self.validateTextBoxes()

        # make sure ymin isn't greater than ymax
        try:
            ymin = float(self.yminBox.GetValue())
            ymax = float(self.ymaxBox.GetValue())
        except ValueError, e:
            pass
        else:
            if ymin > ymax:
                self.showError('ymin is greater than ymax')
                hasError = True

        # make sure the x-values aren't outside the range of the 
        # selected chromosome
        error1 = self.check_x_range(self.xminBox, self.xminChrChoice, 'xmin')
        error2 = self.check_x_range(self.xmaxBox, self.xmaxChrChoice, 'xmax')
        hasError = hasError or error1 or error2
        
        # convert the xmin and xmax values in the textboxes that are currently 
        # relative to the chromosome to values relative to the whole graph
        index = self.xminChrChoice.GetSelection()
        minLine, minGroup = self.mainFrame.lines.getLineAndGroup(index)
        index = self.xmaxChrChoice.GetSelection()
        maxLine, maxGroup= self.mainFrame.lines.getLineAndGroup(index)
        try:
            xmin = minGroup.xmin+float(self.xminBox.GetValue())-minGroup.start
            xmax = maxGroup.xmin+float(self.xmaxBox.GetValue())-maxGroup.start
        except ValueError, e:
            pass
        else:
            # make sure xmin is less than xmax
            if xmin > xmax:
                self.showError('xmin is greater than xmax')
                hasError = True
            else:
                self.info[0] = xmin
                self.info[1] = xmax

        if hasError:
            # keep the dialog box alive
            return
        else:
            event.Skip()

    def validateTextBoxes(self):
        """Checks all text boxes to see if their value is a number"""

        hasError = False
        if self.info[0] == None:
            s = "The xmin value '%s' cannot be converted to a number" % \
                                                        self.xminBox.GetValue()
            self.showError(s)
            hasError = True
        if self.info[1] == None:
            s = "The xmax value '%s' cannot be converted to a number" % \
                                                        self.xmaxBox.GetValue()
            self.showError(s)
            hasError = True
        if self.info[2] == None:
            s = "The ymin value '%s' cannot be converted to a number" % \
                                                        self.yminBox.GetValue()
            self.showError(s)
            hasError = True
        if self.info[3] == None:
            s = "The ymax value '%s' cannot be converted to a number" % \
                                                        self.ymaxBox.GetValue()
            self.showError(s)
            hasError = True

        return hasError

    def check_x_range(self, textbox, choice, label):
        """Check to see if the value in the text box are within the range
        of the selected chromosome in the choice box.
        
        """

        index = choice.GetSelection()
        line = self.mainFrame.lines.getLine(index)

        try:
            x = float(textbox.GetValue())
        except:
            return False
        else:
            if x < line.pos[0] or x > line.pos[-1]:
                s="%s value of %f is outside the chromosome range of %f-%f" % \
                                           (label, x, line.pos[0], line.pos[-1])
                self.showError(s)
                return True
            else:
                return False

    def showError(self, message):
            """Shows a warning message box"""
                
            dlg = wx.MessageDialog(None, message, "Error!", 
                                   wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
