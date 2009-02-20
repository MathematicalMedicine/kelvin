"""Utility functions to help create toolbars"""


import wx


def addToolbarLabel(parent, toolbar, label, bitmap, shortHelp,longHelp,handler):
    """Add one label to the toolbar"""

    if label == '':
        toolbar.AddSeparator()
        return

    id = wx.NewId()
    toolbar.AddLabelTool(id,label,bitmap,shortHelp=shortHelp,longHelp=longHelp)

    if handler is not None:
        parent.Bind(wx.EVT_TOOL, handler, id=id)

def addToolbarLabels(parent, toolbar, toolbarData):
    """Add labels to a toolbar"""

    for eachEntry in toolbarData:
        addToolbarLabel(parent, toolbar, *eachEntry)

