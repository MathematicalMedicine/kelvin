"""Used to generate events"""

import wx

def createEvent(action, widget, integer=None):
    """Given the event action type and widget, generate such an event"""


    evt = wx.CommandEvent(action.evtType[0], widget.GetId())
    evt.SetEventObject(widget)
    if integer is not None:
        evt.SetInt(integer)
    widget.Command(evt)
