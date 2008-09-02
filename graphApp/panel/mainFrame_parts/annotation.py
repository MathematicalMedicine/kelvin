"""Contains functions that deal with annotations and their interactions with the mainFrame"""

import wx

from panel.annotateThresholdDialog import AnnotateThresholdDialog
from transforms import pixelToAxes
from widgets.annotations import my_mpl_annotation

def moveAnnotation(mainFrame, event):
    """Update the position of the currently selected annotation"""

    if mainFrame.navTb.myMode == 'ANNOTATE_MOVE' and \
             mainFrame.move_annotation is not None:
        # move the selected annotation to the mouse position
        x, y = pixelToAxes(mainFrame.axes, event.x, event.y)
        mainFrame.move_annotation.set_position((x, y))
        mainFrame.move_annotation.setArrow()
        mainFrame.draw()
        mainFrame.move_annotation.autoLayout = False

def annotatePoint(mainFrame, x, y):
    """Create an annotation, with x and y the point being annotated 
    in data space.
    
    """

    mainFrame.annotationManager.annotate(x, y)
    mainFrame.annotationManager.layout()
    mainFrame.draw()

def annotateThreshold(mainFrame):
    """Annotate all peaks above or troughs below a certain threshold"""

    # ask user about what to annotate,
    # the variable info is a dictionary of parameters 
    # that are decided from the dialog box
    info = dict(threshold=0.0, annotate_peaks=True)
    dlg = AnnotateThresholdDialog(mainFrame, -1, info)
    dlg.CenterOnParent()
    if dlg.ShowModal() == wx.ID_OK:
        val = info['threshold']
        annotate_peaks = info['annotate_peaks']
        mainFrame.annotationManager.annotateThreshold(val, annotate_peaks)
        mainFrame.draw()
    dlg.Destroy()

def annotateMax(mainFrame, x0, y0, x1, y1):
    """Given the coordinates of a rectangle in data space, annotate
    the highest peak in that rectangle
    
    """

    mainFrame.annotationManager.annotateMax(x0, y0, x1, y1)

def annotateMin(mainFrame, x0, y0, x1, y1):
    """Given the coordinates of a rectangle in data space, annotate
    the lowest valley in that rectangle
    
    """

    mainFrame.annotationManager.annotateMin(x0, y0, x1, y1)

def updateAnnotations(mainFrame):
    """Updates the text in annotations"""

    mainFrame.annotationManager.update()

    # this is needed so the information text area on the right side
    # reflects the latest annotation display
    mainFrame.updateInfoText()

def createAnnPopupMenu(mainFrame):
    """Create a popup menu dealing with annotations"""

    data = (('Edit Annotation', lambda evt: onEditAnnotation(mainFrame, evt)),
            ('Remove Annotation',lambda evt:onRemoveAnnotation(mainFrame, evt)))

    mainFrame.popupmenu = wx.Menu()
    for menuItem in data:
        item = mainFrame.popupmenu.Append(-1, menuItem[0])
        mainFrame.Bind(wx.EVT_MENU, menuItem[1], item)

def onEditAnnotation(mainFrame, event):
    """Event when a user right clicks on an annotation and chooses this"""

    text = mainFrame.last_annotation.get_text()
    dlg = wx.TextEntryDialog(None, 'Please enter text', 'Edit Annotation',
                             text, style=wx.OK|wx.CANCEL)
    dlg.CenterOnParent()
    if dlg.ShowModal() == wx.ID_OK:
        # need to let the Annotation class know the change
        if isinstance(mainFrame.last_annotation, my_mpl_annotation):
            # go up one to get to Annotation class
            mainFrame.last_annotation.parent.set_text(dlg.GetValue())
        else:
            # this is Annotation class
            mainFrame.last_annotation.set_text(dlg.GetValue())
        # need to let sidebar know about changes
        mainFrame.updateInfoText()

    dlg.Destroy()
    mainFrame.draw()

def onRemoveAnnotation(mainFrame, event):
    """Event when a user right clicks on an annotation and chooses this"""

    mainFrame.annotationManager.remove(mainFrame.last_annotation)
    mainFrame.draw()

