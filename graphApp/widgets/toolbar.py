"""Toolbars used in the graph program"""

import Image
import ImageOps
import os
import wx
from matplotlib.backend_bases import NavigationToolbar2
from matplotlib.backends.backend_wx import _load_bitmap
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.backends.backend_wxagg import DEBUG_MSG
from matplotlib.backends.backend_wxagg import NavigationToolbar2WxAgg 

from config import graphInfo
import panel.mainFrame_parts.annotation as ann
from transforms import pixelToData
from transforms import dataToPixel
import utility

def getIcon(name):
    """Given a name of a local icon file, return the bitmap for it"""
    return wx.Bitmap(graphInfo.getImagePath(name))

# i'm create subclasses of all these toolbar classes because in the
# matplotlib source code it assumes that the toolbar's parent is the same
# as the canvas's parent, which isn't the case if you have the canvas
# inside another panel, while you want the toolbar to be present in the frame,
# in which case the parent of the toolbar should be the frame. 
# otherwise, placing the toolbar will result in errors

class MyNavigationToolbar2Wx(NavigationToolbar2Wx):
    """My own version of said class..."""

    def __init__(self, parent, canvas):
        # SOURCE: taken from source code
        # don't use the canvas's parent as the toolbar's parent!
        #wx.ToolBar.__init__(self, canvas.GetParent(), -1)  # <-- old way
        wx.ToolBar.__init__(self, parent, -1)               # <-- new way

        NavigationToolbar2.__init__(self, canvas)
        self.canvas = canvas
        self._idle = True
        self.statbar = None

class MyNavigationToolbar2WxAgg(MyNavigationToolbar2Wx,NavigationToolbar2WxAgg):
    """My own version of said class..."""

    def __init__(self, parent, canvas):
        MyNavigationToolbar2Wx.__init__(self, parent, canvas)

class NavigationToolbar(MyNavigationToolbar2WxAgg):
    """
    My own version of the builtin matplotlib toolbar2, toolbar that
    lets you zoom, pan, etc.
    
    """

    def __init__(self, parent, canvas):
        MyNavigationToolbar2WxAgg.__init__(self, parent, canvas)

        self.parent = parent
        self.myMode = ''
        self.axes = self.canvas.figure.get_axes()[0]
        self.mouse_press = None

        # find ymin and ymax in pixel space
        ymin, ymax = self.axes.get_ylim()
        self.ymin = dataToPixel(self.axes, y=ymin)
        self.ymax = dataToPixel(self.axes, y=ymax)

        self.canvas.mpl_connect('button_press_event', self.onButtonPress)
        self.canvas.mpl_connect('button_release_event', self.onButtonRelease)
        self.canvas.mpl_connect('motion_notify_event', self.onMouseMotion)

    def _update_view(self):
        NavigationToolbar2WxAgg._update_view(self)

        # the parent needs to know!
        self.parent.onUpdateView(None)

    def _init_toolbar(self):
        # SOURCE: this function is taken from the source code, 
        # in the class NavigationToolbar2WxAgg, with some of my
        # own additions...

        DEBUG_MSG("_init_toolbar", 1, self)

        self._parent = self.canvas.GetParent()
        _NTB2_HOME            = wx.NewId()
        self._NTB2_BACK       = wx.NewId()
        self._NTB2_FORWARD    = wx.NewId()
        self._NTB2_PAN        = wx.NewId()
        self._NTB2_ZOOM       = wx.NewId()
        self._NTB2_ZOOM_X     = wx.NewId()
        self._NTB2_ZOOM_Y     = wx.NewId()
        self._NTB2_ANN_POINT  = wx.NewId()
        self._NTB2_ANN_THRESH = wx.NewId()
        self._NTB2_ANN_MAX    = wx.NewId()
        self._NTB2_ANN_MIN    = wx.NewId()
        self._NTB2_ANN_MOVE   = wx.NewId()
        _NTB2_SAVE            = wx.NewId()
        _NTB2_SUBPLOT         = wx.NewId()        

        self.SetToolBitmapSize(wx.Size(24,24))

        self.AddSimpleTool(_NTB2_HOME, _load_bitmap('home.xpm'),
                           'Home', 'Reset original view')
        self.AddSimpleTool(self._NTB2_BACK, _load_bitmap('back.xpm'),
                           'Back', 'Back navigation view')
        self.AddSimpleTool(self._NTB2_FORWARD, _load_bitmap('forward.xpm'),
                           'Forward', 'Forward navigation view')
        self.AddCheckTool(self._NTB2_PAN, _load_bitmap('move.xpm'),
                           shortHelp='Pan', 
                           longHelp='Pan with left, zoom with right')
        self.AddCheckTool(self._NTB2_ZOOM, _load_bitmap('zoom_to_rect.xpm'),
                           shortHelp='Zoom', longHelp='Zoom to rectangle')
        self.AddCheckTool(self._NTB2_ZOOM_X, getIcon('zoom_x.xpm'), 
                          shortHelp='Zoom x-axis',
                          longHelp='Zoom to rectangle: x-axis')
        self.AddCheckTool(self._NTB2_ZOOM_Y, getIcon('zoom_y.xpm'), 
                          shortHelp='Zoom y-axis',
                          longHelp='Zoom to rectangle: y-axis')

        self.AddSeparator()

        ann_point_bmp = getIcon('annotate_point.xpm')
        ann_thresh_bmp = getIcon('annotate_threshold.xpm')
        ann_max_bmp = getIcon('annotate_max.xpm')
        ann_min_bmp = getIcon('annotate_min.xpm')
        ann_move_bmp = getIcon('move_annotation.xpm')
        self.AddCheckTool(self._NTB2_ANN_POINT, ann_point_bmp, 
                          shortHelp='Annotate Point', longHelp='Annotate Point')
        self.AddSimpleTool(self._NTB2_ANN_THRESH, ann_thresh_bmp, 
                     'Annotate Threshold', 'Annotate All Peaks Above Threshold')
        self.AddCheckTool(self._NTB2_ANN_MAX, ann_max_bmp, 
                          shortHelp='Annotate Peak', 
                          longHelp='Annotate Peak of Area')
        self.AddCheckTool(self._NTB2_ANN_MIN, ann_min_bmp, 
                          shortHelp='Annotate Valley', 
                          longHelp='Annotate Valley of Area')
        self.AddCheckTool(self._NTB2_ANN_MOVE, ann_move_bmp, 
                          shortHelp='Move Annotation', 
                          longHelp='Move Annotation Mode')

        self.AddSeparator()

        self.AddSimpleTool(_NTB2_SUBPLOT, _load_bitmap('subplots.xpm'),
                           'Configure plot', 'Configure plot parameters')

        self.Bind(wx.EVT_TOOL, self.home, id=_NTB2_HOME)
        self.Bind(wx.EVT_TOOL, self.forward, id=self._NTB2_FORWARD)
        self.Bind(wx.EVT_TOOL, self.back, id=self._NTB2_BACK)
        self.Bind(wx.EVT_TOOL, self.zoom, id=self._NTB2_ZOOM)
        self.Bind(wx.EVT_TOOL, self.zoom_x, id=self._NTB2_ZOOM_X)
        self.Bind(wx.EVT_TOOL, self.zoom_y, id=self._NTB2_ZOOM_Y)
        self.Bind(wx.EVT_TOOL, self.pan, id=self._NTB2_PAN)
        self.Bind(wx.EVT_TOOL, self.onAnnotatePoint, id=self._NTB2_ANN_POINT)
        self.Bind(wx.EVT_TOOL, self.onAnnotateThreshold,
                                                    id=self._NTB2_ANN_THRESH)
        self.Bind(wx.EVT_TOOL, self.onAnnotateMax, id=self._NTB2_ANN_MAX)
        self.Bind(wx.EVT_TOOL, self.onAnnotateMin, id=self._NTB2_ANN_MIN)
        self.Bind(wx.EVT_TOOL, self.onAnnotateMove, id=self._NTB2_ANN_MOVE)
        self.Bind(wx.EVT_TOOL, self.configure_subplot, id=_NTB2_SUBPLOT)
        self.Bind(wx.EVT_TOOL, self.save, id=_NTB2_SAVE)

        self.Realize()

        # this is a list of all my own toggle button id's
        self.buttonIds = [self._NTB2_ZOOM, self._NTB2_PAN, self._NTB2_ANN_POINT,
                          self._NTB2_ANN_MAX, self._NTB2_ANN_MIN, 
                          self._NTB2_ANN_MOVE, self._NTB2_ZOOM_X, 
                          self._NTB2_ZOOM_Y]

    def zoom(self, *args):
        """Event when the zoom button is pressed"""

        # SOURCE: overwriting builtin-function
        self.toggleButtons(self._NTB2_ZOOM)
        if self.myMode == 'ZOOM':
            self.myMode = ''
        else:
            self.myMode = 'ZOOM'
        MyNavigationToolbar2WxAgg.zoom(self, *args)
    
    def release_zoom(self, event):
        """Event when the mouse button is released in zoom mode"""

        # SOURCE: overwriting a builtin function
        # save annotation positions before zoom happens
        MyNavigationToolbar2WxAgg.release_zoom(self, event)

    def home(self, *args):
        """Event when the home button is pressed"""

        # SOURCE: overwriting builtin-function
        # make sure the correct viewing limits are used when pressing home
        self.parent.annotationManager.push_back()
        self.parent.setToDefaultViewLimits()
        self.parent.onUpdateView(None)
        self.parent.canvas.draw()
        self.push_current()

    def back(self, *args):
        """Event when the back button is pressed"""

        # SOURCE: overwriting a built-in function
        self.parent.annotationManager.push_forward()
        MyNavigationToolbar2WxAgg.back(self, *args)
        self.parent.annotationManager.pop_back()

    def forward(self, *args):
        """Event when the forward button is pressed"""

        # SOURCE: overwriting a built-in function
        self.parent.annotationManager.push_back(False)
        MyNavigationToolbar2WxAgg.forward(self, *args)
        self.parent.annotationManager.pop_forward()

    def pan(self, *args):
        """Event when the pan button is pressed"""

        # SOURCE: overwriting builtin-function
        self.toggleButtons(self._NTB2_PAN)
        if self.myMode == 'PAN':
            self.myMode = ''
        else:
            self.myMode = 'PAN'
        MyNavigationToolbar2WxAgg.pan(self, *args)

    def save(self, evt, filename=None):
        # SOURCE: overwriting builtin-function

        # this may be called using the command-line command -export,
        # if so, then evt is None, and filename is a string.

        if evt == graphInfo.cmdline:
            # called by command-line, filename is supplied
            dirname = os.path.dirname(filename)
            filename = os.path.basename(filename)
        else:
            # need to ask user for filename,
            # Fetch the required filename and file type.
            filetypes = self.canvas._get_imagesave_wildcards()
            dlg =wx.FileDialog(self._parent, "Save to file", "", "", filetypes,
                               wx.SAVE|wx.OVERWRITE_PROMPT|wx.CHANGE_DIR)
            if dlg.ShowModal() == wx.ID_OK:
                dirname  = dlg.GetDirectory()
                filename = dlg.GetFilename()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return

        DEBUG_MSG('Save file dir:%s name:%s' % (dirname, filename), 3, self)
        self.canvas.print_figure(os.path.join(dirname, filename))

    def save_grayscale(self, evt):
        """Save the graph as a grayscale image"""

        # Fetch the required filename and file type.
        filetypes = self.get_grayscale_imagesave_wildcards()
        dlg =wx.FileDialog(self._parent, "Save to file", "", "", filetypes,
                           wx.SAVE|wx.OVERWRITE_PROMPT|wx.CHANGE_DIR)
        if dlg.ShowModal() == wx.ID_OK:
            # check if the filename has the right extension
            filename = dlg.GetPath()
            extension = os.path.splitext(filename)[1]
            if not extension or not extension.lower() == '.png':
                filename = filename + '.png'
            self.canvas.print_figure(filename)
            im = Image.open(filename)
            im = ImageOps.grayscale(im)
            im.save(filename)

    def get_grayscale_imagesave_wildcards(self):
        """return the wildcard string for the grayscale filesave dialog"""
        return "PNG (*.png)|*.png"

    def onAnnotatePoint(self, event):
        """Event when the annotate point button is pressed"""

        self.toggleButtons(self._NTB2_ANN_POINT)

        if self.myMode == 'ANNOTATE_POINT':
            self.myMode = ''
        else:
            self.myMode = 'ANNOTATE_POINT'

    def onAnnotateThreshold(self, event):
        """Event when the annotate threshold button is pressed"""
        ann.annotateThreshold(self.parent)

    def onAnnotateMax(self, event):
        """Event when the annotate max toolbar button is pressed"""

        self.toggleButtons(self._NTB2_ANN_MAX)

        if self.myMode == 'ANNOTATE_MAX':
            self.myMode = ''
        else:
            self.myMode = 'ANNOTATE_MAX'

    def onAnnotateMin(self, event):
        """Event when the annotate min toolbar button is pressed"""

        self.toggleButtons(self._NTB2_ANN_MIN)

        if self.myMode == 'ANNOTATE_MIN':
            self.myMode = ''
        else:
            self.myMode = 'ANNOTATE_MIN'

    def onAnnotateMove(self, event):
        """Event when the annotate min toolbar button is pressed"""

        self.toggleButtons(self._NTB2_ANN_MOVE)

        if self.myMode == 'ANNOTATE_MOVE':
            self.myMode = ''
        else:
            self.myMode = 'ANNOTATE_MOVE'

    def zoom_x(self, event):
        """Event when zooming only the x-axis"""

        self.toggleButtons(self._NTB2_ZOOM_X)

        if self.myMode == 'ZOOM_X':
            self.myMode = ''
        else:
            self.myMode = 'ZOOM_X'

    def zoom_y(self, event):
        """Event when zooming only the y-axis"""

        self.toggleButtons(self._NTB2_ZOOM_Y)

        if self.myMode == 'ZOOM_Y':
            self.myMode = ''
        else:
            self.myMode = 'ZOOM_Y'

    def onButtonPress(self, event):
        """Event when the mouse button is pressed"""

        if event.inaxes and event.button == 1:
            if self.myMode == 'ANNOTATE_POINT':
                self.annotatePoint(event.xdata, event.ydata)
            elif self.myMode == 'ANNOTATE_MAX':
                self.mouse_press = (event.x, event.y)
            elif self.myMode == 'ANNOTATE_MIN':
                self.mouse_press = (event.x, event.y)
            elif self.myMode == 'ZOOM_X':
                ymin, ymax = self.axes.get_ylim()
                ymin = dataToPixel(self.axes, y=ymin)
                self.mouse_press = (event.x, ymin)
            elif self.myMode == 'ZOOM_Y':
                xmin, xmax = self.axes.get_xlim()
                xmin = dataToPixel(self.axes, x=xmin)
                self.mouse_press = (xmin+1, event.y)

    def onButtonRelease(self, event):
        """Event when the mouse button is released"""

        if self.mouse_press == None:
            return

        if self.myMode == 'ANNOTATE_MAX':
            x0, y0 = pixelToData(self.axes, *self.mouse_press)
            x1, y1 = pixelToData(self.axes, event.x, event.y)
            ann.annotateMax(self.parent, x0, y0, x1, y1)
            self.parent.draw()
        elif self.myMode == 'ANNOTATE_MIN':
            x0, y0 = pixelToData(self.axes, *self.mouse_press)
            x1, y1 = pixelToData(self.axes, event.x, event.y)
            ann.annotateMin(self.parent, x0, y0, x1, y1)
            self.parent.draw()
        elif self.myMode == 'ZOOM_X':
            x0, y0 = pixelToData(self.axes, *self.mouse_press)
            x1, y1 = pixelToData(self.axes, event.x, event.y)
            self.parent.zoom_x(min(x0, x1), max(x0, x1))
        elif self.myMode == 'ZOOM_Y':
            x0, y0 = pixelToData(self.axes, *self.mouse_press)
            x1, y1 = pixelToData(self.axes, event.x, event.y)
            self.parent.zoom_y(min(y0, y1), max(y0, y1))

        self.mouse_press = None
        self.release(event)

    def onMouseMotion(self, event):
        """Event when mouse moves"""

        if self.myMode == 'ANNOTATE_MAX' or self.myMode == 'ANNOTATE_MIN':
            if self.mouse_press and self.axes.in_axes(event.x, event.y):
                x0, y0 = self.mouse_press
                self.draw_rubberband(None, x0, y0, event.x, event.y)
        if self.myMode == 'ZOOM_X':
            if self.mouse_press and self.axes.in_axes(event.x, event.y):
                x0, y0 = self.mouse_press
                ymin, ymax = self.axes.get_ylim()
                ymax = dataToPixel(self.axes, y=ymax)
                self.draw_rubberband(None, x0, y0, event.x, ymax-1)
        if self.myMode == 'ZOOM_Y':
            if self.mouse_press and self.axes.in_axes(event.x, event.y):
                x0, y0 = self.mouse_press
                xmin, xmax = self.axes.get_xlim()
                xmax = dataToPixel(self.axes, x=xmax)
                self.draw_rubberband(None, x0, y0, xmax, event.y)

    def annotatePoint(self, x, y):
        ann.annotatePoint(self.parent, x, y)

    def toggleButtons(self, id):
        """Given an id of a toggle button being clicked, unclicks all
        other appropriate buttons
        
        """

        # this checks if some other toggle button other than the two
        # builtin ones (pan and zoom) was clicked while pan and zoom
        # was on.  if so, pretend that the particular button was pressed,
        # thus turning off any kind of internal flags and stuff for them
        if id != self._NTB2_PAN and id != self._NTB2_ZOOM:
            if self.GetToolState(self._NTB2_PAN):
                MyNavigationToolbar2WxAgg.pan(self, None)
            elif self.GetToolState(self._NTB2_ZOOM):
                MyNavigationToolbar2WxAgg.zoom(self, None)

        for eachId in self.buttonIds:
            if eachId != id:
                self.ToggleTool(eachId, False)

class FileToolbar(wx.ToolBar):

    def __init__(self, parent):
        wx.ToolBar.__init__(self, parent, -1, style=wx.TB_HORIZONTAL)

        # create toolbar buttons
        utility.toolbar.addToolbarLabels(parent, self, self.getToolBarData())

        # add a choice box then a button for zooming in on a certain chromosome
        choices = [str(x) for x in range(24)]
        self.chrChoice = wx.Choice(self, -1, choices=choices)
        self.AddControl(self.chrChoice)

        zoom_bmp = getIcon('zoom.xpm')
        zoom_data = [('Zoom', zoom_bmp, 'Zoom', 'Zoom to Chromosome', 
                       parent.onZoomToChr),]
        utility.toolbar.addToolbarLabels(parent, self, zoom_data)
        self.chrChoice.Bind(wx.EVT_CHOICE, self.onChrChoice, self.chrChoice)

        self.Realize()

    def getToolBarData(self):
        """The labels, bitmaps, and handles used for the toolbar"""

        parent = self.GetParent()

        # get some bitmaps
        opengraph_bmp = getIcon('open_graph.xpm')
        savegraph_bmp = getIcon('save_graph.xpm')
        opentext_bmp = getIcon('open_text.xpm')
        saveImage_bmp = getIcon('save_image.xpm')
        print_bmp = wx.ArtProvider.GetBitmap(wx.ART_PRINT, wx.ART_TOOLBAR)
        hline_bmp = getIcon('line.xpm')

        # the labels, bitmaps, and handlers used for the toolbar buttons
        toolbarData = [('Open', opengraph_bmp, 'Open grx file', 
                        'Open grx file', parent.onOpenGRX),
                       ('Import', opentext_bmp, 'Import Data From Text File', 
                        'Import Data From Text File', parent.onOpenTextFile),
                       ('Save', savegraph_bmp, 'Save grx file', 
                        'Save grx file', parent.onSaveGRX),
                       ('Export', saveImage_bmp, 'Export To Image',
                        'Export to Image', parent.navTb.save),
                       ('Print', print_bmp, 'Print', 'Print Current View',
                                                         parent.onPrint),
                       ('', None, None, None, None),
                       ('Insert Horizontal Line', hline_bmp,
                        'Insert Horizontal Line',
                        'Insert a static horizontal line', 
                        parent.insertHLine),
                       ('', None, None, None, None)]

        return toolbarData

    def onChrChoice(self, event):
        """Event when the chromosome choice box is changed"""
        # nothing is done here since nothing happens until the zoom to chr 
        # button is pressed
        pass

    def update(self):
        """Refreshes information in toolbar controls"""

        self.chrChoice.Clear()
        choices = self.GetParent().lines.getNames()
        for eachChoice in choices:
            self.chrChoice.Append(eachChoice)
        self.chrChoice.SetSelection(0)
