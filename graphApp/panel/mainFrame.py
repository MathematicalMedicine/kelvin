"""The main window of the application"""

import pickle
import os
import subprocess
import sys
from tempfile import TemporaryFile
import urllib2
from urlparse import urlparse
import wx

from pylab import arange
import matplotlib
from matplotlib.figure import Figure
from matplotlib.text import Annotation

from colors import hexcolor
from config import misc
from config.config import Config
import config.graphInfo as graphInfo
import io.graph
import io.grx
import io.textFileReader as txtfr
from panel.configAllLinesDialog import ConfigAllLinesDialog
from panel.configDialog import ConfigDialog
from panel.configGeneMarkersDialog import ConfigGeneMarkerDialog
from panel.configLineSetsDialog import ConfigLineSetsDialog
from panel.fileFormatDialog import FileFormatDialog
from panel.figureCanvas import FigureCanvas
from panel.insertGeneMarkerDialog import InsertGeneMarkerDialog
import panel.mainFrame_parts.annotation as ann
import panel.mainFrame_parts.command_line as cl
from panel.overlapDialog import OverlapExistingChrDialog
from panel.scrolledPanel import ScrolledPanel
from panel.viewLimitsDialog import ViewLimitsDialog
from transforms import pixelToAxes
from widgets import graph
from widgets.annotations import my_mpl_annotation
from widgets.geneMarkers import GeneMarkers, GeneMarkersGroup
from widgets.toolbar import NavigationToolbar, FileToolbar
from utility import menuBar
from utility import toolbar
from utility import util

matplotlib.use('WXAgg')

latest_mainFrame = None

def setMainFrame(mainFrame):
    """Set the latest mainFrame object"""
    global latest_mainFrame
    latest_mainFrame = mainFrame

def getMainFrame():
    """Return the latest instance of mainFrame. Handy when un-pickling"""
    return latest_mainFrame

class MainFrame(wx.Frame):
    """The main frame (window) that the application uses"""

    # the filter used for what files to show in the dialog box
    wildcard = "Graph file (*.graph)|*.graph|All files (*.*)|*.*"
    wildcard2 = "Graph file (*.grx)|*.grx|All files (*.*)|*.*"

    def __init__(self):

        wx.Frame.__init__(self, None, -1, title=misc.mainWindowTitle, 
                        size=misc.mainWindowSize, pos=misc.mainWindowPosition)

        # register myself as the latest mainFrame
        setMainFrame(self)

        self.figure = Figure()
        self.axes = self.figure.gca()
        self.legend = self.axes.legend()

        # the current filename that this graph is saved under
        self.filename = None

        # the last directories used when opening and saving files
        self.last_dir = os.getcwd()

        # object that stores all lines, settings, comments, and annotations
        self.graph = graph.Graph(self)

        # all graph settings
        self.config = self.graph.config

        # all the lines in the graph
        self.lines = self.graph.lineCollection

        # a list of floats that specify y-values of horizontal lines
        self.HLines = self.graph.HLines

        # will handle annotations
        self.annotationManager = self.graph.annotationManager

        # list of geneMarkersGroups
        self.geneMarkerGroups = self.graph.geneMarkerGroups

        # a counter that determines what color the next line inserted will be
        self.colorIndex = 0

        # the minimum length the x-axis must span in 
        # order for tickmarks to appear
        self.minTickLength = 2200

        # the tickInterval special value when there are no tickmarks present
        self.noLabelInterval = 10000

        # the tickInterval special value when only tickmarks 
        # are shown at group boundaries
        self.boundaryInterval = 9999

        # the interval value, usually the spacing between tick marks,
        # or may be a special value, (see:noLabelInterval and boundaryInverval)
        self.tickInterval = self.noLabelInterval

        # whether a mouse button is pressed down right now or not
        self.mouseBtn = False

        # the annotation that is being moved
        self.move_annotation = None

        # the last cursor type that was displayed
        self.lastCursor = None

        # the point where zooming starts
        self.zoomPoint = None

        self.legend.set_visible(self.graph.config.showLegend)
        self.create_xlabel()

        # create the window and other doodads
        mainSizer = wx.BoxSizer(wx.VERTICAL)
        self.createWindow(mainSizer)
        self.createToolBar(mainSizer)
        menuBar.createMenuBar(self, self.getMenuData())
        ann.createAnnPopupMenu(self)
        self.CreateStatusBar()

        # bind some events using the matplotlib api
        self.canvas.mpl_connect('motion_notify_event', self.onMotion)
        self.canvas.mpl_connect('button_press_event', self.onButtonPress)
        self.canvas.mpl_connect('button_release_event', self.onButtonRelease)
        self.canvas.mpl_connect('pick_event', self.onPick)

        # set the default size of the graph relative to the window
        self.axes.set_position(graphInfo.defaultFigureSize)

        # a flag determining if the program should be shown or not.
        # it can be not shown by using a command-line command.
        self.show = True

        self.do_command_line()

    def do_command_line(self):
        """Parse and execute any command-line commands"""

        if hasattr(misc, 'openCmdFile'):
            f = open(misc.openCmdFile, 'r')
            lines = f.readlines()
            f.close()
            s = ''
            for line in lines:
                if len(line) > 0 and line[0] != '#':
                    s += line
            s = s.split()
            cl.parse_execute_commands(self, s)
        else:
            cl.parse_execute_commands(self, sys.argv[1:])

    def getMenuData(self):
        return \
          [('File', (
             ('&Open GRX File...\tCtrl+O', 'Open GRX File', self.onOpenGRX),
             ('&Save GRX File...\tCtrl+S', 'Save GRX Graph', self.onSaveGRX),
             ('&Save As GRX File...', 'Save GRX File As', self.onSaveAsGRX),
             ('&Open Graph File...','Open Graph File',self.onOpenGraph),
             ('&Save Graph File...', 'Save Graph', self.onSaveGraph),
             ('&Save As Graph File...','Save Graph File As',self.onSaveAsGraph),
             ('&Import Data from File...', 
              'Import chromosomes from a text file', self.onOpenTextFile),
             ('&Import Data from Remote File...', 
              'Import data from a url location', self.onOpenRemoteTextFile),
             ('&Export Graph To Color Image...\tCtrl+E', 
              'Export graph to PNG or EPS file', self.navTb.save),
             ('&Export Graph To Grayscale Image...', 
              'Export graph to Grayscale PNG', self.navTb.save_grayscale),
             ('Page Setup..', 'Printer Setup', self.canvas.Printer_Setup),
             ('Print Preview...', 'Print Preview',self.canvas.Printer_Preview),
             ('&Print...\tCtrl+P', 'Print', self.onPrint),
             ('', '', ''),
             ('&Quit\tCtrl+Q', 'Quit', self.onClose))),
           ('Edit', (
             ('&Clear Graph', 'Clear Graph', self.onClearGraph),
             ('&Copy Graph\tCtrl+C', 'Copy figure to clipboard', 
                                       self.canvas.Copy_to_Clipboard),
             ('Remove &Chromosome...', 
              'Remove a chromosome from the graph', self.onRemoveLine),
             ('Remove &Horizontal Line...', 
              'Remove a horizontal line from the graph', self.onRemoveHLine))),
           ('Options', (
             ('&Configure Graph...', 'Configure Graph Options', 
                                                self.onConfigure),
             ('&Configure Line Sets...', 'Configure line sets', 
                                                self.onConfigureLineSets),
             ('&Configure All Lines...', 'Configure all lines at once', 
                                                self.onConfigureAllLines),
             ('&Configure Gene Markers...', 'Configure Gene Markers', 
                                                self.onConfigureGeneMarkers),
             ('&Set x-axis, y-axis limits...', 
              'Set x and y axis limits for graph', self.onSetViewLimits))),
           ('Insert', (
             ('&Horizontal Line','Insert Horizontal Line',self.insertHLine),
             ('&Gene Markers','Gene Markers',self.onInsertGeneMarkers))),
           ('Overlap', (
             ('Insert Overlapping Chromosomes from File...', 
              'Insert Overlapping Chromosomes from File', 
               self.insertOverlapLineFromFile),
             ('Overlap existing chromosomes...', 
              'Overlap two chromosomes that are already in the graph...',
               self.onOverlapChr))),
           ('Help', (
             ('graphApp Manual', 'Open Online Help Manual', self.onHelp),))]
    
    def createWindow(self, mainSizer):
        """Create all the panels and subpanels in the mainFrame"""

        splitter = wx.SplitterWindow(self, -1, style=wx.SP_LIVE_UPDATE)
        splitter.SetSashGravity(1.0)
        self.mainPanel = splitter
        self.canvas = FigureCanvas(splitter, -1, self.figure)

        # create the panel on the right with a texbox for comments,
        # and a text area to display information about graph
        sidePanel = wx.Panel(splitter)
        self.commentBox = wx.TextCtrl(sidePanel, -1, style=wx.TE_MULTILINE | 
                                                           wx.TE_WORDWRAP)
        self.commentBox.Bind(wx.EVT_TEXT, self.onCommentText, self.commentBox)
        self.commentBox.MacCheckSpelling(False)
        self.infoPanel = ScrolledPanel(sidePanel)
        self.infoText = wx.StaticText(self.infoPanel, -1)
        font = wx.Font(12, wx.DEFAULT, wx.NORMAL, wx.NORMAL)
        self.infoText.SetFont(font)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.infoText, 1, wx.EXPAND | wx.UP, 5)
        self.infoPanel.SetSizer(sizer)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.commentBox, 2, wx.EXPAND)
        sizer.Add(self.infoPanel, 1, wx.EXPAND)

        sidePanel.SetSizer(sizer)
        self.infoPanel.SetAutoLayout(1)
        self.infoPanel.SetupScrolling()
        splitter.SplitVertically(self.canvas, sidePanel, -150)
        mainSizer.Add(splitter, 1, wx.EXPAND)
        self.SetSizer(mainSizer)

    def createToolBar(self, sizer):
        """Creates the toolbars used"""

        # two toolbars are used, one is the built in one from matplotlib,
        # and the other is a custom one.  there are two because it crashes
        # if the built-in one is used and additional buttons are added to it

        # the built-in toolbar, with some added modifications
        self.navTb = NavigationToolbar(self, self.canvas)

        # a custom toolbar for my own buttons and functions 
        # (mainly file operations)
        self.fileTb = FileToolbar(self)
        self.SetToolBar(self.fileTb)

        # note: fileTb isn't being added to the sizer because in linux 
        # there would be an extra blank toolbar space being shown
        sizer.Insert(0, self.navTb, 0, wx.EXPAND)

    def GetToolBar(self):
        """Needed by wxpython to determine the current toolbar"""
        return self.fileTb

    def onOpenTextFile(self, event):
        """Event when a text file with data wants to be opened"""

        dlg = wx.FileDialog(self, "Open File", self.last_dir, style=wx.OPEN)
                           
        if dlg.ShowModal() == wx.ID_OK:
            self.openTextFile(dlg.GetPath())
            self.last_dir = os.path.dirname(dlg.GetPath())
        dlg.Destroy()

    def onOpenRemoteTextFile(self, event):
        """Open a remote file described by a url"""

        # ask the user what the url is
        title = 'Enter Url'
        message = 'Please enter the url of the remote file.'
        text = 'http://'
        dlg = wx.TextEntryDialog(self, message, title, text)
        if dlg.ShowModal() == wx.ID_OK:
            try:
                file = urllib2.urlopen(dlg.GetValue())
            except urllib2.URLError, e:
                m = "Error accessing file '%s'" % dlg.GetValue()
                dlg = wx.MessageDialog(self, m,"Error!", wx.OK|wx.ICON_ERROR)
                dlg.ShowModal()
                dlg.Destroy()
            else:
                # the file object returned by urlopen does not have all methods
                # of a normal file object, but readTextFile expects a normal
                # file object, so create a temporary file and place contents
                # of remote file into the temporary file
                tempFile = TemporaryFile()
                tempFile.write(file.read())
                tempFile.seek(0)

                # since the name of the temp file is random and should not
                # be set to another value (there will be an error when the 
                # file closes), find the filename from the url and pass it 
                # on as an argument
                path = urlparse(dlg.GetValue()).path
                filename = os.path.basename(path)

                txtfr.readTextFile(self, tempFile, filename, False)
                tempFile.close()
                file.close()

        dlg.Destroy()

    def onPrint(self, event):
        """Event when the print graph occurs"""

        self.canvas.Printer_Print(event=event)

    def onClose(self, event):
        """Event when close is chosen from the menu bar"""

        self.Destroy()

    def onClearGraph(self, event):
        """Event when the clear graph menu option is chosen"""

        # remove all chromosomes
        for eachName in self.lines.getNames():
            self.lines.removeLine(eachName)

        # remove all horizontal lines
        for eachLine in self.HLines[:]:
            self.removeHLine(eachLine, False)

        # remove all annotations
        for eachAnn in self.annotationManager[:]:
            self.annotationManager.remove(eachAnn)

        self.axes.clear()
        self.axes.set_xlim(xmin=0, xmax=1)
        self.axes.set_ylim(ymin=0, ymax=1)
        self.draw()

    def onConfigure(self, event):
        """Event when the user wants to configure the graph"""

        dlg = ConfigDialog(self)
        dlg.ShowModal()

    def onConfigureLineSets(self, event):
        """Event when the user wants to configure line sets"""

        dlg = ConfigLineSetsDialog(self)
        dlg.ShowModal()

    def onConfigureAllLines(self, event):
        """Event when the user want to configure all lines at once"""

        dlg = ConfigAllLinesDialog(self)
        dlg.ShowModal()

    def onConfigureGeneMarkers(self, event):
        """Event when the user wants to configure gene markers"""

        dlg = ConfigGeneMarkerDialog(self)
        dlg.ShowModal()

    def onInsertGeneMarkers(self, event):
        """Event when the insert gene markers menu item is selected"""

        dlg = wx.FileDialog(self, "Open Gene Marker File", 
                            self.last_dir, style=wx.OPEN)
                           
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            dlg = InsertGeneMarkerDialog(self, -1, filename)
            dlg.CenterOnParent()
            if dlg.ShowModal() == wx.ID_OK:
                self.last_dir = os.path.dirname(filename)
        dlg.Destroy()

    def onSetViewLimits(self, event):
        """Event when the user selects to change the graph viewing limits"""

        xmin, xmax = self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()
        info = [xmin, xmax, ymin, ymax]
        dlg = ViewLimitsDialog(self, info)
        if dlg.ShowModal() == wx.ID_OK:

            self.setViewLimits((info[0], info[1]), (info[2], info[3]))

    def setViewLimits(self, xlimits=None, ylimits=None, firstView=False):
        """Change the viewing limits of the graph, arguments
        should be in data space.
        
        """

        # store the positions of the annotations in the old view
        # don't do this when loading from a .grx file, since it will
        # reset annotation autoLayout information
        if not firstView:
            self.annotationManager.push_back()

        # if the viewing stack is empty push the current one 
        # so the home button is established
        if self.navTb._views.empty():
            self.pushCurrentView()

        if xlimits is not None:
            self.axes.set_xlim(xmin=xlimits[0], xmax=xlimits[1])
        if ylimits is not None:
            self.axes.set_ylim(ymin=ylimits[0], ymax=ylimits[1])
        self.onUpdateView(None)

        # push the current view so the home, back, 
        # forward buttons work right
        self.pushCurrentView()

    def onRemoveLine(self, event):
        """Event when the user wants to remove a line from the graph"""

        choices = self.lines.getNames()
        dlg = wx.MultiChoiceDialog(self, 'Choose chromosome(s) to remove', 
                                   'Remove Chromosome', choices)

        if dlg.ShowModal() == wx.ID_OK:
            selections = dlg.GetSelections()

            # reversing this list is important because you want to remove lines
            # starting from the last index, otherwise the index of the other
            # elements will change after removing lines if you start from the 
            # first index
            selections.reverse()

            for eachIndex in selections:
                self.lines.removeLine(eachIndex)
            if len(selections) > 0:
                self.drawLines()

        dlg.Destroy()

    def removeLine(self, index):
        """Given a line index, remove that line"""
        self.lines.removeLine(index)

    def onRemoveHLine(self, event):
        """Event when the user wants to remove a horizontal line,
        a dialog is shown of the y-values of the horizontal lines in the
        graph, and the user chooses which ones to remove
        
        """

        choices = [str(x) for x in self.HLines]
        dlg = wx.MultiChoiceDialog(self, 'Choose a line(s) to remove',
                                   'Remove Horizontal Line', choices)

        if dlg.ShowModal() == wx.ID_OK:
            selections = dlg.GetSelections()
            values = [float(choices[x]) for x in selections]
            for each_y_value in self.HLines:
                for eachValue in values:
                    if each_y_value == eachValue:
                        self.removeHLine(eachValue, False)
            self.drawLines()

    def openTextFile(self, filename, overlap=False):
        """Open a text file containing data needing to be graphed"""

        try:
            f = open(filename, 'r')
        except IOError, e:
            message = "Error opening file '%s' : %s" % (filename, e.strerror)
            dlg = wx.MessageDialog(self, message,"Error!", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return
        else:
            filename = os.path.basename(f.name)
            txtfr.readTextFile(self, f, filename, overlap)
            f.close()

    def onOpenGraph(self, event):
        """Event when a graph is loaded from a file"""

        # ask user for filename
        dlg = wx.FileDialog(self, "Open File", self.last_dir, 
                            style=wx.OPEN, wildcard=self.wildcard)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.last_dir = os.path.dirname(filename)
            io.graph.openGraph(self, filename)

        dlg.Destroy()

    def onSaveGraph(self, event):
        """Event when the graph is saved to a .graph file"""

        # this happens when the graph is saved for the first time
        if self.filename is None:
            self.onSaveAsGraph(event)
            return

        try:
            io.graph.save(self.graph, self.filename)
        except IOError, e:
            s = "Error saving file '%s' : %s" % (self.filename, e.strerror)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return False
        except pickle.PicklingError, e:
            s = "Error saving file '%s' : %s" % (self.filename, e)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return False
        except Exception, e:
            print 'warning!  caught some exception when saving to .graph file!'
            print e
            s = "Error saving file '%s'" % (self.filename)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return False

        return True

    def onSaveAsGraph(self, event, filename=None):
        """Event when the Save As for a graph is selected"""

        # this function is also called when the -save command is used
        # on the command line, have to take that into account

        if event == graphInfo.cmdline:
            # must have been called from command-line, filename is 
            # already supplied
            pass
        else:
            # user picked to save, show save dialog to get filename
            dlg = wx.FileDialog(self, "Save as...", self.last_dir, 
                                style=wx.SAVE | wx.OVERWRITE_PROMPT, 
                                wildcard=self.wildcard)
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return

        # save the file
        self.last_dir = os.path.dirname(filename)

        # if the new filename doesn't have an extension, add it
        if not os.path.splitext(filename)[1]:
            filename = filename + '.graph'

        # try to save it
        oldFilename = self.filename
        self.filename = filename
        if self.onSaveGraph(None):
            # set title to name of file
            baseFilename = os.path.basename(filename)
            self.SetTitle('graph - ' + baseFilename)
        else:
            self.filename = oldFilename

    def onOpenGRX(self, event):
        """Event when a .grx is loaded from a file"""

        # ask user for filename
        dlg = wx.FileDialog(self, "Open File", self.last_dir, 
                            style=wx.OPEN, wildcard=self.wildcard2)
        if dlg.ShowModal() == wx.ID_OK:
            filename = dlg.GetPath()
            self.last_dir = os.path.dirname(filename)
            self.openGRX(filename)

        dlg.Destroy()

    def openGRX(self, filename):
        try:
            newGraph = io.grx.load(self, filename)
        except IOError, e:
            s = "Error opening file '%s' : %s" % (filename,e.strerror)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        except Exception, e:
            print 'Warning! Caught some exception when opening a .grx file!'
            print type(e), e
            s = "Error opening file '%s'" % (filename)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            # set some variables
            self.graph = newGraph
            self.filename = filename
            self.lines = newGraph.lineCollection
            self.annotationManager = newGraph.annotationManager
            self.geneMarkerGroups = newGraph.geneMarkerGroups
            self.HLines = newGraph.HLines
            self.config = newGraph.config
            self.lines.initializeLineSets()

            # update graph
            self.drawLines()
            self.graph.config.update()
            self.annotationManager.update()
            self.setComments(self.graph.comments)
            self.SetTitle(misc.mainWindowTitle + ' - ' + 
                          os.path.basename(filename))

            # need to reset the home, back, and forward buttons
            self.navTb._views.clear()
            xmin, xmax, ymin, ymax = self.graph.config.getViewLimits()

            # there may be a chance that the grx file did not set
            # view limits, so check if they are all still zero
            if xmin == 0 and xmax == 0 and ymin == 0 and ymax == 0:
                self.setToDefaultViewLimits()
                self.graph.config.updateValues()
            else:
                self.setViewLimits((xmin, xmax), (ymin, ymax), firstView=True)

            self.annotationManager.restore_positions()

    def onSaveGRX(self, event):
        """Event when the graph is saved to a .grx file (graph xml)"""

        # this happens when a graph file is saved for the first time
        if self.filename is None or (self.filename.endswith('.graph') and \
                                     event != graphInfo.cmdline):
            self.onSaveAsGRX(event)
            return

        try:
            io.grx.save(self.graph, self.filename)
        except IOError, e:
            s = "Error saving file '%s' : %s" % (self.filename, e.strerror)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return False
        except Exception, e:
            print 'Warning! Caught some exception when saving to .grx file!'
            print type(e), e
            s = "Error saving file '%s'" % (self.filename)
            dlg = wx.MessageDialog(self, s, "Error!", wx.OK | wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return False

        return True

    def onSaveAsGRX(self, event, filename=None):
        """Event when the Save As for a grx is selected"""

        # this function is also called when the -save command is used
        # on the command line, have to take that into account

        if event == graphInfo.cmdline:
            # must have been called from command-line, filename is 
            # already supplied
            pass
        else:
            # user picked to save, show save dialog to get filename
            dlg = wx.FileDialog(self, "Save as...", self.last_dir, 
                                style=wx.SAVE | wx.OVERWRITE_PROMPT, 
                                wildcard=self.wildcard2)
            if dlg.ShowModal() == wx.ID_OK:
                filename = dlg.GetPath()
                dlg.Destroy()
            else:
                dlg.Destroy()
                return

        # save the file
        self.last_dir = os.path.dirname(filename)

        # if the new filename doesn't have an extension, add it
        if not os.path.splitext(filename)[1]:
            filename = filename + '.grx'

        # try to save it
        oldFilename = self.filename
        self.filename = filename
        if self.onSaveGRX(event):
            # set title to name of file
            baseFilename = os.path.basename(filename)
            self.SetTitle('graph - ' + baseFilename)
        else:
            self.filename = oldFilename

    def onButtonPress(self, event):
        """Event when any mouse button is pressed"""

        self.mouseBtn = True
        if self.navTb.myMode == 'ZOOM' and event.inaxes:
            # starting to zoom, keep track of starting point
            self.zoomPoint = (event.x, event.y)

        if self.navTb.myMode == 'PAN' and event.inaxes:
            self.annotationManager.push_back()
            self.annotationManager.beginPan(event)

    def onButtonRelease(self, event):
        """Event when any mouse button is released"""

        # check to see if a zoom is about to happen
        # if so, save annotation positions
        if self.navTb.myMode == 'ZOOM' and self.zoomPoint:
            # check if area is above 5 pixel threshold
            lastx, lasty = self.zoomPoint
            x, y = event.x, event.y
            self.zoomPoint = None
            if abs(x-lastx) >= 5 and abs(y-lasty) >= 5:
                self.annotationManager.push_back()
                self.onUpdateView(None)

        # check if a pan is ending
        if self.isPanning():
            self.onUpdateView(None)

        self.mouseBtn = False
        self.move_annotation = None

    def onMotion(self, event):
        """Event when the mouse moves"""

        self.updateMousePosText(event)
        self.updateMouseCursor(event)
        ann.moveAnnotation(self, event)

        if self.isPanning():
            self.onUpdateView(None)
            self.annotationManager.updatePan(event)

        # only update the position when not drawing a rubberband, otherwise
        # when the canvas is redrawn the rubberband box doesn't show up
        if not self.isDrawingRubberband():
            self.draw()

    def onHelp(self, event):
        """Event when the open manual help menu item is chosen"""

        info = wx.AboutDialogInfo()
        info.Name = "graphApp Manual"
        site = 'http://hodgkin.ccri.net/software/graphApp/graphApp_manual.htm'
        info.WebSite = (site, site)
        wx.AboutBox(info)

    def updateMousePosText(self, event):
        """Update the text that displays the mouse pointer coordinates
        relative to the chromosome.
        
        """

        if hasattr(self, 'mousePosText'):
            if event.inaxes:
                for eachGroup in self.lines.getGroups():
                    if eachGroup.isWithin(event.xdata):
                        x = float(event.xdata)-eachGroup.xmin+eachGroup.start
                        s = 'chr %s (%.4f, %.4f)' % \
                                              (eachGroup.label,x,event.ydata)
                        metadata = self.lines.findClosestMetada(event.xdata, 
                                                                event.ydata)
                        for eachItem in metadata:
                            s += ' ' + eachItem
                        self.mousePosText.set_text(s)
                        break
            else:
                self.mousePosText.set_text(' ')

    def updateMouseCursor(self, event):
        """Change the mouse cursor depending on the current mode"""

        # update the mouse cursor
        cross_cursor_modes = ['ANNOTATE_MAX', 'ANNOTATE_MIN', 'ZOOM_X', 
                              'ZOOM_Y']
        if self.navTb.myMode in cross_cursor_modes and event.inaxes:
            if self.lastCursor != wx.CURSOR_CROSS:
                self.canvas.SetCursor(wx.StockCursor(wx.CURSOR_CROSS))
                self.lastCursor = wx.CURSOR_CROSS        
        elif self.lastCursor != None:
            self.canvas.SetCursor(wx.NullCursor)
            self.lastCursor = None

    def onCommentText(self, event):
        """Event when text is entered in the comments text box"""
        self.setComments(self.commentBox.GetValue())

    def setComments(self, comments):
        """Set the comment box to a certain value, along with graph.comments"""

        self.graph.comments = comments
        self.commentBox.ChangeValue(comments)

    def onUpdateView(self, event):
        """Event when the view of the graph is changed, from zooming or panning
        or hitting the back, forward, or home toolbar buttons.
        
        """

        self.updateTickmarks(False)
        self.updateLegend(False)
        self.drawLabels(False)
        if not self.isPanning():
            self.annotationManager.layout()
        self.updateInfoText()

        self.draw()

    def drawLines(self):
        """Clear the figure and draw all the current lines"""

        # the actual x-position of where tickmarks will go. there will be 
        # one tick at the end of each line (chromosome), so the tickmarks 
        # will be hidden by the vertical lines that separate chromosomes
        ticks = []

        # the accumulating total of where each line is offset
        offset = 0

        # clear the figure and draw the chromosome lines
        axes = self.axes
        axes.clear()
        color = hexcolor('grey60')
        for eachGroup in self.lines.getGroups():
            axes.axvline(offset, color=color)
            eachGroup.plot(offset)
            offset += eachGroup.length
            ticks.append(offset)
            eachGroup.text = axes.text(0, 0, '', horizontalalignment='center',
                                       verticalalignment='top', 
                                       transform=axes.transAxes)
            axes.axvline(offset, color=color)

        # make some text that shows the mouse position
        self.mousePosText = axes.text(0.0, -0.1, '',horizontalalignment='left',
                           verticalalignment='bottom',transform=axes.transAxes)

        self.drawHLines(False)

        # re-add annotations
        self.annotationManager.add_annotations()

        axes.set_xticks(ticks)
        axes.set_xticklabels([])
        self.config.update()
        self.setToDefaultViewLimits()
        self.onUpdateView(None)
        self.fileTb.update()
        
        # needed to properly update the home, back, forward buttons
        self.pushCurrentView()

    def drawLabels(self, refresh=True):
        """Draws the chromosome labels underneath the x-axis"""

        axes = self.figure.gca()
        xmin, xmax = axes.get_xlim()
        ymin, ymax = axes.get_ylim()

        # go through each group, and add the labels for all chromosomes in
        # that group   
        for eachGroup in self.lines.getGroups():
            # blank out the text labels if they aren't shown in the graph
            if xmax <= eachGroup.xmin:
                eachGroup.text.set_text('')
                continue
            if xmin >= eachGroup.xmax:
                eachGroup.text.set_text('')
                continue

            start = max(eachGroup.xmin, xmin)
            end = min(eachGroup.xmax, xmax)
            point = (end - start) / 2.0 + start
            x = (point - xmin)  / (xmax - xmin)
            y = graphInfo.chrlabel_y_pos

            eachGroup.text.set_text(eachGroup.label)
            eachGroup.text.set_x(x)
            eachGroup.text.set_y(y)

        if refresh:
            self.draw()

    def changeLineName(self, line=None, group=None, newName=None):
        """Change the name of a line"""

        self.lines.changeLineName(line, group, newName)

        # the file toolbar needs to know of line name changes
        self.fileTb.update()

        # update the chromosome name labels underneath the graph
        self.drawLabels()

        # legend should know about this
        self.updateLegend(False)

    def isPanning(self):
        """Returns True if panning is going on else False"""

        # need to determine whether panning is occuring or not
        # note that self.navTb._NTB2_PAN is a variable defined
        # in the matplotlib source code, so this could break with
        # later releases
        isPanOn = self.navTb.GetToolState(self.navTb._NTB2_PAN)
        if isPanOn and self.mouseBtn:
            return True
        else:
            return False

    def isZooming(self):
        """Returns True if zooming is going on else False"""

        # SOURCE:
        # need to determine whether zooming is occuring or not
        # note that self.navTb._NTB2_ZOOM is a variable defined
        # in the matplotlib source code, so this could break with
        # later releases
        isZoomOn = self.navTb.GetToolState(self.navTb._NTB2_ZOOM)
        if isZoomOn and self.mouseBtn:
            return True
        else:
            return False

    def isDrawingRubberband(self):
        """Returns True if a rubberband is currently being drawn"""

        # a rubberband is the box that is shown while a user chooses
        # what area to zoom in, or what area to annotate the max or min with
        # you need to know this because if the canvas is drawn, the
        # rubberband will disappear
        if self.mouseBtn:
            mode = self.navTb.myMode
            if mode=='ANNOTATE_MAX' or mode=='ANNOTATE_MIN' or \
                 mode=='ZOOM_X' or mode == 'ZOOM_Y' or self.isZooming():
                return True
        return False

    def setToDefaultViewLimits(self):
        """Set the viewing limits of the graph to default values"""

        axes = self.figure.gca()

        xmin = self.lines.get_x_min()
        xmax = self.lines.get_x_max()
        axes.set_xlim(xmin=xmin, xmax=xmax)
        axes.set_ylim(ymin=0, ymax=1)

    def insertHLine(self, event, value=None):
        """Insert a horizontal line at a specified y-value"""

        if event == graphInfo.cmdline:
            # must have been called from command-line, 
            # value is already supplied
            self.HLines.append(value)
        else:
            # need a dialog box to ask the user where to put horizontal lines
            text1 = "Specify the y-values of the line \n"
            text2 = "(separate multiple values with commas)"
            dlg = wx.TextEntryDialog(self, text1 + text2, 
                                     "Please specify the y-values of the line")
            if dlg.ShowModal() == wx.ID_OK:
                for eachValue in dlg.GetValue().split(', '):
                    try:
                        value = float(eachValue)
                    except ValueError:
                        s = 'Error: ' + eachValue + ' is not a valid number.'
                        dlg2 = wx.MessageDialog(self, s, "Error!", wx.OK | 
                                                wx.ICON_ERROR)
                        dlg2.ShowModal()
                        dlg2.Destroy()
                    else:
                        self.HLines.append(value)
            dlg.Destroy()

        self.HLines.sort()

        # draw the lines
        self.drawHLines()

    def insertOverlapLineFromFile(self, event):
        dlg = wx.FileDialog(self, "Open File", self.last_dir, style=wx.OPEN)
                           
        if dlg.ShowModal() == wx.ID_OK:
            self.openTextFile(dlg.GetPath(), overlap=True)
            self.last_dir = os.path.dirname(dlg.GetPath())

    def onOverlapChr(self, event):
        """Overlap two existing lines on the graph togeher"""

        # first check if there's at least two lines on the graph
        if len(self.lines) < 2:
            s1 = "Error, there are less than two lines in the graph."
            s2 = "Open a file to insert more lines in the graph first."
            message = s1 + '\n' + s2
            dlg = wx.MessageDialog(self, message,"Error!",wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return

        # info will contain what chromosome to move where, the format is
        # note that all values are indices to the chromosomes
        # <list of source chromosomes> <destination chromosome>
        info = [None, None]
        dlg = OverlapExistingChrDialog(self, -1, info, self.lines.getNames())
        if dlg.ShowModal() == wx.ID_OK:
            self.lines.overlap(*info)
            self.drawLines()

        dlg.Destroy()

    def drawHLines(self, refresh=True):
        """Draw horizontal lines specified by user"""

        axes = self.figure.gca()

        # the view boundaries may change when drawing the lines 
        # (for some stupid reason...)
        # so save them and re-apply them later
        xmin, xmax = axes.get_xlim()
        ymin, ymax = axes.get_ylim()

        for each_y_value in self.HLines:
            axes.axhline(each_y_value, color='k')

        axes.set_xlim(xmin=xmin, xmax=xmax)
        axes.set_ylim(ymin=ymin, ymax=ymax)

        if refresh:
            self.draw()

    def removeHLine(self, value, refresh=False):
        """Remove a horizontal line with the value specified"""

        if self.HLines.count(value) > 0:
            self.HLines.remove(value)

        if refresh:
            self.drawLines()

    def updateInfoText(self):
        """Update the text in the sidebar that show graph information"""

        # show all the troughs and peaks
        text = ''
        for eachAnn in self.annotationManager:
            if eachAnn.get_visible():
                text += eachAnn.getInfoText()

        self.infoText.SetLabel(text)
        self.infoPanel.FitInside()

    def updateLegend(self, refresh=True):
        """Update the contents of the legend"""

        self.legend.set_visible(self.graph.config.showLegend)
        if self.graph.config.showLegend is False:
            return

        if self.config.showLineSetsInLegend:
            # show each line set in the legend.
            # first gather representative lines for each set
            lines = []
            labels = []
            for eachSet in self.lines.lineSets:
                lines.append(eachSet.getRepresentativeLine())
                labels.append(eachSet.name)
        else:
            # show each individual line in the graph.
            # check which lines are needed in the legend,
            # only put a line in the legend if it is actually being shown
            xmin, xmax = self.axes.get_xlim()
            lines = []
            labels = []
            inside = False

            for eachGroup in self.lines.getGroups():
                # take advantage that the groups appear in 
                # increasing order along the x-axis
                if eachGroup.xmin >= xmax:
                    break
                if not inside: 
                    if eachGroup.xmin >= xmin and eachGroup.xmax <= xmax:
                        inside = True
                    elif eachGroup.xmin >= xmin and eachGroup.xmin < xmax:
                        inside = True
                    elif eachGroup.xmax > xmin and eachGroup.xmax <= xmax:
                        inside = True
                    elif eachGroup.xmin <= xmin and eachGroup.xmax >= xmax:
                        inside = True
                if inside:
                    labels.extend(eachGroup.getNames())
                    lines.extend(eachGroup.get_mpl_lines())

            if len(lines) == 0:
                # nothing to show, otherwise there will be a blank box in
                # the middle of the graph
                self.legend.set_visible(False)
                return

        self.legend = self.axes.legend(lines, labels, 
                                       loc=self.graph.config.legendPosition)
        self.legend.draw_frame(self.graph.config.showLegendFrame)

        if refresh:
            self.draw()

    def updateTickmarks(self, refresh=True):
        """Determine appropriate spacing between tickmarks and create them"""

        # tickmark spacing is based on the length of the x-axis currently shown
        xmin, xmax = self.axes.get_xlim()
        xdistance = xmax - xmin

        # first check if tickmarks are needed at all
        # tick marks are not needed if zoomed out enough
        if xdistance >= self.minTickLength:
            # this checks if tickmarks are currently cleared out
            if self.tickInterval != self.noLabelInterval:
                # clear out the tick marks
                self.tickInterval = self.noLabelInterval
                ticks=[eachGroup.xmax for eachGroup in self.lines.getGroups()]
                tick_labels = []
                self.applyTickMarks(ticks, tick_labels)
        else:
            # tickmarks are needed, now determine what interval to use

            # the params variable is a list of the different parameters 
            # to use when deciding what spacing of tickmarks to use
            # format is (<min x distance>, <value to set tickInverval to>, 
            #              (<interval amount>, <min distance between ticks>, 
            #               <number of decimal places for tickmarks>))
            # note that each entry in params must be listed in decreasing
            # order of the min x distance
            params = [(1300,self.boundaryInterval,(self.boundaryInterval,0,0)),
                      (900, 100, (100, 40, 0)),
                      (420, 50, (50, 35, 0)),
                      (180, 25, (25, 15, 0)),
                      (100, 10, (10, 5, 0)),
                      (25, 5, (5, 3, 0)),
                      (8, 1, (1, 0.4, 0)),
                      (3, 0.5, (0.5, .4, 1)),
                      (1, 0.25, (0.25, 0.16, 2)),
                      (0.35, 0.1, (0.1, 0.07, 2)),
                      (0.15, 0.05, (0.05, 0.025, 2)),
                      (0.01, 0.01, (0.01, 0.01, 2)),
                      (0, 0.001, (0.001, 0.0007, 3))]

            for eachEntry in params:
                if xdistance > eachEntry[0]: 
                    self.tickInterval = eachEntry[1]
                    ticks,tick_labels = self.createTickMarks(*eachEntry[2])
                    self.applyTickMarks(ticks, tick_labels)
                    break

        if refresh:
            self.draw()

    def createTickMarks(self, interval, minDistance, precision):
        """Create a list of tick positions and labels given the interval of 
        the tickmarks and the minimum distance between each tickmark. The
        precision argument is how many decimal places to show the numbers.
        
        """

        def getStartPoint(group):
            """Given a group, find where tick marks should start"""

            # the left boundary does not count as the first tick mark,
            # the one after that is the first one. 
            # this function returns the offset amount past the left boundary.
            # left boundary can be group boundary or left side of graph

            # check for easy case first
            if group.xmin >= xmin and group.xmin <= xmax:
                return interval
            else:
                x = max(xmin, eachGroup.xmin) - eachGroup.xmin
                r = x % interval
                x = x-r
                if x + group.xmin < xmin:
                    x += interval
                if x == 0:
                    x += interval
                return x

        def numToText(num):
            """Turn a number to text with the correct precision"""
            return util.numToText(num, precision, 'extend')

        # create only tickmarks in the interval that is currently being shown
        xmin, xmax = self.axes.get_xlim()
        ticks = []
        tick_labels = []
        groups = self.lines.getGroups()

        # add the very first tickmark if the leftmost group is visible
        if xmin <= 0 and len(groups) > 0:
            ticks.append(0)
            tick_labels.append(numToText(groups[0].start))

        # loop through each group, adding tickmarks. the tickmarks added by
        # this loop include the first one after the left boundary of the group
        # all the way until the right boundary
        for eachGroup in groups:
            if eachGroup.doesOverlap(xmin, xmax):
                x = getStartPoint(eachGroup)
                offset = eachGroup.start
                while x < eachGroup.length and x + eachGroup.xmin < xmax:
                    ticks.append(x + eachGroup.xmin)
                    tick_labels.append(numToText(x + offset))
                    x += interval
                # remove the last tick mark if it's 
                # too close to the right boundary
                if len(ticks) > 0 and eachGroup.xmax - ticks[-1] < minDistance:
                    ticks.pop()
                    tick_labels.pop()
                # now decide what to put at the right boundary, either 
                # the group's max or zero. if it's the very first tickmark for
                # a group, make it zero, otherwise use the group's max
                if eachGroup.xmax >= xmin and eachGroup.xmax <= xmax:
                    if len(ticks) > 0:
                        ticks.append(eachGroup.xmax)
                        tick_labels.append(numToText(eachGroup.end))
                    else:
                        ticks.append(eachGroup.xmax)
                        tick_labels.append(numToText(0))

        # take care of the special case that if there is only one tickmark 
        # at the end, don't label the end of the group as 0, but as the 
        # group's xmax. this happens when only the rightmost part of the 
        # last group is visible. 
        if len(ticks) == 1 and len(groups) > 0:
            group = self.lines.getGroups()[-1]
            if group.xmax >= xmin and group.xmax <= xmax:
                tick_labels[0] = numToText(group.end)

        return ticks, tick_labels

    def applyTickMarks(self, ticks, tick_labels):
        """Given a list of tick positions and labels, applies them"""

        # need to save the current view since changing 
        # tick marks will change the viewing parameters
        xmin, xmax = self.axes.get_xlim()
        ymin, ymax = self.axes.get_ylim()

        self.axes.set_xticks(ticks)
        self.axes.set_xticklabels(tick_labels)

        self.axes.set_xlim(xmin=xmin, xmax=xmax)
        self.axes.set_ylim(ymin=ymin, ymax=ymax)

    def create_xlabel(self):
        """Create the xlabel for the graph"""

        # because matplotlib is stupid and i can't adjust the y-position
        # of the builtin xlabel, create a text object that will
        # serve as the xlabel
        x1 = graphInfo.defaultFigureSize[0]
        x2 = graphInfo.defaultFigureSize[2] + x1
        x = (x1+x2) / 2.0
        label = self.graph.config.xlabel
        self.xlabel = self.figure.text(x, graphInfo.xlabel_y_pos, label, 
                 horizontalalignment='center', verticalalignment='top',
                 transform=self.figure.transFigure)

    def getNextColor(self):
        """Returns the next color along the list of predefined colors"""

        color = graphInfo.colors[self.colorIndex]
        self.colorIndex = (self.colorIndex + 1) % len(graphInfo.colors)
        return color

    def onZoomToChr(self, event):
        """Event when user clicked the zoom to 
        chromosome button on the toolbar
        
        """

        # if the viewing stack is empty push the current one 
        # so the home button is established
        if self.navTb._views.empty():
            self.pushCurrentView()
        self.annotationManager.push_back()

        index = self.fileTb.chrChoice.GetSelection()
        line = self.lines.getLine(index)
        self.axes.set_xlim(xmin=line.getStart(), xmax=line.getEnd())
        self.onUpdateView(None)

        # push the current view so the home, back, forward buttons work right
        self.pushCurrentView()

    def zoom_x(self, xmin, xmax):
        """Zoom in to the given x-coordinates"""

        self.annotationManager.push_back()

        # if the viewing stack is empty push the current one 
        # so the home button is established
        if self.navTb._views.empty():
            self.pushCurrentView()

        self.axes.set_xlim(xmin=xmin, xmax=xmax)
        self.onUpdateView(None)

        # push the current view so the home, back, forward buttons work right
        self.pushCurrentView()

    def zoom_y(self, ymin, ymax):
        """Zoom in to the given y-coordinates"""

        self.annotationManager.push_back()

        # if the viewing stack is empty push the current one 
        # so the home button is established
        if self.navTb._views.empty():
            self.pushCurrentView()

        self.axes.set_ylim(ymin=ymin, ymax=ymax)
        self.onUpdateView(None)

        # push the current view so the home, back, forward buttons work right
        self.pushCurrentView()

    def onPick(self, event):
        """Event when user clicks on an annotation"""

        #if isinstance(event.artist, Annotation):
        if isinstance(event.artist, my_mpl_annotation):
            if event.mouseevent.button == 1:
                # need the Annotation class, not the mpl_annotation class
                self.move_annotation = event.artist.parent
            elif event.mouseevent.button == 3:
                self.last_annotation = event.artist
                self.onShowPopup(event)

    def onShowPopup(self, event):
        """Event when the right click popup menu is shown"""

        self.mainPanel.PopupMenu(self.popupmenu)

        # this is needed because doing a right-click on an annotation to show
        # this menu will register the button press, but not the button release.
        # this is a problem when, for example, you edit an annotation, then
        # immediately press the pan button. the program will think you are 
        # panning since the button release was not registered by the program
        self.mouseBtn = False

    def pushCurrentView(self):
        """Store the current view for use with the home, 
        forward, and back buttons
        """

        self.navTb.push_current()

    def draw(self):
        self.canvas.draw()

