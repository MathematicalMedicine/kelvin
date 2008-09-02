"""Contains functions that reads and writes a .graph format file. Used primarily in MainFrame class. This format is obsolete, should use .grx format instead"""

import cPickle as pickle
import os
import wx

from config import misc

def save(graph, filename):
    """Save the given graph to the given filename in the .graph format.  The
    graph argument should be a Graph object (see graph.py). The pickle module
    will be used to save the object to a file.
    
    """

    # make sure everyone has the latest values
    graph.config.updateValues()
    graph.annotationManager.save_positions()

    # do the pickling of the object
    try:
        s = pickle.dumps(graph)
    except:
        raise
    else:
        file = open(filename, 'w')
        file.write(s)
        file.close()

def load(mainFrame, filename):
    """Given a filename, load the .graph file from the file 
    and returns the graph object
    
    """

    file = open(filename)
    try:
        graph = pickle.load(file)
    finally:
        file.close()

    return graph

def openGraph(mainFrame, filename, tryAgain=True):
    try:
        newGraph = load(mainFrame, filename)
    except IOError, e:
        s = "Error opening file '%s' : %s" % (filename,e.strerror)
        dlg = wx.MessageDialog(mainFrame, s, "Error!", wx.OK|wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()
    except pickle.PicklingError, e:
        s = "Error opening file '%s' : %s" % (filename, e)
        dlg = wx.MessageDialog(mainFrame, s, "Error!", wx.OK|wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()
    except AttributeError, e:
        if not tryAgain:
            # give up
            s = "Error opening file '%s' : %s" % (filename, e)
            dlg = wx.MessageDialog(mainFrame, s, "Error!", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        elif str(e).count('module') > 0:
            # this means that some old version of .graph was saved before
            # some modules got moved around. have to open the actual file
            # and fix them
            f = open(filename, 'r')
            lines = f.readlines()
            f.close()

            # the LineGroup and LineCollection object got moved to
            # their own modules
            for i in range(len(lines)):
                if lines[i].count('LineCollection') > 0:
                    # put in new module path in previous line
                    lines[i-1] = lines[i-1].replace('widgets.line', 
                                                  'widgets.lineCollection')
                if lines[i].count('LineGroup') > 0:
                    # put in new module path in previous line
                    lines[i-1] = lines[i-1].replace('widgets.line', 
                                                    'widgets.lineGroup')

            # write the line back to the file
            f = open(filename, 'w')
            f.writelines(lines)
            f.close()

            # try opening .graph file again
            openGraph(mainFrame, filename, tryAgain=False)
    except Exception, e:
        print 'warning!  caught some exception when unpickling!'
        print type(e), e
        s = "Error opening file '%s'" % (filename)
        dlg = wx.MessageDialog(mainFrame, s, "Error!", wx.OK|wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()
    else:
        # set some variables
        mainFrame.graph = newGraph
        mainFrame.filename = filename
        mainFrame.lines = newGraph.lineCollection
        mainFrame.annotationManager = newGraph.annotationManager
        mainFrame.geneMarkerGroups = newGraph.geneMarkerGroups
        mainFrame.HLines = newGraph.HLines
        mainFrame.config = newGraph.config

        # check if line sets are used, if not, make them up
        if len(mainFrame.lines) > 0 and len(mainFrame.lines.lineSets) == 0:
            mainFrame.lines.initializeLineSets()


        mainFrame.drawLines()
        mainFrame.graph.config.update()
        mainFrame.setComments(mainFrame.graph.comments)
        mainFrame.SetTitle(misc.mainWindowTitle + ' - ' + 
                      os.path.basename(filename))

        # need to reset the home, back, and forward buttons
        mainFrame.navTb._views.clear()
        xmin, xmax, ymin, ymax = mainFrame.graph.config.getViewLimits()
        mainFrame.setViewLimits((xmin, xmax), (ymin, ymax))

        mainFrame.annotationManager.restore_positions()

        # populate line id's and genemarker id's, if they aren't there.
        # this happens if the .graph file was created before id's were used
        for eachLine in mainFrame.lines:
            if eachLine.id == 0:
                eachLine.id = mainFrame.lines.getNextId()

        for eachGeneMarkerGroup in mainFrame.geneMarkerGroups:
            if eachGeneMarkerGroup.id == 0:
                eachGeneMarkerGroup.id = geneMarkers.getNextGroupId(
                                                     mainFrame.geneMarkerGroups)
                i = 1
                for eachSet in eachGeneMarkerGroup:
                    eachSet.id = i
                    i += 1

