"""This module has classes that deal with gene markers. (The dots along a chromosome that signify different gene positions"""

from bisect import bisect_left
import copy
import os

import colors
from config import graphInfo

def getNextGroupId(geneMarkerGroups):
    """Given a list of gene marker groups, return a gene marker group id that
    is unique.
    
    """

    # basically run through all previous gene marker groups, then return an id
    # number that is one greater than the max found
    maxId = 0
    for group in geneMarkerGroups:
        maxId = max(group, maxId)

    return maxId + 1

def getNewGeneMarkerGroup(mainFrame):
    """Return a new gene marker group"""

    geneMarkersGroup = GeneMarkersGroup()
    geneMarkersGroup.id = getNextGroupId(mainFrame.geneMarkerGroups)
    return geneMarkersGroup

def getNewGeneMarkers(pos, line):
    """Given a the line associated with a GeneMarkers, and the marker
    positions, return a GeneMarkers object.
    
    """

    geneMarkers = GeneMarkers(pos, line)
    return geneMarkers

def loadGeneMarkerFile(mainFrame, info):
    """Given the main frame and a dict of parameters, load a file containing
    gene marker positions.
    
    """

    # try to open the file
    try:
        f = open(info['filename'], 'r')
    except IOError, e:
        s1="Error opening marker file '%s' : %s\n" % (filename,e.strerror)
        s2 = "Marker data not loaded.\n"
        message = s1 + s2
        dlg = wx.MessageDialog(self, message,"Error!", wx.OK|wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()
        return 

    # load marker information into a dictionary
    chr_col = info['chr_col']
    pos_col = info['pos_col']
    markers = dict()
    for eachLine in f.readlines():
        line = eachLine.split()
        if len(line) == 0 or line[0][0] == '#':
            # skip blank lines and comments
            continue
        if line[chr_col].lower() == 'x' or line[chr_col].lower() == 'y':
            # sex chromosomes
            num = 23
        else:
            num = int(line[chr_col])
        if markers.has_key(num):
            markers[num].append(float(line[pos_col]))
        else:
            markers[num] = [float(line[pos_col])]

    # add markers to lines
    geneMarkersGroup = getNewGeneMarkerGroup(mainFrame)
    for eachLine in mainFrame.lines:
        if markers.has_key(eachLine.number):
            pos = markers[eachLine.number]
            geneMarkers = getNewGeneMarkers(pos, eachLine)
            eachLine.addGeneMarkers(geneMarkers)
            geneMarkersGroup.addMarkers(geneMarkers)
    geneMarkersGroup.name = info['name']
    geneMarkersGroup.set_color(info['color'])
    geneMarkersGroup.set_marker(info['symbol'])
    geneMarkersGroup.set_markersize(info['size'])

    mainFrame.geneMarkerGroups.append(geneMarkersGroup)
    mainFrame.drawLines()
    f.close()

class GeneMarkers():
    """A class representing one set of gene markers. This represents gene 
    positions along *one* chromosome.
    
    """

    def __init__(self, pos=None, line=None, real=True):
        """line should be the line associated with the gene markers, and pos 
        should be a list of position values, group should be the
        geneMarkerGroup that this is associated with.
        
        """

        # the positions and ppl values of the gene markers
        if line is not None:
            self.pos = self.get_pos(line, pos)
            self.ppl = self.get_ppl(line)
        else:
            self.pos = pos
            self.ppl = None

        # the id number of the GeneMarkerGroup that this gene marker belongs to
        self.group_id = 0

        # a unique id that identifies the gene marker in the group
        self.id = 0

        # the actual matplotlib line class
        self.line = None

        # the color of the markers 
        self.color = colors.hexcolor(graphInfo.markers_color)

        # the symbol that is shown at gene positions
        self.symbol = graphInfo.markers_symbol

        # the size of the markers
        self.size = graphInfo.markers_size

        if real:
            self.set_default_color(line)

    def plot(self, axes, offset):
        """Given an axes and an offset in data space, plot the markers"""

        # plot positions in the same offset as the line, and make sure
        pos = [x+offset for x in self.pos]
        line = axes.plot(pos, self.ppl, color=self.color, linestyle='None',
                         marker=self.symbol, markersize=self.size)
        self.line = line[0]

    def get_pos(self, line, pos):
        """Given the line assoicated with the gene markers and a list of 
        position values, makes sure its in the right order.
        
        """

        # sort list, remove any duplicates, and remove any markers that are
        # outside the line boundaries
        u = {}
        for x in pos:
            u[x] = 1
        pos = u.keys()
        pos.sort()
        minVal = line.pos[0]
        maxVal = line.pos[-1]
        pos = [x for x in pos if x>=minVal and x<=maxVal]
        return pos

    def get_ppl(self, line):
        """Given the line associated with this gene markers, finds the correct
        ppl values for the gene positions.
        
        """

        ppl = []
        l = len(line.pos)
        maxpos = line.pos[-1]
        for eachPos in self.pos:
            # use binary search to find position of eachPos in pos list
            i = bisect_left(line.pos, eachPos)
            if i == l:
                # this position is outside the range
                continue
            elif line.pos[i] == eachPos:
                # hit an exact position
                ppl.append(line.ppl[i])
            else:
                # the exact position was not found, need to interpolate 
                # between the two closest points
                a = line.ppl[i-1]
                b = line.ppl[i]
                x1 = line.pos[i-1]
                x2 = line.pos[i]
                val = a + float(eachPos-x1) / (x2-x1) * (b-a)
                ppl.append(val)
        return ppl

    def set_color(self, color):
        """Set a color to the markers"""

        self.color = colors.hexcolor(color)
        if self.line is not None:
            self.line.set_color(color)

    def set_default_color(self, line):
        """Sets the gene markers to a default color. Right now it is the 
        compliment of the line color.
        
        """

        if line is None:
            return

        try:
            # see if self.color is a letter, like 'b'
            index = graphInfo.colors.index(line.color)
            color = colors.mpl_color(line.color)
        except ValueError, e:
            # line.color must be hex value
            color = colors.hex2rgb(line.color)
        finally:
            color = colors.getCompliment(color)
            self.set_color(colors.rgb2hex(color))

    def set_marker(self, marker):
        """Set the symbol of the markers"""

        self.symbol = marker
        if self.line is not None:
            self.line.set_marker(marker)

    def set_markersize(self, size):
        """Set the symbol size of the markers"""

        self.size = size
        if self.line is not None:
            self.line.set_markersize(size)

    def get_color(self):
        """Return the current color of the marker as an rgb tuple"""

        color = colors.hexcolor(self.color)
        return colors.hex2rgb(color)

    def setLine(self, line):
        """Given a line object, create ppl values that correspond to the line"""

        # this is used when the GeneMarker object is created without a line,
        # and the line is added later on. this happens when loading a .grx
        # file, for example

        self.pos = self.get_pos(line, self.pos)
        self.ppl = self.get_ppl(line)

    def backup(self):
        """Create a backup of values that might need to be restored later"""

        self.color_backup = self.color
        self.symbol_backup = self.symbol
        self.size_backup = self.size

    def restore(self):
        """Restore backup values"""

        try:
            self.color = self.color_backup
            self.symbol = self.symbol_backup
            self.size = self.size_backup

            # update any changes to values
            self.set_color(self.color)
            self.set_marker(self.symbol)
            self.set_markersize(self.size)
        except Exception, e:
            print 'warning! caught exception while in restore in Line class'
            print e
            pass

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = GeneMarkers(real=False).__dict__
        newKeys = set(newState.keys())
        curKeys = set(curState.keys())

        # class attributes that are in the new state but not the current state
        # need to be added
        for eachKey in newKeys - curKeys:
            curState[eachKey] = copy.copy(newState[eachKey])

        # class attributes that are in the current state but not the new state
        # need to be deleted
        for eachKey in curKeys - newKeys:
            del curState[eachKey]

        return curState

    def __len__(self):
        """Returns how many gene markers there are"""
        return len(self.pos)

    def __getstate__(self):
        """Returns the state used when pickling."""

        # remove all matplotlib references, since they cannot be pickled
        state = self.__dict__.copy()
        state['line'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        self.__dict__ = state

class GeneMarkersGroup():
    """This class represents a group of gene markers, i.e. all GeneMarkers that
    are loaded from one file. For example, if a file was loaded that had gene 
    markers for chr 1, 2, 3, then three GeneMarkers objects and one 
    GeneMarkersGroup object would be created. This class manages the color,
    size, and symbol of all gene markers in the group.
    
    """

    def __init__(self, real=True):
        # the name of the gene marker group
        self.name = ''

        # the list of gene markers that make up this group
        self.geneMarkers = []

        # the color of the markers 
        self.color = graphInfo.markers_color

        # the symbol that is shown at gene positions
        self.symbol = graphInfo.markers_symbol

        # the size of the markers
        self.size = graphInfo.markers_size

        # an id number to uniquely identify this group
        self.id = 0

        # a counter that assigns a unique id to each GeneMarker that is added
        # to this group. must set to to None because there could be a situation
        # where this object is created when loading a .grx file, then new
        # geneMarkers would be added later
        self.marker_id = None

    def addMarkers(self, markers):
        """Add a sets of markers to this group. Argument should be 
        GeneMarkers object.
        
        """

        markers.id = self.getNextId()
        markers.group_id = self.id
        self.geneMarkers.append(markers)

    def getNextId(self):
        """Return the next id to assign a GeneMarker and increment the counter"""

        if self.marker_id is None:
            # look through all the markers i have, and find the max id
            self.marker_id = 0
            for markers in self.geneMarkers:
                self.marker_id = max(self.marker_id, markers.id)
            self.marker_id += 1

        id = self.marker_id
        self.marker_id += 1
        return id

    def set_color(self, color):
        """Set a color to the gene markers"""

        if color == None:
            # this means just set your color to the color of the first
            # genemarker, and don't change any genemarker color. this is used
            # when the default color for genemarkers is used, and so the
            # genemarker's color is set elsewhere, and the genemarkers in the
            # group could be different colors, but just use the first one
            self.color = self.geneMarkers[0].color
        else:
            self.color = color
            for eachMarker in self.geneMarkers:
                eachMarker.set_color(color)

    def set_marker(self, marker):
        """Set the marker (symbol) of the markers"""

        self.symbol = marker
        for eachMarker in self.geneMarkers:
            eachMarker.set_marker(marker)

    def set_group_id(self):
        """Set the group id for all markers"""

        for eachMarker in self.geneMarkers:
            eachMarker.group_id = self.id

    def set_markersize(self, size):
        """Set the symbol size of the markers"""

        size = int(size)
        self.size = size
        for eachMarker in self.geneMarkers:
            eachMarker.set_markersize(size)

    def update(self):
        """Apply all configuration values to the gene markers"""

        # this is used to make sure all gene markers adhere to the group
        # values, expecially when loading from a .grx file
        self.set_color(self.color)
        self.set_marker(self.symbol)
        self.set_markersize(self.size)
        self.set_group_id()

    def backup(self):
        """Create a backup of values that might need to be restored later"""

        self.name_backup = self.name
        self.color_backup = self.color
        self.symbol_backup = self.symbol
        self.size_backup = self.size
        for eachMarker in self.geneMarkers:
            eachMarker.backup()

    def restore(self):
        """Restore backup values"""

        try:
            self.name = self.name_backup
            self.color = self.color_backup
            self.symbol = self.symbol_backup
            self.size = self.size_backup
            for eachMarker in self.geneMarkers:
                eachMarker.restore()
        except Exception, e:
            print 'warning! caught exception while in restore',
            print 'in GeneMarkersGroup'
            print e

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = GeneMarkersGroup(real=False).__dict__
        newKeys = set(newState.keys())
        curKeys = set(curState.keys())

        # class attributes that are in the new state but not the current state
        # need to be added
        for eachKey in newKeys - curKeys:
            curState[eachKey] = copy.copy(newState[eachKey])

        # class attributes that are in the current state but not the new state
        # need to be deleted
        for eachKey in curKeys - newKeys:
            del curState[eachKey]

        return curState

    def __iter__(self):
        """Used so this object can be iterated through with a for construct"""

        for eachMarker in self.geneMarkers:
            yield eachMarker

    def __len__(self):
        """Returns how many items are in the group"""
        return len(self.geneMarkers)

    def __getitem__(self, key):
        """This is used so the class can be accessed like a list"""
        return self.geneMarkers[key]

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        self.__dict__ = state
