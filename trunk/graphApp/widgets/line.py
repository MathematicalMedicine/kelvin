"""The Line class represents one line in the graph"""

import copy

from colors import hexcolor
import config.graphInfo as graphInfo
import panel.mainFrame

class Line():
    """A class representing one line (i.e. one chromosome)"""

    def __init__(self, axes=None, color=graphInfo.colors[0], name=None, num=0, 
                 real=True):

        # the axes that this line will be drawn on
        self.axes = axes

        # an id number that uniquely identifies this line. used primarily
        # to reference the line in the .grx file format
        self.id = 0

        # the id number of the line set that this line belongs to
        self.lineSet_id = 0

        # the actual matplotlib line class
        self.line = None

        # the chromosome number (1-23)
        self.number = num

        # the name given to the chromosome
        if name == None:
            # use the number as the name
            self.name = str(num)
        else:
            self.name = str(name)

        # a list of all positions (should be float)
        self.pos = []

        # a list of all ppl values (should be float)
        self.ppl = []

        # a 2D list of data that is associated with each data point
        self.metadata = []

        # the color this line should be drawn in
        self.color = hexcolor(color)

        # the thickness of the line
        self.thickness = 2

        # the linestyle (solid, dashed, dotted, etc.)
        self.linestyle = '-'

        # since each position of every chromosome starts at 0, but you
        # want the lines drawn next to each other instead of overlapping
        # each other, need to add an offset to each position value
        self.offset = 0

        # the symbol that is shown at data points
        self.marker = 'None'

        # the size of the symbols for markers for the data points
        self.markersize = 6

        # a list of gene markers that is assoicated with this line
        self.geneMarkers = []

        # backup values that are used in the backup and restore functions
        self.name_backup = self.name
        self.color_backup = self.color
        self.thickness_backup = self.thickness
        self.linestyle_backup = self.linestyle
        self.marker_backup = self.marker
        self.markersize_backup = self.markersize

    def plot(self, offset):
        """Given an offset in data space, plot the line on the axes"""

        self.offset = offset
        pos = [x+offset for x in self.pos]
        label = self.name
        self.line = self.axes.plot(pos, self.ppl, label=label, color=self.color,
                                   linewidth=self.thickness, 
                                   linestyle=self.linestyle, marker=self.marker,
                                   markersize=self.markersize)
        self.line = self.line[0]

        # plot any gene markers also
        self.plot_markers()

    def plot_markers(self):
        """Plot the markers on the line"""

        for eachMarker in self.geneMarkers:
            eachMarker.plot(self.axes, self.offset)

    def addMetadata(self, data):
        """Add some metadata to the line, data is a list"""
        self.metadata.append(data)

    def addGeneMarkers(self, markers):
        """Add information about gene markers. Argument markers should be a
        GeneMarkers object. 
        
        """

        self.geneMarkers.append(markers)

    def isWithin(self, x):
        """Given a point on the x-axes, determines if this line occupies
        that point in the x-axes when the offset is taken into account
        
        """

        if self.getStart() <= x and self.getEnd() >= x:
            return True
        else:
            return False

    def getStart(self):
        """Returns the actual starting point of the line in data space"""
        return self.pos[0] + self.offset

    def getEnd(self):
        """Returns the actual ending point of the line in data space"""
        return self.pos[-1] + self.offset

    def set_color(self, color):
        """Set a color to the line"""

        color = hexcolor(color)
        self.color = color
        if self.line is not None:
            self.line.set_color(color)

    def set_linewidth(self, width):
        """Set the thickness of a line"""

        self.thickness = width
        if self.line is not None:
            self.line.set_linewidth(width)

    def set_linestyle(self, linestyle):
        """Set the linestyle of the line"""

        self.linestyle = linestyle
        if self.line is not None:
            self.line.set_linestyle(linestyle)

    def set_marker(self, marker):
        """Set the marker of the line"""

        self.marker = marker
        if self.line is not None:
            self.line.set_marker(marker)

    def set_markersize(self, markersize):
        """Set the markersize of the line"""

        self.markersize = markersize
        if self.line is not None:
            self.line.set_markersize(markersize)

    def backup(self):
        """Create a backup of values that might need to be restored later"""

        # this function is called when configuring graph options
        self.name_backup = self.name
        self.color_backup = self.color
        self.thickness_backup = self.thickness
        self.linestyle_backup = self.linestyle
        self.marker_backup = self.marker
        self.markersize_backup = self.markersize

        for eachMarker in self.geneMarkers:
            eachMarker.backup()

    def restore(self):
        """Restore backup values"""

        try:
            self.name = self.name_backup
            self.color = self.color_backup
            self.thickness = self.thickness_backup
            self.linestyle = self.linestyle_backup
            self.marker = self.marker_backup
            self.markersize = self.markersize_backup

            for eachMarker in self.geneMarkers:
                eachMarker.restore()

            # update any changes to values
            self.set_color(self.color)
            self.set_linewidth(self.thickness)
            self.set_linestyle(self.linestyle)
            self.set_marker(self.marker)
        except Exception, e:
            print 'warning! caught exception while in restore in Line class:', e
            pass

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = Line(real=False).__dict__
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
        """Returns how many values are in the line"""
        return len(self.ppl)

    def __getstate__(self):
        """Returns the state used when pickling."""

        # remove all matplotlib references, since they cannot be pickled
        state = self.__dict__.copy()
        state['axes'] = None
        state['line'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        state['axes'] = panel.mainFrame.getMainFrame().axes
        self.__dict__ = state
