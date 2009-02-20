"""The LineSet class represents a group of lines that are associated together. The lines have the same color, linestyle, marker, etc."""

import copy

from colors import hexcolor
import config.graphInfo as graphInfo
import panel.mainFrame

class LineSet(object):
    """A class representing a group of lines that are associated together. The
    lines have the same color, linestyle, marker, etc.
    
    """


    def __init__(self, id=0, name='', color='b', thickness=2, linestyle='-', 
                 marker='None', markersize=6, lines=None, real=True):

        # a unique id number for this set
        self.id = id

        # the name of the set. this name is shown when configuring line
        # set properties, and also show up when the legend shows line sets
        self.name = name

        # the color this line should be drawn in
        self.color = hexcolor(color)

        # the thickness of the line
        self.thickness = thickness

        # the linestyle (solid, dashed, dotted, etc.)
        self.linestyle = linestyle

        # the symbol that is shown at data points
        self.marker = marker

        # the size of the symbols for markers for the data points
        self.markersize = markersize

        # the lines in this set (list of Line objects)
        if lines is None:
            self.lines = []
        else:
            self.lines = lines

        if len(self.lines) > 0:
            # correspond the class variables with the variables of a line
            self.setValues(lines[0])
            self.setLines()

        # backup values that are used in the backup and restore functions
        self.name_backup = self.name
        self.color_backup = self.color
        self.thickness_backup = self.thickness
        self.linestyle_backup = self.linestyle
        self.marker_backup = self.marker
        self.markersize_backup = self.markersize

    def getRepresentativeLine(self):
        """Return a line in the set that has the same properties as the set"""

        if len(self.lines) == 0:
            return None

        # need to check each line, and return the first line in the set that has
        # the same properties as the set. need to do this because the user
        # can change individual line properties, and it might be possible that
        # the legend could look wrong if that particular line was changed.
        # count how many parameters each line shares with the set
        # Note that this returns the matplotlib line, not Line object.
        count = []
        for eachLine in self.lines:
            i = 0
            if eachLine.color != self.color:
                i += 1
            if eachLine.thickness != self.thickness:
                i += 1
            if eachLine.linestyle != self.linestyle:
                i += 1
            if eachLine.marker != self.marker:
                i += 1
            if eachLine.markersize != self.markersize:
                i += 1

            # line must be exactly the same
            if i == 5:
                return line.line
            else:
                count.append(i)

        # if execution hits here, that means no line exactly matches, find
        # line that matches most closely
        max_i = max(count)
        i = count.index(max_i)
        return self.lines[i].line

    def setValues(self, line):
        """Given a Line object, set all parameters to be the same as the line"""

        self.color = line.color
        self.thickness = line.thickness
        self.linestyle = line.linestyle
        self.marker = line.marker
        self.markersize = line.markersize

    def setLines(self):
        """For all the lines in the set, change their lineSet_id to this 
        set's id.
        
        """

        for eachLine in self.lines:
            eachLine.lineSet_id = self.id

    def set_color(self, color):
        """Set a color to the line"""

        color = hexcolor(color)
        self.color = color
        for eachLine in self.lines:
            eachLine.set_color(color)

    def set_linewidth(self, width):
        """Set the thickness of a line"""

        self.thickness = width
        for eachLine in self.lines:
            eachLine.set_linewidth(width)

    def set_linestyle(self, linestyle):
        """Set the linestyle of the line"""

        self.linestyle = linestyle
        for eachLine in self.lines:
            eachLine.set_linestyle(linestyle)

    def set_marker(self, marker):
        """Set the marker of the line"""

        self.marker = marker
        for eachLine in self.lines:
            eachLine.set_marker(marker)

    def set_markersize(self, markersize):
        """Set the markersize of the line"""

        self.markersize = markersize
        for eachLine in self.lines:
            eachLine.set_markersize(markersize)

    def backup(self):
        """Create a backup of values that might need to be restored later"""

        # this function is called when configuring line set options
        self.name_backup = self.name
        self.color_backup = self.color
        self.thickness_backup = self.thickness
        self.linestyle_backup = self.linestyle
        self.marker_backup = self.marker
        self.markersize_backup = self.markersize

    def restore(self):
        """Restore backup values"""

        self.name = self.name_backup
        self.set_color(self.color_backup)
        self.set_linewidth(self.thickness_backup)
        self.set_linestyle(self.linestyle_backup)
        self.set_marker(self.marker_backup)

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = LineSet(real=False).__dict__
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

    def __setstate__(self, state):
        """Sets the state of the object when unpickled"""

        state = self.updateSelf(state)
        self.__dict__ = state
