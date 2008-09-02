"""LineCollection is a manager of all line-related things for the graph"""

import copy
from math import sqrt
import types

from widgets.line import Line
from widgets.lineGroup import LineGroup
from widgets.lineSet import LineSet

class LineCollection():
    """This class represents a collection of lines, i.e. all lines shown in the 
    graph.  Internally the lines are stored as a list of groups of lines.
    
    """

    def __init__(self, real=True):
        # the list of groups of lines
        self.groups = []

        # list of lineSet objects
        self.lineSets = []

        # a counter that assigns the id number to a line
        self.nextId = None

        # a counter that assigns the id number to a line set
        self.nextLineSetId = None

    def getNextId(self):
        """Return the next id number to be used. Increment the counter."""

        # this happens when the first line is added, and also happens
        # when a grx file is loaded and then a line is added
        if self.nextId is None:
            # look through all the lines i have, and find the max id
            self.nextId = 0
            for eachLine in self:
                self.nextId = max(self.nextId, eachLine.id)
            self.nextId += 1

        id = self.nextId
        self.nextId += 1
        return id

    def getNextLineSetId(self):
        """Return the next id number to be used for a line set. Increment the 
        counter.
        
        """

        # this happens when the first line set is created, and also happens
        # when a grx file is loaded and then another line set is created 
        if self.nextLineSetId is None:
            # look through all the line sets i have, and find the max id
            self.nextLineSetId = 0
            for eachSet in self.lineSets:
                self.nextLineSetId = max(self.nextLineSetId, eachSet.id)
            self.nextLineSetId += 1

        id = self.nextLineSetId
        self.nextLineSetId += 1
        return id

    def addLineSet(self, lineSet, newId=True):
        """Add an existing line set to this line collection. The newId arg
        indicates whether to assign the lineSet a new id or not. A new one is
        not needed when loading a .grx file, for example.
        
        """

        if newId:
            lineSet.id = self.getNextLineSetId()
        self.lineSets.append(lineSet)

    def addNewLineSet(self, name='', lines=None):
        """Create a new line set, add it to this line collection"""

        id = self.getNextLineSetId()
        lineSet = LineSet(id=id, name=name, lines=lines)
        self.addLineSet(lineSet, newId=False)

        return lineSet

    def initializeLineSets(self):
        """Group the lines in the collection into the appropriate line sets. If
        no line sets are present, create sets and try to partition the lines as
        logically as possible. This is called when loading a .grx file.
        
        """

        # Note: this function assumes all line sets and all lines have been
        # already loaded into the line collection

        lines = self.getLines()
        if len(lines) > 0 and len(self.lineSets) == 0:
            # no line sets are present, but there are some lines. this means
            # and older .grx file is being opened, when line sets were not
            # used. try to group lines into sets as best as possible
            self.partitionLinesIntoSets()
        else:
            # line sets are being used. build up their list of lines
            for eachLine in lines:
                for eachSet in self.lineSets:
                    if eachLine.lineSet_id == eachSet.id:
                        eachSet.lines.append(eachLine)

    def partitionLinesIntoSets(self):
        """Partition the lines into logical sets, and create line sets out of
        them.
        
        """

        # this is used for backward compatibility, since line sets were not
        # used for a long time. bascially the way this works is that it will
        # group lines that have the same color, thickness, line style, marker,
        # and marker size into the same line set. if two lines have the same
        # attributes and the same number, then they are placed in different
        # sets, as it is unlikely two lines with the same number would be in
        # the same line set.

        lines = self.getLines()
        if len(lines) == 0:
            # nothing to do
            return

        def sameLine(line1, line2):
            # returns true if both lines have same attributes
            if line1.color != line2.color:
                return False
            if line1.thickness != line2.thickness:
                return False
            if line1.linestyle != line2.linestyle:
                return False
            if line1.marker != line2.marker:
                return False
            if line1.markersize != line2.markersize:
                return False
            return True

        def isPresent(line1, set):
            # returns True if a line with the same number as line1 is present
            # in set, which is a list of lines
            for eachLine in set:
                if eachLine.number == line1.number:
                    return True
            return False

        self.lineSets = []
        sets = []
        while len(lines) > 0:
            curLine = lines.pop(0)
            curColor = curLine.color
            curThickness = curLine.thickness
            curLineStyle = curLine.linestyle
            curMarker = curLine.marker
            curMarkerSize = curLine.markersize
            set = [curLine]
            for eachLine in lines[:]:
                if sameLine(curLine, eachLine) and not isPresent(eachLine, set):
                    set.append(eachLine)
                    lines.remove(eachLine)
            sets.append(set)

        # create line sets based on the sets found earlier
        for eachSet in sets:
            lineSet = self.addNewLineSet(lines=eachSet)

    def addLine(self, line, newId=True):
        """Add a line to the collection. The line argument is a Line object.
        the newId argument specifies if the line should be assigned a new id
        number. You don't want to assign new id numbers when lines are being
        loaded from a .grx file, for example.
        
        """

        # make a new group for the new line
        group = LineGroup()
        self.groups.append(group)
        self.addLineToCollection(line, group, newId)

    def addOverlappingLine(self, newLine, overlapLine, newId=True):
        """Add a line that is to be overlapped with an existing line.  The 
        newLine argument is a Line object, and the overlapLine argument is
        the index or the name of the line to be overlapped with. See addLine
        for an explanation of newId.
        
        """

        line, group = self.getLineAndGroup(overlapLine)
        self.addLineToCollection(newLine, group, newId)

    def addLineToCollection(self, line, group, newId=True):
        """Given a line object and the group, add a line to the collection.
        This should always be called when a line is added to the graph.
        
        """

        if newId:
            line.id = self.getNextId()
        group.append(line)

    def overlap(self, source, dest):
        """Takes the line(s) specified in source and moves them to the same 
        group as the dest line.  The source argument is a list of line indices, 
        and dest is a line index.
        
        """

        # the group the source lines will be moved to
        line, destGroup = self.getLineAndGroup(dest)

        # this will be a list of tuples with format (<line>, <group>)
        sourceLines = []

        # since moving lines around will change the index of them, first
        # get all the source lines and put them in a list
        for i in source:
            sourceLines.append(self.getLineAndGroup(i))

        for eachEntry in sourceLines:
            # remove the line from its current group
            eachEntry[1].remove(eachEntry[0])
            if len(eachEntry[1]) == 0:
                self.groups.remove(eachEntry[1])

            # add line to destination group
            destGroup.append(eachEntry[0])

    def removeLine(self, line):
        """Remove the line from the collection with the given index or name"""

        line, group = self.getLineAndGroup(line)
        group.remove(line)
        if len(group) == 0:
            self.groups.remove(group)

    def getGroups(self):
        """Returns the list of groups"""
        return self.groups

    def getLine(self, line):
        """Returns the line with the given index or name"""
        return self.getLineAndGroup(line)[0]

    def getLines(self):
        """Return all lines in this collection"""

        lines = []
        for eachGroup in self.groups:
            for eachLine in eachGroup:
                lines.append(eachLine)
        return lines

    def getLineById(self, id):
        """Return the line with the given line id"""

        for eachLine in self:
            if eachLine.id == id:
                return eachLine

        return None

    def getGroup(self, line):
        """Returns the group given a line index or name"""
        return self.getLineAndGroup(line)[1]

    def getGroupByX(self, x):
        """Returns the first group that occupies the value on the x-axis"""

        for eachGroup in self.groups:
            if eachGroup.isWithin(x):
                return eachGroup

        # point is not in any group
        return None

    def getLineAndGroup(self, line):
        """Returns the line with the given index or name or Line reference, 
        along with the group its in.
        
        """

        if type(line) == types.IntType:
            # treat the first line in the first group as index 0, the second
            # line in the first group as index 1, etc.
            count = 0
            for eachGroup in self.groups:
                for eachLine in eachGroup:
                    if count == line:
                        return eachLine, eachGroup
                    count += 1

        if type(line) in types.StringTypes:
            # return the first line encountered with the same name
            for eachGroup in self.groups:
                for eachLine in eachGroup:
                    if line == eachLine.name:
                        return eachLine, eachGroup
        
        if isinstance(line, Line):
            for eachGroup in self.groups:
                for eachLine in eachGroup:
                    if line is eachLine:
                        return eachLine, eachGroup

        # line is not found, return None
        return None, None

    def getNames(self):
        """Returns a list of strings of all the line names."""

        names = []
        for eachGroup in self.groups:
            for eachLine in eachGroup:
                names.append(eachLine.name)

        return names

    def get_x_min(self):
        """Returns the minimum x-value that a line starts, also taking into 
        account the offset.
        
        """

        if len(self.groups) < 1:
            return 0

        return min(eachLine.getStart() for eachLine in self.groups[0])

    def get_x_max(self):
        """Returns the maximum x-value of the lines in the last group, 
        also taking into account the offset.
        
        """

        if len(self.groups) < 1:
            return 0

        return max(eachLine.getEnd() for eachLine in self.groups[-1])

    def getMaxDataPoint(self, xmin, xmax, ymin=None, ymax=None):
        """Given a range of x-values and an optional range of y-values, 
        find the maximum data point that lies within that range
        
        """

        maxPoint = None
        for eachGroup in self.groups:
            if eachGroup.doesOverlap(xmin, xmax):
                groupMax = eachGroup.getMaxDataPoint(xmin, xmax, ymin, ymax)
                if groupMax is not None:
                    if maxPoint is None or maxPoint[1] < groupMax[1]:
                        maxPoint = groupMax

        return maxPoint

    def getMinDataPoint(self, xmin, xmax, ymin=None, ymax=None):
        """Given a range of x-values and an optional range of y-values, 
        find the minimum data point that lies within that range
        
        """

        minPoint = None
        for eachGroup in self.groups:
            if eachGroup.doesOverlap(xmin, xmax):
                groupMin = eachGroup.getMinDataPoint(xmin, xmax, ymin, ymax)
                if groupMin is not None:
                    if minPoint is None or groupMin[1] < minPoint[1]:
                        minPoint = groupMin

        return minPoint

    def getPeak(self, xmin, xmax, ymin, ymax):
        """Given a rectangle, find the highest peak in that rectangle"""
        return self.getPeakOrTrough(xmin, xmax, ymin, ymax, True)

    def getTrough(self, xmin, xmax, ymin, ymax):
        """Given a rectangle, find the lowest trough in that rectangle"""
        return self.getPeakOrTrough(xmin, xmax, ymin, ymax, False)

    def getPeakOrTrough(self, xmin, xmax, ymin, ymax, find_peaks):
        """Given a rectangle and a flag indicating peaks or troughs,
        find either the highest peak or lowest trough in the rectangle.  
        Note that all points are in data space.
        
        """

        def isInside(point, xmin, xmax, ymin, ymax):
            """Returns true if the value x is between min and max"""
            x, y = point
            if x >= xmin and x <= xmax and y >= ymin and y <= ymax:
                return True
            else:
                return False

        def checkPoint(newLine, i, point, oldLine):
            """Checks if the i-th point in the new line is better 
            than the current point
            
            """

            x = newLine.pos[i] + newLine.offset
            y = newLine.ppl[i]
            if point is None or compare(y, point[1]):
                return (x, y), newLine
            else:
                return point, oldLine

        def compare_peaks(a, b):
            return a >= b

        def compare_troughs(a, b):
            return a <= b

        if find_peaks:
            compare = compare_peaks
        else:
            compare = compare_troughs

        point = None
        line = None
        for eachGroup in self.groups:
            if eachGroup.doesOverlap(xmin, xmax):
                for eachLine in eachGroup.lines:
                    length = len(eachLine)
                    for i in range(0, length):
                        x = eachLine.pos[i] + eachLine.offset
                        y = eachLine.ppl[i]
                        if isInside((x, y), xmin, xmax, ymin, ymax):
                            # case when it's the first point of the line
                            if i == 0 and compare(eachLine.ppl[0], 
                                                  eachLine.ppl[1]):
                                point, line=checkPoint(eachLine, i, point, line)
                            # case when it's the last point in the line
                            elif i == length-1 and compare(eachLine.ppl[-1], 
                                                           eachLine.ppl[-2]):
                                point, line=checkPoint(eachLine, i, point, line)
                            # case when point is in the middle
                            elif compare(eachLine.ppl[i], eachLine.ppl[i-1]) \
                               and compare(eachLine.ppl[i], eachLine.ppl[i+1]):
                                point, line=checkPoint(eachLine, i, point, line)
        return point, line

    def getAllPeaksOrTroughs(self, threshold, find_peaks):
        """Given a threshold value and a flag indicating peaks or troughs,
        find all points in the graph that fit the criteria.  Note that all
        points are in data space. Points along with their corresponding lines
        are returned.
        
        """

        def addPoint(line, i):
            """Add a point to the result that will be returned"""
            x = line.pos[i] + line.offset
            y = line.ppl[i]
            points.append((x, y))
            lines.append(line)

        def compare_peaks(a, b):
            return a >= b

        def compare_troughs(a, b):
            return a <= b

        if find_peaks:
            compare = compare_peaks
        else:
            compare = compare_troughs

        points = []
        lines = []
        for eachLine in self:
            length = len(eachLine.pos)
            for i in range(0, len(eachLine)):
                if compare(eachLine.ppl[i], threshold):
                    # case when it's the first point of the line
                    if i == 0 and compare(eachLine.ppl[0], eachLine.ppl[1]):
                        addPoint(eachLine, i)
                    # case when it's the last point in the line
                    elif i == length-1 and compare(eachLine.ppl[-1], 
                                                   eachLine.ppl[-2]):
                        addPoint(eachLine, i)
                    # case when point is in the middle
                    elif compare(eachLine.ppl[i], eachLine.ppl[i-1]) and \
                         compare(eachLine.ppl[i], eachLine.ppl[i+1]):
                        addPoint(eachLine, i)

        return points, lines

    def findClosestMetada(self, x, y):
        """Given a point, find the closest data point to that point and
        return that point's metadata
        
        """

        # check if there is actually some metadata to show
        if len(self) == 0 or len(self.groups[0].lines[0].metadata) == 0:
            return []

        def distance(p1, p2):
            """Return distance between two points"""

            # transform everything to pixel space
            xd = p2[0] - p1[0]
            yd = p2[1] - p1[1]
            return sqrt(xd**2 + yd**2)

        # first find the closest point to the data point
        # info keeps track of the closest one so far
        # format is <line> <index> <distance>
        info = None
        for eachGroup in self.groups:
            if eachGroup.isWithin(x):
                for eachLine in eachGroup.lines:
                    for i in range(0, len(eachLine)):
                        x2 = eachLine.pos[i] + eachLine.offset
                        y2 = eachLine.ppl[i]
                        dist = distance((x, y), (x2, y2))
                        if info is None or dist < info[2]:
                            info = (eachLine, i, dist)

        line = info[0]
        i = info[1]
        return line.metadata[i]

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = LineCollection(real=False).__dict__
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

    def changeLineName(self, line=None, group=None, newName=None):
        """Change the name of a line"""

        if line is None:
            # this indicates that any line name might have been changed
            for eachGroup in self.groups:
                eachGroup.redoLabel()
        else:
            # only updating one line
            line.name = newName
            group.redoLabel()

    def __iter__(self):
        """Used so this object can be iterated through with a for construct"""

        for eachGroup in self.groups:
            for eachLine in eachGroup:
                yield eachLine

    def __len__(self):
        """Returns how many lines are in the collection"""
        return sum(len(eachGroup) for eachGroup in self.groups)

    def __getitem__(self, key):
        """This is used so the class can be accessed like a list"""

        return self.groups[key]

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        self.__dict__ = state
