"""The LineGroup class represents all lines in one vertical section of the graph"""

import copy

import config.graphInfo as graphInfo

class LineGroup():
    """This class represents a group of lines, i.e. all lines that overlap each
    other and share a partitioned space on the graph"""

    def __init__(self, line=None, real=True):
        # the list of lines that are part of the same group
        self.lines = []

        # the offset amount that the group is at in data space
        self.offset = 0

        # parameters of the group in data space
        self.xmin = 0
        self.xmax = 0
        self.length = 0

        # these are the x-values at the beginning and end of the group relative
        # to the chromosomes. for example, if the only line in the group had 
        # positions from 50 to 200, then start = 50 and end = 200.
        self.start = 0
        self.end = 0

        # the string used to label the group in the graph, shows up at the
        # bottom of the graph
        self.label = ''

        # a matplotlib text object that is show underneath the graph, usually
        # denoting what number chromosomes are in the group
        self.text = None

        # a flag to determine if duplicate numbers in the label are allowed.
        # for example, the label, '1, 1, 1' becomes '1'
        self.allowDuplicates = graphInfo.allowDuplicates

        if line is not None:
            self.lines.append(line)
            self.label += line.name
            self.setBoundaries()

    def plot(self, offset):
        """Given an offset in data space, plot all lines in the group"""

        # plot all lines such that the line with the minimum starting
        # position is flush left to the group, and all other lines
        # are drawn relative to that line.

        self.offset = offset
        minPos = min(eachLines.pos[0] for eachLines in self.lines)
        for eachLine in self.lines:
            eachLine.plot(offset - minPos)

        # since given a new offset, figure out the rest 
        # of the boundary parameteres
        self.setBoundaries()

    def setBoundaries(self):
        """Calculate the group boundary parameters"""

        self.start = min(eachLine.pos[0] for eachLine in self.lines)
        self.end = max(eachLine.pos[-1] for eachLine in self.lines)
        self.length = self.end - self.start

        self.xmin = self.offset
        self.xmax = self.xmin + self.length

    def append(self, line):
        """Add a line to the end of the group"""

        self.lines.append(line)
        self.setBoundaries()
        self.redoLabel()

    def remove(self, line):
        """Remove a line from the group """

        self.lines.remove(line)
        self.redoLabel()

    def redoLabel(self):
        """Re-creates the label"""

        self.label = ''
        oldNames = []
        if self.allowDuplicates:
            for eachLine in self.lines:
                self.label += eachLine.name + ', '
        else:
            for eachLine in self.lines:
                if not eachLine.name in oldNames:
                    self.label += eachLine.name + ', '
                    oldNames.append(eachLine.name)

        self.label = self.label[:-2]

    def isWithin(self, x):
        """Given a point on the x-axis, determines if this group occupies
        that point
        
        """

        for eachLine in self.lines:
            if eachLine.isWithin(x):
                return True

        return False

    def doesOverlap(self, xmin, xmax):
        """Given an x range, determines if this groups overlaps it"""

        if self.xmin >= xmin and self.xmax <= xmax:
            return True
        elif self.xmin <= xmin and self.xmax >= xmax:
            return True
        elif self.xmin >= xmin and self.xmin <= xmax:
            return True
        elif self.xmax >= xmin and self.xmax <= xmax:
            return True

        return False

    def hasLineName(self, name):
        """Returns true if the group has a line with the same name as name"""

        for eachLine in self.lines:
            if eachLine.name == name:
                return True

        return False

    def getLine(self, name):
        """Returns a line with the given name, or None if no such line exists"""

        for eachLine in self.lines:
            if eachLine.name == name:
                return eachLine

        return None

    def getNames(self):
        """Return all names of all lines in the list"""

        return [eachLine.name for eachLine in self.lines]

    def get_mpl_lines(self):
        """Returns all the matplotlib lines in the group"""
        return [x.line for x in self.lines]

    def getMaxDataPoint(self, xmin, xmax, ymin=None, ymax=None):
        """Given an x range and an optional y range, find the 
        maximum data point lying in that range.
        
        """

        def isWithin(x, y, xmin, xmax, ymin, ymax):
            """Given a point and a range, determine whether the point
            lies inside those ranges
            
            """

            if x >= xmin and x <= xmax:
                if ymin is None or ymax is None:
                    return True
                if y >= ymin and y <= ymax:
                    return True
            return False

        # assume that xmin and xmax take the offset into account
        maxPoint = None
        for eachLine in self.lines:
            for eachPoint, eachPPL in zip(eachLine.pos, eachLine.ppl):
                eachPoint += eachLine.offset
                if isWithin(eachPoint, eachPPL, xmin, xmax, ymin, ymax):
                    if maxPoint is None or eachPPL > maxPoint[1]:
                        maxPoint = (eachPoint, eachPPL)
                if eachPoint > xmax:
                    break;

        return maxPoint

    def getMinDataPoint(self, xmin, xmax, ymin=None, ymax=None):
        """Given an x range and an optional y range, find the 
        mininum data point lying in that range.
        
        """

        def isWithin(x, y, xmin, xmax, ymin, ymax):
            """Given a point and a range, determine whether the point
            lies inside those ranges
            
            """

            if x >= xmin and x <= xmax:
                if ymin is None or ymax is None:
                    return True
                if y >= ymin and y <= ymax:
                    return True
            return False

        # assume that xmin and xmax take the offset into account
        minPoint = None
        for eachLine in self.lines:
            for eachPoint, eachPPL in zip(eachLine.pos, eachLine.ppl):
                eachPoint += eachLine.offset
                if isWithin(eachPoint, eachPPL, xmin, xmax, ymin, ymax):
                    if minPoint is None or eachPPL < minPoint[1]:
                        minPoint = (eachPoint, eachPPL)
                if eachPoint > xmax:
                    break

        return minPoint

    def getLongestLine(self):
        """Returns the line with the longest length in the group"""

        def length(line):
            """Given a line, return its length"""
            return line.pos[-1] - line.pos[0]

        line = None
        for eachLine in self.lines:
            if line is None or length(line) < length(eachLine):
                line = eachLine
        return line

    def backup(self):
        """Create a backup of values that might need to be restored later"""
        self.allowDuplicates_backup = self.allowDuplicates

    def restore(self):
        """Restore backup values"""

        try:
            self.allowDuplicates = self.allowDuplicates_backup
        except:
            print 'warning! caught exception while in restore in LineGroup'
            pass

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = LineGroup(real=False).__dict__
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

        for eachLine in self.lines:
            yield eachLine

    def __len__(self):
        """Returns how many items are in the group"""
        return len(self.lines)

    def __getitem__(self, key):
        """This is used so the class can be accessed like a list"""
        return self.lines[key]

    def __getstate__(self):
        """Returns the state used when pickling."""

        # remove all matplotlib references, since they cannot be pickled
        state = self.__dict__.copy()
        state['text'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        self.__dict__ = state
