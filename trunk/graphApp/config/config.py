"""Holds configurable information about the graph"""

import copy

import graphInfo
import panel.mainFrame

class Config():

    def __init__(self, mainFrame=None, real=True):
        """Holds configuration data for a figure"""

        self.mainFrame = mainFrame
        if real:
            self.figure = mainFrame.figure
            self.axes = self.figure.gca()
        else:
            self.figure = None
            self.axes = None

        self.title = ''
        self.xlabel = ''
        self.ylabel = ''

        if real:
            self.bgcolor = self.axes.get_axis_bgcolor()
        else:
            self.bgcolor = 'w'
        self.showLegend = False
        self.showLegendFrame = True
        self.legendPosition = graphInfo.legendPositions[0]

        # this flag indicates what to show in the legend. the legend can either
        # show every individual line, or simply show the different line sets
        self.showLineSetsInLegend = False

        # the text size of the title and axis labels
        self.titleSize = 14
        self.axisSize = 12

        # the viewing limits of the graph
        self.xmin = 0
        self.xmax = 0
        self.ymin = 0
        self.ymax = 0

        # the precision of the annotation text refering to ppl.
        # if the flag defaultAnnPPLPrec is True, default precision is used.
        # with default precision, two decimal places are shown, except
        # when ppl <.02, in which case three decimal places are used.
        # if the flag is False, the user specified annPPLPrec is used.
        # it needs to be at least 2 because a percent sign is always shown.
        # affects annotation text and the sidebar
        self.defaultAnnPPLPrec = True
        self.annPPLPrec = 2

        # the precision of the annotation text refering to position.
        # shown in the annotation text and the sidebar
        self.annPosPrec = 0

        # flag to determine whether to display all annotations as same color
        # or to have annotations the same color as their lines
        self.isAnnSameColor = True

        # flag to determine whether to allow duplicates in the line group labels
        # for example, '1, 1, 1' becomes '1' if allowDuplicates is false
        self.allowDuplicates = graphInfo.allowDuplicates

    def update(self):
        """Applies all the values to the figure and axes"""

        self.axes.set_axis_bgcolor(self.bgcolor)
        self.axes.set_title(self.title, size=self.titleSize)
        self.mainFrame.xlabel.set_text(self.xlabel)
        self.mainFrame.xlabel.set_size(self.axisSize)
        self.axes.set_ylabel(self.ylabel, size=self.axisSize)
        for eachGroup in self.mainFrame.lines.getGroups():
            eachGroup.text.set_size(self.axisSize)

    def updateValues(self):
        """Updates the values in config"""

        # these variables only need to be collected when needed
        self.xmin, self.xmax = self.axes.get_xlim()
        self.ymin, self.ymax = self.axes.get_ylim()

    def getViewLimits(self):
        """Returns the view limits in a tuple"""
        return (self.xmin, self.xmax, self.ymin, self.ymax)

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = Config(real=False).__dict__
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

    def __getstate__(self):
        """Returns the state used when pickling."""

        # remove all matplotlib references, since they cannot be pickled
        state = self.__dict__.copy()
        state['mainFrame'] = None
        state['figure'] = None
        state['axes'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)

        state['mainFrame'] = panel.mainFrame.getMainFrame()
        state['figure'] = state['mainFrame'].figure
        state['axes'] = state['mainFrame'].axes

        # since elsewhere in the code i make copies of this object,
        # i need a copy of the state so the original and copy don't refernce 
        # the same __dict__ instance
        self.__dict__ = copy.copy(state)
