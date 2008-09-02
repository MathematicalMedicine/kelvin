"""A class that represents all the information in a graph"""

import copy

from config.config import Config
from widgets.annotationManager import AnnotationManager
from widgets.geneMarkers import GeneMarkersGroup
from widgets.lineCollection import LineCollection

class Graph():
    """A class that represents all the information in a graph"""

    def __init__(self, parent=None, real=True):
        # the Config object, contains all the settings of the graph
        self.config = Config(parent, real)

        # a collection of lines
        self.lineCollection = LineCollection()

        # handles all the annotations
        figure = self.config.figure
        axes = self.config.axes
        self.annotationManager = AnnotationManager(parent, axes, real)

        # a list of all the gene marker groups
        self.geneMarkerGroups = []

        # a list of floats that specify y-values of horizontal lines
        self.HLines = []

        # the comments about the graph on the side of the frame
        self.comments = ''

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = Graph(real=False).__dict__
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
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        self.__dict__ = state

