Note: this file contains certain notes about the code.  This file has nothing to do with how to use the program.  See graphApp_manual.htm for that.

---------------------------------------------------------------
ACKNOWLEDGMENTS
---------------------------------------------------------------
Much of the code from the configuration dialog is from the library MPlot

http://cars9.uchicago.edu/~newville/Python/MPlot/

---------------------------------------------------------------
VERSION NUMBERS
---------------------------------------------------------------
For some of the objects used, there are version numbers applied to them.  This is because those particular objects get pickled when a .graph file is saved.  There is a problem when, for instance, a certain version of an object is pickled, then later some new class variables are added to that class.  When the old version of the class gets loaded and unpickled, it won't have those new class variables.  Version numbers are used determine if variables have to manually be added to the object.  Note that every time a new class variable is introduced to a class that has a version number, 

1. The 'latest_version' variable of that object must be increased.
2. Code must be written to the class's __setstate__ function to deal with older versions of the object.  This usually means creating said missing class variables and putting them in the __dict__ variable.

An example can be found in the __setstate__ function in the Config object in config/config.py.  

Note:  If the code involves making shallow copies of the class, then in the __setstate__ function you cannot simply state that 'self.__dict__ = state', you have to make a copy of the state and assign that.  This is because when shallow copies are made, the __setstate__ function is called, which would result in two objects referencing the same __dict__ object, which basically means that they are the same object.  An example of this is the Config class in config/config.py.

current objects with a version number:
--------------------------------------
Config() in config/config.py
Line() in widgets/line.py
LineGroup() in widgets/line.py
LineCollection() in widgets/line.py
AnnotationManager() in widgets/annotationManager.py
Interval() in widgets/annotationManager.py
Annotation() in widgets/annotations.py
Graph() in widgets/graph.py -- ?
GeneMarkers() in widgets/geneMarkers.py
GeneMarkersGroup() in widgets/geneMarkers.py

---------------------------------------------------------------
COORDINATE SPACES
---------------------------------------------------------------
When you read through the code, you may notice that the comments talk about point being in different spaces (i.e. 'this point is in data space...', 'given a point in axes space...'). There are three different type of spaces (coordinate systems) used: pixel space, data space, and axes space.

pixel space -- Position in terms of pixels.  The lower left corner of the window is (0, 0), and upper right corner of the window is width of x in pixels, width of y in pixels.

data space -- Coordinates that correspond to points defined by the graph. Boundaries of this space depends on the current viewing limits and the scale of the x and y axis. The origin of the graph is (0,0). Think of the actual points being plotted being in data space. For example, say there are two chromosomes to plot. The first one plots points from x=0 to x=100. Now the second chromosome needs to be plotted, and its first point is (0, 10). It needs to be plotted after the first chromosome, to the actual point being plotted is (0+100, 10). So the point (100, 10) is in data space.

axes space -- A fractional coordinate space for the graph axes.  x ranges from 0 to 1 and y ranges from 0 to 1. These values only extend to the section that the graph is in, not the whole window. For example, the sidebar on the right and the toolbars on top are not part of axes space. The gray area that immediately surrounds the graph are not part of axes space.

Refer to transforms.py for functions that convert between the different spaces.
