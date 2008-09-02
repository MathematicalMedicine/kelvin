"""A custom version of the builtin FigureCanvasWxAgg class, the panel that 
draws the figure

"""

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg 

from config import misc

class FigureCanvas(FigureCanvasWxAgg):

    def __init__(self, parent, id, figure):
        # SOURCE
        FigureCanvasWxAgg.__init__(self, parent, id, figure)

        # this is needed just in case there is a command-line instruction
        # to save the graph as an image before the window even shows up.
        # this ensures that the resulting image file is close to the proper
        # size. will not be exact becase of other windows in the program.
        # not perfect, but as close as i can get.
        self.SetSize(misc.mainWindowSize)
