"""Settings used for the graph program"""

import os
import wx
import images

# the default boundaries of the figure (the graph)
# format is <left> <bottom> <width> <height>
defaultFigureSize = (0.06, 0.1, 0.90, 0.83)

# this is the y value position of the chromosome labels at the bottom
# it's in axes space, where 0 is the bottom of the axes, 1 is top
chrlabel_y_pos = -0.035

# this is the y value position of the x-axis label
# it's in figure space, where 0 is the bottom of the figure and 1 is top
xlabel_y_pos = 0.04

# a list of all colors used to color lines if multiple colors are wanted
colors = ['b', 'g', 'c', 'r', 'm', 'y']

# the path to the images folder, which is always one folder above 
# the folder that contains this file
#imagePath = os.path.abspath(images.__name__)
imagePath = os.path.dirname(__file__)
imagePath = imagePath + os.path.sep + os.path.pardir + \
                    os.path.sep + images.__name__

def getImagePath(filename):
    """Given a filename to an image file, creates the entie filepath to it """
    return os.path.join(imagePath, filename)

# all the available line styles
# format is <description> <style>
lineStyleEntry = [('solid', '-'), ('dashed', '--'), ('dotted', ':'), 
                  ('dash-dot', '-.'), ('steps', 'steps'), ('no line', 'None')]

lineStyleDescriptions = [x[0] for x in lineStyleEntry]
lineStyleTags = [x[1] for x in lineStyleEntry]
lineStyleDescrToTag = dict(lineStyleEntry)

lineSymbolEntry = [('no symbol', 'None'), ('square', 's'), ('+', '+'),
                   ('o', 'o'), ('x', 'x'), ('diamond', 'D'), 
                   ('thin diamond', 'd'), ('^', '^'), ('v', 'v'), ('>', '>'), 
                   ('<', '<'), ('|', '|'), ('_', '_'), ('hexagon','h'), 
                   ('pentagon','p'), ('tripod 1', '1'), ('tripod 2', '2')]

lineSymbolDescriptions = [x[0] for x in lineSymbolEntry]
lineSymbolTags = [x[1] for x in lineSymbolEntry]
lineSymbolDescrToTag = dict(lineSymbolEntry)

# all possible positions to place the legend
legendPositions = ['upper right' , 'upper left', 'upper center',
                   'lower right',  'lower left', 'lower center',
                   'center left',  'center right', 'right', 'center']

# defaults for gene marker attributes 
markers_size = 4
markers_symbol = 'o'
markers_color = '#FF0000'

# flag to determine whether to allow duplicates in the chromosome labels.
# for example, '1, 1, 1' becomes '1' if allowDuplicates is false
allowDuplicates = False

# this is a string that the command_line module uses to convey to
# other functions it call that it is being called as part of a
# command line instruction
cmdline = 'Command Line'

