"""Functions that deal with transforming from one coordinate space to another"""

# there are three coordinate spaces to deal with, 
# pixel space, data space, and axes space

# pixel space -- position in terms of pixels.  the lower left is (0, 0), and
# upper right is width of x in pixels, width of y in pixels

# data space -- coordinates that corresponds to points defined by the graph, 
# depends on the current viewing limits and the scale of the x and y axis

# axes space -- a fractional coordinate space for the axes.  x ranges from
# 0 to 1 and y ranges from 0 to 1

# the following are helper functions to convert from 
# one coordinate space to another

def axesToPixel(axes, x=None, y=None):
    """Axes space to pixel space"""
    
    if x is None:
        return axes.transAxes.xy_tup((0, y))[1]
    elif y is None:
        return axes.transAxes.xy_tup((x, 0))[0]
    else:
        return axes.transAxes.xy_tup((x, y))

def dataToPixel(axes, x=None, y=None):
    """Data space to pixel space"""

    if x is None:
        return axes.transData.xy_tup((0, y))[1]
    elif y is None:
        return axes.transData.xy_tup((x, 0))[0]
    else:
        return axes.transData.xy_tup((x, y))

def pixelToAxes(axes, x=None, y=None):
    """Pixel space to axes space"""

    if x is None:
        return axes.transAxes.inverse_xy_tup((0, y))[1]
    elif y is None:
        return axes.transAxes.inverse_xy_tup((x, 0))[0]
    else:
        return axes.transAxes.inverse_xy_tup((x, y))

def pixelToData(axes, x=None, y=None):
    """Pixel space to data space"""

    if x is None:
        return axes.transData.inverse_xy_tup((0, y))[1]
    elif y is None:
        return axes.transData.inverse_xy_tup((x, 0))[0]
    else:
        return axes.transData.inverse_xy_tup((x, y))

def axesToData(axes, x=None, y=None):
    """Axes space to data space"""
    
    if x is None:
        y = axesToPixel(axes, y=y)
        return pixelToData(axes, y=y)
    elif y is None:
        x = axesToPixel(axes, x=x)
        return pixelToData(axes, x=x)
    else:
        x, y = axesToPixel(axes, x, y)
        return pixelToData(axes, x, y)

def dataToAxes(axes, x=None, y=None):
    """Data space to axes space"""
    
    if x is None:
        y = dataToPixel(axes, y=y)
        return pixelToAxes(axes, y=y)
    elif y is None:
        x = dataToPixel(axes, x=x)
        print 'x is', x
        return pixelToAxes(axes, x=x)
    else:
        x, y = dataToPixel(axes, x, y)
        return pixelToAxes(axes, x, y)

