"""Class that deals with graph annotations"""

import copy
import math

import matplotlib
from matplotlib.patches import bbox_artist, YAArrow
from matplotlib.transforms import identity_transform
from matplotlib.transforms import lbwh_to_bbox

from hirsch import Hirsch_layout
import panel.mainFrame
from transforms import axesToPixel, dataToPixel, pixelToAxes, pixelToData
from transforms import axesToData, dataToAxes

class my_YAArrow(YAArrow):
    """My own yet another arrow class.  Doesn't draw an arrow, but just a 
    line segment
    
    """

    def get_verts(self):
        """Return a list of vertices to draw"""

        # SOURCE: overwriting a already defined function
        # instead of the vertices of an arrow, just return the vertices
        # to draw a line segment

        # to find the vertices, the getpoints function is used, which
        # returns two points that are perpendicular to the line segment 
        # and intersects the second point and are distance width away
        xb, yb = self.xybase  # x base, y base
        xt, yt = self.xytip   # x tip, y tip

        # find correct thickness of the line
        width = self.width*self.dpi.get()/72./2.

        # get the two points closest to the base and the two closest
        # to the tip of the arrow
        xb1, yb1, xb2, yb2 = self.getpoints(xt, yt, xb, yb, width)
        xt1, yt1, xt2, yt2 = self.getpoints(xb, yb, xt, yt, width)

        x_list = self.convert_xunits([xb1, xb2, xt2, xt1])
        y_list = self.convert_yunits([yb1, yb2, yt2, yt1])

        return zip(x_list, y_list)

class my_mpl_annotation(matplotlib.text.Annotation):
    """My own version of matplotlib's annotation class"""

    # SOURCE: many of the functions use private variables defined in
    # matplotlib.text.Annotation, so future versions of matplotlib
    # may make this break

    def __init__(self, parent, s, xy, xycoords='data', xytext=None, 
                 textcoords=None, arrowprops=None, **kwargs):

        matplotlib.text.Annotation.__init__(self, s, xy, xycoords, xytext, 
                                            textcoords, arrowprops, **kwargs)
        self.parent = parent

    def enableArrow(self):
        """Shows an arrow pointing to the data point"""

        if self.arrowprops is None:
            self.arrowprops = self.arrowprops_backup

    def disableArrow(self):
        """Turns off the display of the arrow pointing"""

        if self.arrowprops:
            self.arrowprops_backup = self.arrowprops
            self.arrowprops = None

        if self.arrow:
            self.arrow.set_visible(False)

    def set_position2(self, (x, y)):
        """x and y are in the coordinate space that this 
        class was initialized with, usually axes space.

        """

        self.xytext = (x, y)

    def get_position_pixel(self):
        """Returns the position of the text in pixel space"""
        return axesToPixel(self.axes, *self.xytext)

    def get_position_axes(self):
        """Returns the position of the text in axes space"""
        return self.xytext

    def update_positions(self, renderer):
        """Figure out position of text and arrow"""

        # SOURCE: overwriting this function from the source code.
        # needed to change the way the arrow was being drawn

        x, y = self.xytext
        self._x, self._y = self._get_xy(x, y, self.textcoords)


        x, y = self.xy
        x, y = self._get_xy(x, y, self.xycoords)

        if self.arrowprops:
            x0, y0 = x, y
            l,b,w,h = self.get_window_extent(renderer).get_bounds()
            dpi = self.figure.dpi.get()
            r = l+w
            t = b+h
            xc = 0.5*(l+r)
            yc = 0.5*(b+t)

            # choose where to place the x value of the base of the arrow
            # if x0 is within a certain percent of xc, place it at xc,
            # otherwise linearly interpolate toward the correct corner
            p = .9
            amt = w*(1-p)/2.0
            xcs = l + amt
            xce = r - amt
            if x0 >= xcs and x0 <= xce:
                x = xc
            elif xc < x0:
                percent = (min(r, x0) - xce) / (r-xce)
                length = percent*(xc-l)
                x = xc + length 
            else:
                percent = (xcs - max(l, x0)) / (xcs-l)
                length = percent*(xc-l)
                x = xc - length

            # pick the y corner of the text bbox closest to point annotated
            dsu = [(abs(val-y0), val) for val in b, t, yc]
            dsu.sort()
            d, y = dsu[0]

            d = self.arrowprops.copy()
            width = d.pop('width', 4)
            headwidth = d.pop('headwidth', 12)
            frac = d.pop('frac', 0.1)
            shrink = d.pop('shrink', 0.0)
            head_max = d.pop('head_max', 0.0)

            theta = math.atan2(y-y0, x-x0)
            r = math.sqrt((y-y0)**2. + (x-x0)**2.)
            shrink = shrink / r
            dx = shrink*r*math.cos(theta)
            dy = shrink*r*math.sin(theta)

            headwidth_max = 10
            if r*frac > headwidth_max:
                frac = headwidth_max / r

            frac = 0
            headwidth = width+1

            self.arrow = my_YAArrow(self.figure.dpi, (x0+dx ,y0+dy), 
                            (x-(dx/2.), y-(dx/2.)),
                            width=width, headwidth=headwidth, frac=frac,
                            **d)
            self.arrow.set_clip_box(self.get_clip_box())

class Annotation():
    """Class that deals with graph annotations"""

    # note: changed definition of shrink, it is now the number of pixels
    # the arrow tip is away from the point. also added head_max, which
    # is the maximum size in pixels of the length of the head
    arrowprops = dict(facecolor='black', width=1.0, headwidth=5, head_max=10, 
                      shrink=7)

    # the number of pixels the text has to be away from the data point
    # in order for the arrow to show up
    arrowLimit = 15
    arrowLimit_sqr = arrowLimit*arrowLimit

    # the default position of the text will be an offset away from the
    # data point. this offset is in pixel units.
    offset = (0, 10)

    # padding for the bounding box of the text
    pad = 2


    def __init__(self, axes=None, x=0, y=0, text='', real=True):
        """x and y is the location of the data point for the 
        annotation in data space.
        
        """

        self.axes = axes
        self.x = x
        self.y = y
        self.text = text

        self.isPeak = False
        self.isTrough = False

        # this is the line (line as in chromosome) that the annotation 
        # is associated with. this is set only when the 
        # annotate max or annotate min is used
        self.line = None

        # the id number of the associated line. this is useful when loading
        # the annotation from a .grx file
        self.lineId = 0

        # a flag that indicates if this annotation should be used in the
        # automatic layout annotation algorithm. right now if an annotation has
        # been manually moved, it should not be automatically layout-ed again
        # for that view
        self.autoLayout = True

        if real:
            self.annotation = self.create_mpl_annotation()
            self.add_annotation()
        else:
            self.annotation = None

    def create_mpl_annotation(self):
        """Create the actual matplotlib annotation object"""

        textPos = self.getTextStartPosition()
        annotation = my_mpl_annotation(self, 
                                       self.text, 
                                       xy=(self.x, self.y), 
                                       xycoords='data', 
                                       xytext=textPos,
                                       textcoords='axes fraction', 
                                       arrowprops=self.arrowprops,
                                       horizontalalignment='center', 
                                       verticalalignment='bottom')
        annotation.set_picker(True)
        annotation.set_transform(identity_transform())
        self.axes._set_artist_props(annotation)
        annotation.disableArrow()

        return annotation

    def set_default_position(self):
        """Set the annotation's position to the default"""

        default_x, default_y = self.getTextStartPosition()
        self.set_position((default_x, default_y))

    def getTextStartPosition(self):
        """Find the position the text will start at"""

        # the start position of the text is based on an offset in pixels,
        # so first find the data point position in pixels, then
        # add the offset, then convert to axes coordinates
        x, y = dataToPixel(self.axes, self.x, self.y)
        x += self.offset[0]
        y += self.offset[1]
        return pixelToAxes(self.axes, x, y)
        
    def enableArrow(self):
        """Shows an arrow pointing to the data point"""
        self.annotation.enableArrow()

    def disableArrow(self):
        """Turns off the display of the arrow pointing"""
        self.annotation.disableArrow()

    def setArrow(self):
        """Considering the current text location, determines if the arrrow
        should show up or not
        
        """

        textx, texty = self.annotation.get_position_pixel()
        datax, datay = self.x, self.y
        datax, datay = dataToPixel(self.axes, datax, datay)
        x = textx - datax
        y = texty - datay
        dist_squared= x*x + y*y
        if dist_squared > self.arrowLimit_sqr:
            self.enableArrow()
        else:
            self.disableArrow()

    def setLine(self, line):
        """Associate this annotation with a given Line object"""

        self.line = line
        self.lineId = line.id

    def set_position(self, (x, y)):
        """Set the position of the text. x and y are in the coordinate space 
        that the annotation was originally created in, usually axes space.
        
        """

        self.annotation.set_position2((x, y))

    def set_visible(self, isVisible):
        """Sets visibility"""
        self.annotation.set_visible(isVisible)

    def set_text(self, text):
        """Change the text of the annotation to a new value"""

        self.text = text
        self.annotation.set_text(text)

    def set_color(self, color):
        """Change the color of the annotation"""
        self.annotation.set_color(color)

    def get_window_extent(self, renderer=None):
        """Return the bounding box of the text, in pixels"""

        l, b, w, h = self.annotation.get_window_extent(renderer).get_bounds()
        return l-self.pad, b-self.pad, w+self.pad, h+self.pad

    def get_position(self):
        """Get the position from the mpl_annotation"""
        return self.annotation.get_position()

    def get_position_axes(self):
        """Return the position of the annotation"""
        return self.annotation.get_position_axes()

    def get_position_pixel(self):
        """Return the position in pixel space"""
        return self.annotation.get_position_pixel()

    def get_text(self):
        """Return text of annotation"""
        return self.annotation.get_text()

    def isVisible(self, xmin, xmax, ymin, ymax):
        """Given a range in data space, determines if the annotations's
        data point is visible in that range
        
        """

        if self.x >= xmin and self.x <= xmax:
            if self.y >= ymin and self.y <= ymax:
                return True

        return False

    def get_visible(self):
        """Returns whether the annotation is visible or not"""
        return self.annotation.get_visible()

    def add_annotation(self):
        """Register the annotation with the axes"""

        if self.annotation is None:
            # need to create the mpl_annotation object
            self.annotation = self.create_mpl_annotation()
        self.axes.add_artist(self.annotation)

    def get_arrowtip(self):
        """Return the position of the tip of the arrow in pixel space"""

        if hasattr(self.annotation, 'arrow'):
            if self.annotation.arrow is not None:
                if self.annotation.arrow.get_visible():
                    return self.annotation.arrow.xytip
        return None, None

    def get_arrowbase(self):
        """Return the position of the base of the arrow in pixel space"""

        if hasattr(self.annotation, 'arrow'):
            if self.annotation.arrow is not None:
                if self.annotation.arrow.get_visible():
                    return self.annotation.arrow.xybase
        return None, None

    def getInfoText(self):
        """Return text that is suitable to be shown in the info panel"""

        text = ''
        if self.isPeak:
            text = 'chr ' + self.line.name + ': p - ' + self.text + '\n'
        elif self.isTrough:
            text = 'chr ' + self.line.name + ': v - ' + self.text + '\n'

        return text

    def draw(self, renderer):
        """Draws the annotation using the renderer"""
        self.annotation.draw(renderer)

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = Annotation(real=False).__dict__
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
        state['annotation'] = None
        state['axes'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)
        state['axes'] = panel.mainFrame.getMainFrame().axes
        self.__dict__ = state
