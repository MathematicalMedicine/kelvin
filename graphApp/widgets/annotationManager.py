"""Manages the layout of annotations and other duties"""

import copy

from annotations import Annotation
from annotations import my_mpl_annotation
import panel.mainFrame
from transforms import pixelToAxes, pixelToData
from transforms import axesToData, dataToAxes
from widgets.hirsch import Hirsch_layout
from utility.util import numToText

class AnnotationManager():
    """Manages the layout of annotations and other duties"""

    # this is the number of intervals (i.e. vertical strips) the screen
    # is divided into to keep track of empty space in them
    numIntervals = 20

    # draw the intervals or not, used for debug purposes
    showIntervals = False

    def __init__(self, mainFrame=None, axes=None, real=True):
        self.mainFrame = mainFrame
        self.axes = axes

        # a list of all the annotations (Annotation class)
        self.annotations = []

        # the screen is split into vertical strips called intervals
        self.intervals = []
        self.createIntervals()

        # lines to show the intervals, used for debugging
        self.interval_hlines = []
        self.interval_vlines = []

        # a stack that keeps track of annotation locations 
        # between different views, locations are in axes space
        # this stack  holds the list for all the views going back
        # format for each item is <annotation>, <x, y location>
        self.back_loc_stack = []

        # a stack that keeps track of annotation locations 
        # between different views, locations are in axes space
        # this stack  holds the list for all the views going forward
        # format for each item is <annotation>, <x, y location>
        self.forward_loc_stack = []

        # a list of all visible annotations and their locations in axes space,
        # used when saving and loading the graph
        self.saved_locations = []

        # a stack that keeps track of the autoLayout flag for previous views
        # (views that happen when you hit the back button)
        self.back_layout_stack = []

        # a stack that keeps track of the autoLayout flag for forward views
        # (views that happen when you hit the forward button)
        self.forward_layout_stack = []

        # an object that implements the hirsh label layout algorithm
        # does the dirty work of how to automatically place annotations
        self.hirsch_layout = Hirsch_layout(self)

    def annotate(self, x, y):
        """Add an annotation with the specified x and y points. 
        x and y are in data coordinates. Creates text of annotation
        automatically.

        """

        text = self.getAnnotationText(x, y)
        return self.createAnnotation(x, y, text)

    def createAnnotation(self, x, y, text):
        """Create and add an annotation with specified attributes"""

        annotation = Annotation(self.axes, x, y, text)
        self.annotations.append(annotation)
        self.sort()

        return annotation

    def getAnnotationText(self, x, y):
        """Given an x, y position, return the text that should be there"""

        # the x, y position should be in axes spaces
        # first change the x-coordinate so it is relative to the group it is in
        x2 = 0
        for eachGroup in self.mainFrame.lines.getGroups():
            if eachGroup.xmin <= x and eachGroup.xmax >= x:
                x2 = x - eachGroup.xmin + eachGroup.start
                break

        # round x to number of decimal places specified in config
        xtext = numToText(x2, self.mainFrame.config.annPosPrec, 'round')

        if self.mainFrame.config.defaultAnnPPLPrec == True:
            # normally round the ppl value to two decimal places, 
            # if ppl is <.02, then use three decimal places
            if y < .02:
                y = round(y, 3)
                ytext = '%.1f%%' % (y*100,)
            else:
                y = round(y, 2)
                ytext = '%i%%' % (y*100,)
        else:
            # round y to number of decimal places specified in config.
            # since it's always shown as a percentage, need to change 
            # y value and precision value
            ytext = numToText(y*100, self.mainFrame.config.annPPLPrec-2,'round')
            ytext = ytext + '%'

        return ytext + '@' + xtext

    def hasDefaultText(self, ann):
        """Given an annotation, see if it contains the default text"""

        # check if there is any substring of the annotation text that contains
        # a '%@', and a number before and after that
        index = ann.text.find('%@')
        if index == -1:
            return False

        if ann.text[index-1].isdigit() and ann.text[index+2].isdigit():
            return True
        else:
            return False

    def getDefaultText(self, ann):
        """Given an annotation, return the default text if there is any."""

        # first check if there is some default text in the annotation
        if not self.hasDefaultText(ann):
            return None

        # parse out the numbers before and after '%@'
        text = ann.text
        index = text.find('%@')
        i = index - 1
        num1 = ''
        while i >= 0 and (text[i].isdigit() or text[i] == '.'):
            num1 = text[i] + num1
            i -= 1
        i = index + 2
        num2 = ''
        length = len(text)
        while i < length and (text[i].isdigit() or text[i] == '.'):
            num2 = num2 + text[i]
            i += 1

        return num1 + '%@' + num2

    def update(self):
        """Update all annotations. Usually called when the precision required 
        for position or ppl changes or color of annotation changes.
        
        """

        for eachAnn in self.annotations:
            # first change precision of position, determine if this annotation 
            # contains the default text anywhere
            if self.hasDefaultText(eachAnn):
                # replace the old text with the new one
                oldText = self.getDefaultText(eachAnn)
                newText = self.getAnnotationText(eachAnn.x, eachAnn.y)
                eachAnn.set_text(eachAnn.text.replace(oldText, newText))

            # change the color of the annotation
            if self.mainFrame.config.isAnnSameColor:
                # set all annotations to black
                eachAnn.set_color('k')
            else:
                # set annotation to line color
                eachAnn.set_color(eachAnn.line.color)

    def annotateThreshold(self, threshold, annotate_peaks):
        """Annotate all peaks above or troughs below a certain threshold.
        annotate_peaks is a boolean indicating if peaks or troughs are used
        
        """

        points, lines = self.mainFrame.lines.getAllPeaksOrTroughs(threshold, 
                                                                 annotate_peaks)
        for eachPoint, eachLine in zip(points, lines):
            annotation = self.annotate(*eachPoint)
            annotation.setLine(eachLine)
            if annotate_peaks:
                annotation.isPeak = True
            else:
                annotation.isTrough = True

        self.layout()
        self.mainFrame.updateInfoText()

    def annotateMax(self, x0, y0, x1, y1):
        """Given the coordinates of a rectangle in data space, annotate
        the highest peak in that rectangle
        
        """

        xmax = max(x0, x1)
        xmin = min(x0, x1)
        ymax = max(y0, y1)
        ymin = min(y0, y1)

        maxPoint, line = self.mainFrame.lines.getPeak(xmin, xmax, ymin, ymax)
        if maxPoint is not None:
            annotation = self.annotate(*maxPoint)
            annotation.isPeak = True
            annotation.line = line
            if not self.mainFrame.config.isAnnSameColor:
                annotation.set_color(line.color)
            self.layout()
            self.mainFrame.updateInfoText()

    def annotateMin(self, x0, y0, x1, y1):
        """Given the coordinates of a rectangle in data space, annotate
        the lowest valley in that rectangle
        
        """

        xmax = max(x0, x1)
        xmin = min(x0, x1)
        ymax = max(y0, y1)
        ymin = min(y0, y1)

        minPoint, line = self.mainFrame.lines.getTrough(xmin, xmax, ymin, ymax)
        if minPoint is not None:
            annotation = self.annotate(*minPoint)
            annotation.isTrough = True
            annotation.line = line
            if not self.mainFrame.config.isAnnSameColor:
                annotation.set_color(line.color)
            self.layout()
            self.mainFrame.updateInfoText()

    def createIntervals(self):
        """Create and initialize the intervals of the screen"""

        # the intervals are defined in axes space, where x and y
        # values range from 0 to 1
        intervalWidth = 1.0 / self.numIntervals

        for i in range(self.numIntervals):
            xmin = i*intervalWidth
            xmax = xmin + intervalWidth
            interval = Interval(i, xmin, xmax, 0, 1)

            # in the beginning it's all empty, so the y-value
            # emtpy space interval is 0 to 1
            interval.append([0, 1])
            self.intervals.append(interval)

    def getMaxDataPoint(self, xmin, xmax):
        """Given an x range, find the maximum data point found in that range"""

        # add some padding, in pixel space
        pad = 4
        y1 = pixelToData(self.axes, y=0)
        y2 = pixelToData(self.axes, y=pad)
        pad = y2 - y1

        # convert limits to data space first
        xmin = axesToData(self.axes, x=xmin)
        xmax = axesToData(self.axes, x=xmax)

        max = self.mainFrame.lines.getMaxDataPoint(xmin, xmax)
        if max is not None:
            x, y = max

            # answer was in data space, convert to axes space
            # add some padding, too
            return dataToAxes(self.axes, x, y+pad)
        else:
            return max

    def layout(self):
        """Layout the visible annotations so they don't overlap"""

        # the fist part of the layout basically determines where vertically the 
        # annotations can be placed without overlapping each other.
        # it does not move annotations left and right

        if len(self.annotations) < 1:
            return

        # get intervals ready
        self.resetIntervals()

        # the annotations list is currently sorted from left to right.
        # create a new list that will place annotations from the edges first,
        # leaving the ones in the middle for last.  this will result in
        # more aesthetically pleasing placements
        newAnn = []
        oldAnn = [x for x in self.annotations[:] if x.autoLayout]
        for i in range(0, len(oldAnn)):
            if i % 2 == 0:
                newAnn.append(oldAnn.pop(0))
            else:
                newAnn.append(oldAnn.pop())

        view_xmin, view_xmax = self.axes.get_xlim()
        view_ymin, view_ymax = self.axes.get_ylim()

        # first place the annotations that shouldn't move.
        # they still have to be inserted into the intervals
        exemptAnn = [x for x in self.annotations if x.autoLayout == False]
        for eachAnn in exemptAnn:
            eachAnn.draw(self.mainFrame.canvas.get_renderer())
            l, b, w, h = eachAnn.get_window_extent()

            # change to axes space
            xmin = pixelToAxes(self.axes, x=l)
            xmax = pixelToAxes(self.axes, x=l+w)
            ymin = pixelToAxes(self.axes, y=b)
            ymax = pixelToAxes(self.axes, y=b+h)

            intervals = self.getOverlappingIntervals(xmin, xmax)
            self.insertIntoIntervals(intervals, ymin, ymax-ymin)

        for eachAnn in newAnn:
            if eachAnn.isVisible(view_xmin, view_xmax, view_ymin, view_ymax):
                # set annotation to initial position
                eachAnn.set_visible(True)
                eachAnn.set_default_position()
                default_x, default_y = eachAnn.getTextStartPosition()

                # get bounding box of text
                # left, bottom, width, height
                eachAnn.draw(self.mainFrame.canvas.get_renderer())
                l, b, w, h = eachAnn.get_window_extent()

                # change to axes space
                xmin = pixelToAxes(self.axes, x=l)
                xmax = pixelToAxes(self.axes, x=l+w)
                ymin = pixelToAxes(self.axes, y=b)
                ymax = pixelToAxes(self.axes, y=b+h)

                intervals = self.getOverlappingIntervals(xmin, xmax)
                y = self.findFirstEmptySpace(intervals, ymin, ymax-ymin)

                eachAnn.set_position((default_x, y))
                eachAnn.setArrow()
            else:
                # annotation is not visible
                eachAnn.set_visible(False)
                eachAnn.disableArrow()

        self.hirsch_layout.layout()

    def getOverlappingIntervals(self, xmin, xmax):
        """Given an x range in axes space, returns a list of all 
        intervals that overlap that range
        
        """

        intervals = []
        for eachInterval in self.intervals:
            if eachInterval.xmin < xmin and eachInterval.xmax > xmin:
                intervals.append(eachInterval)
            elif eachInterval.xmin >= xmin and eachInterval.xmax <= xmax:
                intervals.append(eachInterval)
            elif eachInterval.xmin <= xmin and eachInterval.xmax >= xmax:
                intervals.append(eachInterval)
            elif eachInterval.xmin < xmax and eachInterval.xmax > xmax:
                intervals.append(eachInterval)
        return intervals

    def findFirstEmptySpace(self, intervals, ymin, height):
        """Given a list of intervals and a height, find the lowest point
        that a slot with that height can be placed in
        
        """

        epsilon = 0.001

        # first find the lowest possible emtpy space over all intervals
        new_ymin = True
        while new_ymin:
            new_ymin = False
            for eachInterval in intervals:
                for eachSlot in eachInterval:
                    if eachSlot[0] <= ymin and eachSlot[1]-ymin >= height:
                        # this slot is ok, no need to change ymin
                        break
                    if eachSlot[0] >= ymin and eachSlot[1]-eachSlot[0]>=height:
                        # need to change ymin, and go over all intervals again
                        ymin = eachSlot[0]
                        new_ymin = True
                        break

        # insert the new slot into each interval
        self.insertIntoIntervals(intervals, ymin, height)

        return ymin

    def insertIntoIntervals(self, intervals, ymin, height):
        """Given a list of intervals, and an occupied space, 
        inserts that occupied space into all intervals
        
        """

        # clarification: the name 'slot' applies to pairs of numbers
        # e.g. [y1, y2], that designate stretches of empty space

        epsilon = 0.001

        def checkSlot(interval, slot):
            """If a slot is small enough, remove it from the interval"""
            if slot[1]-slot[0] < epsilon:
                interval.remove(slot)

        # insert the new slot into each interval
        ymax = ymin + height
        for eachInterval in intervals:
            for i in range(len(eachInterval)):
                # need to take into account the cases where the new 
                # occupied space is in the middle of an empty slot, 
                # or overlaps an empty slot and occupied space at the same time
                slot = eachInterval[i]
                if slot[0] <= ymin and slot[1] >= ymax:
                    # since inserting something in the middle of an empty
                    # slot, the old slot is removed and split into two
                    slot1 = (slot[0], ymin)
                    slot2 = (ymax, slot[1])
                    index = eachInterval.index(slot)
                    eachInterval[index:index] = [slot1, slot2]
                    eachInterval.remove(slot)

                    # eliminate empty slots if they are small enough
                    checkSlot(eachInterval, slot1)
                    checkSlot(eachInterval, slot2)
                    break

                if i > 0:
                    prevSlot = eachInterval[i-1]
                    if prevSlot[1] <= ymin and slot[0] >= ymin:
                        if slot[0] <= ymax and slot[1] >= ymax:
                            slot[0] = ymax
                            checkSlot(slot)
                            break

                    if prevSlot[0] <= ymin and prevSlot[1] >= ymin:
                        if prevSlot[1] <= ymax and slot[0] >= ymax:
                            prevSlot[1] = ymin
                            checkSlot(prevSlot)
                            break

                # two special cases, in which there is only one
                # empty slot, and occupied space at the bottom or top
                if len(eachInterval) == 1: 
                    if slot[0] <= ymin and slot[1] >= ymin:
                        if slot[1] <= ymax and eachInterval.ymax >= ymax:
                            slot[1] = ymin
                            checkSlot(slot)
                            break
                    if slot[0] <= ymax and slot[1] >= ymax:
                        if eachInterval.ymin <= ymin and slot[0] >= ymin:
                            slot[0] = ymax
                            checkSlot(slot)
                            break

    def resetIntervals(self):
        """Clear all information in intervals, and set it to one empty space
        that is above the max data point in the interval
        
        """

        for interval in self.intervals:
            maxPoint = self.getMaxDataPoint(interval.xmin, interval.xmax)
            if maxPoint is None:
                maxPoint = (0, 0)
            interval.clear()
            ymin = max(0, maxPoint[1])
            interval.ymin = ymin
            interval.append([ymin, 1])

        if self.showIntervals:
            self.drawIntervals()

    def drawIntervals(self):
        """Draws lines diving the partitions, and draws a horizontal line
        over the max data point in that partition
        
        """

        xmin_view, xmax_view = self.axes.get_xlim()
        ymin_view, ymax_view = self.axes.get_ylim()

        # do this to clear old lines and clearly show new ones
        for eachLine in self.interval_hlines:
            eachLine.set_visible(False)
        for eachLine in self.interval_vlines:
            eachLine.set_visible(False)

        self.axes.set_xlim(xmin=xmin_view, xmax=xmax_view)
        self.axes.set_ylim(ymin=ymin_view, ymax=ymax_view)

        vlines = []
        hlines = []
        for interval in self.intervals:
            x = axesToData(self.axes, x=interval.xmin)
            vlines.append(x)

            y = axesToData(self.axes, y=interval[0][0])
            hlines.append((y, interval.xmin, interval.xmax))

        for i in range(len(vlines)):
            line = self.axes.axvline(vlines[i], color='k')
            self.interval_vlines.append(line)
            y, xmin, xmax = hlines[i]
            line = self.axes.axhline(y=y, xmin=xmin, xmax=xmax, color='k')
            self.interval_hlines.append(line)

        self.axes.set_xlim(xmin=xmin_view, xmax=xmax_view)
        self.axes.set_ylim(ymin=ymin_view, ymax=ymax_view)

    def remove(self, annotation):
        """Given an annotation, remove it"""
        
        if isinstance(annotation, my_mpl_annotation):
            # argument is an mpl_annotation class,
            # but need Annotation class instance
            annotation = annotation.parent
        annotation.set_visible(False)
        annotation.disableArrow()
        self.annotations.remove(annotation)

    def sort(self):
        """Sort the annotations by x values, then y values"""

        # need to define a compare method for the annotation objects
        def annotation_compare(a, b):
            if a.x > b.x:
                return 1
            elif a.x < b.x:
                return -1
            else:
                # a.x == b.x, so sort by y
                if a.y > b.y:
                    return 1
                elif a.y < b.y:
                    return -1
                else:
                    # both x and y are equal
                    return 0

        self.annotations.sort(annotation_compare)

    def printOut(self):
        """print out all annotations in order"""

        print '-'*30
        for eachAnn in self.annotations:
            print '(%.2f, %.2f)' % (eachAnn.x, eachAnn.y)
        print '-'*30

    def add_annotations(self):
        """Register all annotations the manager has with the axes"""

        for eachAnn in self.annotations:
            eachAnn.add_annotation()

    def save_positions(self):
        """Save all annotations and their positions"""

        # this is used when saving the graph to a .graph file
        self.saved_locations = []
        for eachAnn in self.annotations:
            if eachAnn.get_visible():
                item = (eachAnn, eachAnn.get_position_axes())
                self.saved_locations.append(item)

    def restore_positions(self):
        """Restore all previous annotation positions"""

        # this is used when the graph is loaded from a .graph file
        for eachAnn, eachPos in self.saved_locations:
            eachAnn.set_visible(True)
            eachAnn.set_position(eachPos)
            eachAnn.setArrow()

    def push_back(self, newView=True):
        """Push all annotation locations on the back_loc_stack and
        store autoLayout flags in the back_layout_stack.
        
        """

        self.pushLocation(self.back_loc_stack)
        self.pushLayout(self.back_layout_stack)

        # a new view, so reset autoLayout flags
        for eachAnn in self.annotations:
            eachAnn.autoLayout = True

        # as per how the back and forward buttons work, when a new view
        # is pushed in the back stack, all the views in the forward
        # stack are removed
        if newView:
            self.forward_loc_stack = []
            self.forward_layout_stack = []

    def push_forward(self):
        """Push all annotation locations on the forward_loc_stack and
        store autoLayout flags in the forward_layout_stack.
        
        """

        self.pushLocation(self.forward_loc_stack)
        self.pushLayout(self.forward_layout_stack)

    def pop_back(self):
        """Pop all annotation locations on the back_loc_stack and
        pop all flags in the back_layout_stack
        
        """

        self.popLocation(self.back_loc_stack)
        self.popLayout(self.back_layout_stack)

    def pop_forward(self):
        """Pop all annotation locations on the forward_loc_stack and
        pop all flags in the forward_layout_stack
        
        """

        self.popLocation(self.forward_loc_stack)
        self.popLayout(self.forward_layout_stack)

    def pushLocation(self, stack):
        """Push all visible annotation locations on the stack"""

        info = []
        for eachAnn in self.annotations:
            if eachAnn.get_visible():
                item = (eachAnn, eachAnn.get_position_axes())
                info.append(item)
        stack.append(info)

    def popLocation(self, stack):
        """Pop the top entry of the stack and set the locations"""

        if len(stack) > 0:
            info = stack.pop()
            for eachAnn, eachPos in info:
                if eachAnn in self.annotations:
                    eachAnn.set_position(eachPos)
                    eachAnn.setArrow()
            self.mainFrame.draw()

    def pushLayout(self, stack):
        """Push all autoLayout flags on the stack"""

        info = []
        for eachAnn in self.annotations:
            item = (eachAnn, eachAnn.autoLayout)
            info.append(item)
        stack.append(info)

    def popLayout(self, stack):
        """Pop the top entry of the stack and apply the flags"""

        if len(stack) > 0:
            info = stack.pop()
            for eachAnn, eachFlag in info:
                if eachAnn in self.annotations:
                    eachAnn.autoLayout = eachFlag

    def beginPan(self, event):
        """Event when a pan begins"""

        self.panPoint = (event.x, event.y)
        self.panStartLocs = []
        for eachAnn in self.annotations:
            if eachAnn.get_visible():
                item = (eachAnn, eachAnn.get_position_pixel())
                self.panStartLocs.append(item)

    def updatePan(self, event):
        """Event when in the middle of a pan"""

        # need to adjust the positions of all annotations so that they move
        # the exact same amount as the pan

        if len(self.annotations) < 1:
            return

        # calculate the amount to move
        view_xmin, view_xmax = self.axes.get_xlim()
        view_ymin, view_ymax = self.axes.get_ylim()
        offset = (event.x-self.panPoint[0], event.y-self.panPoint[1])

        # check to see if user is using special panning modes
        if event.key == 'x':
            # only use the x offset
            offset = (offset[0], 0)
        elif event.key == 'y':
            # only use the y offset
            offset = (0, offset[1])

        # go through each annotation and update the position as necessary
        old_anns = [x[0] for x in self.panStartLocs]
        old_panStartLocs = self.panStartLocs[:]
        for eachAnn in self.annotations:
            if eachAnn.isVisible(view_xmin, view_xmax, view_ymin, view_ymax):
                # this annotation is visible, check if it was visible before
                # or it newly became visible because of panning
                if eachAnn in old_anns:
                    # already keeping track of this annotation, update position
                    for eachItem in old_panStartLocs:
                        if eachItem[0] == eachAnn:
                            x, y = eachItem[1]
                            x += offset[0]
                            y += offset[1]
                            x, y = pixelToAxes(self.axes, x, y)
                            eachAnn.set_position((x, y))
                            break
                else:
                    # a new annotation that became visible, set it to
                    # a default position
                    eachAnn.set_visible(True)
                    eachAnn.set_default_position()
                    eachAnn.setArrow()
                    pos = eachAnn.get_position_pixel()
                    x = pos[0] - offset[0]
                    y = pos[1] - offset[1]
                    item = (eachAnn, (x, y))
                    self.panStartLocs.append(item)
            elif eachAnn in old_anns:
                # this is an annnotation that was visible but has moved
                # out of the screen because of panning
                eachAnn.set_visible(False)
                eachAnn.disableArrow()
                for eachItem in self.panStartLocs:
                    if eachItem[0] == eachAnn:
                        self.panStartLocs.remove(eachItem)
                        break

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = AnnotationManager(real=False).__dict__
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

        for eachAnn in self.annotations:
            yield eachAnn

    def __getitem__(self, key):
        """This is used so the class can be accessed like a list"""

        return self.annotations[key]

    def __getstate__(self):
        """Returns the state used when pickling."""

        # remove all matplotlib references, since they cannot be pickled
        state = self.__dict__.copy()
        state['mainFrame'] = None
        state['axes'] = None
        return state

    def __setstate__(self, state):
        """Sets the state of the object when unpickled.  Uses the version
        number to determine if any conversion is needed.
        
        """

        state = self.updateSelf(state)

        state['mainFrame'] = panel.mainFrame.getMainFrame()
        state['axes'] = state['mainFrame'].axes
        state['forward_loc_stack'] = []
        self.__dict__ = state

class Interval(list):
    """Keeps track of available empty space in a vertical strip of the screen"""

    def __init__(self, number=0, xmin=0, xmax=0, ymin=0, ymax=0, real=True,
                 *args, **kwargs):
        # basically this class is just a regular list with two extra
        # member variables, to keep track of the x values this
        # interval represents.  then tuples of y-values of where the empty
        # space is appended to this list

        # the boundaries (in axes space) of the interval
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

        # the number of this interval, with interval 0 being the leftmost
        # (mainly used for debugging)
        self.number = number

        list.__init__(self, *args, **kwargs)

    def clear(self):
        """Remove all items from this list"""

        for i in range(len(self)):
            self.pop()

    def updateSelf(self, curState):
        """Given a dict of the current state of an object, update it to make 
        sure it contains all the current class attributes.
        
        """

        # this is used when a .graph file is loaded, and an object of this 
        # class is unpickled. an old version of the class could be loaded, 
        # so need to update it.
        # make a temp object of the same type, and compare dicts
        newState = Interval(real=False).__dict__
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

