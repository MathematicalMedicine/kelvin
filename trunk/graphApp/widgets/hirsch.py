"""This module contains code to implment Hirsch's algorithm of label placement"""

# basically the gist is to go through each annotation, and find the overlaps
# resulting from it and other annotations, annotation arrows, and the lines
# in the graph.  then from the position and size of the overlap get displacement
# vectors that the annotation will be moved to so the overlaps go away.
# this process is iterated a number of times until no overlaps result or
# the max number of iterations is reached.

from bisect import bisect_left, bisect_right

from transforms import axesToPixel
from transforms import dataToPixel
from transforms import pixelToAxes
from transforms import pixelToData

class bbox():
    """Bounding box class"""

    def __init__(self, xmin, xmax, ymin, ymax):
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def getOverlap(self, box):
        """Given another bounding box, find the overlap, if any.  The return
        value is either None or another bounding box
        
        """

        if self.xmin >= box.xmax or self.xmax <= box.xmin:
            return None
        elif self.ymin >= box.ymax or self.ymax <= box.ymin:
            return None
        xmin = max(self.xmin, box.xmin)
        xmax = min(self.xmax, box.xmax)
        ymin = max(self.ymin, box.ymin)
        ymax = min(self.ymax, box.ymax)
        return bbox(xmin, xmax, ymin, ymax)

class Hirsch_layout():
    """Perform layout according to Hirsch algorithm"""

    def __init__(self, annotationManager):
        self.annManager = annotationManager

    def layout(self):
        # find the boundaries of the graph in pixel space
        self.xmin, self.ymin = axesToPixel(self.annManager.axes, 0, 0)
        self.xmax, self.ymax = axesToPixel(self.annManager.axes, 1, 1)


        # take out the annotations that aren't visible
        annotations = []
        for eachAnn in self.annManager:
            if eachAnn.get_visible():
                annotations.append(eachAnn)

        # the annotations list is currently sorted from left to right.
        # create a new list that will place annotations from the edges first,
        # leaving the ones in the middle for last.  this will result in
        # more aesthetically pleasing placements
        newAnn = []
        for i in range(0, len(annotations)):
            if i % 2 == 0:
                newAnn.append(annotations.pop(0))
            else:
                newAnn.append(annotations.pop())
        annotations = newAnn

        numIterations = 5
        for j in range(numIterations):
            numOverlaps = []
            for eachAnn in [x for x in annotations[:] if x.autoLayout]:
                box = self.getBBox(eachAnn)
                overlaps = self.findOverlaps(eachAnn, box)
                numOverlaps.append(len(overlaps))
                vectors = self.getVectors2(overlaps, eachAnn)
                sumx, sumy = self.sumVectors(vectors, box)
                if j >= numIterations-2 and len(overlaps) > 0:
                    # jolt it a bit
                    sumy += 20
                    sumx *= 1.3
                else:
                    sumy *= 1.2
                    sumx *= 1.2
                self.adjust_position(eachAnn, sumx, sumy)
                if len(overlaps) == 0:
                    annotations.remove(eachAnn)
            if sum(numOverlaps) == 0:
                break

        # make sure annotations aren't outside the graph
        for eachAnn in annotations:
            box = self.getBBox(eachAnn)
            sumx, sumy = 0, 0
            pad = 5
            if box.xmin < self.xmin:
                sumx = self.xmin - box.xmin + pad
            if box.xmax > self.xmax:
                sumx = self.xmax - box.xmax - pad
            if box.ymin < self.ymin:
                sumy = self.ymin - box.ymin + pad
            if box.ymax > self.ymax:
                sumy = self.ymax - box.ymax - pad

            self.adjust_position(eachAnn, sumx, sumy)
            
    def findOverlaps(self, annotation, box):
        """Given an annotation and its bounding box, find the overlaps with 
        other annotations, other arrows, the graph boundaries, 
        and the graph lines.
        
        """ 

        overlaps = []
        #overlaps.extend(self.findAnnotationOverlaps(annotation, box))
        overlaps.extend(self.findGraphBoundaryOverlaps(annotation, box))
        overlaps.extend(self.findArrowOverlaps(annotation, box))
        #overlaps.extend(self.findGraphOverlaps(annotation, box))
        return overlaps

    def findAnnotationOverlaps(self, annotation, box):
        """Given an annotation and its bounding box, find all overlaps 
        with other annotations
        
        """

        overlaps = []
        for eachAnn in self.annManager:
            if eachAnn is not annotation:
                box2 = self.getBBox(eachAnn)
                overlap = box.getOverlap(box2)
                if overlap is not None: 
                    overlaps.append(overlap)

        return overlaps

    def findGraphBoundaryOverlaps(self, annotation, box):
        """Given an annotation and its bounding box, find all overlaps 
        with the graph boundary
        
        """

        pad = 3
        overlaps = []
        if box.xmin < self.xmin:
            overlap = bbox(box.xmin, self.xmin+pad, box.ymin, box.ymax)
            overlaps.append(overlap)
        if box.xmax > self.xmax:
            overlap = bbox(self.xmax-pad, box.xmax, box.ymin, box.ymax)
            overlaps.append(overlap)
        if box.ymin < self.ymin:
            x0 = max(box.xmin, self.xmin)
            x1 = min(box.xmax, self.xmax)
            overlap = bbox(x0, x1, box.ymin, self.ymin+pad)
            overlaps.append(overlap)
        if box.ymax > self.ymax:
            x0 = max(box.xmin, self.xmin)
            x1 = min(box.xmax, self.xmax)
            overlap = bbox(x0, x1, self.ymax-pad, box.ymax)
            overlaps.append(overlap)

        return overlaps

    def findArrowOverlaps(self, annotation, box):
        """Given an annotation and its bounding box, find all overlaps
        from other annotations' arrows
        
        """

        overlaps = []
        for eachAnn in self.annManager:
            if eachAnn is not annotation:
                renderer = self.annManager.mainFrame.canvas.get_renderer()
                eachAnn.draw(renderer)
                x1, y1, = eachAnn.get_arrowbase()
                x2, y2, = eachAnn.get_arrowtip()
                if x1 is not None:
                    xmin = min(x1, x2)
                    xmax = max(x1, x2)
                    ymin = min(y1, y2)
                    ymax = max(y1, y2)
                    if abs(xmax-xmin) <= 2:
                        xmax += 3
                        xmin -= 3
                    if abs(ymax-ymin) <= 2:
                        ymax += 3
                        ymin -= 3
                    arrow_box = bbox(xmin, xmax, ymin, ymax)
                    overlap = box.getOverlap(arrow_box)
                    if overlap is not None:
                        overlaps.append(overlap)

        return overlaps

    def findGraphOverlaps(self, annotation, box):
        """Given an annotation and its bounding box, find all overlaps
        from collisions with the graph.
        
        """

        axes = self.annManager.axes
        overlaps = []
        xmin = pixelToData(axes, x=box.xmin)
        xmax = pixelToData(axes, x=box.xmax)
        for eachGroup in self.annManager.mainFrame.lines.getGroups():
            if eachGroup.doesOverlap(xmin, xmax):
                offset = eachGroup.xmin
                for eachLine in eachGroup:
                    start = bisect_left(eachLine.pos, xmin-offset) - 1
                    start = max(0, start)
                    end = bisect_right(eachLine.pos, xmax-offset) 
                    end = min(len(eachLine.pos) - 1, end)
                    for i in range(start, end):
                        left = eachLine.pos[i] + offset
                        right = eachLine.pos[i+1] + offset
                        up = max(eachLine.ppl[i], eachLine.ppl[i+1]) 
                        down = min(eachLine.ppl[i], eachLine.ppl[i+1])
                        left, down = dataToPixel(axes, left, down)
                        right, up = dataToPixel(axes, right, up)
                        box2 = bbox(left, right, down, up)
                        overlap = box.getOverlap(box2)
                        if overlap is not None:
                            overlaps.append(overlap)
        return overlaps

    def getVectors(self, overlaps, annotation):
        """Given a list of overlaps, return a list of vectors
        that are formed by such overlaps, the vectors represent
        the direction needed to move the annotation in order to
        eliminate the overlap.
        
        """

        # this is the original version, where each overlap is 
        # treated separately

        threshold = 3
        vectors = []
        for eachOverlap in overlaps:
            box = self.getBBox(annotation)
            xmid = (box.xmin + box.xmax) / 2.0
            if xmid > eachOverlap.xmax:
                # needs to be moved to the right
                x = eachOverlap.xmax - box.xmin
            elif xmid < eachOverlap.xmin:
                # needs to be moved to the left
                x = eachOverlap.xmin - box.xmax
            else:
                # overlap must straddle the midpoint, figure out
                # direction by length from midpoint
                left = xmid - eachOverlap.xmin 
                right = eachOverlap.xmax - xmid
                if abs(left - right) < threshold:
                    x = 0
                else:
                    if left > right:
                        # needs to be moved to the right
                        x = eachOverlap.xmax - box.xmin
                    else:
                        # needs to be moved to the left
                        x = eachOverlap.xmin - box.xmax

            ymid = (box.ymin + box.ymax) / 2.0
            if ymid > eachOverlap.ymax:
                # needs to be moved up
                y = eachOverlap.ymax - box.ymin
            elif ymid < eachOverlap.ymin:
                # needs to be moved down
                y = eachOverlap.ymin - box.ymax
            else:
                # overlap must straddle the midpoint, figure out
                # direction by length from midpoint
                down = ymid - eachOverlap.ymin 
                up = eachOverlap.ymax - ymid
                if abs(up-down) < threshold:
                    y = 0
                else:
                    if down > up:
                        # needs to be moved up
                        y = eachOverlap.ymax - box.ymin
                    else:
                        # needs to be moved down
                        y = eachOverlap.ymin - box.ymax

            vectors.append((x, y))

        return vectors

    def getVectors2(self, overlaps, annotation):
        """Given a list of overlaps, return a list of vectors
        that are formed by such overlaps, the vectors represent
        the direction needed to move the annotation in order to
        eliminate the overlap.
        
        """

        # in this version, each overlap is grouped according to what direction
        # the annotation needs to be moved to, then the max/mins of those
        # groups are used as one vector.  this is so multiple overlaps
        # in the same direction won't make the vector too large

        threshold = 1
        vectors = []
        right_list = []
        left_list = []
        up_list = []
        down_list = []

        # group the overlaps into groups based on direction
        box = self.getBBox(annotation)
        xmid = (box.xmin + box.xmax) / 2.0
        ymid = (box.ymin + box.ymax) / 2.0
        for eachOverlap in overlaps:
            if xmid > eachOverlap.xmax:
                # needs to be moved to the right
                right_list.append(eachOverlap)
            elif xmid < eachOverlap.xmin:
                # needs to be moved to the left
                left_list.append(eachOverlap)
            else:
                # overlap must straddle the midpoint, figure out
                # direction by length from midpoint
                left = xmid - eachOverlap.xmin 
                right = eachOverlap.xmax - xmid
                if abs(left - right) > threshold:
                    if left > right:
                        # needs to be moved to the right
                        right_list.append(eachOverlap)
                    else:
                        # needs to be moved to the left
                        left_list.append(eachOverlap)

            if ymid > eachOverlap.ymax:
                # needs to be moved up
                up_list.append(eachOverlap)
            elif ymid < eachOverlap.ymin:
                # needs to be moved down
                down_list.append(eachOverlap)
            else:
                # overlap must straddle the midpoint, figure out
                # direction by length from midpoint
                down = ymid - eachOverlap.ymin 
                up = eachOverlap.ymax - ymid
                if abs(up-down) > threshold:
                    if down > up:
                        # needs to be moved up
                        up_list.append(eachOverlap)
                    else:
                        # needs to be moved down
                        down_list.append(eachOverlap)

        # create only one vector for each direction
        if len(right_list) > 0:
            xmax = max([overlap.xmax for overlap in right_list])
            x = xmax - box.xmin
            vectors.append((x, 0))
        if len(left_list) > 0:
            xmin = min([overlap.xmin for overlap in left_list])
            x = xmin - box.xmax
            vectors.append((x, 0))
        if len(up_list) > 0:
            ymax = max([overlap.ymax for overlap in up_list])
            y = ymax - box.ymin
            vectors.append((0, y))
        if len(down_list) > 0:
            ymin = min([overlap.ymin for overlap in down_list])
            y = ymin - box.ymax
            vectors.append((0, y))

        return vectors

    def sumVectors(self, vectors, box):
        """Given a list of vectors, return the sum of them"""

        sumx = sum([v[0] for v in vectors])
        sumy = sum([v[1] for v in vectors])

        # add some more displacement to y for more padding.
        # the amount of padding is the height of the box
        #if sumy != 0.0:
        #    if sumy > 0:
        #        sumy += (box.ymax - box.ymin)
        #    else:
        #        sumy -= (box.ymax - box.ymin)

        return sumx, sumy

    def getBBox(self, annotation):
        """Return the bounding box of the annotation in pixel space"""

        renderer = self.annManager.mainFrame.canvas.get_renderer()
        annotation.draw(renderer)
        l, b, w, h = annotation.get_window_extent()

        return bbox(l, l+w, b, b+h)

    def adjust_position(self, annotation, xoffset, yoffset):
        """Move the position of the annotation by the offset (x, y)"""

        x, y = annotation.get_position_pixel()
        x += xoffset
        y += yoffset
        annotation.set_position(pixelToAxes(annotation.axes, x, y))
        annotation.setArrow()
