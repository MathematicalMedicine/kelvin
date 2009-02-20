#!/usr/bin/env python

"""
This is a graphing program designed to take data files of chromosome position 
and ppl values and graph them against each other.  Either one file or multiple 
files can be opened at once.

"""

import math
import matplotlib
import wx

from config import misc
from panel.mainFrame import MainFrame

# SOURCE:
# I came across a bug (with matplotlib 0.91) that was crashing this program
# when it was using annotations.  In order to fix it, I have to overwrite
# this particular function.  Later versions of matplotlib will 
# hopefully fix this.
def new_get_verts(self):
    # the base vertices
    x1, y1 = self.xytip
    x2, y2 = self.xybase
    k1 = self.width*self.dpi.get()/72./2.
    k2 = self.headwidth*self.dpi.get()/72./2.
    xb1, yb1, xb2, yb2 = self.getpoints(x1, y1, x2, y2, k1)

    # a point on the segment 20% of the distance from the tip to the base
    theta = math.atan2(y2-y1, x2-x1)
    #r = math.sqrt((y2-y1)**2. + (x2-x1)*2.)     # <--- bug is here!
    r = math.sqrt((y2-y1)**2. + (x2-x1)**2.)     # <--- correct version
    xm = x1 + self.frac * r * math.cos(theta)
    ym = y1 + self.frac * r * math.sin(theta)
    xc1, yc1, xc2, yc2 = self.getpoints(x1, y1, xm, ym, k1)
    xd1, yd1, xd2, yd2 = self.getpoints(x1, y1, xm, ym, k2)


    xs = self.convert_xunits([xb1, xb2, xc2, xd2, x1, xd1, xc1])
    ys = self.convert_yunits([yb1, yb2, yc2, yd2, y1, yd1, yc1])
    return zip(xs, ys)

matplotlib.patches.YAArrow.get_verts = new_get_verts

# SOURCE:
# need to overwrite another function that I was getting errors from.
# the original didn't take care of vertical lines, so I was getting
# divide by zero errors when the slope was calculated
oldFunc = matplotlib.patches.YAArrow.getpoints
def new_getpoints(self, x1, y1, x2, y2, k):
    """
    for line segment defined by x1,y1 and x2,y2, return the points on
    the line that is perpendicular to the line and intersects x2,y2
    and the distance from x2,y2 of the returned points is k
    """
    epsilon = 0.000001
    x1, y1, x2, y2, k = map(float, (x1,y1,x2,y2,k))
    if abs(x2 - x1) <= epsilon:
        # this is pretty much a vertical line, so easy to figure out
        y3a = y3b = y2
        x3a = x2 + k
        x3b = x2 - k
        return x3a, y3a, x3b, y3b
    else:
        try:
            return oldFunc(self, x1, y1, x2, y2, k)
        except ZeroDivisionError, e:
            # also a vertical line
            y3a = y3b = y2
            x3a = x2 + k
            x3b = x2 - k
            return x3a, y3a, x3b, y3b

matplotlib.patches.YAArrow.getpoints = new_getpoints

class GraphApp(wx.App):
    """The class that is the main application"""

    def OnInit(self):
        self.frame = MainFrame()
        self.show = self.frame.show
        if self.frame.show:
            self.frame.CenterOnScreen()
            self.frame.Show()
            self.SetTopWindow(self.frame)
        return True

def main():
    app = GraphApp(False)
    if app.show:
        app.MainLoop()
        pass
    else:
        print misc.mainWindowTitle, 'has exited'

if __name__ == "__main__":
    main()

