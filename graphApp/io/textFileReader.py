"""Contains functions that read text input files in. Used primarily in MainFrame class"""

import os
import wx

from config import misc
from panel.fileFormatDialog import FileFormatDialog
from widgets.line import Line
from widgets.lineSet import LineSet

class FileFormatInfo(dict):
    # class that stores information about the format of a text file

    # this is a dictionary describing the parameters needed to
    # open a data file.  the dialog will gather these parameter values.
    # --- parameters ---
    # pos_col - the column index that has the position
    # ppl_col - the column index that has the ppl values
    # multiple_chr - boolean that tell whether the file contains
    #                data for one chromosome or many chromosomes
    # chr_num - if the file contains only one chromosome, the number
    #           of the chromosome being imported
    # chr_col - if multiple chromosomes are being imported, the column
    #           index that indicates the number of the chromosome, if
    #           no column in the file indicates chromosome number, then
    #           this value should be 0
    # import_all - a boolean to indicate if all chromosomes in the file
    #              should be imported
    # import_list - if import_all is False, a list of chromosome numbers
    #               that would be imported
    # overlap - a boolean indicating if the chromosomes being imported
    #           should overlap with existing ones in the graph
    # overlap_index - the line index (see line.py) of where to begin 
    #                 the overlap
    # header_lines - how many lines at the top of the file are considered 
    #                part of a header. those lines will be ignored when 
    #                parsing the file
    # same_color - a flag indicating if all chromosomes in the file
    #              should be the same color or different colors
    # color - if same_color is true, color that all lines will be set to
    # line_set_name - the name of the line set associated with the lines
    # prefix - text that should be a prefix to chromosome names
    # postfix - text that should be a postfix to chromosome names
    # metadata - columns that indicate metadata associated with data

    def __init__(self):
        self['pos_col'] = 0
        self['ppl_col'] = 0
        self['multiple_chr'] = True
        self['chr_num'] = 0
        self['chr_col'] = 0
        self['import_all'] = True
        self['import_list'] = []
        self['overlap'] = False
        self['overlap_index'] = 0
        self['header_lines'] = 0
        self['same_color'] = True
        self['color'] = 'k'
        self['line_set_name'] = ''
        self['prefix'] = ''
        self['postfix'] = ''
        self['metadata'] = []

def readTextFile(mainFrame, file, filename, overlap=False):
    """Read in a text file that contains data that needs to be graphed"""

    # note: the file argument is a file object.
    # also this function does not close the file

    info = getFileInfo(mainFrame, file, filename, overlap)
    if info is not False:
        parseTextFile(mainFrame, file, filename, overlap, info)
        mainFrame.drawLines()

def getFileInfo(mainFrame, file, filename, overlap):
    """Given stuff, get information about file format from user and return it"""

    # Note: Only when a position is less then the previous position 
    # is when a new chromosome is thought to start in the file. This
    # only applies when the chr_col is not specified

    # info is a dictionary describing the parameters needed to
    # open a data file.  the dialog will gather these parameter values.
    # see FileFormatInfo() for more information.
    info = FileFormatInfo()
    info['overlap'] = overlap
    info['color'] = mainFrame.getNextColor()

    # ask user to specifiy info about file
    dlg = FileFormatDialog(mainFrame, -1, info, file, filename, 
                           mainFrame.lines.getNames())
    dlg.CenterOnParent()
    if dlg.ShowModal() != wx.ID_OK:
        # nothing to do
        return False

    return info

def parseTextFile(mainFrame, file, filename, overlap, info):
    """Given a file and its format information, parse it and created the needed
    data structures.
    
    """

    # some helper functions
    def showWarnings():
        # if specific chromosomes were specified and they weren't found
        # in the file, warn the user about it
        if len(info['import_list']) > 0:
            chr = ''
            for eachValue in info['import_list']:
                chr = chr + str(eachValue) + ', '
            m='Chromosome(s) ' + chr[:-2] + ' were not found in the file.'
            dlg = wx.MessageDialog(mainFrame, m, "Warning", wx.OK | 
                                    wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()

    def getNumber(s):
        """Given a string, extracts the first number from that string."""

        index = -1
        for i in range(len(s)):
            if s[i].isdigit():
                index = i
                break;

        if index == -1:
            # no digit was found
            return 0

        num = ''
        l = len(s)
        while index < l and s[index].isdigit():
            num += s[index]
            index += 1

        return int(num)

    def getMetadata(line, columns):
        s = []
        for eachNum in columns:
            s.append(line[eachNum])
        return s

    # need column indexes that start at 0, so subtract 1
    pos_col = info['pos_col'] - 1
    ppl_col = info['ppl_col'] - 1
    multiple_chr = info['multiple_chr']
    import_all = info['import_all']
    chr_col = info['chr_col'] - 1
    chr_num = info['chr_num']
    same_color = info['same_color']
    color = info['color']
    line_set_name = info['line_set_name']
    prefix = info['prefix']
    postfix = info['postfix']
    use_metadata = len(info['metadata']) > 0
    header_lines = info['header_lines']
    if use_metadata:
        info['metadata'] = [x-1 for x in info['metadata']]


    # do some error checking, first check if file is totally blank
    if len(file.readlines()) == 0:
        # file is blank
        m = "Error opening file '%s', file is blank." % filename
        dlg = wx.MessageDialog(mainFrame, m,"Error!", wx.OK|wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()
        return
    else:
        file.seek(0)

    # read all lines from the file, and gather Lines in a list
    axes = mainFrame.axes
    chr = 0
    lines = []
    prevPos = None
    curPos = None
    curNum = None
    for eachLine in file.readlines()[header_lines:]:
        prevPos = curPos
        curLine = eachLine.split()
        if len(curLine) == 0 or curLine[0][0] == '#':
            # a blank line or a comment, skip the line
            continue
        curPos = float(curLine[pos_col])
        curPPL = float(curLine[ppl_col])
        metadata = getMetadata(curLine, info['metadata'])

        if (prevPos is None) or (multiple_chr and prevPos>curPos and \
           chr_col < 0) or (chr_col >= 0 and curNum != curLine[chr_col]):
            # must be start of new chromosome, start a new line
            # try to figure out correct name, number, and create Line object
            if multiple_chr and chr_col < 0:
                chr += 1
                line = Line(axes, name=prefix+str(chr)+postfix, num=chr)
            elif multiple_chr and chr_col >= 0:
                num = getNumber(curLine[chr_col])
                name = prefix+str(num)+postfix
                line = Line(axes, name=name, num=num)
            else:
                line = Line(axes, name=prefix+str(chr_num)+postfix, 
                            num=chr_num)

            if chr_col >= 0:
                curNum = curLine[chr_col]

            # set properties of new line
            if not same_color:
                color = mainFrame.getNextColor()
            line.color = color
            lines.append(line)

        line.pos.append(curPos)
        line.ppl.append(curPPL)
        if use_metadata:
            line.addMetadata(metadata)

    # now there is a list with all lines in the file, may have to
    # take some out. if specific chromosomes were indicated, create a list 
    # of those lines in the order they were specified
    if multiple_chr and not import_all:
        oldLines = lines
        newLines = []
        for eachItem in info['import_list'][:]:
            for eachLine in oldLines:
                if eachLine.number == eachItem:
                    newLines.append(eachLine)
                    info['import_list'].remove(eachItem)
        lines = newLines

    # now have list of all lines to be added, add lines to line manager
    if overlap:
        # earlier it was specified which line to start the overlap in, 
        # find the group where the line resides,
        # then add lines to successive groups
        overlap_index = info['overlap_index'] 
        group = mainFrame.lines.getGroup(overlap_index)
        groups = mainFrame.lines.getGroups()
        i = groups.index(group)
        n = len(groups)
        for eachLine in lines:
            if i < n:
                mainFrame.lines.addLineToCollection(eachLine, groups[i]) 
                i += 1
            else:
                # simply append lines to the end
                mainFrame.lines.addLine(eachLine)
    else:
        # no overlapping, append lines to the end
        for eachLine in lines:
            mainFrame.lines.addLine(eachLine)

    # create a new line set for these lines, and add set to line manager
    mainFrame.lines.addNewLineSet(name=line_set_name, lines=lines)

    showWarnings()
