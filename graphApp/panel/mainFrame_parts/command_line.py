"""Contains functions dealing with the command-line commands, e.g. parsing and executing command-line commands."""

import os

from colors import hexcolor
from config import graphInfo
from panel.insertGeneMarkerDialog import getDefaultParameters \
                                                as getGeneMarkerParams
import io.graph
from io.textFileReader import FileFormatInfo, parseTextFile
import panel.mainFrame_parts.annotation as ann
from widgets.geneMarkers import loadGeneMarkerFile

# a list of all commands that can be specified multiple times
multiple_cmds = ['-save', '-save_graph', '-export', '-hline', '-chr_name',
                 '-chr_color', '-chr_linestyle', '-chr_symbol','-chr_thickness',
                 '-gene_markers', '-gene_markers_size', '-gene_markers_color',
                 '-gene_markers_symbol', '-ann_threshold', '-ann_peak', 
                 '-ann_valley', '-overlap', '-overlap_chr_num', 
                 '-overlap_chr_col', '-overlap_chr_select', 
                 '-overlap_chr_prefix', '-overlap_chr_postfix',
                 '-overlap_extra_data_col', '-overlap_chr_same_color',
                 '-overlap_chr_color_all', '-lineset_change_name', 
                 '-lineset_color', '-lineset_linestyle', '-lineset_thickness',
                 '-lineset_symbol', '-lineset_symbol_size'
                ]

# a custom error class
class CmdError(Exception):

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message

def parse_execute_commands(mainFrame, commands):
    """Given a mainFrame and the arguments given to the program, parse and
    execute all specified commands.
    
    """

    cmds = parse(commands)
    run_commands(mainFrame, cmds)

def parse(commands): 
    """Given a list of commands and arguments, parse them into a list. Also
    gets rid of duplicate commands when applicable (the first instance of a
    duplicate is kept). Returns a list of entries, where each entry is one
    command with its arguments.
    
    """

    # parse out all the commands with their arguments
    newCommands = []
    curCommand = []
    for arg in commands:
        if len(arg) > 1 and arg[0] == '-' and arg[1].isalpha():
            if len(curCommand) > 0:
                newCommands.append(tuple(curCommand))
            curCommand = [arg]
        else:
            curCommand.append(arg)

    if len(curCommand) > 0:
        newCommands.append(tuple(curCommand))

    # eliminate duplicate commands that don't allow duplicates, keep track of
    # already seen ones with a dictionary
    cmds = {}
    for entry in newCommands[:]:
        cmd = entry[0].lower()
        if cmds.has_key(cmd) and cmd not in multiple_cmds:
            # duplicate found, remove it
            newCommands.remove(entry)
        else:
            # add commands to dict to keep track of it
            # key is command, the value is its arguments
            cmds[cmd] = entry[1:]

    # rearrange some commands for special cases
    # if there is an open command, move it up to the front
    if cmds.has_key('-open'):
        for entry in newCommands:
            if entry[0].lower() == '-open':
                break
        newCommands.remove(entry)
        newCommands.insert(0, entry)

    # move all save and export commands to the end
    toMove = []
    for entry in newCommands:
        if entry[0].lower() == '-save' or entry[0].lower() == '-export':
            toMove.append(entry)
    for entry in toMove:
        newCommands.remove(entry)
    for entry in toMove:
        newCommands.append(entry)

    return newCommands

def run_commands(mainFrame, cmds):
    """Execute all commands in cmds, given the mainFrame object. cmds should be
    the ordered dictionary returned from the parse function.
    
    """

    for cmd in cmds[:]:
        # need to check if cmd is still valid, since other commands will 
        # depend on other commands, and use them up
        try:
            if cmd in cmds:
                run_command(cmd, cmds, mainFrame)
        except IndexError, e:
            # this happens when arguments to a command are missing
            print 'Warning! Error with', cmd[0]
            print 'Not enough arguments supplied'

def run_command(cmd, cmds, mainFrame):
    """Execute a command. Given is the command, the dict of all commands, and
    the mainFrame.
    
    """

    name = cmd[0].lower()
    args = cmd[1:]

    if name == '-open':
        type = args[0]
        if type == 'graph':
            io.graph.openGraph(mainFrame, args[1])
        elif type == 'grx':
            mainFrame.openGRX(args[1])
        elif type == 'text':
            openTextFile(args[1], mainFrame, cmds)

    elif name == '-save':
        mainFrame.onSaveAsGRX(graphInfo.cmdline, args[0])

    elif name == '-save_graph':
        mainFrame.onSaveAsGraph(graphInfo.cmdline, args[0])

    elif name == '-export':
        mainFrame.navTb.save(graphInfo.cmdline, args[0])

    elif name == '-noshow':
        mainFrame.show = False

    elif name == '-comment':
        comment = ' '.join(args)
        comment = comment.replace('\\n', '\n')
        comment.strip()
        mainFrame.setComments(comment)

    elif name == '-hline':
        insertHLine(mainFrame, args)

    elif name == '-title':
        title = ' '.join(args)
        title = title.replace('\\n', '\n')
        mainFrame.config.title = title
        mainFrame.config.update()

    elif name == '-xlabel':
        label = ' '.join(args)
        label = label.replace('\\n', '\n')
        mainFrame.config.xlabel = label
        mainFrame.config.update()

    elif name == '-ylabel':
        label = ' '.join(args)
        mainFrame.config.ylabel = label
        mainFrame.config.update()

    elif name == '-bgcolor':
        mainFrame.config.bgcolor = hexcolor(args[0])
        mainFrame.axes.set_axis_bgcolor(args[0])
        mainFrame.draw()

    elif name == '-title_text_size':
        setTitleTextSize(mainFrame, args[0])

    elif name == '-axis_text_size':
        setAxisTextSize(mainFrame, args[0])

    elif name == '-show_legend':
        mainFrame.config.showLegend = True
        mainFrame.updateLegend(False)

    elif name == '-show_legend_frame':
        setLegendFrame(mainFrame, args[0])

    elif name == '-legend_show_linesets':
        setLegendShowLineSets(mainFrame, args[0])

    elif name == '-legend_location':
        location = ' '.join(args)
        setLegendLocation(mainFrame, location)

    elif name == '-ann_ppl_precision':
        setPPLPrecision(mainFrame, args[0])

    elif name == '-ann_pos_precision':
        setPosPrecision(mainFrame, args[0])

    elif name == '-ann_same_color':
        mainFrame.config.isAnnSameColor = False
        ann.updateAnnotations(mainFrame)

    elif name == '-no_duplicate_numbers':
        setNoDuplicateNumbers(mainFrame)

    elif name == '-lineset_change_name':
        setLineSetName(mainFrame, args[0], args[1])

    elif name == '-lineset_color':
        setLineSetColor(mainFrame, args[0], args[1])

    elif name == '-lineset_linestyle':
        setLineSetStyle(mainFrame, args[0], args[1])

    elif name == '-lineset_thickness':
        setLineSetThickness(mainFrame, args[0], args[1])

    elif name == '-lineset_symbol':
        setLineSetSymbol(mainFrame, args[0], args[1])

    elif name == '-lineset_symbol_size':
        setLineSetSymbolSize(mainFrame, args[0], args[1])

    elif name == '-chr_name':
        setLineName(mainFrame, args[0], args[1])

    elif name == '-chr_color':
        setLineColor(mainFrame, args[0], args[1])

    elif name == '-chr_linestyle':
        style = ' '.join(args[1:])
        setLineStyle(mainFrame, args[0], style)

    elif name == '-chr_thickness':
        setLineThickness(mainFrame, args[0], args[1])

    elif name == '-chr_symbol':
        symbol = ' '.join(args[1:])
        setLineSymbol(mainFrame, args[0], symbol)

    elif name == '-chr_symbol_size':
        setLineSymbolSize(mainFrame, args[0], args[1])

    elif name == '-view_limits':
        setViewingLimits(mainFrame, args[0], args[1], args[2], args[3], 
                         args[4], args[5])

    elif name == '-xlimits':
        setXViewingLimits(mainFrame, args[0], args[1], args[2], args[3])

    elif name == '-ylimits':
        setYViewingLimits(mainFrame, args[0], args[1])

    elif name == '-zoom_chr':
        zoomToChr(mainFrame, args[0])

    elif name == '-gene_markers':
        insertGeneMarkers(mainFrame, args)

    elif name == '-gene_markers_size':
        setGeneMarkerSize(mainFrame, args[0], args[1])

    elif name == '-gene_markers_color':
        setGeneMarkerColor(mainFrame, args[0], args[1])

    elif name == '-gene_markers_symbol':
        setGeneMarkerSymbol(mainFrame, args[0], args[1])

    elif name == '-overlap':
        openOverlapFile(mainFrame, args, cmd, cmds)

    elif name == '-ann_threshold':
        annotateThreshold(mainFrame, args[0], args[1])

    elif name == '-ann_peak':
        annotatePeak(mainFrame, args[0], args[1], args[2], args[3]) 

    elif name == '-ann_valley':
        annotateValley(mainFrame, args[0], args[1], args[2], args[3]) 

    elif name == '-graph2grx':
        convertGraph2GRX(mainFrame, args)

    elif name in ['-pos_col', '-ppl_col', '-chr_num', '-chr_col', '-chr_select',
                  '-chr_prefix', '-chr_postfix', '-extra_data_col',
                  '-chr_same_color', '-chr_color_all', '-overlap_chr_num',
                  '-overlap_chr_col', '-overlap_chr_select',
                  '-overlap_chr_prefix', '-overlap_chr_postfix',
                  '-overlap_extra_data_col', '-overlap_chr_same_color',
                  '-overlap_chr_color_all', '-overlap_lineset_name', 
                  '-lineset_name', '-header_lines']:
        # commands that help configure other commands, 
        # they do nothing themeselves
        pass
    else:
        print 'Warning! Command \'%s\' not known command! Will be ignored.'%name

def openTextFile(filename, mainFrame, cmds):
    """Open a text file containing chromosome data"""

    # first test to see if the file exists
    try:
        f = open(filename, 'r')
    except IOError, e:
        print "Error opening file '%s' : %s" % (filename, e.strerror)
        print 'File was not opened'
        return
    else:
        info = getFileFormat(mainFrame, cmds)
        filename = os.path.basename(f.name)

        # parse the file, and create the right data structures
        parseTextFile(mainFrame, f, filename, False, info)
        mainFrame.drawLines()

        f.close()

def getFileFormat(mainFrame, cmds):
    """Given the mainFrame and the dict of commands, return a FileFormatInfo
    object that has all the file format information needed.
    
    """

    def setRequired(infoName, commandName, cmds, info, convert):
        # Given a command, the dict of commands, the info object, and a 
        # conversion function to convert the function, set a required value.
        # will complain if not present.
        found = False
        for entry in cmds:
            if entry[0] == commandName:
                info[infoName] = convert(entry[1:])
                found = True
        if not found:
            print 'Error! \'%s\' is not specified!' % commandName

    def setOptional(infoName, commandName, cmds, info, convert):
        # Given a command, the dict of commands, the info object, and a 
        # conversion function to convert the function, set an optional value.
        # if command is not found, will not complain.
        for entry in cmds:
            if entry[0] == commandName:
                info[infoName] = convert(entry[1:])

    def toInt(x):
        return int(x[0])

    def toBool(x):
        if x[0] == 'on':
            return True
        else:
            return False

    def toStr(x):
        return str(x[0])

    def toList(x):
        return [int(a) for a in x]

    def toColor(x):
        return hexcolor(x[0])

    info = FileFormatInfo()
    info['color'] = mainFrame.getNextColor()

    # required information about file format from other commands
    setRequired('pos_col', '-pos_col', cmds, info, toInt)
    setRequired('ppl_col', '-ppl_col', cmds, info, toInt)

    # get optional information
    setOptional('header_lines', '-header_lines', cmds, info, toInt)
    setOptional('chr_num', '-chr_num', cmds, info, toInt)
    setOptional('chr_col', '-chr_col', cmds, info, toInt)
    setOptional('import_list', '-chr_select', cmds, info, toList)
    setOptional('prefix', '-chr_prefix', cmds, info, toStr)
    setOptional('postfix', '-chr_postfix', cmds, info, toStr)
    setOptional('metadata', '-extra_data_col', cmds, info, toList)
    setOptional('same_color', '-chr_same_color', cmds, info, toBool)
    setOptional('color', '-chr_color_all', cmds, info, toColor)
    setOptional('line_set_name', '-lineset_name', cmds, info, toStr)

    # set other file format data that do not translate into commands
    if info['chr_num'] > 0:
        info['multiple_chr'] = False

    if len(info['import_list']) > 0:
        info['import_all'] = False

    return info

def insertHLine(mainFrame, values):
    """Insert horizontal lines at the specified y-values"""

    for eachValue in values:
        try:
            value = float(eachValue)
        except ValueError:
            print 'Warning! Error with -hline'
            print 'hline value of \'%s\' is not a valid number.' % eachValue
        else:
            mainFrame.insertHLine(graphInfo.cmdline, value)

def setTitleTextSize(mainFrame, size):
    """Set the size of the title text"""

    try:
        size2 = int(size)
    except ValueError, e:
        print 'Warning! Error with -title_text_size'
        print 'Value of \'%s\' is not a valid integer' % size
    else:
        mainFrame.config.titleSize = size2
        mainFrame.config.update()

def setAxisTextSize(mainFrame, size):
    """Set the size of the title text"""

    try:
        size2 = int(size)
    except ValueError, e:
        print 'Warning! Error with -axis_text_size'
        print 'Value of \'%s\' is not a valid integer' % size
    else:
        mainFrame.config.axisSize = size2
        mainFrame.config.update()

def setLegendFrame(mainFrame, show):
    """Sets the legend frame visible or not"""

    show2 = show.lower()
    if show2 == 'on':
        show3 = True
    elif show2 == 'off':
        show3 = False
    else:
        print 'Warning! Error with -show_legend_frame'
        print 'Value of \'%s\' is not correct argument.' % show
        print 'Use \'on\' or \'off\''
        return

    mainFrame.config.showLegendFrame = show3
    mainFrame.updateLegend(False)

def setLegendShowLineSets(mainFrame, show):
    """Sets whether the legend lists line sets or lists all lines"""

    show2 = show.lower()
    if show2 == 'on':
        show3 = True
    elif show2 == 'off':
        show3 = False
    else:
        print 'Warning! Error with -legend_show_linesets'
        print 'Value of \'%s\' is not correct argument.' % show
        print 'Use \'on\' or \'off\''
        return

    mainFrame.config.showLineSetsInLegend = show3
    mainFrame.updateLegend(False)

def setLegendLocation(mainFrame, location):
    """Sets where the legend is shown"""

    all_loc = ['upper right' , 'upper left', 'upper center', 'lower right',  
                'lower left', 'lower center', 'center left',  'center right', 
                'right', 'center']
    location2 = location.lower()
    if location2 in all_loc:
        mainFrame.config.legendPosition = location2
        mainFrame.updateLegend(False)
    else:
        loc = ', '.join(all_loc)
        print 'Warning! Error with -legend_location'
        print 'Value of \'%s\' is not a valid location.' % location
        print 'Valid locations are: %s' % loc

def setPPLPrecision(mainFrame, precision):
    """Set the number of decimal places to display for the ppl text in
    annotations.
    
    """

    if precision.lower() == 'default':
        # set it to use the default precision
        mainFrame.config.defaultAnnPPLPrec = True
    else:
        try:
            precision2 = int(precision)
        except ValueError, e:
            print 'Warning! Error with -ann_ppl_precision'
            print 'Value of \'%s\' is not a valid integer' % precision
        else:
            mainFrame.config.defaultAnnPPLPrec = False
            mainFrame.config.annPPLPrec = precision2
            ann.updateAnnotations(mainFrame)

def setPosPrecision(mainFrame, precision):
    """Set the number of decimal places to display for the pos text in
    annotations.
    
    """

    try:
        precision2 = int(precision)
    except ValueError, e:
        print 'Warning! Error with -ann_pos_precision'
        print 'Value of \'%s\' is not a valid integer' % precision
    else:
        mainFrame.config.annPosPrec = precision2
        ann.updateAnnotations(mainFrame)

def setNoDuplicateNumbers(mainFrame):
    """Set it so no duplicate numbers are allowed in the chr labels underneath
    the graph.
    
    """

    mainFrame.config.allowDuplicates = False
    for eachGroup in mainFrame.lines.getGroups():
        eachGroup.allowDuplicates = False
        eachGroup.redoLabel()
    mainFrame.drawLabels()

def setLineSetName(mainFrame, oldname, newname):
    """Change the name of a line set"""

    for eachSet in mainFrame.lines.lineSets:
        if oldname == eachSet.name:
            eachSet.name = newname
            mainFrame.updateLegend(False)
            return

    # line set with oldname was not found
    print 'Warning! Error with -lineset_change_name'
    print 'Line set with name \'%s\' was not found' % oldname

def setLineSetColor(mainFrame, name, color):
    """Change the color of all lines in a line set"""

    for eachSet in mainFrame.lines.lineSets:
        if eachSet.name == name:
            eachSet.set_color(hexcolor(color))
            mainFrame.updateLegend(False)
            return

    # line set with name was not found
    print 'Warning! Error with -lineset_color'
    print 'Line set with name \'%s\' was not found' % name

def setLineSetStyle(mainFrame, name, style):
    """Change the line style of all lines in a line set"""

    for eachSet in mainFrame.lines.lineSets:
        if eachSet.name == name:
            if graphInfo.lineStyleDescrToTag.has_key(style):
                tag = graphInfo.lineStyleDescrToTag[style]
                eachSet.set_linestyle(tag)
                mainFrame.updateLegend(False)
            else:
                print 'Warning! Error with -lineset_linestyle'
                print 'Style: \'%s\' is not a valid style' % style
            return

    # line set with name was not found
    print 'Warning! Error with -lineset_linestyle'
    print 'Line set with name \'%s\' was not found' % name

def setLineSetThickness(mainFrame, name, thickness):
    """Change the color of a line"""

    for eachSet in mainFrame.lines.lineSets:
        if eachSet.name == name:
            try:
                thickness2 = int(thickness)
            except ValueError, e:
                print 'Warning! Error with -lineset_thickness'
                print 'Value \'%s\' is not a valid integer' % thickness
            else:
                eachSet.set_linewidth(thickness2)
                mainFrame.updateLegend(False)
            return

    # line set with name was not found
    print 'Warning! Error with -lineset_thickness'
    print 'Line set with name \'%s\' was not found' % name

def setLineSetSymbol(mainFrame, name, symbol):
    """Change the line symbol of all lines in a line set"""

    for eachSet in mainFrame.lines.lineSets:
        if eachSet.name == name:
            if graphInfo.lineSymbolDescrToTag.has_key(symbol):
                tag = graphInfo.lineSymbolDescrToTag[symbol]
                eachSet.set_marker(tag)
                mainFrame.updateLegend(False)
            else:
                print 'Warning! Error with -lineset_symbol'
                print 'Symbol: \'%s\' is not a valid symbol' % symbol
            return

    print 'Warning! Error with -lineset_symbol'
    print 'Line with name \'%s\' was not found' % name

def setLineSetSymbolSize(mainFrame, name, size):
    """Change the size of the line symbol of all lines in a line set"""

    for eachSet in mainFrame.lines.lineSets:
        if eachSet.name == name:
            try:
                size2 = int(size)
            except ValueError, e:
                print 'Warning! Error with -lineset_symbol_size'
                print 'Value \'%s\' is not a valid integer' % size
            else:
                eachSet.set_markersize(size2)
                mainFrame.updateLegend(False)
            return

    print 'Warning! Error with -lineset_symbol_size'
    print 'Line with name \'%s\' was not found' % name

def setLineName(mainFrame, oldname, newname):
    """Change the name of a line"""

    line, group = mainFrame.lines.getLineAndGroup(oldname)
    if line is not None:
        mainFrame.changeLineName(line, group, newname)
    else:
        print 'Warning! Error with -chr_name'
        print 'Line with name \'%s\' was not found' % oldname

def setLineColor(mainFrame, name, color):
    """Change the color of a line"""

    color = hexcolor(color)
    line, group = mainFrame.lines.getLineAndGroup(name)
    if line is not None:
        line.set_color(color)
        mainFrame.updateLegend(False)
    else:
        print 'Warning! Error with -chr_color'
        print 'Line with name \'%s\' was not found' % name

def setLineStyle(mainFrame, name, style):
    """Change the line style of a line"""

    line, group = mainFrame.lines.getLineAndGroup(name)
    if line is not None:
        if graphInfo.lineStyleDescrToTag.has_key(style):
            tag = graphInfo.lineStyleDescrToTag[style]
            line.set_linestyle(tag)
            mainFrame.updateLegend(False)
        else:
            print 'Warning! Error with -chr_linestyle'
            print 'Style: \'%s\' is not a valid style' % style
    else:
        print 'Warning! Error with -chr_linestyle'
        print 'Line with name \'%s\' was not found' % name

def setLineThickness(mainFrame, name, thickness):
    """Change the thickness of a line"""

    line, group = mainFrame.lines.getLineAndGroup(name)
    if line is not None:
        try:
            thickness2 = int(thickness)
        except ValueError, e:
            print 'Warning! Error with -chr_thickness'
            print 'Value \'%s\' is not a valid integer' % thickness
        else:
            line.set_linewidth(thickness2)
            mainFrame.updateLegend(False)
    else:
        print 'Warning! Error with -chr_thickness'
        print 'Line with name \'%s\' was not found' % name

def setLineSymbol(mainFrame, name, symbol):
    """Change the line symbol of a line"""

    line, group = mainFrame.lines.getLineAndGroup(name)
    if line is not None:
        if graphInfo.lineSymbolDescrToTag.has_key(symbol):
            tag = graphInfo.lineSymbolDescrToTag[symbol]
            line.set_marker(tag)
            mainFrame.updateLegend(False)
        else:
            print 'Warning! Error with -chr_symbol'
            print 'Symbol: \'%s\' is not a valid symbol' % symbol
    else:
        print 'Warning! Error with -chr_symbol'
        print 'Line with name \'%s\' was not found' % name

def setLineSymbolSize(mainFrame, name, size):
    """Change the size of the line symbol of a line"""

    line, group = mainFrame.lines.getLineAndGroup(name)
    if line is not None:
        try:
            size2 = int(size)
        except ValueError, e:
            print 'Warning! Error with -chr_symbol_size'
            print 'Value \'%s\' is not a valid integer' % size
        else:
            line.set_markersize(size2)
            mainFrame.updateLegend(False)
    else:
        print 'Warning! Error with -chr_symbol_size'
        print 'Line with name \'%s\' was not found' % name

def setViewingLimits(mainFrame, chr1, pos1, chr2, pos2, ymin, ymax):
    """Set the viewing limits along the x-axis for the graph"""

    try:
        pos1 = float(pos1)
        pos2 = float(pos2)
        ymin = float(ymin)
        ymax = float(ymax)
    except ValueError, e:
        print 'Warning! Error with -view_limits'
        print 'Value \'%s\' is not a valid number.' % e.message.split()[-1]
        return

    xmin = getCoordinate(mainFrame, chr1, float(pos1))
    xmax = getCoordinate(mainFrame, chr2, float(pos2))

    if xmin is None or xmax is None:
        print 'Warning! Error with -view_limits'
        return

    if xmin > xmax:
        temp = xmin
        xmin = xmax
        xmax = temp

    if ymin > ymax:
        temp = ymin
        ymin = ymax
        ymax = temp

    mainFrame.setViewLimits((xmin, xmax), (ymin, ymax))

def setXViewingLimits(mainFrame, chr1, pos1, chr2, pos2):
    """Set the viewing limits along the x-axis for the graph"""

    try:
        pos1 = float(pos1)
        pos2 = float(pos2)
    except ValueError, e:
        print 'Warning! Error with -xlimits'
        print 'Value \'%s\' is not a valid number.' % e.message.split()[-1]
        return

    xmin = getCoordinate(mainFrame, chr1, float(pos1))
    xmax = getCoordinate(mainFrame, chr2, float(pos2))

    if xmin is None or xmax is None:
        print 'Warning! Error with -xlimits'
        return

    if xmin > xmax:
        temp = xmin
        xmin = xmax
        xmax = temp

    mainFrame.setViewLimits((xmin, xmax), None)

def getCoordinate(mainFrame, chrname, pos):
    """Given a chromosome name and a position, return the x-coordinate value
    of the position in data space.
    
    """

    line, group = mainFrame.lines.getLineAndGroup(chrname)

    # error checking
    if line is None:
        print 'Warning! No chromosome found with name \'%s\'' % chrname
        return None
    if pos > line.pos[-1] or pos < line.pos[0]:
        print 'Warning! Position: %f is outside the range of chromosome %s' % \
                                                          (float(pos), chrname)
        return None
    
    x = group.xmin + pos - group.start
    return x

def setYViewingLimits(mainFrame, ymin, ymax):
    """Set the viewing limits along the y-axis for the graph"""

    try:
        ymin2 = float(ymin)
        ymax2 = float(ymax)
    except ValueError, e:
        print 'Warning! Error with -ylimits'
        print 'Value \'%s\' is not a valid number' % e.message.split()[-1]
    else:
        # make sure min is smaller than max
        if ymin2 > ymax2:
            temp = ymin2
            ymin2 = ymax2
            ymax2 = temp
        mainFrame.setViewLimits(None, (ymin2, ymax2))

def zoomToChr(mainFrame, chrname):
    """Given a name of a chromosome, set the x-limits such that the chromosome
    covers the entire view.
    
    """

    line, group = mainFrame.lines.getLineAndGroup(chrname)

    if line is None:
        print 'Warning! Error with -zoom_chr'
        print 'No chromosome found with name \'%s\'.' % chrname
        return

    xmin = getCoordinate(mainFrame, chrname, line.pos[0])
    xmax = getCoordinate(mainFrame, chrname, line.pos[-1])

    mainFrame.setViewLimits((xmin, xmax), None)

def insertGeneMarkers(mainFrame, args):
    """Load a set of gene markers from a file"""

    # assume that args are in the order: <filename> <chr col> <pos col> 
    # <size> <color> <symbol> <name>

    try:
        length = len(args)

        # the first two args should be the chr column and pos column, check
        # to see if they are there and are actually numbers
        if length < 3:
            s = 'Command requires at least three arguments,' + \
                ' <filename> <chr col> <pos col>'
            raise CmdError(s)
        else:
            # gather required arguments. note that column indices are zero based
            filename = args[0]
            chr_col = int(args[1]) - 1
            pos_col = int(args[2]) - 1

            info = getGeneMarkerParams(filename)
            info['chr_col'] = chr_col
            info['pos_col'] = pos_col

            # read in optional arguments
            if length > 3:
                info['size'] = float(args[3])
            if length > 4:
                info['color'] = hexcolor(args[4])
            if length > 5:
                symbol = args[5]
                if graphInfo.lineSymbolDescrToTag.has_key(symbol):
                    info['symbol'] = graphInfo.lineSymbolDescrToTag[symbol]
                else:
                    raise CmdError("Symbol: '%s' is not a valid symbol" \
                                                                    % symbol)
            if length > 6:
                info['name'] = args[6]

    except ValueError, e:
        print 'Warning! Error with -gene_markers'
        print 'Value \'%s\' is not a valid integer' % e.message.split()[-1]
    except CmdError, e:
        print 'Warning! Error with -gene_markers'
        print e

    try:
        # load the file
        loadGeneMarkerFile(mainFrame, info)
    except Exception, e:
        print 'Warning! Error with -gene_markers'
        m = "Error loading gene markers: %s" % str(e)
        print m

def getGeneMarkerGroup(mainFrame, name):
    """Given a gene marker group name, return the marker group"""

    for eachGroup in mainFrame.geneMarkerGroups:
        if eachGroup.name == name:
            return eachGroup

    # group was not found
    return None

def setGeneMarkerSize(mainFrame, name, size):
    """Given a name of a group of gene markers, sets the size of the group"""

    try:
        markers = getGeneMarkerGroup(mainFrame, name)
        if markers == None:
            raise CmdError("No set of markers found with the name '%s'"  % name)
        markers.set_markersize(int(size))
    except ValueError, e:
        print 'Warning! Error with -gene_markers_size'
        print 'Value \'%s\' is not a valid integer' % e.message.split()[-1]
    except CmdError, e:
        print 'Warning! Error with -gene_markers_size'
        print e

def setGeneMarkerColor(mainFrame, name, color):
    """Given a name of a group of gene markers, sets the color of the group"""

    try:
        markers = getGeneMarkerGroup(mainFrame, name)
        if markers == None:
            raise CmdError("No set of markers found with the name '%s'" % name)
        markers.set_color(hexcolor(color))
    except CmdError, e:
        print 'Warning! Error with -gene_markers_color'
        print e

def setGeneMarkerSymbol(mainFrame, name, symbol):
    """Given a name of a group of gene markers, sets the symbol of the group"""

    try:
        markers = getGeneMarkerGroup(mainFrame, name)
        if markers == None:
            raise CmdError("No set of markers found with the name '%s'" % name)
        if graphInfo.lineSymbolDescrToTag.has_key(symbol):
            tag = graphInfo.lineSymbolDescrToTag[symbol]
            markers.set_marker(tag)
        else:
            raise CmdError("Symbol: '%s' is not a valid symbol"  % symbol)
    except ValueError, e:
        print 'Warning! Error with -gene_markers_symbol'
        print 'Value \'%s\' is not a valid integer' % e.message.split()[-1]
    except CmdError, e:
        print 'Warning! Error with -gene_markers_symbol'
        print e

def openOverlapFile(mainFrame, args, cmd, cmds):
    """Open a text file with chromosome data, and overlap them with 
    existing chromosomes.
    
    """

    if len(args) < 3:
        print 'Warning! Error with -overlap'
        print 'Command needs at least 3 arguments'
        return

    # first test to see if the file exists
    try:
        filename = args[0]
        f = open(filename, 'r')
    except IOError, e:
        print 'Warning! Error with -overlap'
        print "Error opening file '%s' : %s" % (filename, e.strerror)
        print 'File was not opened'
        return

    try:
        info = getOverlapFileFormat(mainFrame, args, cmd, cmds)
    except ValueError, e:
        print 'Warning! Error with -overlap'
        print 'Value \'%s\' is not a valid integer' % e.message.split()[-1]
        return
    except CmdError, e:
        print 'Warning! Error with -overlap'
        print e
        return

    filename = os.path.basename(f.name)

    # parse the file, and create the right data structures
    parseTextFile(mainFrame, f, filename, True, info)
    mainFrame.drawLines()

    f.close()

def getOverlapFileFormat(mainFrame, args, cmd, cmds):
    """Get the parameters for a file being opened that will overlap chr's"""


    def setOptional(infoName, commandName, cmd, cmds, info, convert):
        # Given a command, the list of commands, the info object, and a 
        # conversion function to convert the function, set an optional value.
        # if command is not found, will not complain. Will search the cmd list
        # in order beginning at the -overlap command, and will stop at
        # the next -overlap command
        i = cmds.index(cmd)
        i = min(len(cmds), i+1)
        for entry in cmds[i:]:
            if entry[0] == '-overlap':
                break
            if entry[0] == commandName:
                info[infoName] = convert(commandName, entry[1:])
                break

    def toInt(commandName, x):
        try:
            return int(x[0])
        except ValueError, e:
            print 'Warning! Error with', commandName
            print "Value '%s' is not a valid integer", x[0]

    def toStr(commandName, x):
        return str(x[0])

    def toList(commandName, x):
        return [int(a) for a in x]

    def toColor(commandName, x):
        return hexcolor(x[0])

    def toBool(x):
        if x[0] == 'on':
            return True
        else:
            return False

    info = FileFormatInfo()
    info['overlap'] = True
    info['color'] = mainFrame.getNextColor()

    # get required information about file format from other commands
    info['pos_col'] = int(args[1])
    info['ppl_col'] = int(args[2])
    if len(args) > 3:
        info['overlap_index'] = getLineIndex(mainFrame, args[3])
        if info['overlap_index'] is None:
            raise CmdError("No chromosome found with name '%s'" % name)
    else:
        info['overlap_index'] = 0

    # get optional information
    setOptional('chr_num', '-overlap_chr_num', cmd, cmds, info, toInt)
    setOptional('chr_col', '-overlap_chr_col', cmd, cmds, info, toInt)
    setOptional('import_list', '-overlap_chr_select', cmd, cmds, info, toList)
    setOptional('prefix', '-overlap_chr_prefix', cmd, cmds, info, toStr)
    setOptional('postfix', '-overlap_chr_postfix', cmd, cmds, info, toStr)
    setOptional('metadata', '-overlap_extra_data_col', cmd, cmds, info, toList)
    setOptional('same_color', '-overlap_chr_same_color', cmd,cmds, info, toBool)
    setOptional('color', '-overlap_chr_color_all', cmd, cmds, info, toColor)
    setOptional('line_set_name', '-overlap_lineset_name', cmd, cmds, info,toStr)

    # set other file format data that do not translate into commands
    if info['chr_num'] > 0:
        info['multiple_chr'] = False

    if len(info['import_list']) > 0:
        info['import_all'] = False

    return info

def getLineIndex(mainFrame, name):
    """Given a line name, return the line's index (see line.py for 
    explanation about index number).
    
    """

    names = mainFrame.lines.getNames()

    if name in names:
        return names.index(name)
    else:
        return None

def annotateThreshold(mainFrame, type, threshold):
    """Given what type of thresholding and the threshold value,
    create appropriate annotations that meet the threshold.
    
    """

    try:
        threshold = float(threshold)
    except ValueError, e:
        print 'Warning! Error with -ann_threshold'
        print "Value '%s' is not a valid number" % threshold
        return

    if type.lower() == 'peak':
        ann_peaks = True
    elif type.lower() == 'valley':
        ann_peaks = False
    else:
        # not a valid argument
        print 'Warning! Error with -ann_threshold'
        print "Value '%s' is not a valid type. Use 'peak' or 'valley'" % type
        return

    mainFrame.annotationManager.annotateThreshold(threshold, ann_peaks)
    mainFrame.draw()

def annotatePeak(mainFrame, name1, pos1, name2, pos2):
    """Given a region, annotate the highest peak in that region"""

    # get the x-values specified by the region
    try:
        pos1 = float(pos1)
        pos2 = float(pos2)
    except ValueError, e:
        print 'Warning! Error with -ann_peak'
        print 'Value \'%s\' is not a valid number' % e.message.split()[-1]
        return

    x1 = getCoordinate(mainFrame, name1, pos1)
    x2 = getCoordinate(mainFrame, name2, pos2)

    if x1 is None or x2 is None:
        print 'Warning! Error with -ann_peak'
        return

    # need to get y-values of box, use highest and lowest ppl values
    y1 = mainFrame.lines.getLine(name1).ppl[0]
    y2 = mainFrame.lines.getLine(name1).ppl[0]
    ppl = []
    inRange = False
    for eachLine in mainFrame.lines:
        if eachLine.name == name1:
            inRange = True
        if inRange:
            ppl.extend(eachLine.ppl)
        if eachLine.name == name2:
            break

    y1 = min(ppl)
    y2 = max(ppl)

    mainFrame.annotationManager.annotateMax(x1, y1, x2, y2)

def annotateValley(mainFrame, name1, pos1, name2, pos2):
    """Given a region, annotate the highest peak in that region"""

    # get the x-values specified by the region
    try:
        pos1 = float(pos1)
        pos2 = float(pos2)
    except ValueError, e:
        print 'Warning! Error with -ann_valley'
        print 'Value \'%s\' is not a valid number' % e.message.split()[-1]
        return

    x1 = getCoordinate(mainFrame, name1, pos1)
    x2 = getCoordinate(mainFrame, name2, pos2)

    if x1 is None or x2 is None:
        print 'Warning! Error with -ann_valley'
        return

    # need to get y-values of box, use highest and lowest ppl values
    y1 = mainFrame.lines.getLine(name1).ppl[0]
    y2 = mainFrame.lines.getLine(name1).ppl[0]
    ppl = []
    inRange = False
    for eachLine in mainFrame.lines:
        if eachLine.name == name1:
            inRange = True
        if inRange:
            ppl.extend(eachLine.ppl)
        if eachLine.name == name2:
            break

    y1 = min(ppl)
    y2 = max(ppl)

    mainFrame.annotationManager.annotateMin(x1, y1, x2, y2)

def convertGraph2GRX(mainFrame, args):
    """Convert .graph files to .grx files"""

    # error check
    if len(args) < 1:
        print 'Warning! Error with -graph2grx'
        print 'Command needs at least 1 argument'
        return

    if not os.path.exists(args[0]):
        print 'Warning! Error with -graph2grx'
        print "File or directory '%s' does not exist" % args[0]
        return

    # get arguments
    convert_place = args[0]
    if len(args) > 1:
        if args[1].lower() == 'on':
            isRecursive = True
        elif args[1].lower() == 'off':
            isRecursive = False
        else:
            isRecursive = False
    else:
        isRecursive = False

    # gather all filenames, if needed
    filenames = []
    if os.path.isdir(convert_place):
        filenames.extend(getGraphFiles(convert_place, isRecursive))
    else:
        # must be a filename
        filenames.append(convert_place)

    for eachFilename in filenames:
        # change the file extension to .grx
        i = eachFilename.lower().rfind('.graph')
        if i > -1:
            newFilename = eachFilename.replace('.graph', '.grx', i)
        else:
            newFilename = eachFilename + '.grx'

        io.graph.openGraph(mainFrame, eachFilename)
        mainFrame.onSaveAsGRX(graphInfo.cmdline, newFilename)
        print 'converted %s to %s' % (eachFilename, newFilename)

    mainFrame.onClearGraph(None)
    print 'done with conversion'

def getGraphFiles(directory, isRecursive):
    """Given a directory, return all files in the directory that end in .graph.
    If isRecursive is True, then also return all files in subdirectories.
    
    """

    graphFiles = []

    def getFiles(arg, dirname, filenames):
        # callback function for the os.path.walk function. this function 
        # will be called at every directory and subdirectoy
        for filename in filenames:
            if not os.path.isdir(filename) and filename.rfind('.graph') != -1:
                graphFiles.append(os.path.join(dirname, filename))

    if isRecursive:
        os.path.walk(directory, getFiles, None)
    else:
        filenames = os.listdir(directory)
        getFiles(None, directory, filenames)

    return graphFiles








