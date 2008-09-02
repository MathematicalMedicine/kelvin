"""Contains functions that loads a .grx format file. Uses xml to store the graph information."""

import types

from config import graphInfo
from config import misc
from io import names
from widgets.geneMarkers import GeneMarkers
from widgets.geneMarkers import GeneMarkersGroup
from widgets.line import Line
from widgets.lineSet import LineSet

def toBool(x, default=None):
    """Given a string, convert it to a boolean value"""

    if type(x) in types.StringTypes:
        if x.lower() == 'true':
            return True
        elif x.lower() == 'false':
            return False
    else:
        return default

def getChild(parent, name):
    """Return the first child with a certain name. Returns None if it does not
    exist.
    
    """

    children = parent.getElementsByTagName(name)
    if len(children) > 0:
        return children[0]
    else:
        return None

def getText(parent, nodeName):
    """Given a parent node and a name of a child node containing text,
    return the text in that node.

    """

    # first search for the element node with that name
    for node in parent.childNodes:
        if node.nodeName == nodeName:
            # search for the text node attached to the element node
            for inode in node.childNodes:
                if inode.nodeType == inode.TEXT_NODE:
                    return inode.nodeValue

    # no text was found
    return None

def setValue(root, textNodeName, object, attrName, convert=lambda x: x, default=None):
    """Sets a certain value of an object based on text in an xml node. Text is 
    taken from a child text node with name of textNodeName that is a child of 
    root. Then that text is converted using the convert function (default is 
    no conversion). Finally the the object's attribute is set to the result of 
    the conversion. If no text was found, set the attribute to the default
    value, if one was supplied.
    
    """

    text = getText(root, textNodeName)
    if text is not None:
        value = convert(text)
        setattr(object, attrName, value)
    else:
        # set attribute to default value
        if default is not None:
            setattr(object, attrName, default)

def loadHLine(graph, root):
    """Given a graph object and the root of an xml dom tree, load
    information about horizontal lines.
    
    """

    hlines = getChild(root, names.hlines)
    if hlines is not None:
        # create horizontal lines
        for node in hlines.childNodes:
            if node.nodeName == names.hlineValue:
                value = float(node.childNodes[0].nodeValue)
                graph.HLines.append(value)

def loadComments(graph, root):
    """Given a graph object and the root of an xml dom tree, load
    comments that are in the side panel.
    
    """

    text = getText(root, names.comments)
    if text is not None:
        graph.comments = text

def loadConfig(graph, root):
    """Given a a graph object and the root of an xml dom tree, load
    the configuration data.
    
    """

    def toBoolInverse(x):
        return not toBool(x)
    
    config = graph.config
    config_element = getChild(root, names.config)

    if config_element is None:
        return

    setValue(config_element, names.config_title, config, 'title')
    setValue(config_element, names.config_xlabel, config, 'xlabel')
    setValue(config_element, names.config_ylabel, config, 'ylabel')
    setValue(config_element, names.config_bgcolor, config, 'bgcolor')
    setValue(config_element, names.config_showLegend, config, 'showLegend', 
                                                                        toBool)
    setValue(config_element, names.config_showLineSetsInLegend, config, 
                                                'showLineSetsInLegend', toBool)
    setValue(config_element, names.config_showLegendFrame, config, 
                                                     'showLegendFrame', toBool)
    setValue(config_element, names.config_legendPosition, config, 
                                                              'legendPosition')
    setValue(config_element, names.config_titleSize, config, 'titleSize', int)
    setValue(config_element, names.config_axisSize, config, 'axisSize', int)
    setValue(config_element, names.config_xmin, config, 'xmin', float)
    setValue(config_element, names.config_xmax, config, 'xmax', float)
    setValue(config_element, names.config_ymin, config, 'ymin', float)
    setValue(config_element, names.config_ymax, config, 'ymax', float)
    setValue(config_element, names.config_annPosPrec, config, 'annPosPrec', int)
    setValue(config_element, names.config_isAnnSameColor, config, 
                                              'isAnnSameColor', toBoolInverse)
    setValue(config_element, names.config_allowDuplicates, config, 
                                                      'allowDuplicates', toBool)

    # annPPLPrec is special case
    text = getText(config_element, names.config_annPPLPrec)
    if text is not None:
        if text.lower() == 'default':
            config.defaultAnnPPLPrec = True
        else:
            # text must be a number
            config.defaultAnnPPLPrec = False
            config.annPPLPrec = int(text)

def loadLines(mainFrame, graph, root):
    """Given a mainFrame, a graph object, and the root of an xml dom tree, load
    the lines.
    
    """

    lines_element = getChild(root, names.lines)

    if lines_element is None:
        return

    # load up all line sets
    lineSets = lines_element.getElementsByTagName(names.lineset)
    for eachSet in lineSets:
        lineSet = loadLineSet(mainFrame, graph, eachSet)
        graph.lineCollection.addLineSet(lineSet, newId=False)

    # load up all line groups
    lines_groups = lines_element.getElementsByTagName(names.line_group)
    for eachGroup in lines_groups:
        loadLineGroup(mainFrame, graph, eachGroup)

def loadLineSet(mainFrame, graph, root):
    """Given the mainFrame, a graph object and the root of a line set xml tree,
    load the line set.
    
    """

    lineSet = LineSet()
    setValue(root, names.lineset_id, lineSet, 'id', int)
    setValue(root, names.lineset_name, lineSet, 'name')
    setValue(root, names.lineset_color, lineSet, 'color')
    setValue(root, names.lineset_thickness, lineSet, 'thickness', int)
    setValue(root, names.lineset_style, lineSet, 'linestyle')
    setValue(root, names.lineset_marker, lineSet, 'marker')
    setValue(root, names.lineset_marker_size, lineSet, 'markersize', int)

    return lineSet

def loadLineGroup(mainFrame, graph, root):
    """Given the mainFrame, a graph object, and the root of a line group xml
    tree, load the line group.
    
    """

    # get all the child nodes that represent a line
    lines = root.getElementsByTagName(names.line)

    if lines is None:
        return

    if len(lines) > 0:
        # for the first line, just add it normally so a new group is formed
        firstLine = loadLine(mainFrame, graph, lines[0])
        if firstLine.id == 0:
            graph.lineCollection.addLine(firstLine, newId=True)
        else:
            graph.lineCollection.addLine(firstLine, newId=False)

        # for all subsequent lines in the same group, add them as overlapping
        # lines so they get inserted into the same group
        for lineNode in lines[1:]:
            line = loadLine(mainFrame, graph, lineNode)
            if line.id == 0:
                graph.lineCollection.addOverlappingLine(line, firstLine,
                                                                    newId=True)
            else:
                graph.lineCollection.addOverlappingLine(line, firstLine,
                                                                    newId=False)

def loadLine(mainFrame, graph, root):
    """Given the mainFrame, a graph object and the root of a line xml tree, load
    the line.
    
    """

    line = Line(mainFrame.axes)
    setValue(root, names.line_id, line, 'id', int)
    setValue(root, names.line_set_id, line, 'lineSet_id', int)
    setValue(root, names.line_name, line, 'name')
    setValue(root, names.line_number, line, 'number', int)
    setValue(root, names.line_color, line, 'color')
    setValue(root, names.line_thickness, line, 'thickness', int)
    setValue(root, names.line_style, line, 'linestyle')
    setValue(root, names.line_marker, line, 'marker')
    setValue(root, names.line_marker_size, line, 'markersize', int)

    # load chr position values
    pos = getText(root, names.line_pos)
    pos = [float(x) for x in pos.split()]
    line.pos = pos

    # load chr ppl values
    ppl = getText(root, names.line_ppl)
    ppl = [float(x) for x in ppl.split()]
    line.ppl = ppl

    # load metadata values
    metadata = []
    metadataNode = getChild(root, names.line_metadata)
    if metadataNode is not None:
        for eachChild in metadataNode.childNodes:
            if eachChild.nodeName == names.line_metadata_set:
                set = []
                for eachChild2 in eachChild.childNodes:
                    if eachChild2.nodeName == names.line_metadata_value:
                        value = eachChild2.childNodes[0].nodeValue
                        set.append(value)
                metadata.append(set)
        line.metadata = metadata

    # load gene marker info
    geneMarkers = getChild(root, names.line_gene_markers)
    if geneMarkers is not None:
        for eachChild in geneMarkers.childNodes:
            if eachChild.nodeName == names.line_gene_marker:
                # load one gene markers set
                # get the id numbers
                group_id = int(getText(eachChild, names.line_gene_marker_group))
                set_id = int(getText(eachChild, names.line_gene_marker_id))

                # now find the actual GeneMarker object with correct ids
                for eachGroup in graph.geneMarkerGroups:
                    if eachGroup.id == group_id:
                        for eachSet in eachGroup.geneMarkers:
                            if eachSet.id == set_id:
                                eachSet.setLine(line)
                                line.geneMarkers.append(eachSet)
    return line

def loadAnnotations(graph, root):
    """Given a a graph object and the root of an xml dom tree, load
    the annotations.
    
    """

    anns_element = getChild(root, names.anns)
    if anns_element is None:
        return

    # get all child nodes that are annotations
    anns_children = anns_element.getElementsByTagName(names.ann)

    for ann in anns_children:
        loadAnnotation(graph, ann)

def loadAnnotation(graph, root):
    """Given a a graph object and the root of an annotation xml
    tree, load an annotation.
    
    """

    text = getText(root, names.ann_text)
    x = float(getText(root, names.ann_x))
    y = float(getText(root, names.ann_y))
    isPeak = toBool(getText(root, names.ann_isPeak), False)
    isTrough = toBool(getText(root, names.ann_isTrough), False)
    autoLayout = toBool(getText(root, names.ann_autoLayout), True)

    # recover position of annotation if possible
    xytext = getText(root, names.ann_text_pos)
    if xytext is not None:
        xytext = xytext.split()
        xytext = float(xytext[0]), float(xytext[1])

    annotation = graph.annotationManager.createAnnotation(x, y, text)
    annotation.autoLayout = autoLayout
    if xytext is not None:
        annotation.set_position(xytext)

    # if the annotation records a peak or trough, recover the associated line
    annotation.isPeak = isPeak
    annotation.isTrough = isTrough
    if isPeak or isTrough:
        id = int(getText(root, names.ann_lineId))
        if id is not None:
            line = graph.lineCollection.getLineById(id)
            annotation.line = line

        if id is None or line is None:
            # id or line with id was not found, can't be peak or trough
            annotation.isPeak = False
            annotation.isTrough = False

def loadGeneMarkers(graph, root):
    """Given a graph object, and the root of an xml dom tree, load
    gene markers.
    
    """

    geneMarkerGroups = []

    # get all child nodes that are gene marker groups
    groupNodes = root.getElementsByTagName(names.gene_markers_group)

    for groupNode in groupNodes:
        group = loadGeneMarkersGroup(graph, groupNode)
        geneMarkerGroups.append(group)

    graph.geneMarkerGroups = geneMarkerGroups

def loadGeneMarkersGroup(graph, root):
    """Given a graph object and the root of a gene marker group
    xml tree, load the group.
    
    """

    geneMarkersGroup = GeneMarkersGroup()
    setValue(root, names.gene_markers_group_id, geneMarkersGroup, 'id', int)
    setValue(root, names.gene_markers_group_name, geneMarkersGroup, 'name')
    setValue(root, names.gene_markers_group_color, geneMarkersGroup, 'color')
    setValue(root, names.gene_markers_group_symbol, geneMarkersGroup, 'symbol')
    setValue(root, names.gene_markers_group_size, geneMarkersGroup, 'size', int)

    # load each set of gene markers
    markerSetsNode = root.getElementsByTagName(names.gene_markers)
    if markerSetsNode is not None:
        for eachMarkerSetNode in markerSetsNode:
            markerSet = loadGeneMarkerSet(graph, eachMarkerSetNode)
            geneMarkersGroup.geneMarkers.append(markerSet)

    # enforce all group config values onto the gene marker sets
    geneMarkersGroup.update()
    return geneMarkersGroup

def loadGeneMarkerSet(graph, root):
    """Given a graph object, and the root of a gene marker set
    tree, load the gene marker set.
    
    """

    # get the gene markers id
    id = getText(root, names.gene_markers_id)
    if id is None:
        # what to do here?, if there is no id, then they're just screwed
        id = 0
    else:
        id = int(id)

    # load the positions
    pos = getText(root, names.gene_markers_pos)
    pos = [float(x) for x in pos.split()]

    geneMarkers = GeneMarkers(pos=pos, line=None)
    geneMarkers.id = id

    return geneMarkers
