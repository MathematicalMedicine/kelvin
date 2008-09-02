"""Contains functions that saves and writes a .grx format file. Used primarily in MainFrame class. This format uses xml to store the graph information."""

from io import names

def addText(doc, node, text):
    """Given the document node, add a text node to a certain node with the
    specified text.
    
    """

    text = str(text)
    if text == '':
        return
    textNode = doc.createTextNode(text)
    node.appendChild(textNode)

def addTextNode(doc, parent, newNodeName, text):
    """Add a new element node to parent with the given name, and add a text
    node containing the given text to the new element node.
    
    """

    newElementNode = doc.createElement(newNodeName)
    addText(doc, newElementNode, text)
    parent.appendChild(newElementNode)

def saveHLine(graph, doc):
    """Create the element node that holds information about horizontal lines"""

    hlines = doc.createElement(names.hlines)

    for n in graph.HLines:
        addTextNode(doc, hlines, names.hlineValue, n)

    return hlines

def saveComments(graph, doc):
    """Create the element node that holds information about comments"""

    # the comments mentioned are the ones entered in the 
    # side panel, next to the graph
    comments = doc.createElement(names.comments)
    addText(doc, comments, graph.comments)
    return comments

def saveConfig(graph, doc):
    """Create the element node that holds information about configuration
    information. Basically all the important stuff in the Config object
    (see config.py)
    
    """

    # create the new node, and get the config object
    configNode = doc.createElement(names.config)
    config = graph.config

    # need to do this to get current viewing limits
    config.updateValues()

    # add in text xml nodes about configuration data
    addTextNode(doc, configNode, names.config_title, config.title)
    addTextNode(doc, configNode, names.config_xlabel, config.xlabel)
    addTextNode(doc, configNode, names.config_ylabel, config.ylabel)
    addTextNode(doc, configNode, names.config_bgcolor, config.bgcolor)
    addTextNode(doc, configNode, names.config_showLegend, config.showLegend)
    addTextNode(doc, configNode, names.config_showLineSetsInLegend, 
                                                   config.showLineSetsInLegend)
    addTextNode(doc, configNode, names.config_showLegendFrame, 
                                                        config.showLegendFrame)
    addTextNode(doc, configNode, names.config_legendPosition, 
                                                         config.legendPosition)
    addTextNode(doc, configNode, names.config_titleSize, config.titleSize)
    addTextNode(doc, configNode, names.config_axisSize, config.axisSize)
    addTextNode(doc, configNode, names.config_xmin, config.xmin)
    addTextNode(doc, configNode, names.config_xmax, config.xmax)
    addTextNode(doc, configNode, names.config_ymin, config.ymin)
    addTextNode(doc, configNode, names.config_ymax, config.ymax)
    addTextNode(doc, configNode, names.config_allowDuplicates, 
                                                   config.allowDuplicates)

    if config.defaultAnnPPLPrec:
        addTextNode(doc, configNode, names.config_annPPLPrec, 'default')
    else:
        addTextNode(doc, configNode, names.config_annPPLPrec, config.annPPLPrec)
    addTextNode(doc, configNode, names.config_annPosPrec, config.annPosPrec)
    addTextNode(doc, configNode, names.config_isAnnSameColor, 
                                                     not config.isAnnSameColor)
    return configNode

def saveAnnotations(graph, doc):
    """Create the element node that holds information about annotations"""

    annotationsNode = doc.createElement(names.anns)
    for eachAnn in graph.annotationManager:
        annotationNode = saveAnnotation(doc, eachAnn)
        annotationsNode.appendChild(annotationNode)

    return annotationsNode

def saveAnnotation(doc, annotation):
    """Create the element node that holds information about one annotation"""

    annotationNode = doc.createElement(names.ann)

    addTextNode(doc, annotationNode, names.ann_x, annotation.x)
    addTextNode(doc, annotationNode, names.ann_y, annotation.y)
    addTextNode(doc, annotationNode, names.ann_text, annotation.text)
    pos = ' '.join([str(x) for x in annotation.get_position_axes()])
    addTextNode(doc, annotationNode, names.ann_text_pos, pos)
    addTextNode(doc, annotationNode, names.ann_isPeak, annotation.isPeak)
    addTextNode(doc, annotationNode, names.ann_isTrough, annotation.isTrough)
    addTextNode(doc, annotationNode, names.ann_autoLayout,annotation.autoLayout)

    if annotation.isPeak or annotation.isTrough:
        addTextNode(doc, annotationNode, names.ann_lineId, annotation.line.id)

    return annotationNode

def saveLines(graph, doc):
    """Create the element node that holds information about lines"""

    linesNode = doc.createElement(names.lines)

    # save information about line sets
    for eachSet in graph.lineCollection.lineSets:
        setNode = saveLineSet(doc, eachSet)
        linesNode.appendChild(setNode)

    # save information about line groups
    for eachGroup in graph.lineCollection.getGroups():
        groupNode = saveLineGroup(doc, eachGroup)
        linesNode.appendChild(groupNode)

    return linesNode

def saveLineSet(doc, lineSet):
    """Create the element node that holds information about one line set"""

    setNode = doc.createElement(names.lineset)
    addTextNode(doc, setNode, names.lineset_id, lineSet.id)
    addTextNode(doc, setNode, names.lineset_name, lineSet.name)
    addTextNode(doc, setNode, names.lineset_color, lineSet.color)
    addTextNode(doc, setNode, names.lineset_thickness, lineSet.thickness)
    addTextNode(doc, setNode, names.lineset_style, lineSet.linestyle)
    addTextNode(doc, setNode, names.lineset_marker, lineSet.marker)
    addTextNode(doc, setNode, names.lineset_marker_size, lineSet.markersize)
    return setNode

def saveLineGroup(doc, group):
    """Create the element node that holds information about one line group"""

    groupNode = doc.createElement(names.line_group)
    for eachLine in group.lines:
        lineNode = saveLine(doc, eachLine)
        groupNode.appendChild(lineNode)

    return groupNode

def saveLine(doc, line):
    """Create the element node that holds information about one line"""

    lineNode = doc.createElement(names.line)
    addTextNode(doc, lineNode, names.line_id, line.id)
    addTextNode(doc, lineNode, names.line_set_id, line.lineSet_id)
    addTextNode(doc, lineNode, names.line_name, line.name)
    addTextNode(doc, lineNode, names.line_number, line.number)
    addTextNode(doc, lineNode, names.line_color, line.color)
    addTextNode(doc, lineNode, names.line_thickness, line.thickness)
    addTextNode(doc, lineNode, names.line_style, line.linestyle)
    addTextNode(doc, lineNode, names.line_marker, line.marker)
    addTextNode(doc, lineNode, names.line_marker_size, line.markersize)

    # add gene marker information
    geneMarkersNode = doc.createElement(names.line_gene_markers)
    for eachMarkers in line.geneMarkers:
        geneMarkerNode = doc.createElement(names.line_gene_marker)
        addTextNode(doc, geneMarkerNode, names.line_gene_marker_group, 
                                                         eachMarkers.group_id)
        addTextNode(doc, geneMarkerNode, names.line_gene_marker_id, 
                                                               eachMarkers.id)
        geneMarkersNode.appendChild(geneMarkerNode)
    lineNode.appendChild(geneMarkersNode)

    # add in position data
    positions = [str(x) for x in line.pos]
    positions = ' '.join(positions)
    addTextNode(doc, lineNode, names.line_pos, positions)

    # add in ppl data
    ppls = [str(x) for x in line.ppl]
    ppls = ' '.join(ppls)
    addTextNode(doc, lineNode, names.line_ppl, ppls)

    # add in metadata
    metadataNode = doc.createElement(names.line_metadata)
    for eachSet in line.metadata:
        metadataSetNode = doc.createElement(names.line_metadata_set)
        for eachValue in eachSet:
            addTextNode(doc,metadataSetNode,names.line_metadata_value,eachValue)
        metadataNode.appendChild(metadataSetNode)
    lineNode.appendChild(metadataNode)

    return lineNode

def saveGeneMarkers(graph, doc):
    """Create the element node that holds information about gene markers"""

    geneMarkersNode = doc.createElement(names.gene_markers_groups)

    for eachGeneMarkerGroup in graph.geneMarkerGroups:
        geneMarkersGroupNode = saveGeneMarkersGroup(doc, eachGeneMarkerGroup)
        geneMarkersNode.appendChild(geneMarkersGroupNode)

    return geneMarkersNode

def saveGeneMarkersGroup(doc, group):
    """Given a group of gene markers, create the element node that holds
    information about it.
    
    """

    groupNode = doc.createElement(names.gene_markers_group)
    addTextNode(doc, groupNode, names.gene_markers_group_id, group.id)
    addTextNode(doc, groupNode, names.gene_markers_group_name, group.name)
    addTextNode(doc, groupNode, names.gene_markers_group_color, group.color)
    addTextNode(doc, groupNode, names.gene_markers_group_symbol, group.symbol)
    addTextNode(doc, groupNode, names.gene_markers_group_size, group.size)

    # add nodes about each individual gene markers (a set of gene markers for
    # one chromosome)
    for eachMarkerSet in group:
        markersNode = saveGeneMarkersSet(doc, eachMarkerSet)
        groupNode.appendChild(markersNode)

    return groupNode

def saveGeneMarkersSet(doc, markers):
    """Given a set of gene markers (gene markers that apply to one chromosome),
    return an element node that holds information about it.
    
    """

    markersNode = doc.createElement(names.gene_markers)
    addTextNode(doc, markersNode, names.gene_markers_id, markers.id)

    # add position information
    pos = [str(x) for x in markers.pos]
    pos = ' '.join(pos)
    addTextNode(doc, markersNode, names.gene_markers_pos, pos)

    return markersNode
