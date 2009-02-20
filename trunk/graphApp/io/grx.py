"""Contains functions that reads and writes a .grx format file. Used primarily in MainFrame class. This format uses xml to store the graph information."""

import time
import xml.dom.ext
import xml.dom.minidom

from config import graphInfo
from config import misc
from io.grx_load import loadAnnotations
from io.grx_load import loadComments
from io.grx_load import loadConfig
from io.grx_load import loadGeneMarkers
from io.grx_load import loadHLine
from io.grx_load import loadLines
from io.grx_save import saveAnnotations
from io.grx_save import saveComments
from io.grx_save import saveConfig
from io.grx_save import saveGeneMarkers
from io.grx_save import saveHLine
from io.grx_save import saveLines
from widgets.graph import Graph

def save(graph, filename):
    """Save the given graph to the given filename in the .grx format.  The
    graph argument should be a Graph object (see graph.py). 
    
    """

    # create the document node and the graph element
    # the graph element is the root node of the tree
    doc = xml.dom.minidom.Document()
    graph_element = doc.createElement('graph')
    doc.appendChild(graph_element)

    # add information about horizontal lines
    hline_element = saveHLine(graph, doc)
    graph_element.appendChild(hline_element)

    # add information about comments
    comment_element = saveComments(graph, doc)
    graph_element.appendChild(comment_element)

    # add configuration information
    config_element = saveConfig(graph, doc)
    graph_element.appendChild(config_element)

    # add information about lines in the graph
    lines_element = saveLines(graph, doc)
    graph_element.appendChild(lines_element)

    # add information about annotations
    ann_element = saveAnnotations(graph, doc)
    graph_element.appendChild(ann_element)

    # add information about gene markers
    geneMarkers_element = saveGeneMarkers(graph, doc)
    graph_element.appendChild(geneMarkers_element)

    # write the xml tree to a file
    f = open(filename, 'w')
    xml.dom.ext.PrettyPrint(doc, f)
    f.close()

def load(mainFrame, filename):
    """Given a filename, load the .grx file from the 
    filename and return the graph object.  
    
    """

    # read the input file
    f = open(filename, 'r')
    doc = xml.dom.minidom.parse(f).documentElement
    f.close()

    # the new graph object
    graph = Graph(mainFrame)

    # the first child of doc is a DocumentType node, the second is the
    # graph node, which is the root of the xml tree
    #root = doc.childNodes[1]
    root = doc

    # load the main components
    # load horizontal lines
    loadHLine(graph, root)

    # load comments
    loadComments(graph, root)

    # load configuration information
    loadConfig(graph, root)

    # load gene marker information
    # Note: this has to be done before lines are loaded
    loadGeneMarkers(graph, root)

    # load line information
    # Note: must be loaded after gene markers are loaded, but before
    # annotations are loaded
    loadLines(mainFrame, graph, root)

    # load annotations
    # Note: this must be done after lines are loaded
    loadAnnotations(graph, root)

    return graph
