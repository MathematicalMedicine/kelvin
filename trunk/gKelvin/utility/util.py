"""Miscellaneous functions that are used here and there"""

import time
import wx

def createEvent(action, widget, integer=None):
    """Given the event action type and widget, generate such an event"""

    evt = wx.CommandEvent(action.evtType[0], widget.GetId())
    evt.SetEventObject(widget)
    if integer is not None:
        evt.SetInt(integer)
    widget.Command(evt)

def callAncestorFunction(obj, functionName, *parameters):
    """Given an object, go up the wxPyton parent tree to find the 
    first ancestor with a function of the given name and call that 
    function with the given parameters
    
    """

    parent = obj.GetParent()
    while parent is not None and not hasattr(parent, functionName):
        parent = parent.GetParent()

    if parent is not None and hasattr(parent, functionName):
        function = getattr(parent, functionName)
        return function(*parameters)

def callAncestorTypeFunction(obj, ancestorType, functionName, *parameters):
    """Given an object, go up the wxPyton parent tree to find the first ancestor
    that is the given type and has a function of the given name. Call that 
    function with the given parameters.
    
    """

    # Note that the checking of the object type is not totally specific.
    # isinstance will return true if the object type is the specified type
    # and will also return true if the object type is a subclass of the type

    parent = obj.GetParent()
    while parent is not None and (not isinstance(parent, ancestorType) or \
                                  not hasattr(parent, functionName)):
        parent = parent.GetParent()

    if parent is not None and isinstance(parent, ancestorType) and \
       hasattr(parent, functionName):
        function = getattr(parent, functionName)
        function(*parameters)

def bindTraitTypeChange(obj, function):
    """Given an object, and a function that needs to be called when the trait
    type changes, go up the parent tree to register that function with that 
    ancestor
    
    """

    registerFunctionName = 'bindTraitTypeChange'
    callAncestorFunction(obj, registerFunctionName, function)

def bindAnalysisTypeChange(obj, function):
    """Given an object, and a function that needs to be called when the analysis
    type changes, go up the parent tree to register that function with that 
    ancestor
    
    """

    registerFunctionName = 'bindAnalysisTypeChange'
    callAncestorFunction(obj, registerFunctionName, function)

def bindLiabilityClassChange(obj, function):
    """Given an object, and a function that needs to be called when the number 
    of liability class changes changes, go up the parent tree to register that 
    function with that ancestor
    
    """

    registerFunctionName = 'bindLiabilityClassChange'
    callAncestorFunction(obj, registerFunctionName, function)

def bindValuesChange(obj, function):
    """Given an object, and a function that needs to be called when something 
    in the values panel changes, go up the parent tree to register that 
    function with that ancestor.
    
    """

    registerFunctionName = 'bindValuesChange'
    callAncestorFunction(obj, registerFunctionName, function)

def showTimings(description, function, *parameters):
    """Given a function and its parameters, time how long the function
    takes to run. Print out result using the description.
    
    """

    t1 = time.time()
    function(*parameters)
    t2 = time.time()
    print '%s took %0.3f ms' % (description, (t2-t1)*1000.0)


