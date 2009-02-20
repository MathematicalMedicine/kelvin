"""A class that holds information about a token"""

class Token(object):
    """A class that holds information about a token"""
    def __init__(self, type="", str="", row=0, column=0):
        self.type = type
        self.str = str
        self.row = row
        self.column = column
        self.isConstraint = False
        self.isVariableNumber = False

    def printInfo(self):
        print "type:", self.type
        print "str:", self.str
        print "row:", self.row
        print "column:", self.column

    def update(self, newValue):
        self.str = str(newValue)
