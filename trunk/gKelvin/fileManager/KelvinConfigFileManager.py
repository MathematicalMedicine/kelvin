"""Parse and manage a Kelvin config file"""

import os

from config import KelvinInfo
from config import misc
from fileManager.KelvinConfigParser import KelvinConfigParser
from fileManager import types
from fileManager.token import Token

class KelvinConfigFileManager(object):
    """Parse and manage a Kelvin config file"""

    def __init__(self, fileName, isNew, defaultsManager=None):
        self.priorityTags = [x[0] for x in KelvinInfo.priorityTagList]
        if isNew:
            # start a new file from scratch, using a template
            self.fileName = ''
            templatefilename = os.path.join(misc.newConfigFilePath, fileName)
            parser = KelvinConfigParser(templatefilename)
            self.lines = parser.parse()
        else:
            # open a file that's on the disk
            self.fileName = fileName
            parser = KelvinConfigParser(fileName)
            self.lines = parser.parse()

        # an object that keeps track of default values of a kelvin config file
        self.kelvinDefaults = defaultsManager

        # a flag to determine if the file has been modified or not since
        # the last save
        self.isModified = False

        # a flag that indicates whether to strip the file, aka take out
        # all comments and blanks in the file before saving it
        self.insertComments = self.checkCommentState()

        # this is a flag to indicate whether the file is using the
        # threshold DT/QT trait type. since in the file the tag is simply
        # QT, other panels may need to know if using threshold trait or not
        if self.isPresent('QT') and self.isPresent('TT'):
            self.trait_threshold = True
        else:
            self.trait_threshold = False

    def save(self):
        if misc.savingEnabled:
            try:
                file = open(self.fileName, 'w')
            except IOError, e:
                print 'Error in writing file %s : %s',(self.fileName,e.strerror)
            else:
                self.writeLines(file)
                file.close
                self.isModified = False

    def saveAs(self, fileName):
        self.fileName = fileName
        self.save()

    def writeLines(self, file):
        """Write all the lines in the parsetree to a file"""


        if not self.insertComments:
            lines = self.getStrippedVersion()
        else:
            lines = self.lines

        for eachLine in lines:
            if eachLine[0].type == types.NEWLINE:
                file.write('\n')
                continue

            s = ""
            for eachToken in eachLine:
                s += eachToken.str + ' '

            # take out whitespace before and after line
            s = s.strip()

            # add a newline
            s += '\n'

            # so it doesn't look weird, replace " ;" with just ";"
            # so GF 0.2 ; 0.3 ; 0.4 becomes GF 0.2; 0.3; 0.4
            s = s.replace(' ;', ';')

            file.write(s)

    def isPresent(self, tag):
        """Given a tag, returns true if the tag appears at least
           once in the parse tree.
        
        """

        l = [x for x in self.lines if x[0].str == tag]
        return len(l) > 0

    def isValuePresent(self, tag):
        """Given a tag, returns true if the tag appears at least
           once in the parse tree, and the line is a value,
           e.g. Th 0.1; 0.2 (point value) or Tm 0.0 0.5 0.01 (range value).

        """

        # first gather all lines that begin with the tag, then test
        # if the second item in the line is a valid number
        lines = [x for x in self.lines if x[0].str == tag]
        for eachLine in lines:
            if len(eachLine) > 1:
                # if the line is point value, the rest of the line after 
                # the tag is one token, so need to parse it
                s = eachLine[1].str
                s = s.partition(';')[0]
                s = s.strip()
                try:
                    num = float(s)
                except ValueError, e:
                    pass
                else:
                    return True
        return False

    def isPointValuePresent(self, tag):
        """Given a tag, returns true if the tag appears at least
           once in the parse tree, and the line is a point value,
           e.g. Th 0.1 or GF 0.1; 0.2; 0.3.

        """

        # Note: this assumes that point and range values are not
        # mixed into the same line

        # first gather all lines that begin with the tag
        lines = [x for x in self.lines if x[0].str == tag]
        for eachLine in lines:
            # a point value line always has two tokens because one is the tag,
            # while everything after the token is moved to one token
            if len(eachLine) == 2:
                # if the line is point value, the rest of the line after 
                # the tag is one token, so need to parse it,
                # check if the first item after the tag is a number, and
                # that the third item is not a number
                s = eachLine[1].str
                head, sep, tail = s.partition(';')
                head = head.strip()
                try:
                    num = float(head)
                except ValueError, e:
                    pass
                else:
                    return True
        return False

    def isRangeValuePresent(self, tag):
        """Given a tag, returns true if the tag appears at least
           once in the parse tree, and the line is a range value,
           e.g. Tm 0.0 0.5 0.01 or Tf 0.0 0.5 0.01; 0.4 0.8 0.05.

        """

        # first gather all lines that begin with the tag, then test
        # if the second item in the line is a valid number
        lines = [x for x in self.lines if x[0].str == tag]
        for eachLine in lines:
            if len(eachLine) % 4 == 0:
                # check if the first two tokens after the tag are numbers
                s1 = eachLine[1].str.strip()
                s2 = eachLine[2].str.strip()
                try:
                    num1 = float(s1)
                    num2 = float(s2)
                except ValueError, e:
                    pass
                else:
                    return True
        return False

    def hasLine(self, line):
        """Returns True if the line is already in the file, otherwise false"""

        # go through each line, checking if the number of tokens is the same,
        # and if so check if the strings are the same

        num = len(line)
        for eachLine in self.lines:
            if len(eachLine) == num:
                bool = True
                for i in range(num):
                    if eachLine[i].str != line[i].str:
                        bool = False
                        break
                if bool:
                    return True

        return False

    def getLines(self, tag):
        """Given a tag, returns a list of the lines in which the
           tag is the first item in that line.
           
        """

        lines = [line for line in self.lines if len(line) > 0 and \
                                                line[0].str == tag]
        return lines

    def createLine(self, *args):
        """Given an argument of strings, treat each string as one token in a 
        new line.  Does not insert the new line into the file.

        """

        newLine = []

        for eachString in args:
            string = str(eachString)
            newToken = Token(type=self.getType(string), str=string)
            newLine.append(newToken)

        return newLine

    def insert(self, *args):
        """Given an argument of strings, treat each string as one token in a 
        new line, and insert that line in the parsetree

        """

        newLine = self.createLine(*args)
        self.addLine(newLine)
        self.isModified = True

        return newLine

    def addLine(self, line):
        """Given a line (a list of tokens), add it to the parsetree. Inserts 
        the line, a comment for the line, and appropriate blank lines.
        
        """

        # first find the correct place to insert the line
        # insertExtras is a flag indicating whether to insert comments or blank
        # lines or just to insert the line in that spot
        index, insertExtras = self.findInsertionPoint(line)

        if insertExtras == False:
            # don't add any blank lines or comments
            self.lines.insert(index, line)
            return

        linesToAdd = [line]

        # add the comment in
        i = self.priorityTags.index(line[0].str)
        commentLine = self.createLine('# ' + KelvinInfo.priorityTagList[i][1])
        linesToAdd.insert(0, commentLine)

        # figure out where to put blank lines
        length = len(self.lines)

        # if the line is at the top and there isn't a blank line underneath,
        # insert a blank line underneath
        if index == 0 and length > 0:
            if self.lines[index][0].type != types.NEWLINE:
                linesToAdd.append(self.createLine('\n'))
        # else if inserting a line at the bottom and no blank line above,
        # insert blank line above
        elif index == length and length > 0:
            if self.lines[index-1][0].type != types.NEWLINE:
                linesToAdd.insert(0, self.createLine('\n'))
        # must be in middle of file, insert lines at the top and 
        # bottom if they aren't present
        elif length > 0:
            if self.lines[index-1][0].type != types.NEWLINE:
                linesToAdd.insert(0, self.createLine('\n'))
            if self.lines[index][0].type != types.NEWLINE:
                linesToAdd.append(self.createLine('\n'))

        # insert chosen lines into parsetree
        for i in range(len(linesToAdd)-1, -1, -1):
            self.lines.insert(index, linesToAdd[i])

    def findInsertionPoint(self, line, useGenotypeFunction=True):
        """Given a line find the correct place to insert it in the parsetree"""

        # the 'useGenotypeFunction' argument lets the function know if it can
        # call findInsertionPointGenotype() or not. this is needed because
        # findInsertionPointGenotype() will call this function, and an
        # infinite loop will occur

        # an inner function that is used frequently in this function
        def getLines(tag, isConstraint, differGenotypes=True):
            """Given a tag and whether to look for constraint lines or not, 
            look for lines that begin with the same tag.  Returning a list 
            of those lines. If nothing is found, an empty list 
            is returned. Also takes into account if its a constraint.
            The differGenotypes flag indicates whether or not if looking
            for lines that start with a genotype, to treat all genotypes as
            a block of the same thing, or to differentiate between them 
            (i.e. DD, Dd, dd).
            
            """

            isGenotype = KelvinInfo.genotypes.count(tag) > 0

            if isConstraint:
                # make sure to look for files that are not just the same tag, 
                # but same type of constraint (don't worry about std
                # constraints since they start with a unique tag)
                if isGenotype:
                    lines = [x for x in self.lines \
                              if KelvinInfo.genotypes.count(x[0].str) > 0]
                elif KelvinInfo.thetaTags.count(tag) > 0:
                    lines = [x for x in self.lines \
                              if KelvinInfo.thetaTags.count(x[0].str) > 0]
                else:
                    lines = [x for x in self.lines if x[0].str == tag]
            else:
                if not differGenotypes and isGenotype:
                    lines = [x for x in self.lines \
                              if KelvinInfo.genotypes.count(x[0].str) > 0]
                else:
                    lines = [x for x in self.lines if x[0].str == tag]

            return [x for x in lines if x[0].isConstraint == isConstraint]

        # check if it is a genotype line, if so use a special function for it
        isGenotype = KelvinInfo.genotypes.count(line[0].str)
        if isGenotype and useGenotypeFunction:
            return self.findInsertionPointGenotype(line)

        # first try to find a line with the exact same tag.
        # if so, add the new line after the last line with the same tag
        lines = getLines(line[0].str, line[0].isConstraint)

        if len(lines) > 0:
            # same line type found, add it to the end, don't add
            # any comments or blank lines
            index = self.lines.index(lines[-1]) + 1
            return index, False
        else:
            # warning: *extremely* ad-hock
            # use the priority list to find out where to insert it
            # first find where the line's tag is in the priority list,
            # look for lines with the tag in the priority list above its tag, 
            # and if any are found, insert the line after those lines.  
            # if not, look at the tag after this line's tag, and look for 
            # lines there, keep going back and forth until something is found

            # first find where in the priority list the line is
            lines = [x for x in KelvinInfo.priorityTagList \
                      if x[0]==line[0].str and x[2]==line[0].isConstraint]
            index = KelvinInfo.priorityTagList.index(lines[0])

            top = index-1
            bottom = index+1
            length = len(self.priorityTags)
            found = False
            if isGenotype and not line[0].isConstraint:
                differGenotypes = True
            else:
                differGenotypes = False
            while not found:
                if bottom < length:
                    tag = KelvinInfo.priorityTagList[bottom][0]
                    isConstraint = KelvinInfo.priorityTagList[bottom][2]
                    lines = getLines(tag, isConstraint, differGenotypes)
                    if len(lines) > 0:
                        found = True
                        insertionPoint = self.lines.index(lines[0]) 
                        break
                    else:
                        bottom += 1
                if top >= 0:
                    tag = KelvinInfo.priorityTagList[top][0]
                    isConstraint = KelvinInfo.priorityTagList[top][2]
                    lines = getLines(tag, isConstraint, differGenotypes)
                    if len(lines) > 0:
                        found = True
                        insertionPoint = self.lines.index(lines[-1]) + 1
                        break
                    else:
                        top -= 1

            # there may be comments in the way, keep going up until its clear
            while insertionPoint != 0 and \
                    self.lines[insertionPoint-1][0].type == types.COMMENT:
                insertionPoint -= 1

            # if the current line is a genotype and the line below or above
            # the insertion point is a genotype, don't add the extra comment
            if isGenotype:
                tag = self.lines[insertionPoint][0].str
                if KelvinInfo.genotypes.count(tag) > 0:
                    return insertionPoint, False
                if insertionPoint != 0:
                    tag = self.lines[insertionPoint-1][0].str
                    if KelvinInfo.genotypes.count(tag) > 0:
                        return insertionPoint, False

            return insertionPoint, True

    def findInsertionPointGenotype(self, line):
        """Given a line that is a genotype, find the correct place to insert
        it in the file.
        
        """

        # since genotypes are really special cases and don't really fit into
        # the framework of normally finding the insertion point, this function
        # does all the work 

        # an inner function that is used frequently in this function
        def getLines(tag, isConstraint, differGenotypes=True):
            """Given a tag and whether to look for constraint lines or not, 
            look for lines that begin with the same tag.  Returning a list 
            of those lines. If nothing is found, an empty list 
            is returned. Also takes into account if its a constraint.
            The differGenotypes flag indicates whether or not if looking
            for lines that start with a genotype, to treat all genotypes as
            a block of the same thing, or to differentiate between them 
            (i.e. DD, Dd, dd).
            
            """

            isGenotype = KelvinInfo.genotypes.count(tag) > 0

            if isConstraint:
                # make sure to look for files that are not just the same tag, 
                # but same type of constraint (don't worry about std
                # constraints since they start with a unique tag)
                if isGenotype:
                    lines = [x for x in self.lines \
                              if KelvinInfo.genotypes.count(x[0].str) > 0]
                elif KelvinInfo.thetaTags.count(tag) > 0:
                    lines = [x for x in self.lines \
                              if KelvinInfo.thetaTags.count(x[0].str) > 0]
                else:
                    lines = [x for x in self.lines if x[0].str == tag]
            else:
                if not differGenotypes and isGenotype:
                    lines = [x for x in self.lines \
                              if KelvinInfo.genotypes.count(x[0].str) > 0]
                else:
                    lines = [x for x in self.lines if x[0].str == tag]

            return [x for x in lines if x[0].isConstraint == isConstraint]

        # first try to find the same type of line as the one being inserted
        isConstraint = line[0].isConstraint
        tag = line[0].str
        lines = getLines(tag, isConstraint)
        if len(lines) > 0:
            # same line type found, add it to the end. 
            # don't add any comments or blank lines
            index = self.lines.index(lines[-1]) + 1
            return index, False

        # check if there are any other type of genotype lines out there.
        # go up-and-down the genotype list to figure out which tag to check for
        index = KelvinInfo.genotypes.index(tag)
        top = index - 1
        bottom = index + 1
        length = len(KelvinInfo.genotypes)
        found = False
        while not found and (top >= 0 or bottom < length):
            if top >= 0:
                tag = KelvinInfo.genotypes[top]
                lines = getLines(tag, isConstraint)
                if len(lines) > 0:
                    found = True
                    insertionPoint = self.lines.index(lines[-1]) + 1
                    break
                else:
                    top -= 1
            if bottom < length:
                tag = KelvinInfo.genotypes[bottom]
                lines = getLines(tag, isConstraint)
                if len(lines) > 0:
                    found = True
                    insertionPoint = self.lines.index(lines[0]) 
                    break
                else:
                    bottom += 1

        if not found:
            # no other genotype lines found, use the normal 
            # findInsertionPoint() function
            return self.findInsertionPoint(line, False)
        else:
            # make sure no other comments get inserted
            return insertionPoint, False

    def remove(self, line):
        """Remove the given line if it is present"""

        linesToRemove = [line]

        index = self.lines.index(line)
        top = index
        bottom = index
        length = len(self.lines)

        # some ad-hock rules about removing lines

        # if there is exactly one comment above the line, then
        # a blank line above the comment, and a blank line underneath the
        # line, remove the comment
        if index != 0:
            if self.lines[index-1][0].type == types.COMMENT:
                if index == 1:
                    if index < length-1:
                        if self.lines[index+1][0].type == types.NEWLINE:
                            linesToRemove.append(self.lines[index-1])
                            top -= 1
                    else:
                        linesToRemove.append(self.lines[index-1])
                        top -= 1
                elif self.lines[index-2][0].type == types.NEWLINE:
                    if index < length-1:
                        if self.lines[index+1][0].type == types.NEWLINE:
                            linesToRemove.append(self.lines[index-1])
                            top -= 1
                    else:
                        linesToRemove.append(self.lines[index-1])
                        top -= 1

        # if there is a blank above and below the region to be deleted,
        # remove the bottom blank line
        if top != 0 and bottom+1 != length:
            if self.lines[top-1][0].type == types.NEWLINE:
                if self.lines[bottom+1][0].type == types.NEWLINE:
                    linesToRemove.append(self.lines[bottom+1])
                    bottom += 1

        # if the deleted region is at the top and right below it is a 
        # blank line, remove the blank line
        if top == 0 and bottom+1 != length:
            if self.lines[bottom+1][0].type == types.NEWLINE:
                linesToRemove.append(self.lines[bottom+1])
                bottom += 1

        # if the deleted region is at the bottom and there is a blank line
        # right above it, remove the blank line
        if bottom+1 == length and top != 0:
            if self.lines[top-1][0].type == types.NEWLINE:
                linesToRemove.append(self.lines[top-1])
                top -= 1

        if len(linesToRemove) > 0:
            self.isModified = True

        for eachLine in linesToRemove:
            self.lines.remove(eachLine)

    def append(self, line, value):
        """Given a line and a string value, appends the string value onto 
        that line, and returns the token created for value.
        
        """

        # create a new token using value as the token string value
        token = Token(str=str(value))

        # append the new token onto the line
        line.append(token)
        self.isModified = True

        return token

    def getVariableNumberLines(self):
        """Gets all the lines that deal with a variable number tag"""

        return [line for line in self.lines if line[0].isVariableNumber]

    def getConstraintLines(self):
        """Gets all the lines that specifies a constraint"""

        return [line for line in self.lines if line[0].isConstraint]

    def getType(self, str):
        """Given a string, returns the type of what that string is, for
        example, given '1', the function will return types.NUMBER
        
        """


        # check if it's an operator
        if KelvinInfo.operatorTags.count(str) > 0:
            return types.OPERATOR
        elif str.isdigit():
            return types.NUMBER
        elif str.isalnum():
            return types.ID
        elif str == ';':
            return types.SEMICOLON
        elif len(str) > 0 and str[0] == '#':
            return types.COMMENT
        elif str == '\n':
            return types.NEWLINE


        # not a type that should be there
        return types.UNKNOWN

    def getFirstConstraint(self, line):
        """Given a line from the file, gets the first constraint"""

        # there may be multiple constraints in one line joined by a semicolon,
        # and need to get the first one, handy for trying to figure out
        # what type of constraint a line is
        newLine = []
        for eachToken in line:
            if eachToken.type == types.SEMICOLON:
                break
            newLine.append(eachToken)

        return newLine

    def getStrippedVersion(self):
        """Returns a version of the parsetree that doesn't have any
        newlines or comments in it.
        
        """

        return [x for x in self.lines if x[0].type != types.NEWLINE and \
                                         x[0].type != types.COMMENT]

    def checkCommentState(self):
        """Returns True if there are any comments in the file, else False"""

        for eachLine in self.lines:
            if eachLine[0].type == types.NEWLINE:
                return True

        return False

    def getDefaults(self, tag, isConstraint):
        """Given a tag and whether a constraint is wanted or not,
        return a list of lines with the tag that are the default values
        for that tag.
        
        """

        # consult the defaults manager
        if self.kelvinDefaults is not None:
            return self.kelvinDefaults.getDefaults(self, tag, isConstraint)
        else:
            return []

    def getValueDefaults(self, analysis=None, trait=None):
        """Return a list of lines that are the default values considering 
        the current file condition.
        
        """

        if self.kelvinDefaults is not None:
            return self.kelvinDefaults.getValueDefaults(self, analysis, trait)
        else:
            return []
