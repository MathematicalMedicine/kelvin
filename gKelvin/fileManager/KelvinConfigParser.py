"""A parser to determine that a kelvin config file has the correct syntax"""

from fileManager import types
from fileManager.parser import Parser
from config import KelvinInfo
from fileManager.parseError import ParseError

class KelvinConfigParser(Parser):
    """A parser to determine that a kelvin config file has the correct syntax"""

    def __init__(self, fileName):
        Parser.__init__(self, fileName)

        # gather some information 

        # tags that have a set number of items after them
        self.setNumberTags = self.getSetNumberTags()

        # tags that have a varialbe number of items after them
        self.variableNumberTags = self.getVariableNumberTags()

    def parse(self):
        """The function to parse the file.  Will return a list where each \n\
           item represents a line, and in each item is a list of token \n\
           objects that was in that line"""

        # goto the first token in the file
        self.advance()

        # keep getting tokens until reach end of input
        while not self.checkType(types.ENDOFINPUT):
            # check for every possible type of line
            if self.checkType(types.COMMENT):
                self.addLine([self.match(types.COMMENT)])
                self.match(types.NEWLINE)
            elif self.checkType(types.NEWLINE):
                # keep blank lines
                self.addLine([self.match(types.NEWLINE)])
            elif self.checkType(types.ID):
                self.checkTagLines()
            elif self.isInteger():
                # the only way an integer starts a line is if it's a constraint
                self.addLine(self.getConstraintList())
                error = "Extra items after the constraint, is a ';' needed?"
                self.match(types.NEWLINE, error)
            else:
                # raise an error exception
                s1 = "Error opening file '%s'!\n" % self.fileName
                s2 = "'%s' is an illegal way to begin a line\n" % \
                     self.currentToken.str
                s3 = "Error at line:%i column:%i\n" % \
                     (self.currentToken.row, self.currentToken.column)
                s = s1 + s2 + s3 
                raise ParseError(s)

        return self.lines

    def checkTagLines(self):
        # figured out that the first token is an id, so it must be some tag
        # now check what kind of tag it is

        # check for tags that fit a standard pattern
        # note that constraints need to be checked first, to differentiate
        # between a line like "DD >= dd" versus "DD 0.2"
        if self.checkConstraintList(): 
            pass
        elif self.checkFileTags():
            pass
        elif self.checkSetNumberTag():
            pass
        elif self.checkVariableNumberTag():
            pass

        # check for other miscellaneous tags that don't fit previous patterns
        elif self.checkToken("QT"):
            self.getQuantitativeTag()
        else:
            # unknown tag, assume it's some parameter that is variable number
            line = [self.match(types.ID)]
            line.extend(self.getNumberValuesList())
            error = 'Extra items on the line'
            self.match(types.NEWLINE, error)
            line[0].isVariableNumber = True
            self.addLine(line)

    def checkFileTags(self):
        # checks for all possible tags that specify a file, 
        # format is <tag> <filename>
        filetags = KelvinInfo.getFileNameTags()

        error1 = 'Error with filename specification'
        error2 = 'Extra items after filename specification'
        for eachfiletag in filetags:
            if self.checkToken(eachfiletag):
                self.addLine([self.match(types.ID),self.match(types.ID,error1)])
                self.match(types.NEWLINE, error2)
                return True

        # no matches were found, must be a different tag
        return False

    def getSetNumberTags(self):
        # tags that have a standard number after them
        # the list is that the first item is the tag, then
        # the second item is the number of numbers following the tag
        # example: ('AS', 3) means the AS tag is followed by three numbers
        tags = []

        # add dichotomous trait type
        dichotomous = KelvinInfo.traitTypesEntry[0]
        tags.append((dichotomous[0], dichotomous[2]))

        # add analysis type tags
        analysisTags = KelvinInfo.getAnalysisTypeTags()
        analysisNumbers = KelvinInfo.getAnalysisTypeNumbers()
        tags.extend(zip(analysisTags, analysisNumbers))

        # add unknown id tag
        unknownId = KelvinInfo.unknownPersonId
        tags.append((unknownId[0], unknownId[1]))

        # add liability class tag
        tags.append(KelvinInfo.liabilityClassEntry)

        # add affection status tag
        tags.append(KelvinInfo.affectionStatusEntry)

        # add disease allele tag
        diseaseAllele = KelvinInfo.diseaseAlleles
        tags.append((diseaseAllele[0], diseaseAllele[1]))

        # add polynomial evaluation tag
        tags.append(KelvinInfo.polynomialEvaluationEntry)

        # add other tags
        tags.append(KelvinInfo.TMEntry)
        tags.append(KelvinInfo.MMEntry)
        tags.append(KelvinInfo.AMEntry)
        tags.append(KelvinInfo.chromosomeEntry)
        tags.append(KelvinInfo.MINEntry)
        tags.append(KelvinInfo.MAXEntry)
        tags.append(KelvinInfo.TMINEntry)
        tags.append(KelvinInfo.TMAXEntry)

        return tags

    def checkSetNumberTag(self):
        # checks tags that have a standard number of numbers after it
        error1 = 'Looking for exactly %i number(s) in the %s specification'
        error2 = 'Extra items after the %s specification'
        for eachTag, eachNumber in self.setNumberTags:
            if self.checkToken(eachTag):
                line = [self.match(types.ID)]
                for i in range(eachNumber):
                    line.append(self.match(types.NUMBER, 
                                error1%(eachNumber-i, eachTag)))
                self.match(types.NEWLINE, error2 % eachTag)
                self.addLine(line)
                return True

        # no tag matches found
        return False

    def getVariableNumberTags(self):
        # these tags are followed by either one number or three numbers
        # also numbers can be separated by semicolons
        # example: GF 0.2
        # example: GF 0 1 0.01
        # example: GF 0.2; 0 1 0.01; 0.4
        info = KelvinInfo.variableNumberInfo
        return [x.tag for x in info]

    def checkVariableNumberTag(self):
        # these tags are followed by either one number or three numbers
        # also numbers can be separated by semicolons

        error = "Extra items after the %s specification, is a ';' needed?"
        for eachTag in self.variableNumberTags:
            if self.checkToken(eachTag):
                line = [self.match(types.ID)]
                line.extend(self.getNumberValuesList())
                self.match(types.NEWLINE, error % eachTag)
                self.addLine(line)

                # set a flag
                line[0].isVariableNumber = True

                return True

        # no tag matches were found
        return False

    def getNumberValuesList(self):
        # gets a number values list, which is either one or three numbers
        # that are separated by semicolons

        # there is at least one number
        error = "Looking for a number to set a parameter to a value or range." 
        l = [self.match(types.NUMBER, error)]
        if self.checkType(types.NUMBER):
            # another number means it's three numbers to specify a range
            error = "Looking for a %s number to set a parameter to a range."
            l.append(self.match(types.NUMBER, error % 'second'))
            l.append(self.match(types.NUMBER, error % 'third'))
        if self.checkType(types.SEMICOLON):
            # another set of numbers coming
            l.append(self.match(types.SEMICOLON))
            l.extend(self.getNumberValuesList())
        return l

    def checkConstraintList(self):
        """check if the line has one or more constraints"""

        # constraint starts with at least one tag or number
        if self.checkType(types.ID):
            line = [self.match(types.ID)]
        else:
            line = [self.match(types.NUMBER)]

        # could be more tags or numbers after the first one
        while self.checkType(types.ID) or self.checkType(types.NUMBER):
            if self.checkType(types.ID):
                line.append(self.match(types.ID))
            else:
                line.append(self.match(types.NUMBER))

        # the id's and numbers should be followed by an operator
        if not self.checkType(types.OPERATOR):
            isConstraint = False
        else:
            isConstraint = True

        # push all tokens that have been read so far
        line.append(self.currentToken)
        self.saveTokens(line)
        self.advance()

        # if it is a constraint, read it in
        if isConstraint:
            self.addLine(self.getConstraintList())
            self.match(types.NEWLINE);

        return isConstraint

    def getConstraintList(self):
        # five types of constraints
        line = []

        # one type involves only the sex-specific theta (Tm or Tf)
        # for example, Tm >= Tf
        if self.checkToken('Tf') or self.checkToken('Tm'):
            line = self.getThetaConstraint()

        # two types involves the first token being a genotypes,
        # either just constraining the genotype or a constraint
        # of a genotype and a liability class
        # for example, DD > Dd  or 00 >= 01, DD 1 > Dd 2
        if self.isGenotype():
            line = self.getGenotypeConstraint()

        # the last two types of constraints start with a parameter
        # and can have the liability class or not
        # for example, P1 DD > P1 Dd, P1 Dd 1 > P1 DD 2
        if self.isParameter():
            line = self.getParameterConstraint()

        # there could be a semicolon joining more constraints
        if(self.checkType(types.SEMICOLON)):
            line.append(self.match(types.SEMICOLON))
            line.extend(self.getConstraintList())

        # set some flags
        line[0].isConstraint = True

        return line

    def getGenotypeConstraint(self):
        """Match a constraint that starts with a genotype, there are two cases,
        just a genotype constraint, like DD > dd or with a liability class,
        for example DD 1 > dd 2"""

        line = []
        # remember exactly what type of genotype (id or int) it was
        if self.isGenotype() and self.isInteger():
            currentType = types.NUMBER
        else:
            currentType = types.ID
        
        # the first token is always a genotype (id or int)
        line.append(self.match(currentType))

        # there could be a liability class for the second token
        if self.checkType(types.NUMBER):
            liability = True
            line.append(self.match(types.NUMBER))
        else:
            liability = False

        line.append(self.match(types.OPERATOR))

        # match the same type of genotype as the previous one
        error = "Expecting a valid genotype as the same type as before"
        line.append(self.match(currentType, error))

        # match a liability class if there was one before
        if liability:
            error = "Expecting a liability class"
            line.append(self.match(types.NUMBER, error))

        return line

    def getThetaConstraint(self):
        """Parse a constraint that deals with theta, for example, Tf > Tm"""

        line = []
        error1 = "Expecting either 'Tf' or 'Tm'"
        error2 = "Expecting an operator in a theta constraint"

        line.append(self.match(types.ID, error1))
        line.append(self.match(types.OPERATOR, error2))
        if self.checkToken('Tf'):
            line.append(self.match(types.ID, str='Tf'))
        elif self.checkToken('Tm'):
            line.append(self.match(types.ID, str='Tm'))
        else:
            # there is an error, use these statements to
            # generate a logical sounding error message
            line.append(self.match(types.ID, error1, str="'Tm' or 'Tf'"))

        return line

    def getParameterConstraint(self):
        """Match the types of constraints that start with a parameter
        and may have an optional liability class,
        for example, P1 DD > P1 Dd, P1 Dd 1 > P1 DD 2

        """

        line = []

        # assume the first token was already 
        # checked to see if it was a parameter
        line.append(self.match(types.ID))

        # next token is always a genotype tag
        if self.isGenotype():
            line.append(self.match(types.ID))
        else:
            # statements to create a logical error message
            error = "Expecting either 'DD', 'Dd', or 'dd'"
            line.append(self.match(types.ID,error,str="'DD' or 'Dd' or 'dd'"))

        #optional liability class
        if self.isInteger():
            line.append(self.match(types.NUMBER))
            liability = True
        else:
            liability = False

        line.append(self.match(types.OPERATOR))

        if self.isParameter():
            line.append(self.match(types.ID))
        else:
            # statements to create a logical error message
            error = "Expecting a parameter, e.g., 'P1'"
            line.append(self.match(types.ID, error, str="P1"))

        # another genotype tag
        if self.isGenotype():
            line.append(self.match(types.ID))
        else:
            # statements to create a logical error message
            error = "Expecting either 'DD', 'Dd', or 'dd'"
            line.append(self.match(types.ID,error,str="'DD' or 'Dd' or 'dd'"))

        # a matching liability class if there was one before
        if liability:
            error = "Expecting a liability class"
            line.append(self.match(types.NUMBER, error))

        return line

    def isGenotype(self):
        """Checks if the current token is a valid genotype, of the form 
        'DD', 'Dd', 'dd', or an integer that represents the genotype

        """

        # check for valid tags
        for eachtag in KelvinInfo.genotypes:
            if self.checkToken(eachtag):
                return True

        # a valid genotype is also an integer
        if self.isInteger():
            return True

        # nothing matches, must not be genotype
        return False

    def isInteger(self):
        """Checks if the current token is an integer"""
        return self.checkType(types.NUMBER) and \
               self.currentToken.str.find('.') == -1

    def isParameter(self):
        """Check if the current token represents a parameter"""

        # it is a parameter if the token is of length 2, and it's
        # of the form P?, where ? is any integer
        str = self.currentToken.str
        if self.checkType(types.ID) and len(str) == 2 and \
           str[0] == 'P' and str[1].isdigit():
            return True
        else:
            return False

    def getQuantitativeTag(self):
        """Parses lines with the QT directive"""

        error1 = 'Expecting a word describing the type of distribution'
        error2 = 'Expecting a parameter describing the distribution'
        error3 = 'Extra items after the QT specification'
        
        line = []

        line.append(self.match(types.ID, str='QT'))

        dist = self.currentToken.str

        # format depends on type of distribution
        foundit = False
        for eachString,eachLabel,eachNumParams in KelvinInfo.distributionEntry:
            if self.checkToken(eachString):
                line.append(self.match(types.ID))
                foundit = True
                for i in range(0, eachNumParams):
                    line.append(self.match(types.NUMBER, error2))

        if not foundit:
            # not a known distribution type, make a logical error message
            error = "Expecting a distribution type"
            dist = [x[0] for x in KelvinInfo.distributionEntry]
            dist = ' or '.join(dist)
            self.match(types.ID, error, str=dist)

        self.addLine(line)
        self.match(types.NEWLINE, error3)
