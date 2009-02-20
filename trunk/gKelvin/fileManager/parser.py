"""A base class for a parser"""

from token import Token
from lexer import Lexer
from fileManager.parseError import ParseError
from fileManager import types

class Parser(object):
    """A parser to determine that a kelvin config file has the correct syntax"""
    def __init__(self, fileName):
        self.fileName = fileName
        self.currentToken = ""
        self.lexMachine = Lexer(self.fileName)

        # a list of previously lexed tokens, used when needed
        # to go backwards in the file
        self.tokenList = []

        # the list of lines of tokens after the file is parsed
        self.lines = []

    def parse(self):
        """The function that parses the file and returns a list of tokens"""
        # implement parse function here

        return self.lines

    def match(self, type, errorMsg='', str=''):
        currentStr = self.currentToken.str
        if currentStr == '\n':
            currentStr = types.NEWLINE

        if not self.checkType(type):
            # wrong syntax!
            #self.raiseWrongTypeError(type, currentStr, errorMsg)
            # the a is added for it to make sense
            self.raiseParseError('a ' + type, currentStr, errorMsg)
        if str != '' and not self.checkToken(str):
            # right type, but wrong word!
            #self.raiseWrongTokenError(str, errorMsg)
            # add quotes around what was expected
            self.raiseParseError("'"+str+"'", currentStr, errorMsg)
        token = self.currentToken
        self.advance()
        return token

    def checkType(self, type):
        """Checks whether the current token is the same as the type passed in"""
        return self.currentToken.type == type

    def checkToken(self, str):
        """Checks whether the current token has the same token string as
           the string passed in"""
        # note that this is case sensitive
        return self.currentToken.str == str

    def advance(self):
        if len(self.tokenList) == 0:
            try:
                self.currentToken = self.lexMachine.lex()
            except ParseError:
                self.lexMachine.closeFile()
                raise
        else:
            self.currentToken = self.tokenList.pop(0)
        return self.currentToken

    def saveTokens(self, tokens):
        """Save one token or a list of tokens on the queue"""
        if type(tokens) == list:
            # must be a list of tokens
            self.tokenList.extend(tokens)
        else:
            # argument was simply one token
            self.tokenList.append(tokens)

    def addLine(self, line):
        """Adds an entry to the lines list"""
        self.lines.append(line)

    def raiseParseError(self, expecting, actuallyFound, errorMsg=''):
        """Raise an error when the wrong type of token is encountered"""

        s1 = "Error opening file %s!\n" % self.fileName
        s2 = "Was expecting %s, but found '%s' instead\n" % \
              (expecting, actuallyFound)
        s3 = "Error location -- line:%i column:%i\n" % \
              (self.currentToken.row, self.currentToken.column)

        if errorMsg == '':
            s = s1 + s2 + s3
        else:
            s = s1 + s2 + errorMsg + '\n' + s3 

        raise ParseError(s)
