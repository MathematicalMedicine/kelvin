"""This class analysis text in a file, and returns token objects \
   based on the file contents."""

from token import Token
from fileManager import types

class Lexer(object):
    """This class analysis text in a file, and returns token objects based on 
    the file contents.
    
    """

    def __init__(self, fileName):
        # the current position of where the scanning is
        self.row = 1
        self.column = 0

        self.fileName = fileName
        self.file = open(self.fileName, 'r')

    def lex(self):
        """recognizes the next token in the file, and returns it"""

        self.skipWhiteSpace()

        ch = self.readChar()

        if ch == "":
            # must be at the end of the file, so close it
            self.closeFile()
            return Token(types.ENDOFINPUT, types.ENDOFINPUT, 
                         self.row, self.column)
        elif ch == '\n':
            # token is a newline
            column = self.column
            row = self.row
            self.column = 0
            self.row += 1
            return Token(types.NEWLINE, ch, row, column)
        elif ch == ';':
            return Token(types.SEMICOLON, ';', self.row,self.column)
        elif ch == '#':
            # line is a comment
            row = self.row
            column = self.column
            self.pushChar()

            s = self.readline()

            # readline keeps the newline at the end, 
            # so take it out if it is present
            if(s[-1] == '\n'):
                s = s[:-1]
                self.pushChar()
            return Token(types.COMMENT, s, row, column)
        elif ch == '>':
            # check if just greater than or greater and equal
            s = ch
            row = self.row
            column = self.column
            ch = self.readChar()
            if ch == '=':
                s += ch
                return Token(types.OPERATOR, s, row, column)
            else:
                self.pushChar()
                return Token(types.OPERATOR, s, row, column)
        elif ch == '!':
            # a not equal operator
            s = ch
            row = self.row
            column = self.column
            ch = self.readChar()
            if ch == '=':
                s += ch
                return Token(types.OPERATOR, s, row, column)
            else:
                self.pushChar()
                return Token(types.UNKNOWN, s, row, column)
        elif ch == '=':
            # an equal operator
            s = ch
            row = self.row
            column = self.column
            ch = self.readChar()
            if ch == '=':
                s += ch
                return Token(types.OPERATOR, s, row, column)
            else:
                self.pushChar()
                return Token(types.UNKNOWN, s, row, column)
        elif ch.isalpha() or ch == '/' or ch == '\\' or ch == '.':
            # token is an id
            s = ""
            row = self.row
            column = self.column

            # have to allow for any character that can be in a filename
            while not ch.isspace() and ch != ';':
                s += ch
                ch = self.readChar()
            self.pushChar()
            return Token(types.ID, s, row, column)
        elif ch.isdigit() or ch == '-':
            # token is a number
            self.pushChar()
            return self.lexNumber(False)
        else:
            # character is not recognized at all
            return Token(types.UNKNOWN, ch, self.row, self.column)

    def readChar(self):
        """reads the next character in the file"""
        if self.file.closed:
            return ""
        else:
            self.column += 1
            return self.file.read(1)

    def pushChar(self):
        """goes backwards once in the file"""
        self.file.seek(-1, 1)
        self.column -= 1

    def readline(self):
        """reads the rest of the line in the file"""
        s = self.file.readline()
        self.column += len(s)
        return s

    def skipWhiteSpace(self):
        # skip all white spaces except for newlines
        ch = self.readChar()
        while ch.isspace() and ch != '\n':
            ch = self.readChar()

        # don't push the character back if at the end of the file
        if ch != "":
            self.pushChar()

    def lexNumber(self, foundDot):
        ch = self.readChar()
        if ch.isdigit() or ch == '-':
            s = ch
        else:
            self.pushChar()
            return
        ch = self.readChar()
        row = self.row
        column = self.column
        while ch.isdigit():
            s += ch
            ch = self.readChar()
        if ch == '.' and not foundDot:
            # a decimal point, keep scanning
            s += ch
            s += self.lexNumber(True).str
        else:
            self.pushChar()

        return Token(types.NUMBER, s, row, column)

    def closeFile(self):
        if not self.file.closed:
            self.file.close()
