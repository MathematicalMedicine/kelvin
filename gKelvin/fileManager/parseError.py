"""This module defines a parse error class that is thrown when a file being 
parsed contains illegal symbols or the structure of the file is incorrect"""

class ParseError(Exception):
    """Exception raised when a file contains illegal symbols or the 
    structure of the file is incorrect"""

    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message
