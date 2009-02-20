#!/usr/bin/env python

"""gKelvin is a GUI for Kelvin.  Run this module as the main module, or
simply call main().

"""

import wx
import panel.main
import fileManager.types
from fileManager import lexer
from fileManager.KelvinConfigParser import KelvinConfigParser
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager

class gKelvinApp(wx.App):
    """The class that is the main application"""

    def OnInit(self):
        self.frame = panel.main.MainFrame()
        self.frame.CenterOnScreen()
        self.frame.Show()
        self.SetTopWindow(self.frame)
        return True

def main():
    app = gKelvinApp(False)
    app.MainLoop()

def testLexer(fileName):
    # print out all tokens the lexer lexes
    lexMachine = lexer.Lexer(fileName)
    a = lexMachine.lex()
    while a.type != fileManager.types.ENDOFINPUT:
       a.printInfo()
       print 
       a = lexMachine.lex()
    a.printInfo()

def testParser(fileName):
    # prints out the entire parse tree, with indentions for displaying the level
    def printParseTree(tree, level):
        indent = 4
        if type(tree[0]) is not list:
            #must be a token
            s = ' ' * (level-1) * indent
            for eachToken in tree:
                if eachToken.type == fileManager.types.NEWLINE:
                    break
                s += eachToken.str + ' '
            print s
            return
        for eachItem in tree:
            printParseTree(eachItem, level+1)

    parser = KelvinConfigParser(fileName)
    a = parser.parse()
    printParseTree(a, 0)

if __name__ == "__main__":
    main()
