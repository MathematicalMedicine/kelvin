"""set of functions to create a menubar"""

# to use, give createMenuBar a frame where the menu will be created, 
# and a list of the menuData.  menu name, then status bar text, then handler.
# Here is an example of a menuData:
#[("File", (
#     ("&New", (
#              ("&Kelvin Config File", "Kelvin", None),
#              ("&Pedigree File", "Pedigree", None),
#              ("&Marker File", "Marker", None))),
#     ("&Open", "Open File", None),
#     ("&Save", "Save File", None),
#     ("Save &As", "Save File As", None),
#     ("", "", ""),
#     ("&Quit", "Quit", self.OnCloseWindow))),
# ("Run", (
#         # notice the tuple of one item needs a trailing comma
#         ("&Run Kelvin", "Run Kelvin with Config File", None),))]

import wx

def createMenuBar(frame, menuData):
        menuBar = wx.MenuBar()
        for eachMenuData in menuData:
            menuLabel = eachMenuData[0]
            menuItems = eachMenuData[1]
            menuBar.Append(createMenu(frame, menuItems), menuLabel)
        frame.SetMenuBar(menuBar)
        return menuBar

def createMenu(frame, menuData):
    menu = wx.Menu()
    for eachItem in menuData:
        if len(eachItem) == 2:
            label = eachItem[0]
            subMenu = createMenu(frame, eachItem[1])
            menu.AppendMenu(-1, label, subMenu)
        else:
            createMenuItem(frame, menu, *eachItem)
    return menu

def createMenuItem(frame, menu, label, status, handler, kind=wx.ITEM_NORMAL):
    if not label:
        menu.AppendSeparator()
        return
    menuItem = menu.Append(-1, label, status, kind)
    frame.Bind(wx.EVT_MENU, handler, menuItem)
