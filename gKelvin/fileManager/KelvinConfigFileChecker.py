"""Checks a Kelvin config file for errors, and shows prompts about them"""

import wx
import os
import sys
from config import KelvinInfo
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager

class KelvinConfigFileChecker(object):
    """Checks a Kelvin config file for errors, and shows errors about them"""

    def __init__(self, kelvinConfigPanel, fileManager, fileName):
        self.kelvinConfigPanel = kelvinConfigPanel
        self.fm = fileManager
        self.fileName = fileName

    def checkFileForErrors(self):
        """Does validation checks on the file for errors in the semantics"""

        # if there were no errors at all, returns False, if there was an
        # error found, return True

        functions = [self.checkMinMax, self.checkTL, self.checkTT]
        result = False
        for eachFunc in functions:
            result = result or eachFunc()
        return result

    def checkFileForWarnings(self):
        """Checks the file for errors that result in warnings"""

        # if there were no errors at all, returns False, if there was an
        # error found, return True

        functions = [self.checkIfFileNamesExist, self.checkTheta]
        result = False
        for eachFunc in functions:
            result = result or eachFunc()
        return result

    def checkMinMax(self):
        """If the trait type is quantitative and trait threshold is specified,
        then TMIN and TMAX must be specified.
        
        """

        # note: there is an identical check in traitTypePanel.py
        tmin = KelvinInfo.TMINEntry[0]
        tmax = KelvinInfo.TMAXEntry[0]
        if self.fm.isPresent('QT') and self.fm.isPresent('TT'):
            if not self.fm.isPresent(tmin) or not self.fm.isPresent(tmax):
                s1 = "The file is configured for a quantitative trait and "
                s2 = "has specified a trait threshold, but is missing Min "
                s3 = "and/or Max values. Please input these values in the "
                s4 = "Settings panel before running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                self.showError(message)
                return True

        return False

    def checkTL(self):
        """If multi-point is specified, then TL must be specified"""

        # note: you can also specify TM instead of TL, but it's really rare,
        # so don't mention it in the error message
        SS = self.fm.isPresent('SS')
        SA = self.fm.isPresent('SA')
        TL = self.fm.isPresent('TL')
        TL = self.fm.isPresent('TL')
        TM = self.fm.isPresent('TM')
        if (SS or SA) and not (TL or TM):
            s1 = "The file is configured for a multi-point analysis, but "
            s2 = "a trait location is not specified. Please specify this value "
            s3 = "in the Grid Specification before running Kelvin.\n"
            message = s1 + s2 + s3
            self.showWarning(message)
            return True
        return False

    def checkTT(self):
        """If the trait type selected is threshold dt/qt, then 
        TT (trait threshold) must be specified.
        
        """

        # first check if threshold dt/qt is specified, which for right now
        # means that the choice box has the third option selected
        choice = self.kelvinConfigPanel.traitTypePanel.traitChoice
        if choice.GetSelection() == 2:
            if not self.fm.isPresent('TT'):
                s1 = "The trait type is Threshold DT/QT, which requires a "
                s2 = "trait threshold be specified. Please specify a trait "
                s3 = "threshold in the Grid Specification before "
                s4 = "running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                self.showWarning(message)
                return True
        return False

    def checkTheta(self):
        """Check if the file correctly specified theta values"""

        TP = self.fm.isPresent('TP')
        SS = self.fm.isPresent('SS')
        SA = self.fm.isPresent('SA')
        Th = self.fm.isValuePresent('Th')
        Tm = self.fm.isValuePresent('Tm')
        Tf = self.fm.isValuePresent('Tf')
        check1 = self.checkThetaTwoPoint(TP, Th, Tm, Tf)
        check2 = self.checkThetaMultiPoint(SS, SA, Th, Tm, Tf)

        return check1 or check2

    def checkThetaTwoPoint(self, TP, Th, Tm, Tf):
        """If the file is twopoint, theta must be expressed"""

        if TP:
            if not Th and not (Tm or Tf):
                s1 = "The file is configured for a two point analysis, but "
                s2 = "theta is not specified.  Please specify values for "
                s3 = "theta in the Grid Specification before running Kelvin.\n"
                message = s1 + s2 + s3
                self.showWarning(message)
                return True
            elif (Tm and not Tf) or (not Tm and Tf):
                s1 = "The file is configured for a two point analysis, but "
                s2 = "either theta-male or theta-female is not specified. "
                s3 = "Please specify values for theta-male and theta-female "
                s4 = "in the Grid Specification before running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                self.showWarning(message)
                return True
            elif Th and Tm and Tf:
                s1 = "The file is configured for a two point analysis, and "
                s2 = "both sex average and sex specific theta values are "
                s3 = "specified. Please remove unnecessary theta values "
                s4 = "in the Grid Specification before running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                self.showWarning(message)
                return True
        return False

    def checkThetaMultiPoint(self, SS, SA, Th, Tm, Tf):
        """If the file is multi-point, theta should not be in the config file"""

        if SS or SA:
            if Th or Tm or Tf:
                s1 = "The file is configured for a multi-point analysis, but "
                s2 = "either theta or theta-male or theta-female is specified. "
                s3 = "Please remove values for theta "
                s4 = "in the Grid Specification before running Kelvin.\n"
                message = s1 + s2 + s3 + s4
                self.showWarning(message)
                return True
        return False

    def checkIfFileNamesExist(self):
        """Checks whether the filenames specified in the file exist or not"""

        # a filename must be present in order for this test to make sense
        if not self.fileName:
            return False

        # a list of all the filenames specified that weren't found
        badFileNames = []

        for eachEntry in KelvinInfo.fileNameInfo:

            # don't check for files that aren't required
            if eachEntry.isRequired is False:
                continue

            # get the line where the tag is specified,
            # assume there's a max of one line for each file tag
            lines = self.fm.getLines(eachEntry.tag)
            
            if len(lines) > 0:
                # a filename for this tag was specified
                fileName = lines[0][1].str
            else:
                # was not specified in file, use default
                fileName = eachEntry.default

            fileName = self.getRelativeFilePath(fileName)

            # check if each filename exists
            if not os.path.exists(fileName):
                badFileNames.append((eachEntry.description, fileName))

        if len(badFileNames) > 0:
            # show an error message
            files = ''
            indent = 0
            for eachFile in badFileNames:
                files += ' '*indent + eachFile[0] + ': ' + eachFile[1] + '\n'

            s1 = "The following files specified in '%s'" % self.fm.fileName
            s2 = " do not exist or cannot be opened:\n\n" 
            s3 = "Please make sure all these files exist before running Kelvin."

            message = s1 + s2 + files + '\n' + s3
            self.showWarning(message)
            return True

        return False

    def showWarning(self, message):
        """Shows a warning message box"""
            
        print message
        dlg = wx.MessageDialog(None, message, "Warning!", 
                               wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def showError(self, message):
        """Shows a warning message box"""
            
        print message
        dlg = wx.MessageDialog(None, message, "Error!", 
                               wx.OK | wx.ICON_ERROR)
        dlg.ShowModal()
        dlg.Destroy()

    def getRelativeFilePath(self, fileName):
        """Given a filename, make sure the path is correct relative to the file 
        being open.  This is used because otherwise the filename will be checked
        relative to the current working directory of the python interpreter.
        
        """

        # if fileName is an absolute address, nothing to do
        if os.path.isabs(fileName):
            return fileName

        # fileName is a relative address, make it an absolute one by joining
        # it with the path of the file being opened
        path = os.path.dirname(self.fileName)
        return os.path.join(path, fileName)
