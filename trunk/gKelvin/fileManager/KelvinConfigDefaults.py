"""An object to keep track of default values for a kelvin config file"""

from config import KelvinInfo
import copy
from fileManager.KelvinConfigFileManager import KelvinConfigFileManager

class KelvinConfigDefaults(object):
    """An object to keep track of default values for a kelvin config file"""

    # the filenames of the files holding default values
    TP_DT = KelvinInfo.newConfigFile[0][1]
    TP_QT = 'twopoint_qt.conf'
    MP_DT = KelvinInfo.newConfigFile[1][1]
    MP_QT = 'multipoint_qt.conf'

    def __init__(self):
        # this class has four filemanagers, with different values of
        # two-point/multi-point and DT/QT
        self.tp_dt_fm = KelvinConfigFileManager(self.TP_DT, True)
        self.tp_qt_fm = KelvinConfigFileManager(self.TP_QT, True)
        self.mp_dt_fm = KelvinConfigFileManager(self.MP_DT, True)
        self.mp_qt_fm = KelvinConfigFileManager(self.MP_QT, True)

        if self.tp_dt_fm.isPresent('AF'):
            print 'tp_dt_fm is present!'
        if self.tp_qt_fm.isPresent('AF'):
            print 'tp_qt_fm is present!'
        if self.mp_dt_fm.isPresent('AF'):
            print 'mp_dt_fm is present!'
        if self.mp_qt_fm.isPresent('AF'):
            print 'mp_qt_fm is present!'

    def getDefaults(self, fm, tag, isConstraint):
        """Returns a list of lines that are the default values for
        the given tag. The argument isConstraint specifies whether to 
        return constraints or not. The argument fm is the filemanager 
        that is looking for this default value.

        """

        myfm = None
        if fm.isPresent('TP') and fm.isPresent('DT'):
            myfm = self.tp_dt_fm
        elif fm.isPresent('TP') and fm.isPresent('QT'):
            myfm = self.tp_qt_fm
        elif (fm.isPresent('SS') or fm.isPresent('SA')) and fm.isPresent('DT'):
            myfm = self.mp_dt_fm
        elif (fm.isPresent('SS') or fm.isPresent('SA')) and fm.isPresent('QT'):
            myfm = self.mp_qt_fm

        newLines = []
        if myfm is not None:
            lines = myfm.getLines(tag)
            lines = [x for x in lines if x[0].isConstraint == isConstraint]

            # create a copy of the newlines so the original parsetree 
            # cannot be modified by anyting
            for eachLine in lines:
                newLines.append(copy.deepcopy(eachLine))

        return newLines

    def getValueDefaults(self, fm, analysis=None, trait=None):
        """Returns a list of lines that are the default values.  The argument 
        fm is the filemanager that is looking for the default values. The
        analysis argument is what analysis type is wanted, and the trait
        argument is what trait is wanted.

        """

        # possible values of the argument analysis are 'tp' for two-point and
        # 'mp' for multi-point. possible values of the argument trait are
        # 'dt' for dichtomous, 'qt' for quantative. if any of these are none,
        # then use the current values in the file to determine which defaults
        # to get

        # figure out which filemanager to use based on current values
        myfm = None

        if analysis == None or trait == None:
            # get defaults based on current values in the file
            tp = fm.isPresent('TP')
            ss = fm.isPresent('SS')
            sa = fm.isPresent('SA')
            dt = fm.isPresent('DT')
            qt = fm.isPresent('QT')

            if tp and dt:
                myfm = self.tp_dt_fm
            elif tp and qt:
                myfm = self.tp_qt_fm
            elif (ss or sa) and dt:
                myfm = self.mp_dt_fm
            elif (ss or sa) and qt:
                myfm = self.mp_qt_fm
        else:
            # get defaults based on arguments
            if analysis == 'tp' and trait == 'dt':
                myfm = self.tp_dt_fm
            elif analysis == 'tp' and trait == 'qt':
                myfm = self.tp_qt_fm
            elif analysis == 'mp' and trait == 'dt':
                myfm = self.mp_dt_fm
            elif analysis == 'mp' and trait == 'qt':
                myfm = self.mp_qt_fm

        newLines = []
        if myfm is not None:
            lines = myfm.getVariableNumberLines()

            # create a copy of the newlines so the original parsetree 
            # cannot be modified by anything
            for eachLine in lines:
                newLines.append(copy.deepcopy(eachLine))

        return newLines
