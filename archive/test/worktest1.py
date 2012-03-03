#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Unittest for 15 days of output of the qtcm package.

Compares with output from a benchmark run.

Masked and unmasked Numeric, numarray, and numpy output data are tested
in this module.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: package_version.py,v 1.9 2004/06/04 00:43:47 jlin Exp $
#
# Modification History:
# - 14 Jun 2007:  Original by Johnny Lin, North Park University.
#   Passed passably reasonable tests.@@@
#
# Notes:
# - Written for Python 2.4.
# - See import/reload statements throughout module for dependencies.
#
# Copyright (c) 2004-2008 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
#=======================================================================




#---------------- Module General Import and Declarations ---------------

#- If you're importing this module in testing mode, or you're running 
#  pydoc on this module via the command line, import user-specific 
#  settings to make sure any non-standard libraries are found:

import os
import sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user


#- Import other modules at module level:

import __main__
import copy
import unittest
from utilities import read_benchmark
from utilities import read_output




#--------------------------- Module Functions --------------------------

def reload_select_mods():
    """Reload selected modules in sys.modules.
    
    Reloads modules with names beginning with "seaice", "gemath", and
    "modelutil".  This function is used to make sure that all modules 
    that can be compiled for either Numeric or numarray arrays consis-
    tently use the same array type.
    @@@need to revise or maybe remove since don't want gemath to be
    required?
    """
    #- Make list of modules to redo:

    #@@@redo_mods = ['qtcm', 'gemath']
    redo_mods = ['gemath', ]


    #- Cycle through all modules in sys.modules and reload if their
    #  values are not None and have names whose prefix are found in
    #  the redo_mods list:

    for amodname in sys.modules.keys():
        result = [amodname.startswith(i) for i in redo_mods]
        if (result.count(True) != 0) and (sys.modules[amodname] != None):
            reload(sys.modules[amodname])




#------ Parent Classes Setting Up Numeric/numarray etc. Selection ------

class VariablesTestCase(unittest.TestCase):
    """Data initialization common to all array packages
    """
    def init_data(self):
        """Initialize common data.

        Note in this test script, different from many of the other unit-
        tests, the model is NOT instantiated here, but is left to be
        instantiated in the test methods.

        Attributes created:
        * array_type:  String, "Numeric", "numarray", "numpy".  Default
          is None.
        * var_list:  List of netCDF field names (strings).
        """
        self.var_list = [ 'WD', 'GMs1', 'u1', 'cdn', 'S0', 'Ts', 'Prec' \
                        , 'top', 'advT1', 'dfsq1', 'GMq1', 'div0', 'q1' \
                        , 'Evap', 'lon', 'u0', 'OLR', 'div1', 'wet' \
                        , 'advq1', 'FLWds', 'FLWus', 'Runs', 'dfsT1', 'T1' \
                        , 'v0', 'v1', 'vort0', 'Runf', 'lat', 'FTs', 'Evapi' \
                        , 'FSWut', 'FSWus', 'us', 'vs', 'psi0', 'stype' \
                        , 'time', 'taux', 'tauy', 'cl1', 'FSWds' ]
        self.array_type = None


class NumericVariablesTestCase(VariablesTestCase):
    def setUp(self):
        import Numeric
        import MA
        self.N = Numeric
        self.MA = MA
        reload_select_mods()
        self.init_data()
        self.array_type = "Numeric"

    def test_AreWeUsingNumeric(self):
        self.failUnlessEqual(self.N.__name__, 'Numeric')
        self.failUnlessEqual(self.MA.__name__, 'MA')


class NumarrayVariablesTestCase(VariablesTestCase):
    def setUp(self):
        import numarray
        import numarray.ma.MA
        self.N = numarray
        self.MA = numarray.ma.MA
        reload_select_mods()
        self.init_data()
        self.array_type = 'numarray'

    def test_AreWeUsingNumarray(self):
        self.failUnlessEqual(self.N.__name__, 'numarray')
        self.failUnlessEqual(self.MA.__name__, 'numarray.ma.MA')


class NumPyVariablesTestCase(VariablesTestCase):
    def setUp(self):
        import numpy
        self.N = numpy
        self.MA = numpy.ma
        reload_select_mods()
        self.init_data()
        self.array_type = 'numpy'

    def test_AreWeUsingNumpy(self):
        self.failUnlessEqual(self.N.__name__, 'numpy')
        self.failUnlessEqual(self.MA.__name__, 'numpy.core.ma')




#-------- Parent Class With Tests Common to All Array Packages ---------

class Tests(object):
    """Tests to conduct on qtcm.
    """
    def test_exceptions(self):
        """Test some exceptions.
        """
        #@@@- Check places with snow cannot be places of open water:

        #setcopy = self.set.copy()
        #setcopy['t0'].add( \
        #    Variable( A.array([0.0, 0.2, 0.25, 0.25, 0.25, 0.25]) \
        #         , id='hI', units='m' ) )
        #self.failUnlessRaises(ValueError, Semtner0, setcopy)
        pass




#------------------- Class Tests:  Test Single Runs --------------------

    def test_all_fields_unmasked_full_aquaplanet(self):
        """Test all aquaplanet output fields match the benchmark.

        Test uses the assumes self.compiled_form is 'full'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'aquaplanet' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'full_365_aqua' \
                                  , array_type=self.array_type )

            self.failUnless( self.N.allclose(data[0], myoutput[0]) )
            self.failUnless( self.N.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_full_landon(self):
        """Test all landon output fields match the benchmark.

        Test uses the assumes self.compiled_form is 'full'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'landon' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'full_365_landon' \
                                  , array_type=self.array_type )

            self.failUnless( self.N.allclose(data[0], myoutput[0]) )
            self.failUnless( self.N.allclose(data[1], myoutput[1]) )




#---------------------------- Main Program -----------------------------
#
# Set variable run_verbose to True or False, depending on what I want
# to do (i.e. verbose testing or not verbose testing).

if __name__ == "__main__":
    #- Run model via the OS (Unix).  I have to do this because using
    #  execfile to execute the scripts, for some strange reason, does
    #  not work on my OS.  The settings persist from one run to the
    #  other:

    print 'Running model to obtain output (this may take a while) ...'
    os.system('python2.4 create_full_aqua.py > create.log')
    os.system('python2.4 create_parts_aqua.py >> create.log')
    print 'Running tests on output ...'


    run_verbose = False          #- Set this to control verbosity
    suite = unittest.TestSuite()

    #- Add each array type test case if that array package exists:
    #  @@@need comment back in Numeric and numarray tests once i run 
    #     the model using those array packages (so far in the run_short
    #     line above, i run the model using only numpy, so i only test 
    #     for numpy in this test script.
    #
    #try:
    #    import Numeric
    #    class NumericTests(  NumericVariablesTestCase , Tests ):  pass
    #    suite.addTest(unittest.makeSuite(NumericTests))
    #except: pass
    #
    #try:
    #    import numarray
    #    class NumarrayTests( NumarrayVariablesTestCase, Tests ):  pass
    #    suite.addTest(unittest.makeSuite(NumarrayTests))
    #except: pass

    try:
        import numpy
        class NumPyTests(    NumPyVariablesTestCase,    Tests ):  pass
        suite.addTest(unittest.makeSuite(NumPyTests))
    except: pass

    if run_verbose:
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.TextTestRunner(verbosity=1).run(suite)




# ===== end file =====
