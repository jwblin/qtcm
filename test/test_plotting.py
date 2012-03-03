#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Unittest for the plot module and Qtcm plot* methods.

This module is run using the following command line:

>>> python test_plotting.py

Masked and unmasked numpy data are tested in this module.  Note the
numpy function isscalar is used for all types of data, as it appears
to work, returning False for arguments of any array type, no matter
what the rank.

Besides executing a unittest, the module also generates PNGs of
plots.  These PNGs are in ./rundir/plotm_test and can be visually
compared with files of the same name in ./benchmarks/plotm_test.
All plots with an "a" after the numeral in the name should be the
same as the plot without the numeral, except filled contour (for
contour plots).  plotm1.png and plotm1a.png should be identical,
because it is a line plot (plotm1a.png is created to test that the
filled keyword has no effect for line plots).  plotm6*.png and
plotm9*.png should be the same, with the difference being plotm9*
specifies contour levels by the levels keyword.  likewise for
plotm5*.png vs. plotm8*.png.  The m6*, m9* series tests the situation
for lat vs. lon contour plots, and the m5* vs. m8* series tests the
situation for non-lat vs. lon contour plots.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: test_plotting.py 39 2008-07-15 00:18:10Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, Physics Department, North 
#   Park University.  Passed passably reasonable tests.
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
import qtcm
import qtcm.defaults as defaults
from qtcm.plot import nice_levels
from qtcm import Qtcm
from qtcm import Field
import shutil
import tempfile
import unittest
import utilities




#--------------------------- Module Functions --------------------------

def reload_select_mods():
    """Reload selected modules in sys.modules.
    
    Reloads modules with names beginning with "seaice", "gemath", and
    "modelutil".  This function is used to make sure that all modules 
    that can be compiled for either Numeric or numarray arrays consis-
    tently use the same array type (if that feature is desired).
    """
    #- Make list of modules to redo:

    #@@@currently set so no modules are redone:
    #redo_mods = ['qtcm', 'gemath']
    #redo_mods = ['gemath',]
    redo_mods = []


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
        var_list:  List of netCDF field names (strings).
        array_type:  String, "Numeric", "numarray", "numpy".  Default is
           None.
        """
        self.var_list = None
        self.array_type = None


class NumPyVariablesTestCase(VariablesTestCase):
    def setUp(self):
        import numpy
        self.N = numpy
        self.MA = numpy.ma
        from numpy import isscalar
        self.isscalar = isscalar
        reload_select_mods()
        self.init_data()
        self.array_type = 'numpy'

    def test_AreWeUsingNumPy(self):
        self.failUnlessEqual(self.N.__name__, 'numpy')
        self.failUnless(self.MA.__name__ in ['numpy.core.ma','numpy.ma'])




#--------- Parent Class With Tests Common to All Array Packages --------

class Tests(object):
    """Tests to conduct on qtcm.
    """
    def test_exceptions(self):
        """Test some exceptions for plotm.
        """
        #- Run a 5 day test case for use to test method against:

        model = qtcm.Qtcm(compiled_form='parts', ntout=1, ntouti=1)
        rundirname = 'test'
        dirbasepath = tempfile.mkdtemp()
        model.outdir.value = dirbasepath
        model.runname.value = rundirname
        model.bnddir.value = \
            os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
        model.SSTdir.value = \
            os.path.join( os.getcwd(), 'bnddir', 'r64x42', 'SST_Reynolds' )
        model.lastday.value = 5
        model.run_session()

        self.failUnlessRaises(ValueError, model.plotm,
            'Qc', lat=0.0, lon=[20,100], time=[0,4])

        if os.path.exists(dirbasepath):  shutil.rmtree(dirbasepath)
        if os.path.exists('qtcm_00010101.restart'):
            os.remove('qtcm_00010101.restart')




#--------- Private Methods to Do Things For Manual Comparison  ---------

def _gen_plots():
    """Generate plots for manual comparison.

    This method only works with numpy.
    """
    import numpy as N
    rundirname = 'plotm_test'
    dirbasepath = utilities.prepare_outdir(rundirname)

    inputs = {}
    inputs['dt'] = 1200.
    inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
    inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
    inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42',
                                     'SST_Reynolds' )
    inputs['outdir'] = dirbasepath
    inputs['runname'] = rundirname
    inputs['landon'] = 0
    inputs['year0'] = 1
    inputs['month0'] = 11
    inputs['day0'] = 1
    inputs['lastday'] = 60
    inputs['ntout'] = 1
    inputs['ntouti'] = 1
    inputs['noout'] = 0
    inputs['mrestart'] = 0
    inputs['compiled_form'] = 'full'

    model = Qtcm(**inputs)
    model.run_session()

    model.plotm('Qc',lat=[-40,40], time=58, lon=[0.,],
        fn=os.path.join(dirbasepath, 'plotm1.png'))
    model.plotm('FLWut',lat=[-40,40], time=52,
        fn=os.path.join(dirbasepath, 'plotm2.png'))
    model.plotm('Ts',lon=[20,300], time=50, approx_nlev=15,
        fn=os.path.join(dirbasepath, 'plotm3.png'))
    model.plotm('Qc',lat=[-20, 50], lon=247.5, time=[30,50],
        fn=os.path.join(dirbasepath, 'plotm4.png'))
    model.plotm('us', time=52,
        fn=os.path.join(dirbasepath, 'plotm5.png'))
    model.plotm('Qc',lat=5.625, lon=[50,200], time=[20,50],
        fn=os.path.join(dirbasepath, 'plotm6.png'))
    model.plotm('q1',lat=5.625, lon=112.5, time=[20,50],
        fn=os.path.join(dirbasepath, 'plotm7.png'))
    model.plotm('vs', time=50, levels=N.arange(-10,12,2),
        fn=os.path.join(dirbasepath, 'plotm8.png'))
    model.plotm('Qc',lat=5.625, lon=[50,200], time=[20,50], 
        levels=N.arange(-20, 320, 20),
        fn=os.path.join(dirbasepath, 'plotm9.png'))

    model.plotm('Qc',lat=[-40,40], time=58, lon=[0.,], filled=1,
        fn=os.path.join(dirbasepath, 'plotm1a.png'))
    model.plotm('FLWut',lat=[-40,40], time=52, filled=1,
        fn=os.path.join(dirbasepath, 'plotm2a.png'))
    model.plotm('Ts',lon=[20,300], time=50, approx_nlev=15, filled=1,
        fn=os.path.join(dirbasepath, 'plotm3a.png'))
    model.plotm('Qc',lat=[-20, 50], lon=247.5, time=[30,50], filled=1,
        fn=os.path.join(dirbasepath, 'plotm4a.png'))
    model.plotm('us', time=52, filled=1,
        fn=os.path.join(dirbasepath, 'plotm5a.png'))
    model.plotm('Qc',lat=5.625, lon=[50,200], time=[20,50], filled=True,
        fn=os.path.join(dirbasepath, 'plotm6a.png'))
    model.plotm('vs', time=50, levels=N.arange(-10,12,2), filled=True,
        fn=os.path.join(dirbasepath, 'plotm8a.png'))
    model.plotm('Qc',lat=5.625, lon=[50,200], time=[20,50],
        levels=N.arange(-20, 320, 20), filled=True,
        fn=os.path.join(dirbasepath, 'plotm9a.png'))

    if os.path.exists('qtcm_00011115.restart'):  
        os.remove('qtcm_00011115.restart')




#-------------------------- Create Test Case ---------------------------
#
# Create test suite and add each array type test case if that array 
# package exists.  Also add the doctest unit tests of the docstrings
# of qtcm.plot

suite = unittest.TestSuite()

try:
    import numpy
    class NumPyTests(    NumPyVariablesTestCase,    Tests ):  pass
    suite.addTest(unittest.makeSuite(NumPyTests))
    import qtcm.plot, doctest
    suite.addTest(doctest.DocTestSuite(qtcm.plot))
except: pass



#---------------------------- Main Program -----------------------------
#
# Set variable run_verbose to True or False, depending on what I want
# to do (i.e. verbose testing or not verbose testing).

if __name__ == "__main__":

    #- Generate plots for manual comparison:

    _gen_plots()


    #- Run unittest tests:

    run_verbose = False          #- Set this to control verbosity
    if run_verbose:
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.TextTestRunner(verbosity=1).run(suite)




# ====== end file ======
