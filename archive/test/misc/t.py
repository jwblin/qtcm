#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Unittest for the plot module and Qtcm plot* methods.

This module is run using the following command line:

>>> python test_plot.py

Masked and unmasked numpy data are tested in this module.  Note the
numpy function isscalar is used for all types of data, as it appears
to work, returning False for arguments of any array type, no matter
what the rank.

Besides executing a unittest, the module also generates PNGs of
plots.  These PNGs are in ./rundir/plotm_test and can be visually
compared with files of the same name in ./benchmarks/plotm_test.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: t.py 6 2008-06-25 19:34:56Z jlin $
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
        self.failUnlessEqual(self.MA.__name__, 'numpy.core.ma')




#--------- Parent Class With Tests Common to All Array Packages --------

class Tests(object):
    """Tests to conduct on qtcm.
    """
    def test_exceptions(self):
        """Test some exceptions.
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


    def test_nice_levels(self):
        """Test nice_levels function.
        """
        zdata = self.N.array([[4., 128.],[32., 287.]])
        self.failUnless( self.N.allclose(self.N.shape(zdata), (2,2)) )
        self.failUnless( self.N.allclose(nice_levels(zdata),
                                         self.N.arange(0,300,20)) )




#-------------------------- Create Test Case ---------------------------
#
# Create test suite and add each array type test case if that array 
# package exists.

suite = unittest.TestSuite()

try:
    import numpy
    class NumPyTests(    NumPyVariablesTestCase,    Tests ):  pass
    suite.addTest(unittest.makeSuite(NumPyTests))
except: pass




#---------------------------- Main Program -----------------------------
#
# Set variable run_verbose to True or False, depending on what I want
# to do (i.e. verbose testing or not verbose testing).

if __name__ == "__main__":

    #- Generate plots for manual comparison:

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
    inputs['lastday'] = 15
    inputs['ntout'] = 1
    inputs['ntouti'] = 1
    inputs['noout'] = 0
    inputs['mrestart'] = 0
    inputs['compiled_form'] = 'full'

    model = Qtcm(**inputs)
    model.run_session()

#    model.plotm('Qc',lat=[-40,40], time=58, lon=[0.,],
#        fn=os.path.join(dirbasepath, 'plotm1.png'))
    model.plotm('FLWut',lat=[-40,40], time=12, filled=0,
        fn=os.path.join(dirbasepath, 'plotm2.png'))
#    model.plotm('Ts',lon=[20,300], time=50, approx_nlev=5,
#        fn=os.path.join(dirbasepath, 'plotm3.png'))
#   model.plotm('Qc',time=4, tmppreview=True)
#   model.plotm('Qc',lon=[20,100], time=[0,4])  #@@@should give error
#   model.plotm('Qc',lat=5.625, lon=[20,100], time=[0,4])
#    model.plotm('Qc',lat=[-20, 50], lon=247.5, time=[30,50],
#        fn=os.path.join(dirbasepath, 'plotm4.png'))
#    model.plotm('us', time=52,
#        fn=os.path.join(dirbasepath, 'plotm5.png'))

    if os.path.exists('qtcm_00011230.restart'):  
        os.remove('qtcm_00011230.restart')


    #- Run unittest tests:

#    run_verbose = False          #- Set this to control verbosity
#    if run_verbose:
#        unittest.TextTestRunner(verbosity=2).run(suite)
#    else:
#        unittest.TextTestRunner(verbosity=1).run(suite)




# ====== end file ======
