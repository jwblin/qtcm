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
del os


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
    def allclose(self, x, y, rtol=1000., atol=10000.):
        return self.N.allclose(x, y, rtol=rtol, atol=atol)


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




#----------------- Class Tests:  Test Plain Benchmarks -----------------

    def test_all_fields_benchmarks_are_all_different(self):
        """Test all benchmark outputs are different for sub-var_list.

	Only select variables are compared.  This is because the
	slight differences in parameters between aquaplanet and
	aquaplanet_params2 does not propagate to all fields in the
	short time (~15 days) the model is run for these benchmark
	tests.
        """
        sub_var_list = [ 'u1', 'Prec', 'advT1', 'dfsq1', 'GMq1', 'q1' \
                       , 'Evap', 'u0', 'OLR', 'div1', 'advq1', 'FLWds' \
                       , 'dfsT1', 'T1', 'v0', 'v1', 'vort0', 'FTs' \
                       , 'FSWut', 'FSWus', 'us', 'vs', 'psi0' \
                       , 'taux', 'tauy', 'cl1', 'FSWds' ]

        for ivar in sub_var_list:
            data1 = read_benchmark( ivar, 'aquaplanet' \
                                  , array_type=self.array_type )
            data2 = read_benchmark( ivar, 'landon' \
                                  , array_type=self.array_type )
            data3 = read_benchmark( ivar, 'aquaplanet_params2' \
                                  , array_type=self.array_type )

            self.failUnless( not self.allclose(data1[0], data2[0]) )
            self.failUnless( not self.allclose(data1[1], data2[1]) )
            self.failUnless( not self.allclose(data1[0], data3[0]) )
            self.failUnless( not self.allclose(data1[1], data3[1]) )
            self.failUnless( not self.allclose(data2[0], data3[0]) )
            self.failUnless( not self.allclose(data2[1], data3[1]) )




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

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_full_landon(self):
        """Test all landon output fields match the benchmark.

        Test uses the assumes self.compiled_form is 'full'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'landon' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'full_365_landon' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_full_aqua_multi_same(self):
        """Test aqua output match the benchmark for multi same case.

	These are for runs where the results should be identical
	to the benchmark for a default run.  Test uses the assumes
	self.compiled_form is 'full'.
        """
        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'aquaplanet' \
                                      , array_type=self.array_type )
            myoutput1 = read_output( ivar, 'full_365_aqua_multi1' \
                                   , array_type=self.array_type )
            myoutput2 = read_output( ivar, 'full_365_aqua_multi2' \
                                   , array_type=self.array_type )
            myoutput4 = read_output( ivar, 'full_365_aqua_multi4' \
                                   , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput1[0]) )
            self.failUnless( self.allclose(data[1], myoutput1[1]) )
            self.failUnless( self.allclose(data[0], myoutput2[0]) )
            self.failUnless( self.allclose(data[1], myoutput2[1]) )
            self.failUnless( self.allclose(data[0], myoutput4[0]) )
            self.failUnless( self.allclose(data[1], myoutput4[1]) )


    def test_all_fields_unmasked_full_aqua_multi_params2(self):
        """Test aqua output match the benchmark for multi params2 case.

	These are for runs where the results should be identical
	to the benchmark for an aquaplanet_params2 run.  Test uses
	the assumes self.compiled_form is 'full'.
        """
        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'aquaplanet_params2' \
                                      , array_type=self.array_type )
            myoutput3 = read_output( ivar, 'full_365_aqua_multi3' \
                                   , array_type=self.array_type )
            self.failUnless( self.allclose(data[0], myoutput3[0]) )
            self.failUnless( self.allclose(data[1], myoutput3[1]) )


    def test_all_fields_unmasked_parts_aquaplanet(self):
        """Test all aquaplanet output fields match the benchmark.

        Test uses the assumes self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'aquaplanet' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_aqua' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_parts_aqua_inst(self):
        """Test aqua output match the benchmark for inst. case.

        The "inst. case" means testing that running the model with
        the init_with_instance_state keyword set to True will work
        properly.  Test assumes self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'aquaplanet' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_aqua_inst1' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_parts_landon(self):
        """Test all landon output fields match the benchmark.

        Test uses the assumes self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data     = read_benchmark( ivar, 'landon' \
                                     , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_landon' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )


    def test_all_fields_unmasked_parts_aqua_multi_same(self):
        """Test aqua output match the benchmark for multi same case.

	These are for runs where the results should be identical
	to the benchmark for a default run.  Test uses the assumes
	self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'aquaplanet' \
                                      , array_type=self.array_type )
            myoutput1 = read_output( ivar, 'parts_365_aqua_multi1' \
                                   , array_type=self.array_type )
            myoutput2 = read_output( ivar, 'parts_365_aqua_multi2' \
                                   , array_type=self.array_type )
            myoutput4 = read_output( ivar, 'parts_365_aqua_multi4' \
                                   , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput1[0]) )
            self.failUnless( self.allclose(data[1], myoutput1[1]) )
            self.failUnless( self.allclose(data[0], myoutput2[0]) )
            self.failUnless( self.allclose(data[1], myoutput2[1]) )
            self.failUnless( self.allclose(data[0], myoutput4[0]) )
            self.failUnless( self.allclose(data[1], myoutput4[1]) )


    def test_all_fields_unmasked_parts_aqua_multi_params2(self):
        """Test aqua output match the benchmark for multi params2 case.

	These are for runs where the results should be identical
	to the benchmark for an aquaplanet_params2 run.  Test uses
	the assumes self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'aquaplanet_params2' \
                                      , array_type=self.array_type )
            myoutput3 = read_output( ivar, 'parts_365_aqua_multi3' \
                                   , array_type=self.array_type )
            self.failUnless( self.allclose(data[0], myoutput3[0]) )
            self.failUnless( self.allclose(data[1], myoutput3[1]) )


    def test_all_fields_unmasked_parts_landon_multi(self):
        """Test all landon output match the benchmark for multi case.

        Test uses the assumes self.compiled_form is 'parts'.
        """
        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'landon' \
                                      , array_type=self.array_type )
            myoutput1 = read_output( ivar, 'parts_365_landon_multi1' \
                                   , array_type=self.array_type )
            myoutput2 = read_output( ivar, 'parts_365_landon_multi2' \
                                   , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput1[0]) )
            self.failUnless( self.allclose(data[1], myoutput1[1]) )
            self.failUnless( self.allclose(data[0], myoutput2[0]) )
            self.failUnless( self.allclose(data[1], myoutput2[1]) )




#---------------- Class Tests:  Test Continuation Runs -----------------

    def test_all_fields_aqua_cont_benchmarks_different(self):
        """Test aqua output for continuation benchmarks are different.

	The restart system should give results different than doing
	a full run without stopping.  This test checks whether this
	is the case, but only select variables are compared.
        """
        sub_var_list = [ 'u1', 'Prec', 'advT1', 'dfsq1', 'GMq1', 'q1' \
                       , 'Evap', 'u0', 'OLR', 'div1', 'advq1', 'FLWds' \
                       , 'dfsT1', 'T1', 'v0', 'v1', 'vort0', 'FTs' \
                       , 'FSWut', 'FSWus', 'us', 'vs', 'psi0' \
                       , 'taux', 'tauy', 'cl1', 'FSWds' ]

        for ivar in sub_var_list:
            data1 = read_benchmark( ivar, 'aquaplanet_cont40' \
                                  , array_type=self.array_type )
            data2 = read_benchmark( ivar, 'aquaplanet_cont10-30_10' \
                                  , array_type=self.array_type )
            data3 = read_benchmark( ivar, 'aquaplanet_cont10-30_30' \
                                  , array_type=self.array_type )

            myoutput0 = self.N.concatenate([data2[0], data3[0]], axis=0)
            myoutput1 = self.N.concatenate([data2[1], data3[1]], axis=0)
            self.failUnless( not self.allclose(data1[0], myoutput0) )
            self.failUnless( not self.allclose(data1[1], myoutput1) )


    def test_all_fields_aqua_restart_continuation(self):
        """Test aqua output for "restart" continuation runs.

	These are for runs where the results should be identical
	to the benchmark for a default run.
        """
        for ivar in self.var_list:
            data1a = read_benchmark( ivar, 'aquaplanet_cont10-30_10' \
                                   , array_type=self.array_type )
            data1b = read_benchmark( ivar, 'aquaplanet_cont10-30_30' \
                                   , array_type=self.array_type )
            myoutput1a = read_output( ivar, 'full_365_aqua_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput1b = read_output( ivar, 'full_365_aqua_cont10-30_30' \
                                    , array_type=self.array_type )
            myoutput2a = read_output( ivar, 'parts_365_aqua_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput2b = read_output( ivar, 'parts_365_aqua_cont10-30_30' \
                                    , array_type=self.array_type )
            myoutput3a = read_output( ivar, 'parts_365_aqua_nr_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput3b = read_output( ivar, 'parts_365_aqua_nr_cont10-30_30' \
                                    , array_type=self.array_type )

            self.failUnless( self.allclose(data1a[0], myoutput1a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput1a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput1b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput1b[1]) )

            self.failUnless( self.allclose(data1a[0], myoutput2a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput2a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput2b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput2b[1]) )

            self.failUnless( self.allclose(data1a[0], myoutput3a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput3a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput3b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput3b[1]) )


    def test_all_fields_aqua_no_restart_continuation(self):
        """Test aqua output for no restart continuation runs.

	These are for runs where the results should be identical
	to the benchmark for a default run.  A case where nothing
        is changed between the continuations and where u1 is changed
        between continuations is done.  This test also tests to make
        sure for a subset of variables that the case where u1 is
        changed is different from the case where it is not changed,
        prior to the continuation.
        """
        #- Check Qtcm runs match benchmarks:

        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'aquaplanet_cont40' \
                                      , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_aqua_nr1_cont40' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )

            datau1    = read_benchmark( ivar, 'aquaplanet_cont40u1' \
                                      , array_type=self.array_type )
            myoutputu1 = read_output( ivar, 'parts_365_aqua_nr2_cont40' \
                                   , array_type=self.array_type )

            self.failUnless( self.allclose(datau1[0], myoutputu1[0]) )
            self.failUnless( self.allclose(datau1[1], myoutputu1[1]) )


        #- Check continue without change and continue with change are
        #  different:

        sub_var_list = [ 'u1', 'Prec', 'advT1', 'dfsq1', 'GMq1', 'q1' \
                       , 'Evap', 'u0', 'OLR', 'div1', 'advq1', 'FLWds' \
                       , 'dfsT1', 'T1', 'v0', 'v1', 'vort0', 'FTs' \
                       , 'FSWut', 'FSWus', 'us', 'vs', 'psi0' \
                       , 'taux', 'tauy', 'cl1', 'FSWds' ]

        for ivar in sub_var_list:
            data      = read_benchmark( ivar, 'aquaplanet_cont40' \
                                      , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_aqua_nr1_cont40' \
                                  , array_type=self.array_type )
            datau1    = read_benchmark( ivar, 'aquaplanet_cont40u1' \
                                      , array_type=self.array_type )
            myoutputu1 = read_output( ivar, 'parts_365_aqua_nr2_cont40' \
                                   , array_type=self.array_type )

            self.failUnless( not self.allclose(datau1[0], data[0]) )
            self.failUnless( not self.allclose(datau1[1], data[1]) )
            self.failUnless( not self.allclose(myoutputu1[0], myoutput[0]) )
            self.failUnless( not self.allclose(myoutputu1[1], myoutput[1]) )


    def test_all_fields_landon_cont_benchmarks_different(self):
        """Test landon output for continuation benchmarks are different.

	The restart system should give results different than doing
	a full run without stopping.  This test checks whether this
	is the case, but only select variables are compared.
        """
        sub_var_list = [ 'u1', 'Prec', 'advT1', 'dfsq1', 'GMq1', 'q1' \
                       , 'Evap', 'u0', 'OLR', 'div1', 'advq1', 'FLWds' \
                       , 'dfsT1', 'T1', 'v0', 'v1', 'vort0', 'FTs' \
                       , 'FSWut', 'FSWus', 'us', 'vs', 'psi0' \
                       , 'taux', 'tauy', 'cl1', 'FSWds' ]

        for ivar in sub_var_list:
            data1 = read_benchmark( ivar, 'landon_cont40' \
                                  , array_type=self.array_type )
            data2 = read_benchmark( ivar, 'landon_cont10-30_10' \
                                  , array_type=self.array_type )
            data3 = read_benchmark( ivar, 'landon_cont10-30_30' \
                                  , array_type=self.array_type )

            myoutput0 = self.N.concatenate([data2[0], data3[0]], axis=0)
            myoutput1 = self.N.concatenate([data2[1], data3[1]], axis=0)
            self.failUnless( not self.allclose(data1[0], myoutput0) )
            self.failUnless( not self.allclose(data1[1], myoutput1) )


    def test_all_fields_landon_restart_continuation(self):
        """Test landon output for "restart" continuation runs.

	These are for runs where the results should be identical
	to the benchmark for a default run.
        """
        for ivar in self.var_list:
            data1a = read_benchmark( ivar, 'landon_cont10-30_10' \
                                   , array_type=self.array_type )
            data1b = read_benchmark( ivar, 'landon_cont10-30_30' \
                                   , array_type=self.array_type )
            myoutput1a = read_output( ivar, 'full_365_landon_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput1b = read_output( ivar, 'full_365_landon_cont10-30_30' \
                                    , array_type=self.array_type )
            myoutput2a = read_output( ivar, 'parts_365_landon_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput2b = read_output( ivar, 'parts_365_landon_cont10-30_30' \
                                    , array_type=self.array_type )
            myoutput3a = read_output( ivar, 'parts_365_landon_nr_cont10-30_10' \
                                    , array_type=self.array_type )
            myoutput3b = read_output( ivar, 'parts_365_landon_nr_cont10-30_30' \
                                    , array_type=self.array_type )

            self.failUnless( self.allclose(data1a[0], myoutput1a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput1a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput1b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput1b[1]) )

            self.failUnless( self.allclose(data1a[0], myoutput2a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput2a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput2b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput2b[1]) )

            self.failUnless( self.allclose(data1a[0], myoutput3a[0]) )
            self.failUnless( self.allclose(data1a[1], myoutput3a[1]) )
            self.failUnless( self.allclose(data1b[0], myoutput3b[0]) )
            self.failUnless( self.allclose(data1b[1], myoutput3b[1]) )


    def test_all_fields_landon_no_restart_continuation(self):
        """Test landon output for no restart continuation runs.

	These are for runs where the results should be identical
	to the benchmark for a default run.  A case where nothing
        is changed between the continuations and where u1 is changed
        between continuations is done.  This test also tests to make
        sure for a subset of variables that the case where u1 is
        changed is different from the case where it is not changed,
        prior to the continuation.
        """
        #- Check Qtcm runs match benchmarks:

        for ivar in self.var_list:
            data      = read_benchmark( ivar, 'landon_cont40' \
                                      , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_landon_nr1_cont40' \
                                  , array_type=self.array_type )

            self.failUnless( self.allclose(data[0], myoutput[0]) )
            self.failUnless( self.allclose(data[1], myoutput[1]) )

            datau1    = read_benchmark( ivar, 'landon_cont40u1' \
                                      , array_type=self.array_type )
            myoutputu1 = read_output( ivar, 'parts_365_landon_nr2_cont40' \
                                   , array_type=self.array_type )

            self.failUnless( self.allclose(datau1[0], myoutputu1[0]) )
            self.failUnless( self.allclose(datau1[1], myoutputu1[1]) )


        #- Check continue without change and continue with change are
        #  different:

        sub_var_list = [ 'u1', 'Prec', 'advT1', 'dfsq1', 'GMq1', 'q1' \
                       , 'Evap', 'u0', 'OLR', 'div1', 'advq1', 'FLWds' \
                       , 'dfsT1', 'T1', 'v0', 'v1', 'vort0', 'FTs' \
                       , 'FSWut', 'FSWus', 'us', 'vs', 'psi0' \
                       , 'taux', 'tauy', 'cl1', 'FSWds' ]

        for ivar in sub_var_list:
            data      = read_benchmark( ivar, 'landon_cont40' \
                                      , array_type=self.array_type )
            myoutput = read_output( ivar, 'parts_365_landon_nr1_cont40' \
                                  , array_type=self.array_type )
            datau1    = read_benchmark( ivar, 'landon_cont40u1' \
                                      , array_type=self.array_type )
            myoutputu1 = read_output( ivar, 'parts_365_landon_nr2_cont40' \
                                   , array_type=self.array_type )

            self.failUnless( not self.allclose(datau1[0], data[0]) )
            self.failUnless( not self.allclose(datau1[1], data[1]) )
            self.failUnless( not self.allclose(myoutputu1[0], myoutput[0]) )
            self.failUnless( not self.allclose(myoutputu1[1], myoutput[1]) )




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
    #@@@execfile('run_short.py')
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
