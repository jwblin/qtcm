#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Unittest for the default Qtcm object properties.

This module is run using the following command line:

>>> python test_obj_defaults.py

Masked and unmasked Numeric, numarray, and numpy data are tested in 
this module.  Note the numpy function isscalar is used for all types
of data, as it appears to work, returning False for arguments of
any array type, no matter what the rank.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: test_obj_defaults.py 55 2008-08-05 19:03:29Z jlin $
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
import unittest
import qtcm
import qtcm.defaults as defaults
from qtcm import Qtcm
from qtcm.field import Field
import shutil
import utilities




#--------------------------- Module Functions --------------------------

def reload_select_mods():
    """Reload selected modules in sys.modules.
    
    Reloads modules with names beginning with "seaice", "gemath", and
    "modelutil".  This function is used to make sure that all modules 
    that can be compiled for either Numeric or numarray arrays consis-
    tently use the same array type.
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


class NumericVariablesTestCase(VariablesTestCase):
    def setUp(self):
        import Numeric
        import MA
        self.N = Numeric
        self.MA = MA
        from numpy import isscalar
        self.isscalar = isscalar
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
        from numpy import isscalar
        self.isscalar = isscalar
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
        from numpy import isscalar
        self.isscalar = isscalar
        reload_select_mods()
        self.init_data()
        self.array_type = 'numpy'

    def test_AreWeUsingNumarray(self):
        self.failUnlessEqual(self.N.__name__, 'numpy')
        self.failUnless(self.MA.__name__ in ['numpy.core.ma','numpy.ma'])




#--------- Parent Class With Tests Common to All Array Packages --------

class Tests(object):
    """Tests to conduct on qtcm.
    """
    def test_exceptions(self):
        """Test some exceptions.
        """
	#- Test the case where compiled_form is missing raises an
	#  error:

        model = Qtcm(compiled_form='parts')
        delattr(model, 'compiled_form')
        self.failUnlessRaises(ValueError, model.__init__)


	#- Test that varinit cannot be called in certain cases:

        model = Qtcm(compiled_form='full')
        self.failUnlessRaises(AttributeError, model.varinit)

        model = Qtcm(compiled_form='parts')
        model._cont = True
        self.failUnlessRaises(qtcm.qtcm.FieldNotReadableFromCompiledModel, 
            model.varinit)


        #- Test that arrays are not initialized on instantiation, but
        #  rather in a run session (do case with no timesteps):

        model = Qtcm(compiled_form='parts')
        for ikey in defaults.qtcm_fields_ids:
            if not self.isscalar(defaults.qtcm_fields[ikey]['value']):
                self.failUnlessRaises(AttributeError, getattr, model, ikey)

        rundirname = 'test'
        dirbasepath = utilities.prepare_outdir(rundirname)
        model.outdir.value = dirbasepath
        model.runname.value = rundirname
        model.bnddir.value = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
        model.SSTdir.value = \
            os.path.join( os.getcwd(), 'bnddir', 'r64x42', 'SST_Reynolds' )
        model.lastday.value = 0
        model.run_session()
        for ikey in defaults.qtcm_fields_ids:
            if not self.isscalar(defaults.qtcm_fields[ikey]['value']):
                self.failUnless( hasattr(model, ikey) )

        if os.path.exists(dirbasepath):  shutil.rmtree(dirbasepath)
        if os.path.exists('qtcm_00010101.restart'):  
            os.remove('qtcm_00010101.restart')


    def test_Qtcm_object_init_keywords_read(self):
        """Test Qtcm object reads instantiation keywords correctly.

        Test done only for a few cases.  Set compiled_form to 'full' 
        and 'parts'.
        """
        for iform in ['full', 'parts']:
            inputs = {}
            inputs['dt'] = 600.
            inputs['eps_c'] = 0.15
            inputs['compiled_form'] = iform
            model = Qtcm(**inputs)
            self.failUnlessEqual( model.dt.id, 'dt' )
            self.failUnlessEqual( model.dt.value, 600. )
            self.failUnlessEqual( model.eps_c.id, 'eps_c' )
            self.failUnlessEqual( model.eps_c.value, 0.15 )


    def test_Qtcm_object_field_default_on_init(self):
        """Test Qtcm object field attributes on init.

	This test looks a number of Field attributes for a host of
	QTCM fields initialized when the Qtcm instance is initialized.
	It does not test all those default attributes, but for the
        ones it tests, it tests all the Field attributes that describe
        that QTCM field.  In the next test, the fields that are not in
        this test are tested, but only for the correct value.  Set
        compiled_form to 'full' and 'parts'.
        """
        for iform in ['full', 'parts']:
            model = Qtcm(compiled_form=iform)

            self.failUnless( self.N.allclose(model.dt.value, 1200.) )
            self.failUnlessEqual( model.dt.id, 'dt' )
            self.failUnlessEqual( model.dt.units, 's' )
            self.failUnlessEqual( model.dt.long_name, 'time step' )

            self.failUnless( self.N.allclose(model.eps_c.value, 
                0.13888889E-03) )
            self.failUnlessEqual( model.eps_c.id, 'eps_c' )
            self.failUnlessEqual( model.eps_c.units, '1/s' )
            self.failUnlessEqual( model.eps_c.long_name, '1/tau_c NZ (5.7)' )

            self.failUnlessEqual( model.title.value, 'QTCM default title' )
            self.failUnlessEqual( model.title.id, 'title' )
            self.failUnlessEqual( model.title.units, '' )
            self.failUnlessEqual( model.title.long_name, 
                'a descriptive title' )

            self.failUnlessEqual( model.bnddir.value, 
                os.path.join(os.pardir, 'bnddata') )
            self.failUnlessEqual( model.bnddir.id, 'bnddir' )
            self.failUnlessEqual( model.bnddir.units, '' )
            self.failUnlessEqual( model.bnddir.long_name,
                'boundary data other than SST' )

            self.failUnlessEqual( model.SSTdir.value, 
                os.path.join(os.pardir, 'bnddata', 'SST_Reynolds') )
            self.failUnlessEqual( model.SSTdir.id, 'SSTdir' )
            self.failUnlessEqual( model.SSTdir.units, '' )
            self.failUnlessEqual( model.SSTdir.long_name, 
                'where SST files are' )

            self.failUnlessEqual( model.outdir.value,
                os.path.join(os.pardir, 'proc', 'qtcm_output') )
            self.failUnlessEqual( model.outdir.id, 'outdir' )
            self.failUnlessEqual( model.outdir.units, '' )
            self.failUnlessEqual( model.outdir.long_name, 
                'where output goes to' )

            self.failUnlessEqual( model.runname.value, 'runname' )
            self.failUnlessEqual( model.runname.id, 'runname' )
            self.failUnlessEqual( model.runname.units, '' )
            self.failUnlessEqual( model.runname.long_name,
                'string for an output filename' )

            self.failUnlessEqual( model.landon.value, 1 )
            self.failUnlessEqual( model.landon.id, 'landon' )
            self.failUnlessEqual( model.landon.units, '' )
            self.failUnlessEqual( model.landon.long_name,
                'if not 1: land = ocean with fake SST' )

            self.failUnlessEqual( model.SSTmode.value, 'seasonal' )
            self.failUnlessEqual( model.SSTmode.id, 'SSTmode' )
            self.failUnlessEqual( model.SSTmode.units, '' )
            self.failUnlessEqual( model.SSTmode.long_name,
                'decide what kind of SST to use' )

            self.failUnlessEqual( model.year0.value, 0 )
            self.failUnlessEqual( model.year0.id, 'year0' )
            self.failUnlessEqual( model.year0.units, 'yr' )
            self.failUnlessEqual( model.year0.long_name,
                'starting year; if < 0 use year in restart' )

            self.failUnlessEqual( model.month0.value, -1 )
            self.failUnlessEqual( model.month0.id, 'month0' )
            self.failUnlessEqual( model.month0.units, 'mo' )
            self.failUnlessEqual( model.month0.long_name,
                'starting month; if < 0 use mo in restart' )

            self.failUnlessEqual( model.day0.value, -1 )
            self.failUnlessEqual( model.day0.id, 'day0' )
            self.failUnlessEqual( model.day0.units, 'dy' )
            self.failUnlessEqual( model.day0.long_name,
                'starting day; if < 0 use day in restart' )

            self.failUnlessEqual( model.lastday.value, 365 )
            self.failUnlessEqual( model.lastday.id, 'lastday' )
            self.failUnlessEqual( model.lastday.units, 'dy' )
            self.failUnlessEqual( model.lastday.long_name,
                'last day of integration' )

            self.failUnlessEqual( model.dateofmodel.value, 0 )
            self.failUnlessEqual( model.dateofmodel.id, 'dateofmodel' )
            self.failUnlessEqual( model.dateofmodel.units, '' )
            self.failUnlessEqual( model.dateofmodel.long_name,
                'date of model coded as an integer as yyyymmdd' ) 


    def test_Qtcm_object_field_default_values_on_init(self):
	"""Test Qtcm object field values on initialization for full.

	This test looks at Field attributes not tested above to
	make sure that their values are initialized correctly.  But
	on the value in Field is tested, none of the other attributes
	of Field.  Set compiled_form to 'full' and 'parts'.
        """
        for iform in ['full', 'parts']:
            model = Qtcm(compiled_form=iform)
            self.failUnless( self.N.allclose(model.it.value, 1) )
            self.failUnless( self.N.allclose(model.interval.value, 1) )
            self.failUnless( self.N.allclose(model.noout.value, 0) )
            self.failUnless( self.N.allclose(model.nooutr.value, 0) )
            self.failUnless( self.N.allclose(model.ntout.value, -30) )
            self.failUnless( self.N.allclose(model.ntouti.value, 0) )
            self.failUnless( self.N.allclose(model.ntoutr.value, 0) )
            self.failUnless( self.N.allclose(model.mrestart.value, 0) )
            self.failUnless( self.N.allclose(model.mt0.value, 1) )
            self.failUnless( self.N.allclose(model.ziml.value, 500) )
            self.failUnless( self.N.allclose(model.weml.value, 0.01) )
            self.failUnless( self.N.allclose(model.VVsmin.value, 4.5) )
            self.failUnless( self.N.allclose(model.viscxu0.value, 7e5) )
            self.failUnless( self.N.allclose(model.viscyu0.value, 7e5) )
            self.failUnless( self.N.allclose(model.visc4x.value, 7e5) )
            self.failUnless( self.N.allclose(model.visc4y.value, 7e5) )
            self.failUnless( self.N.allclose(model.viscxu1.value, 7e5) )
            self.failUnless( self.N.allclose(model.viscyu1.value, 7e5) )
            self.failUnless( self.N.allclose(model.viscxT.value, 12e5) )
            self.failUnless( self.N.allclose(model.viscyT.value, 12e5) )
            self.failUnless( self.N.allclose(model.viscxq.value, 12e5) )
            self.failUnless( self.N.allclose(model.viscyq.value, 12e5) )


    def test_Qtcm_object_set_and_get_items_from_compiled(self):
        """Test Qtcm set_qtcm_item and get_qtcm_item methods together.

        Also tests the set_qtcm1_item and get_qtcm1_item methods, which
        should behave identically to set_qtcm_item and get_qtcm_item.
        Only scalar parameters are tested.
        """
        for iform in ['full', 'parts']:
            model = Qtcm(compiled_form=iform)
            model.set_qtcm_item('dt', 2400.)
            model.set_qtcm_item('bnddir', 'ooga booga')
            model.set_qtcm_item('ntout', 120)
        
            self.failUnless( self.N.allclose(model.get_qtcm_item('dt'), 2400.) )
            self.failUnlessEqual( model.get_qtcm_item('bnddir'), 'ooga booga' )
            self.failUnlessEqual( model.get_qtcm_item('ntout'), 120 )
        
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'dt', 400)
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'dt', 'hi')
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'bnddir', 20)
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'bnddir', 20.)
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'ntout', 30.)
            self.failUnlessRaises(TypeError, model.set_qtcm_item, 'ntout', 'hi')
            del model

            model = Qtcm(compiled_form=iform)
            model.set_qtcm1_item('dt', 2400.)
            model.set_qtcm1_item('bnddir', 'ooga booga')
            model.set_qtcm1_item('ntout', 120)
        
            self.failUnless( self.N.allclose(model.get_qtcm1_item('dt'), 
                             2400.) )
            self.failUnlessEqual( model.get_qtcm1_item('bnddir'), 'ooga booga' )
            self.failUnlessEqual( model.get_qtcm1_item('ntout'), 120 )
        
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'dt', 400 )
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'dt', 'hi' )
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'bnddir',
                                   20 )
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'bnddir', 
                                   20. )
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'ntout', 
                                   30. )
            self.failUnlessRaises( TypeError, model.set_qtcm1_item, 'ntout', 
                                   'hi' )
            del model




#-------------------------- Create Test Case ---------------------------
#
# Create test suite and add each array type test case if that array 
# package exists.

suite = unittest.TestSuite()

try:
    import Numeric
    class NumericTests(  NumericVariablesTestCase , Tests ):  pass
    suite.addTest(unittest.makeSuite(NumericTests))
except: pass

try:
    import numarray
    class NumarrayTests( NumarrayVariablesTestCase, Tests ):  pass
    suite.addTest(unittest.makeSuite(NumarrayTests))
except: pass

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

    #- Run tests:

    run_verbose = False          #- Set this to control verbosity
    if run_verbose:
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.TextTestRunner(verbosity=1).run(suite)




# ===== end file =====
