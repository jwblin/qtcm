#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Run all unittests for the qtcm package.

This module is run using the following command line:

>>> python test_all.py

It assumes that all unit testing scripts (named test_*.py) and
lengthier model creation scripts (named run_*.py) are in the same
directory as this script.

This file also contains execution of additional doctest tests of
qtcm modules that aren't called in the other test_*.py scripts.
Note that all tests occur with NumPy, but not all tests occur with
numarray and Numeric.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: test_all.py 6 2008-06-25 19:34:56Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, Physics Department, North 
#   Park University.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.4.
# - See import/reload statements throughout module for dependencies.
#
# Copyright (c) 2008 by Johnny Lin.  For licensing, distribution 
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


#- Other imports:

import __main__
import copy
import unittest




#---------------------------- Main Program -----------------------------
#
# Set variable run_verbose to True or False, depending on what I want
# to do (i.e. verbose testing or not verbose testing).

if __name__ == "__main__":
    #- Do all lengthy model runs:

    print 'Running model to obtain output (this may take a while) ...'
    execfile('run_short.py')
    execfile('run_long.py')
    print 'Running tests on output (this may take a while) ...'


    #- Load all subtests in this directory:

    import test_benchmark_short as bench_short
    import test_benchmark_long as bench_long
    import test_obj_defaults as obj
    import test_plotting as plot


    #- Set this to control verbosity:

    run_verbose = False
    suite = unittest.TestSuite()


    #- Add each array type test case if that array package exists:

    try:
        import Numeric
        class BenchShortNumeric(bench_short.NumericTests): pass
        class BenchLongNumeric(bench_long.NumericTests): pass
        class ObjNumeric(obj.NumericTests):  pass
        suite.addTest(unittest.makeSuite(BenchShortNumeric))
        suite.addTest(unittest.makeSuite(BenchLongNumeric))
        suite.addTest(unittest.makeSuite(ObjNumeric))
    except: pass

    try:
        import numarray
        class BenchShortNumarray(bench_short.NumarrayTests): pass
        class BenchLongNumarray(bench_long.NumarrayTests): pass
        class ObjNumarray(obj.NumarrayTests):  pass
        suite.addTest(unittest.makeSuite(BenchShortNumarray))
        suite.addTest(unittest.makeSuite(BenchLongNumarray))
        suite.addTest(unittest.makeSuite(ObjNumarray))
    except: pass

    try:
        import numpy
        class BenchShortNumPy(bench_short.NumPyTests): pass
        class BenchLongNumPy(bench_long.NumPyTests): pass
        class ObjNumPy(obj.NumPyTests):  pass
        class PlotNumPy(plot.NumPyTests):  pass
        suite.addTest(unittest.makeSuite(BenchShortNumPy))
        suite.addTest(unittest.makeSuite(BenchLongNumPy))
        suite.addTest(unittest.makeSuite(ObjNumPy))
        suite.addTest(unittest.makeSuite(PlotNumPy))

        import qtcm.where_close, doctest
        suite.addTest(doctest.DocTestSuite(qtcm.where_close))
    except: pass

    if run_verbose:
        unittest.TextTestRunner(verbosity=2).run(suite)
    else:
        unittest.TextTestRunner(verbosity=1).run(suite)




# ===== end file =====
