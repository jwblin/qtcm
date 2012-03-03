#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Utilities to support creating run output for tests.

Functions include:
* prepare_outdir:  Prepare output directory for run.

* read_benchmark:  Read output in a specific run in the ./benchmarks 
  directory.

* read_output:  Read output data for unit tests in the ./rundir
  directory.  The initialization routine will initialize masked or 
  unmasked Numeric, numarray, or numpy data, depending on the input 
  flag.

* py_can_use:  Function to check if the current version of Python can
  be successfully used.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: utilities.py 37 2008-07-14 23:13:37Z jlin $
#
# Modification History:
# - 21 May 2008:  Original by Johnny Lin, Physics Department, North 
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

import copy
import Scientific.IO.NetCDF as S




#--------------------------- Module Variables --------------------------
#
# These module variables define default paths for prepare_outdir,
# read_benchmark, and read_data.

prepare_outdir_default_path = os.path.join(os.getcwd(), 'rundir')
read_benchmark_default_path = os.path.join(os.getcwd(), "benchmarks")
read_output_default_path = os.path.join(os.getcwd(), "rundir")

#prepare_outdir_default_path = os.path.join('/scra1', 'testing', 'rundir')
#read_benchmark_default_path = os.path.join("/scra1", "testing", "benchmarks")
#read_output_default_path = os.path.join("/scra1", "testing", "rundir")




#------------------- General Function prepare_outdir -------------------

def prepare_outdir(rundirname, path=prepare_outdir_default_path):
    """Prepare output directory for run.

    The output directory is assumed to be os.path.join(path,
    rundirname).  All directories in the path that do not exist are
    created.  If there are instantaneous and mean netCDF output
    files with the rundirname name in the output directory, those
    files are deleted.  A file named stdout.log in the current
    working directory, if it exists, is deleted.

    The default of path is set to the utilities module module
    variable prepare_outdir_default_path, which unless overloaded
    equals os.path.join(os.getcwd(), 'rundir').

    Input Parameter:
    * rundirname:  String describing the run.  This will be the
      name of the output directory in path.

    Output:
    * os.path.join(path, rundirname) is returned.

    Example:
    >>> a = prepare_outdir('test')
    """
    dirbasepath = os.path.join(path, rundirname)
    if not os.path.exists(dirbasepath):  os.makedirs(dirbasepath)

    qi_file = os.path.join( dirbasepath, 'qi_'+rundirname+'.nc' )
    qm_file = os.path.join( dirbasepath, 'qm_'+rundirname+'.nc' )
    log_file = os.path.join( os.getcwd(), 'stdout.log' )
    if os.path.exists(qi_file):   os.remove(qi_file)
    if os.path.exists(qm_file):   os.remove(qm_file)
    if os.path.exists(log_file):  os.remove(log_file)

    return dirbasepath




#------------------ General Function:  read_benchmark ------------------

def read_benchmark(field, suffix, array_type='numpy',
                   path=read_benchmark_default_path):
    """Read instantaneous and mean benchmark data.

    Returns a 2-element tuple of the instantaneous and mean benchmark
    data, as arrays of array_type, for files
    "[qi,qm]_benchmark_"+suffix+".nc" in directory path+"/"+suffix.

    The default of path is set to the utilities module module
    variable read_benchmark_default_path, which unless overloaded
    equals os.path.join(os.getcwd(), 'benchmarks').

    Positional Input:
    * field:  QTCM field name, as given in the netCDF file.  String.
    * suffix:  Descriptor.  String.

    Keyword Input:
    * array_type:  Type of array:  'Numeric', 'numarray', or 'numpy'.
    * path:  Path to the directory the suffix directories holding
      benchmark run output are in.

    Output:
    * 2-element tuple where the first element is the instantaneous
      data and the second element is the mean data.
    """
    if array_type=="Numeric":
        import Numeric as N
    elif array_type=="numarray":
        import numarray as N
    elif array_type=="numpy":
        import numpy as N
    else:
        raise ValueError, "No array package is defined"
        
    fn_inst = os.path.join( path, suffix, "qi_benchmark_" + suffix + ".nc" )
    fn_mean = os.path.join( path, suffix, "qm_benchmark_" + suffix + ".nc" )

    file_inst = S.NetCDFFile(fn_inst, mode='r')
    file_mean = S.NetCDFFile(fn_mean, mode='r')

    output = ( N.array(file_inst.variables[field].getValue())
             , N.array(file_mean.variables[field].getValue()) )

    file_inst.close()
    file_mean.close()

    return output




#-------------------- General Function:  read_output -------------------

def read_output(field, suffix, array_type='numpy',
                path=read_output_default_path):
    """Read instantaneous and mean output data.

    Returns a 2-element tuple of the instantaneous and mean output
    data, as arrays of array_type, for files
    "[qi,qm]_"+suffix+".nc" in directory path+"/"+suffix.

    The default of path is set to the utilities module module
    variable read_output_default_path, which unless overloaded
    equals os.path.join(os.getcwd(), 'rundir').

    Positional Input:
    * field:  QTCM field name, as given in the netCDF file.  String.
    * suffix:  Descriptor.  String.

    Keyword Input:
    * array_type:  Type of array:  'Numeric', 'numarray', or 'numpy'.
    * path:  Path to the directory the suffix directories holding
      run output are in.

    Output:
    * 2-element tuple where the first element is the instantaneous
      data and the second element is the mean data.
    """
    if array_type=="Numeric":
        import Numeric as N
    elif array_type=="numarray":
        import numarray as N
    elif array_type=="numpy":
        import numpy as N
    else:
        raise ValueError, "No array package is defined"
        
    fn_inst = os.path.join( path, suffix, "qi_" + suffix + ".nc" )
    fn_mean = os.path.join( path, suffix, "qm_" + suffix + ".nc" )

    file_inst = S.NetCDFFile(fn_inst, mode='r')
    file_mean = S.NetCDFFile(fn_mean, mode='r')

    output = ( N.array(file_inst.variables[field].getValue())
             , N.array(file_mean.variables[field].getValue()) )

    file_inst.close()
    file_mean.close()

    return output




#-------------------- General Function:  py_can_use --------------------

def py_can_use():
    """Returns True if current version of Python can be used.
    
    The basis of this test is whether or not the current Python
    interpreter can import numpy.  If so, assume that the rest of
    the needed packages exist in that Python and you can use that
    version of Python.
    """
    try:
        import numpy
        return True
    except:
        return False




#-------------------------- Main:  Test Module -------------------------
#
# Simple tests are done here using doctest.  More complex tests are
# done using unittest and are found in the seaice/test/semtner0 dir-
# ectory.

#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module documentation strings.

    Note:  To help ensure that module testing of this file works, the 
    parent directory to the current directory is added to sys.path.
    """
    import doctest
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
