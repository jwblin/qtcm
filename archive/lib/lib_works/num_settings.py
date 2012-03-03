#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Choose and import using numpy, Numeric, or numarray modules.

   Sets the following module variables accordingly:

      N    = All numpy/Numeric/numarray functions, usage e.g. 
             N.arange(4).
      MA   = All numpy/Numeric/numarray masked array functions, usage e.g.
             MA.masked_array(data).
      MLab = All numpy/Numeric/numarray Matlab compatibility functions,
             usage e.g. MLab.squeeze(data).
      typecode = Function to emulate Numeric's typecode() method for
             Numeric objects.

   Thus, by importing using this package you will be able to code to 
   a generic naming convention, regardless of what array package is
   used.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: num_settings.py,v 1.3 2004/08/28 19:54:44 jlin Exp $
#
# Modification History:
# - 28 Aug 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.2.
# - See import statements throughout module for dependencies.
#
# Copyright (c) 2004-2007 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
#=======================================================================




#---------------- Module General Import and Declarations ---------------

#- If you're importing this module in testing mode, or you're running
#  pydoc on this module via the command line, import user-specific
#  settings to make sure any non-standard libraries are found:

import os, sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user
del os, sys


#- Import package version and set module version to package version:

import package_version as _package_version
__version__ = _package_version.version
__author__  = _package_version.author
__date__    = _package_version.date
__credits__ = _package_version.credits


#--------------------------- Module Variables --------------------------
#
# Import numpy/Numeric/numarray as module variables.  Set the module
# functions for each array type:  typecode.
# 

try:
    import numpy as N
    def typecode(arg):
        """Return typecode of arg if arg was a numpy array.

        Input argument:
        * arg:  Any argument that can be read converted by numpy into
          an array.
        """
        return N.array(arg).dtype.char

except ImportError:
    try:
        import numarray
        import numarray.ma.MA as MA
        import numarray.linear_algebra.mlab as MLab
        import numarray as N
        del numarray

    except ImportError:
        try:
            import Numeric as N
            import MA
            import MLab
            def typecode(arg):
                """Return typecode of arg if arg was a Numeric array.

                Input argument:
                * arg:  Any argument that can be read converted by Numeric
                  into an array.
                """
                return N.array(arg).typecode()

        except:
            raise ImportError, 'No array packages found'

#@@@in order to make this work for Numeric and numarray, I need to
#   define an attribute dtype.char of those array types that does
#   the functionality of typecode, since I write this package assuming
#   numpy and use numpy.dtype.char.  also need add MA and MLab for numpy.
#   i haven't programmed in a numarray typecode function yet either.



#--------------------------- Module Functions --------------------------
#
# These functions should work regardless of which type of array package
# is imported as N.





#-------------------------- Main:  Test Module -------------------------

__test__ = {'Additional Tests and Examples':
"""
"""}


if __name__ == "__main__":
    """Test the module documentation strings and __test__.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
