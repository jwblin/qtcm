#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Choose and import using numpy, Numeric, or numarray modules.

   Sets the following module variables accordingly:

      N    = All numpy/Numeric/numarray functions, usage e.g. 
             N.arange(4).
      MA   = All numpy/Numeric/numarray masked array functions, 
             usage e.g.  MA.masked_array(data).
      MLab = All numpy/Numeric/numarray Matlab compatibility 
             functions, usage e.g. MLab.squeeze(data).  For numpy,
             (it appears) many of these function are in the regular
             core library, so MLab for numpy is just numpy.
      isscalar = Function to emulate numpy's isscalar function.
      typecode = Function to emulate Numeric's typecode() method
             for Numeric objects.
      typecodes = Dictionary, with additions, that builds off of a
             copy of the N.typecodes dictionary.

   Thus, by importing using this package you will be able to code
   to a generic naming convention, regardless of what array package
   is used.  The module tries importing numpy first, then numarray,
   and finally Numeric.  The later packages are imported only if
   the earlier ones return an ImportError.  Currently, if you have
   numpy installed, MA and MLab will not be defined, since aliases
   for the numpy versions of those packages have not yet been
   programmed into this module.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: num_settings.py 35 2008-07-14 22:27:39Z jlin $
#
# Modification History:
# - 28 Aug 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed passably reasonable visual tests.
#
# Notes:
# - Written for Python 2.2.
# - See import statements throughout module for dependencies.
#
# Copyright (c) 2004-2008 by Johnny Lin.  For licensing, distribution 
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

try:
    import numpy as N
    try:
        import numpy.core.ma as MA
    except ImportError:
        import numpy.ma as MA
    import numpy as MLab
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
        def typecode(arg):
            """Return typecode of arg if arg was a numarray array.

            Input argument:
            * arg:  Any argument that can be read converted by Numeric
              into an array.
            """
            return N.array(arg).typecode()
        def isscalar(arg):
            raise ValueError, 'isscalar not yet implemented for numarray'

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
            def isscalar(arg):
                raise ValueError, 'isscalar not yet implemented for Numeric'

        except:
            raise ImportError, 'No array packages found'


#- Set typecodes module variable dictionary and make additions:

typecodes = N.typecodes
if 'S' not in typecodes['Character']:
    typecodes['Character'] = typecodes['Character'] + 'S'
if 'c' not in typecodes['Character']:
    typecodes['Character'] = typecodes['Character'] + 'c'




#-------------------------- Main:  Test Module -------------------------

if __name__ == "__main__":
    """Test the module documentation strings."""
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
