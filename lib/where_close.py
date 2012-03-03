#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Single public function module.

   See function docstring for description.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: where_close.py 16 2008-07-09 17:46:32Z jlin $
#
# Modification History:
# - 19 Mar 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed reasonable tests.
# - 11 Oct 2004:  Added capability to handle MA/ma arrays.  Passed
#   passably reasonable tests.
# - 31 May 2008:  Made a part of the qtcm package.
#
# Notes:
# - Written for Python 2.2.2.
# - Function is based on code from the MA module by Paul F. Dubois.
#   Some snippets of code in this function are copied directly from 
#   lines in that module.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python where_close.py".  More complete and complex
#   testing uses the unittest package, found in the test sub-direc-
#   tory at the package level.
# - See import statements throughout for packages/modules required.
#
# Copyright (c) 2004-2008 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
#=======================================================================




#----------------------- Overall Module Imports ------------------------

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


#- Import Numpy or numarray packages:

import num_settings as num
from num_settings import N
from num_settings import MA




#--------------------------- Public Function ---------------------------

def where_close(x, y, rtol=1.e-5, atol=1.e-8):
    """Mask of where x and y are element-wise "equal" to each other.

    Returns an int integer array with elements equal to 1 where x
    and y are "equal", and 0 otherwise.  If x or y are floating
    point, "equal" means where abs(x-y) <= atol + rtol * abs(y).
    This is essentially the same algorithm used in the
    numpy/Numeric/numarray function allclose.  If x and y are
    integer, "equal" means strict equality.  Shape and size of
    output is the same as x and y; if one is an array and the other
    is scalar, shape and size of the output is the same as the
    array.

    Output is an MA/ma masked array, unless both inputs are not
    MA/ma objects, in which case output is a numpy/Numeric/numarray
    array.  If inputs are both unmasked scalars the output is a
    Python integer scalar.

    Positional Input Arguments:
    * x:  Scalar, numpy/Numeric/numarray array, MA/ma array, Python
      list/tuple of any size and shape.  Floating or integer type.
    * y:  Scalar, numpy/Numeric/numarray array, MA/ma array, Python
      list/tuple of any size and shape.  Floating or integer type.

    Keyword Input Arguments:
    * rtol:   "Relative" tolerance.  Default is 1.e-5.  Used in the
      comparison between x and y only if the two are floating point.
    * atol:   "Absolute" tolerance.  Default is 1.e-8.  Used in the
      comparison between x and y only if the two are floating point.

    If either of the inputs are MA/ma masked objects, this function
    uses the MA/ma default algorithm for comparison, i.e., masked
    values are always considered equal.

    Examples:
    >>> from where_close import where_close
    >>> x = [20.,  -32., -1., 2.  , 5., 29.]
    >>> y = [20.1, -31., -1., 2.01, 3., 28.99]
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['0', '0', '1', '0', '0', '0']

    >>> from where_close import where_close
    >>> x = [20.,  -32., -1., 2.            , 5., 29.]
    >>> y = [20.1, -31., -1., 2.000000000001, 3., 28.99]
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['0', '0', '1', '1', '0', '0']

    >>> x = N.array([1,  5,  7, -2, 10])
    >>> y = N.array([1, -5, 17, -2,  0])
    >>> ind = where_close(x, y)
    >>> ['%.1g' % ind[i] for i in range(len(ind))]
    ['1', '0', '0', '1', '0']
    """
    if (not MA.isMA(x)) and (not MA.isMA(y)):
        return _where_close_unmasked(x, y, rtol=rtol, atol=atol)
    else:
        return _where_close_masked(x, y, rtol=rtol, atol=atol)




#------------------ Private Function:  For MA/ma Case ------------------

def _where_close_masked(x, y, rtol=1.e-5, atol=1.e-8):
    """Version of where_close for masked arrays.

    See docstring for where_close for details regarding parameters and
    function.
    """
    abs = MA.absolute


    #- Make sure input is MA/ma type:

    xM = MA.masked_array(x)
    yM = MA.masked_array(y)


    #- Safe compare if floating.  Strict compare if integer.  Any other
    #  type returns an error:

    if (xM.typecode() in MA.typecodes['Float']) or \
       (yM.typecode() in MA.typecodes['Float']):
        output_mask = MA.less_equal(abs(xM-yM), atol+rtol*abs(yM))

    elif (xM.typecode() in MA.typecodes['Integer']) and \
         (yM.typecode() in MA.typecodes['Integer']):
        output_mask = MA.equal(xM, yM)

    else:
        raise ValueError, "where_close:  Inputs must be Float or Integer"


    #- Return output_mask:

    if isinstance(output_mask, int):
        return output_mask
    else:
        return output_mask.astype(MA.Int)




#-------- Private Function:  For numarray/Numeric/numarray Case --------

def _where_close_unmasked(x, y, rtol=1.e-5, atol=1.e-8):
    """Version of where_close for unmasked arrays.

    See docstring for where_close for details regarding parameters and
    function.
    """
    abs = N.absolute


    #- Make sure input is numpy/Numeric/numarray type:

    xN = N.array(x)
    yN = N.array(y)


    #- Safe compare if floating.  Strict compare if integer.  Any other
    #  type returns an error:

    if (num.typecode(xN) in N.typecodes['Float']) or \
       (num.typecode(yN) in N.typecodes['Float']):
        output_mask = N.less_equal(abs(xN-yN), atol+rtol*abs(yN))

    elif (num.typecode(xN) in N.typecodes['Integer']) and \
         (num.typecode(yN) in N.typecodes['Integer']):
        output_mask = N.equal(xN, yN)

    else:
        raise ValueError, "where_close:  Inputs must be Float or Integer"


    #- Return output_mask:

    if isinstance(output_mask, int):
        return output_mask
    else:
        return output_mask.astype(int)




#---------------------- Addition doctest Examples ----------------------

#- Define additional examples for doctest to use:

__test__ = { 'Additional Examples':
    """
    """ }


#-------------------------- Main:  Test Module -------------------------

#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module.

    Note:  To help ensure that module testing of this file works, the
    parent directory to the current directory is added to sys.path.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
