#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Module for class defining QTCM fields.

A field is a model parameters or variables, and includes diagnotic
and prognostic variables, run parameters, coefficients, etc.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# Package Name:
#   qtcm
#
# RCS Revision Code:
#   $Id: __init__.py,v 1.7 2004/04/29 22:47:46 jlin Exp $
#
# Modification History:
# - 12 Feb 2004:  Original by Johnny Lin, Computation Institute,
#   University of Chicago.  Passed passably reasonable tests.@@@
#
# Notes:
# - Written for Python 2.4.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python __init__.py".
# - See import statements throughout for non-"built-in" packages and
#   modules required.
#
# Copyright (c) 2007 by Johnny Lin.  For licensing, distribution 
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


#- Import numpy/Numeric/numarray as appropriate:

import num_settings as num
from num_settings import N


#- Other module imports:

from defaults import qtcm_fields as qtcm_defaults




#---------------------------- Class:  Field ----------------------------

class Field(object):
    """Class for QTCM fields.

    QTCM fields are model parameters or variables, and includes
    diagnotic and prognostic variables, run parameters, coefficients,
    etc.  They can be scalars (numeric or string) or arrays.  The
    values of these QTCM parameter objects can be changed in the
    model by a call at the Python level, though the types of their
    compiled model counterparts cannot be changed (without recompiling
    the compiled model, of course).

    The default values of QTCM fields that are defined in the
    defaults module specify the type of the variables, as well as
    the rank.  This information is used by the Qtcm class to properly
    interface with the compiled model.  Thus, fields that are not
    specified in the defaults module will not properly interface
    with the compiled model.  However, some fields only need to be
    defined at the Python level; those fields do not have to be
    listed in defaults.
    """
    def __init__(self, *args, **kwds):
        """Initialize Field object.

	The Field object is instantiated either with one or two
	positional input parameters, and up to two optional keyword
	input parameters.

        Positional Input Parameters:
	* One argument:  The argument is a string that specifies
	  the name of the field.  The name must match a key in
	  defaults.qtcm_fields.  The value of the Field instance
	  is set to the default value given in defaults.qtcm_fields.

	* Two arguments:  The first argument is a string specifying
	  the name of the field (as in the one argument case).  The
	  second argument is the value that the Field instance is
	  set to.  If the second argument is a Field object, the
	  value of that Field object is extracted as the value for
	  creating the current Field object.

        Keyword Input Parameters:
        * units:  String specifying the units of the field.

        * long_name:  String specifying the long name of the field.
        
        Examples:
        >>> a = Field('dt')
        >>> print a.id
        dt
        >>> print a.value
        1200.0

        >>> a = Field('dt', 1200.)
        >>> print a.id
        dt
        >>> print a.value
        1200.0
        >>> print a.units
        s
        >>> print a.long_name
        time step

        >>> a = Field('dt', Field('dt'), units='h')
        >>> print a.id
        dt
        >>> print a.value
        1200.0
        >>> print a.units
        h 

        >>> a = Field(23)
        Traceback (most recent call last):
            ...
        TypeError: Field id must be string
        """
        if type(args[0]) != type('a'):
            raise TypeError, self.__class__.__name__ + ' id must be string'

        if len(args) == 1:
            self.id = args[0]
            self.value = qtcm_defaults[self.id]['value']
        elif len(args) == 2:
            self.id = args[0]
            if type(args[1]) == type(self):
                if args[1].id != self.id:
                    raise ValueError, 'id mismatch'
                self.value = args[1].value
            else:
                self.value = args[1]
        else:
            raise ValueError, '1-2 required arguments for ' \
                            + self.__class__.__name__

        if kwds.has_key('units'):  
            self.units = kwds['units']
        else:
            self.units = qtcm_defaults[self.id]['units']

        if kwds.has_key('long_name'):
            self.long_name = kwds['long_name']
        else:
            self.long_name = qtcm_defaults[self.id]['long_name']


    def rank(self):
        """Return the rank of self.value.
        
        If self.value does not exist, returns None.
        """
        if hasattr(self, 'value'):
            return N.rank(self.value)
        else:
            return None


    def typecode(self):
	"""Return the typecode of self.value.

	The typecode is determined by first converting self.value
	into an array, and then returning the dtype.char (in numpy).
	This is defined in the module num_settings, and is a function
	of what type of array package you're using.  As a result,
	you shouldn't assume this method is very precise (e.g.,
	don't use it to distinguish between single and double
	precision float), but rather, use it to distinguish between
	different categories of types (e.g., float vs. int).  If
	self.value does not exist, returns None.
        """
        if hasattr(self, 'value'):
            return num.typecode(self.value)
        else:
            return None




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
