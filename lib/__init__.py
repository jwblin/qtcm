#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""qtcm Package.

Some useful online help commands for the package:
* help(qtcm):  Help for the package.  A list of all modules in this
  package is found in the "Package Contents" section of the help
  output.
* help(qtcm.M):  Details of each module "M", where "M" is the module's 
  name.  
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# SVN Revision Code:
#   $Id: __init__.py 3 2008-06-25 01:02:46Z jlin $
#
# Modification History:
# - 30 May 2008:  Original by Johnny Lin, Physics Department, North
#   Park University.  Passed passably reasonable visual tests.
#
# Notes:
# - Written for Python 2.4.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python __init__.py".
# - See import statements throughout for non-"built-in" packages and
#   modules required.
#
# Copyright (c) 2007-2008 by Johnny Lin.  For licensing, distribution 
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


#- Import package version and set module version to package version.
#  Note that in qtcm.__init__.py if you put the statement "del 
#  package_version" then it will be impossible to use the statement
#  "import qtcm.package_version" or "from qtcm import package_ver-
#  sion":

import package_version
__version__ = package_version.version
__author__  = package_version.author
__date__    = package_version.date
__credits__ = package_version.credits


#- List of modules in package:

__all__ = [ "defaults",
            "field",
            "num_settings",
            "package_version",
            "plot",
            "qtcm",
            "where_close" ]




#-------------- Import To Make Available At Package Level --------------

from field import Field
from num_settings import N   #- numpy/Numeric/numarray as appropriate
from qtcm import Qtcm




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
