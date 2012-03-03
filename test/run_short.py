#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Run all not create*_long.py scripts of the create*.py scripts.

This module is run using the following command line:

>>> python run_short.py

This currently works for Unix systems.  Think about how to do it for
Windows systems.  I use os.system in order to send stdout to a file.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: run_short.py 37 2008-07-14 23:13:37Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, North Park University.
#   Passed passably reasonable visual tests.
#
# Notes:
# - Written for Python 2.4.
# - See import/reload statements throughout module for dependencies.
#
# Copyright (c) 2008 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
#=======================================================================


#- Import modules:

import os
import sys
import utilities


#- Get list of all create_* files in the current working directory:

dirlist = []
for ifile in os.listdir(os.getcwd()):
    if ifile.startswith('create_'):
        dirlist.append(ifile)


#- Keep all files that do not end in _long.py:

proglist = []
for ifile in dirlist:
    if not ifile.endswith('_long.py'):
        proglist.append(ifile)
del dirlist


#- Run all proglist scripts:

for i in xrange(len(proglist)):
    print "Running", proglist[i]
    if i == 0:
        if utilities.py_can_use():
            os.system(sys.executable + ' ' + proglist[i] + ' > create.log')
        else:
            os.system('python ' + proglist[i] + ' > create.log')
    else:
        if utilities.py_can_use():
            os.system(sys.executable + ' ' + proglist[i] + ' >> create.log')
        else:
            os.system('python ' + proglist[i] + ' >> create.log')


#- Remove log file:

if os.path.exists('create.log'):  os.remove('create.log')


#====== end of file ======
