#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Create benchmark output files for testing.

This script assumes you are using a Unix system, that the executable
for Python is named "python", that the Python executable is on the
standard system path, and that this script is being called in the
test/benchmarks/create subdirectory of the qtcm package distribution.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: create_benchmarks.py 35 2008-07-14 22:27:39Z jlin $
#
# Modification History:
# - 19 Jun 2008:  Original by Johnny Lin, Physics Department, North
#   Park University.  Passed passably reasonable visual tests.
#
# Notes:
# - Written for Python 2.4.
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


#- Other imports:

import numpy as N
import copy
import shutil




#------------------------ Create the Benchmarks ------------------------

#- Check if we are in create and create variables describing some 
#  directories:

cwddir = os.getcwd()
if os.path.split(cwddir)[-1] != 'create':
    raise ValueError, "not in the create subdirectory"

workdir = os.path.join(cwddir, 'work')
procdir = os.path.join(cwddir, 'proc')


#- Go to work directory and run all shell scripts in the work 
#  directory (assumes scripts all have executable permission):

os.chdir(workdir)
for iname in os.listdir(workdir):
    print "Running", iname
    os.system(os.path.join(os.curdir, iname))


#- Go to proc, delete the executable of qtcm in each output directory.
#  Move the output directory one level up from create and change the
#  name of the output directory via copy and remove:

os.chdir(procdir)
for iname in os.listdir(procdir):
    olddirname = os.path.splitext(iname)[0]
    newdirname = olddirname.replace('benchmark_','')
    inameqtcm = os.path.join(procdir, olddirname, 'qtcm')
    if os.path.exists(inameqtcm):  os.remove(inameqtcm)
    shutil.copytree( os.path.join(procdir, olddirname),
                     os.path.join(procdir, os.pardir, os.pardir, newdirname) )
    shutil.rmtree( os.path.join(procdir, olddirname) )




# ====== end file ======
