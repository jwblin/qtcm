#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Defaults for QTCM fields in the QTCM package.

QTCM fields are model parameters or variables, and includes diagnotic
and prognostic variables, run parameters, coefficients, etc.  They
can be scalars (numeric or string) or arrays.  The values of these
QTCM parameter objects can be changed in the model by a call at the
Python level, though the types of their compiled model counterparts
cannot be changed (without recompiling the compiled model, of
course).

The default values of QTCM fields that are defined in the defaults
module specify the type of the variables, as well as the rank.  This
information is used by the Qtcm class to properly interface with
the compiled model.  Thus, fields that are not specified in the
defaults module will not properly interface with the compiled model.

Module Parameters Defined:
* init_prognostic_dict:  Dictionary giving the default initial
  values of each prognostic variable and right-hand side (as defined
  by the restart file specification).
* qtcm_fields_ids:  List of all the ids in qtcm_fields.
* qtcm_fields:  Dictionary of parameters.

Examples:
>>> print qtcm_fields['dt']['long_name']
time step
>>> print qtcm_fields['dt']['units']
s
>>> print int(qtcm_fields['dt']['value'])
1200
>>> print int(init_prognostic_dict['q1'])
-50
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: defaults.py 4 2008-06-25 01:03:28Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, Physics Department, North
#   Park University.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.4.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python defaults.py".
# - See import statements throughout for non-"built-in" packages and
#   modules required.
#
# Copyright (c) 2008 by Johnny Lin.  For licensing, distribution 
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


#- Import package version and set module version to package version:

import package_version as _package_version
__version__ = _package_version.version
__author__  = _package_version.author
__date__    = _package_version.date
__credits__ = _package_version.credits


#- Import numpy/Numeric/numarray as appropriate and other modules:

from num_settings import N
import copy




#-------------------------- Module Variables ---------------------------

#- Very private variables:
#
#  __viscT, __viscQ, and __viscU are used *only in this module* to set 
#  viscxu0, viscyu0, visc4x, visc4y, viscxu1, viscyu1, viscxT, viscyT, 
#  viscxq, and viscyq.  The settings will not have any effect outside
#  this module because all three variables are set to literals.

__viscT = 12.0e5    #- temperature diffusion parameter [m^2/s]
__viscQ = 12.0e5    #- humidity diffusion parameter [m^2/s]
__viscU = 7.0e5     #- viscocity parameter [m^2/s]


#- init_prognostic_dict:  Dictionary giving the default initial
#  values of each prognostic variable and right-hand side (as defined
#  by the restart file specification).

init_prognostic_dict = {'u1':0.0, 'v1':0.0, 'T1':-100.0, 'q1':-50.0,
                        'u0':0.0, 'v0':0.0, 'vort0':0.0, 'rhsvort0':0.0,
                        'u0bar':0.0, 'rhsu0bar':0.0, 'Ts':295.0, 'psi0':0.0,
                        'WD':0.0}


#- qtcm_fields:  Dictionary of all QTCM parameters that can be
#  changed in the model by a call at the Python level:

qtcm_fields = {}        #- Initialize dictionary


#+ Scalars:

qtcm_fields['title'] = { 'value' : 'QTCM default title' \
                       , 'units' : '' \
                       , 'long_name' : 'a descriptive title' }

qtcm_fields['bnddir'] = { 'value' : os.path.join(os.pardir, 'bnddata') \
                        , 'units' : '' \
                        , 'long_name' : 'boundary data other than SST' }

qtcm_fields['SSTdir'] = { 'value' : \
                          os.path.join(os.pardir, 'bnddata', 'SST_Reynolds') \
                        , 'units' : '' \
                        , 'long_name' : 'where SST files are' }

qtcm_fields['outdir'] = { 'value' : \
                          os.path.join(os.pardir, 'proc', 'qtcm_output') \
                        , 'units' : '' \
                        , 'long_name' : 'where output goes to' }

qtcm_fields['runname'] = { 'value' : 'runname' \
                         , 'units' : '' \
                         , 'long_name' : 'string for an output filename' }

qtcm_fields['landon'] = { 'value' : 1 \
                        , 'units' : '' \
                        , 'long_name' : \
                          'if not 1: land = ocean with fake SST' }

qtcm_fields['SSTmode'] = { 'value' : 'seasonal' \
                         , 'units' : '' \
                         , 'long_name' : 'decide what kind of SST to use' }

qtcm_fields['year0']  = { 'value' : 0 \
                        , 'units' : 'yr' \
                        , 'long_name' : \
                          'starting year; if < 0 use year in restart' }

qtcm_fields['month0']  = { 'value' : -1 \
                         , 'units' : 'mo' \
                         , 'long_name' : \
                           'starting month; if < 0 use mo in restart' }

qtcm_fields['day0']  = { 'value' : -1 \
                       , 'units' : 'dy' \
                       , 'long_name' : \
                         'starting day; if < 0 use day in restart' }

qtcm_fields['lastday']  = { 'value' : 365 \
                          , 'units' : 'dy' \
                          , 'long_name' : 'last day of integration' }

qtcm_fields['dateofmodel']  = { 'value' : 0 \
                              , 'units' : '' \
                              , 'long_name' : \
    'date of model coded as an integer as yyyymmdd' }

qtcm_fields['it']  = { 'value' : 1 \
                     , 'units' : '' \
                     , 'long_name' : 'time of day in time steps' }

qtcm_fields['interval']  = { 'value' : 1 \
                           , 'units' : 'dy' \
                           , 'long_name' : \
                             'atmosphere-ocean coupling interval' }

qtcm_fields['noout'] = { 'value' : 0 \
                       , 'units' : 'dy' \
                       , 'long_name' : 'no output for the first noout days' }

qtcm_fields['nooutr'] = { 'value' : 0 \
                        , 'units' : 'dy' \
                        , 'long_name' : \
                          'no restart file for the first nooutr days' }

qtcm_fields['ntout'] = { 'value' : -30 \
                       , 'units' : 'dy' \
                       , 'long_name' : 'monthly mean output' }

qtcm_fields['ntouti'] = { 'value' : 0 \
                        , 'units' : 'dy' \
                        , 'long_name' : 'monthly instantaneous data output' }

qtcm_fields['ntoutr'] = { 'value' : 0 \
                        , 'units' : 'dy' \
                        , 'long_name' : \
                          'restart file only at end of model run' }

qtcm_fields['nastep'] = { 'value' : 1 \
                        , 'units' : '' \
                        , 'long_name' : \
                          'number of atmosphere time steps within one ' \
                          +'air-sea coupling interval' }

qtcm_fields['mrestart'] = { 'value' : 0 \
                          , 'units' : '' \
                          , 'long_name' : '=1: restart using qtcm.restart' }

qtcm_fields['dt'] = { 'value' : 1200. \
                    , 'long_name' : 'time step' \
                    , 'units' : 's' }

qtcm_fields['eps_c'] = { 'value' : 0.13888889E-03 \
                       , 'long_name' : '1/tau_c NZ (5.7)' \
                       , 'units' : '1/s' }

qtcm_fields['mt0'] = { 'value' : 1 \
                     , 'units' : '' \
                     , 'long_name' : 'barotropic timestep every mt0 timesteps' }

qtcm_fields['ziml'] = { 'value' : 500. \
                      , 'units' : 'm' \
                      , 'long_name' : \
                        'atmosphere mixed layer depth ~ cloud base'}

qtcm_fields['weml'] = { 'value' : 0.01 \
                      , 'units' : 'm/s' \
                      , 'long_name' : 'mixed layer entrainment velocity'}

qtcm_fields['VVsmin'] = { 'value' : 4.5 \
                        , 'units' : 'm/s' \
                        , 'long_name' : 'minimum wind speed for fluxes' }

qtcm_fields['viscxu0'] = { 'value' : __viscU \
                         , 'units' : 'm^2/s' \
                         , 'long_name' : 'viscocity parameter for u0 in x' }

qtcm_fields['viscyu0'] = { 'value' : __viscU \
                         , 'units' : 'm^2/s' \
                         , 'long_name' : 'viscocity parameter for u0 in y' }

qtcm_fields['visc4x'] = { 'value' : __viscU \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : 'del 4 viscocity parameter in x' }

qtcm_fields['visc4y'] = { 'value' : __viscU \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : 'del 4 viscocity parameter in y' }

qtcm_fields['viscxu1'] = { 'value' : __viscU \
                         , 'units' : 'm^2/s' \
                         , 'long_name' : 'viscocity parameter for u1 in x' }

qtcm_fields['viscyu1'] = { 'value' : __viscU \
                         , 'units' : 'm^2/s' \
                         , 'long_name' : 'viscocity parameter for u1 in y' }

qtcm_fields['viscxT'] = { 'value' : __viscT \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : \
                          'temperature diffusion parameter in x' }

qtcm_fields['viscyT'] = { 'value' : __viscT \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : \
                          'temperature diffusion parameter in y' }

qtcm_fields['viscxq'] = { 'value' : __viscQ \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : 'humidity diffusion parameter in x' }

qtcm_fields['viscyq'] = { 'value' : __viscQ \
                        , 'units' : 'm^2/s' \
                        , 'long_name' : 'humidity diffusion parameter in y' }

qtcm_fields['u0bar'] = { 'value' : 0.0 \
                     , 'units' : '' \
                     , 'long_name' : '' }


#+ Arrays:

qtcm_fields['Qc'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : 'K' \
                    , 'long_name' : 'precipitation' }

qtcm_fields['STYPE'] = { 'value' : N.zeros([1,1], dtype=float) \
                     , 'units' : '' \
                     , 'long_name' : \
                       'surface type; ocean or vegetation type over land' }

qtcm_fields['FLWds'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : '' }

qtcm_fields['FLWus'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : '' }

qtcm_fields['FSWds'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : '' }

qtcm_fields['FSWus'] = { 'value' : N.zeros([1,1], dtype=float) \
                     , 'units' : '' \
                     , 'long_name' : '' }

qtcm_fields['Evap'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : '' }

qtcm_fields['FTs'] = { 'value' : N.zeros([1,1], dtype=float) \
                     , 'units' : '' \
                     , 'long_name' : '' }

qtcm_fields['taux'] = { 'value' : N.zeros([1,1], dtype=float) \
                      , 'units' : '' \
                      , 'long_name' : '' }

qtcm_fields['tauy'] = { 'value' : N.zeros([1,1], dtype=float) \
                      , 'units' : '' \
                      , 'long_name' : '' }

qtcm_fields['FLWut'] = { 'value' : N.zeros([1,1], dtype=float) \
                       , 'units' : '' \
                       , 'long_name' : '' }

qtcm_fields['FLW'] = { 'value' : N.zeros([1,1], dtype=float) \
                     , 'units' : '' \
                     , 'long_name' : '' }

qtcm_fields['S0'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : '' }

qtcm_fields['FSWut'] = { 'value' : N.zeros([1,1], dtype=float) \
                     , 'units' : '' \
                     , 'long_name' : '' }

qtcm_fields['FSW'] = { 'value' : N.zeros([1,1], dtype=float) \
                   , 'units' : '' \
                   , 'long_name' : '' }

qtcm_fields['u1'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'm/s' \
                  , 'long_name' : 'current time step baroclinic zonal wind' }

qtcm_fields['v1'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'm/s' \
                  , 'long_name' : '' }

qtcm_fields['T1'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'K' \
                  , 'long_name' : '' }

qtcm_fields['q1'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'K' \
                  , 'long_name' : '' }

qtcm_fields['u0'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'm/s' \
                  , 'long_name' : 'barotropic zonal wind' }

qtcm_fields['v0'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'm/s' \
                  , 'long_name' : 'barotropic meridional wind' }

qtcm_fields['vort0'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : '' \
                  , 'long_name' : '' }

qtcm_fields['psi0'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : '' \
                  , 'long_name' : '' }

qtcm_fields['rhsvort0'] = { 'value' : N.zeros([1,1,1], dtype=float) \
                        , 'units' : '' \
                        , 'long_name' : '' }

qtcm_fields['rhsu0bar'] = { 'value' : N.zeros([1,], dtype=float) \
                        , 'units' : '' \
                        , 'long_name' : '' }

qtcm_fields['Ts'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : 'K' \
                  , 'long_name' : 'surface temperature' }

qtcm_fields['WD'] = { 'value' : N.zeros([1,1], dtype=float) \
                  , 'units' : '' \
                  , 'long_name' : '' }

qtcm_fields['WD0'] = { 'value' : N.zeros([1,], dtype=float) \
                  , 'units' : '' \
                  , 'long_name' : 'field capacity SIB2/CSU (approximately)' }

qtcm_fields['arr1'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 1' }

qtcm_fields['arr2'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 2' }

qtcm_fields['arr3'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 3' }

qtcm_fields['arr4'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 4' }

qtcm_fields['arr5'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 5' }

qtcm_fields['arr6'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 6' }

qtcm_fields['arr7'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 7' }

qtcm_fields['arr8'] = { 'value' : N.zeros([1,1], dtype=float) \
                    , 'units' : '' \
                    , 'long_name' : \
                      'auxiliary optional output array 8' }


#- qtcm_fields_ids:  List of all the ids in qtcm_fields:

qtcm_fields_ids = qtcm_fields.keys()




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
