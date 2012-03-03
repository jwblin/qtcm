#-------------------------------------------------------------------------
# Create model output for multiple runs of multiple instances of Qtcm,
# for compiled_form = 'full'.
#
# The following output should be the same:
# * benchmark aquaplanet
# * full_365_aqua_multi1
# * full_365_aqua_multi2
# * full_365_aqua_multi4
#
# The following output should be the same:
# * benchmark aquaplanet_params2
# * full_365_aqua_multi3
#
# The following output should be the same:
# * benchmark landon
# * full_365_landon_multi1
# * full_365_landon_multi2
#-------------------------------------------------------------------------


#- Imports:

import numpy as N
import sys
import os
import user
import copy

from qtcm import Qtcm
from utilities import prepare_outdir


#- Aqua Run 1:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_aqua_multi1'
dirbasepath = prepare_outdir(rundirname)

inputs = {}
inputs['dt'] = 1200.
inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
                               , 'SST_Reynolds' )
inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname
inputs['landon'] = 0
inputs['year0'] = 1
inputs['month0'] = 11
inputs['day0'] = 1
inputs['lastday'] = 15
inputs['ntout'] = 1
inputs['ntouti'] = 1
inputs['noout'] = 0
inputs['mrestart'] = 0
inputs['compiled_form'] = 'parts'

model1 = Qtcm(**inputs)
model1.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Aqua Run 2:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_aqua_multi2'
dirbasepath = prepare_outdir(rundirname)

inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname

model2 = Qtcm(**inputs)
model2.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Aqua Run 3:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_aqua_multi3'
dirbasepath = prepare_outdir(rundirname)

inputs3 = copy.deepcopy(inputs)
inputs3['outdir'] = dirbasepath
inputs3['runname'] = rundirname
inputs3['ziml'] = 100.

_viscT = 12.0e4
_viscQ = 12.0e4
_viscU = 7.0e4
inputs3['viscxu0'] = _viscU                   
inputs3['viscyu0'] = _viscU                   
inputs3['visc4x'] = _viscU                    
inputs3['visc4y'] = _viscU                    
inputs3['viscxu1'] = _viscU                   
inputs3['viscyu1'] = _viscU                   
inputs3['viscxT'] = _viscT
inputs3['viscyT'] = _viscT
inputs3['viscxq'] = _viscQ
inputs3['viscyq'] = _viscQ
del _viscT, _viscQ, _viscU

model3 = Qtcm(**inputs3)
model3.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Aqua Run 4:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_aqua_multi4'
dirbasepath = prepare_outdir(rundirname)

inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname

model4 = Qtcm(**inputs)
model4.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Landon Run 1:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_landon_multi1'
dirbasepath = prepare_outdir(rundirname)

inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname
inputs['landon'] = 1

model5 = Qtcm(**inputs)
model5.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Landon Run 2:  Make test run directory if it doesn't exist.  Delete
#  old proc and log files as needed, set inputs dictionary, and run:

rundirname = 'full_365_landon_multi2'
dirbasepath = prepare_outdir(rundirname)

inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname
inputs['landon'] = 1

model6 = Qtcm(**inputs)
model6.run_session()
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')




# ====== end file ======
