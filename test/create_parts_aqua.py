#-------------------------------------------------------------------------
# Create model output for case of using the Fortran varinit and the Python
# varinit.
#
# The following output should be the same:
# * benchmark aquaplanet
# * parts_365_aqua_
# * parts_365_aqua_inst1
#-------------------------------------------------------------------------


#- Imports:

import numpy as N
import copy
import sys
import os
import user

from qtcm import Qtcm
from utilities import prepare_outdir


#- Case 1:  Regular aquaplanet run using compiled_form 'parts', and
#  instantiating from scratch using the Fortran varinit:

rundirname = 'parts_365_aqua'
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
inputs['init_with_instance_state'] = False
inputs['compiled_form'] = 'parts'

model = Qtcm(**inputs)
model.run_session()
del model

if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Case 2:  Regular aquaplanet run using compiled_form 'parts', and
#  instantiating from scratch using the Python varinit:

rundirname = 'parts_365_aqua_inst1'
dirbasepath = prepare_outdir(rundirname)

inputs2 = copy.copy(inputs)
inputs2['outdir'] = dirbasepath
inputs2['runname'] = rundirname
inputs2['dateofmodel'] = 0
inputs2['init_with_instance_state'] = True

model = Qtcm(**inputs2)
model.run_session()
del model

if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')




# ====== end file ======
