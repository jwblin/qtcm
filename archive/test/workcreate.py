#-------------------------------------------------------------------------
# Create aquaplanet model output for case of making a run done in multiple
# pieces, with and without using restart.  The total number of days of the
# aquaplanet_cont run is 40.
#
# The following output should be the same:
# * benchmark aquaplanet_cont10-30_10
# * full_365_aqua_cont10-30_10
# * parts_365_aqua_cont10-30_10
# * parts_365_aqua_nr_cont10-30_10
#
# The following output should be the same:
# * benchmark aquaplanet_cont10-30_30
# * full_365_aqua_cont10-30_30
# * parts_365_aqua_cont10-30_30
# * parts_365_aqua_nr_cont10-30_30
#-------------------------------------------------------------------------


#- Imports:

import numpy as N
import copy
import sys
import os
import user

from qtcm import Qtcm
from utilities import prepare_outdir


#- Initialize inputs dictionary:

inputs = {}
inputs['dt'] = 1200.
inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
                               , 'SST_Reynolds' )
inputs['landon'] = 0
inputs['ntout'] = 1
inputs['ntouti'] = 1
inputs['noout'] = 0


#- Case 4:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but same model instance object.  The both passes use the
#  instance state to initialize:

rundirname = 'parts_365_aqua_nr1_cont'
dirbasepath = prepare_outdir(rundirname)

inputs1 = copy.deepcopy(inputs)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = 1
inputs1['month0'] = 11
inputs1['day0'] = 1
inputs1['lastday'] = 10
inputs1['mrestart'] = 0
inputs1['init_with_instance_state'] = True
inputs1['compiled_form'] = 'parts'

model = Qtcm(**inputs1)
model.run_session()
snapshot = model.snapshot

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')

model.run_session(cont=30)
del model

#why is this incorrect?  why when you do a continuation run for another
#  run session do you not get the same results as the 30 day run
#  from case 3?  even if you don't write the same file.  the remove thing
#  didn't work.  what if i try with snapshot once more?  it doesn't work
#  because i hvaen't assigned a new qtcm compiled module?  put that
#  into run_session?  @@@add this note re compiled model to manual

if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')




# ====== end file ======
