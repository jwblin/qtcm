#-------------------------------------------------------------------------
# Create aquaplanet model output for case of making a run done in multiple
# pieces, with and without using restart.  The total number of days of the
# aquaplanet_cont40 run is 40.
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
#
# The following output should be the same:
# * benchmark aquaplanet_cont40
# * parts_365_aqua_nr1_cont40
#
# The following output should be the same:
# * benchmark aquaplanet_cont40u1
# * parts_365_aqua_nr2_cont40
#
# The following output should be different:
# * benchmark aquaplanet_cont40u1
# * full_365_aqua_nr2_cont40
# * parts_365_aqua_nr3_cont40
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


#- Case 1:  Regular aquaplanet run using compiled_form 'full', run
#  in two passes (10 and 30 days), using restart:

rundirname = 'full_365_aqua_cont10-30_10'
dirbasepath = prepare_outdir(rundirname)

inputs1 = copy.deepcopy(inputs)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = 1
inputs1['month0'] = 11
inputs1['day0'] = 1
inputs1['lastday'] = 10
inputs1['mrestart'] = 0
inputs1['compiled_form'] = 'full'

model = Qtcm(**inputs1)
model.run_session()
del model

os.rename( os.path.join(os.getcwd(), 'qtcm_00011110.restart'),
           os.path.join(os.getcwd(), 'qtcm.restart') )

rundirname = 'full_365_aqua_cont10-30_30'
dirbasepath = prepare_outdir(rundirname)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = -1
inputs1['month0'] = -1
inputs1['day0'] = -1
inputs1['lastday'] = 30
inputs1['mrestart'] = 1

model = Qtcm(**inputs1)
model.run_session()
del model

if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')
if os.path.exists('qtcm.restart'):  os.remove('qtcm.restart')


#- Case 2:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), using restart:

rundirname = 'parts_365_aqua_cont10-30_10'
dirbasepath = prepare_outdir(rundirname)

inputs1 = copy.deepcopy(inputs)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = 1
inputs1['month0'] = 11
inputs1['day0'] = 1
inputs1['lastday'] = 10
inputs1['mrestart'] = 0
inputs1['init_with_instance_state'] = False
inputs1['compiled_form'] = 'parts'

model = Qtcm(**inputs1)
model.run_session()
del model

os.rename( os.path.join(os.getcwd(), 'qtcm_00011110.restart'),
           os.path.join(os.getcwd(), 'qtcm.restart') )

rundirname = 'parts_365_aqua_cont10-30_30'
dirbasepath = prepare_outdir(rundirname)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = -1
inputs1['month0'] = -1
inputs1['day0'] = -1
inputs1['lastday'] = 30
inputs1['mrestart'] = 1

model = Qtcm(**inputs1)
model.run_session()
del model

if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')
if os.path.exists('qtcm.restart'):  os.remove('qtcm.restart')


#- Case 3:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but different model instances.  The first pass, do not
#  use instance state to initialize; in the second pass, use the
#  instance state to initialize:

rundirname = 'parts_365_aqua_nr_cont10-30_10'
dirbasepath = prepare_outdir(rundirname)

inputs1 = copy.deepcopy(inputs)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['year0'] = 1
inputs1['month0'] = 11
inputs1['day0'] = 1
inputs1['lastday'] = 10
inputs1['mrestart'] = 0
inputs1['init_with_instance_state'] = False
inputs1['compiled_form'] = 'parts'

model = Qtcm(**inputs1)
model.run_session()
snapshot = model.snapshot
del model

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')

rundirname = 'parts_365_aqua_nr_cont10-30_30'
dirbasepath = prepare_outdir(rundirname)

inputs1 = copy.deepcopy(inputs)
inputs1['outdir'] = dirbasepath
inputs1['runname'] = rundirname
inputs1['lastday'] = 30
inputs1['init_with_instance_state'] = True
inputs1['compiled_form'] = 'parts'

model = Qtcm(**inputs1)
model.sync_set_py_values_to_snapshot(snapshot=snapshot)
model.run_session()
del model

if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')


#- Case 4:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but same model instance object.  Both passes use the instance 
#  state to initialize:

rundirname = 'parts_365_aqua_nr1_cont40'
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
model.run_session(cont=30)
del model

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')
if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')


#- Case 5:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but same model instance object.  Both passes use the instance 
#  state to initialize.  After the first pass, set u1 to twice the 
#  value in certain areas:

rundirname = 'parts_365_aqua_nr2_cont40'
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
model.u1.value[52:55,13:16] = model.u1.value[52:55,13:16] * 2.0
model.run_session(cont=30)
del model

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')
if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')


#- Case 6:  Regular aquaplanet run using compiled_form 'parts', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but same model instance object.  Both passes use the instance 
#  state to initialize.  After the first pass, set u1 to twice the 
#  value in select areas.  Before the second pass, 
#  init_with_instance_state is turned off to see what happens when a 
#  False init_with_instance_state meets a model._cont is True attribute 
#  (this should return output different than parts_365_aqua_nr2_cont40 
#  because varinit decides on whether or not to use the Fortran version 
#  of varinit on the basis of the value of init_with_instance_state 
#  alone):

rundirname = 'parts_365_aqua_nr3_cont40'
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
model.u1.value[52:55,13:16] = model.u1.value[52:55,13:16] * 2.0
model.init_with_instance_state = False
model.run_session(cont=30)
del model

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')
if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')


#- Case 7:  Regular aquaplanet run using compiled_form 'full', run
#  in two passes (10 and 30 days), not using restart, in two output
#  files but same model instance object.  Both passes use the instance 
#  state to initialize.  After the first pass, set u1 to twice the 
#  value for select areas.  Check if cont keyword works for 
#  compiled_form 'full'; it shouldn't since this version doesn't call 
#  the varinit method in the Python level, calling the Fortran init 
#  routines which do not reference the Python field states (or even 
#  existing Fortran field states):

rundirname = 'full_365_aqua_nr2_cont40'
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
inputs1['compiled_form'] = 'full'

model = Qtcm(**inputs1)
model.run_session()
model.u1.value[52:55,13:16] = model.u1.value[52:55,13:16] * 2.0
model.run_session(cont=30)
del model

if os.path.exists('qtcm_00011110.restart'):  os.remove('qtcm_00011110.restart')
if os.path.exists('qtcm_00011210.restart'):  os.remove('qtcm_00011210.restart')


#- Clean-up possibly remaining restarts:

for ifn in ['qtcm_00010105.restart', 'qtcm_00011130.restart',
            'qtcm_00011230.restart', 'qtcm_00020430.restart']:
    if os.path.exists(ifn):  os.remove(ifn)




# ====== end file ======
