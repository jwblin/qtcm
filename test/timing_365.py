#-------------------------------------------------------------------------
# Execute and calculate timings for 365 aquaplanet runs for compiled_form
# full and parts.
#-------------------------------------------------------------------------


#- Imports:

import numpy as N
import copy
import sys
import os
import shutil
import time
import user

from qtcm import Qtcm
from utilities import prepare_outdir


#- Base inputs dictionary:

inputs = {}
inputs['dt'] = 1200.
inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
                               , 'SST_Reynolds' )
inputs['outdir'] = 'foo'
inputs['runname'] = 'foo'
inputs['landon'] = 0
inputs['year0'] = 1
inputs['month0'] = 11
inputs['day0'] = 1
inputs['lastday'] = 365
inputs['ntout'] = 1
inputs['ntouti'] = 1
inputs['noout'] = 0
inputs['mrestart'] = 0
inputs['init_with_instance_state'] = False


#- Timing runs for compiled_form "full":

average_times = []
list_forms = ['full', 'parts']
for iform in list_forms:
    inputs['compiled_form'] = iform
    all_times_form = []
    for i in xrange(3):
        rundirname = 'timing_365-' + inputs['compiled_form'] + str(i+1)
        dirbasepath = prepare_outdir(rundirname)

        inputs2 = copy.copy(inputs)
        inputs2['outdir'] = dirbasepath
        inputs2['runname'] = rundirname
        inputs2['init_with_instance_state'] = False
        model = Qtcm(**inputs2)
        time0 = time.time()
        model.run_session()
        all_times_form.append(time.time()-time0)

        if os.path.exists('qtcm_00021031.restart'):  
            os.remove('qtcm_00021031.restart')
        shutil.rmtree(dirbasepath)
        del model

    average_times.append( N.mean(N.array(all_times_form)) )

print 'Average wall clock time over 3 runs for', list_forms, average_times




# ====== end file ======
