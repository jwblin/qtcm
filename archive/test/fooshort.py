#- Imports:

import numpy as N
import time
import sys
import os
import user
import shutil
import copy

from qtcm import Qtcm
from utilities import prepare_outdir


#- Make test run directory if it doesn't exist.  Delete old proc
#  and log files as needed:

rundirname = 'test'
dirbasepath = prepare_outdir(rundirname)


#- Run:

inputs = {}
inputs['dt'] = 1200.
inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
                               , 'SST_Reynolds' )
inputs['outdir'] = dirbasepath
inputs['runname'] = rundirname
inputs['landon'] = 1
inputs['year0'] = 11
inputs['month0'] = 3
inputs['day0'] = 14
inputs['lastday'] = 15
inputs['ntout'] = 1
inputs['ntouti'] = 1
inputs['noout'] = 0
inputs['mrestart'] = 0
inputs['dateofmodel'] = 341023
inputs['init_with_instance_state'] = False
inputs['init_with_instance_state'] = True
inputs['compiled_form'] = 'full'
inputs['compiled_form'] = 'parts'

model = Qtcm(**inputs)
#a = model.get_qtcm_item('Qc')
#print '&&&', a
#model.sync_all_py_values_to_qtcm_items()
#execfile('/Users/jlin/.stop.py')  #@@@e
#b = model.get_qtcm_item('u1')
model.run_session()
#b = model.get_qtcm_item('Qc')
#print '@@@', b #[[  2.6744113    8.97750568
#model.run_session()
#print '###', b
#print '***', model.get_qtcm_item('Qc')  #[  3.13557339   9.62195015  
#c = model.get_qtcm_item('u1')
#print N.shape(a)
#print N.shape(b)




# ====== end file ======
