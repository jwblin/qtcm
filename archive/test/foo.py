#- Imports:

import numpy as N
import time
import sys
import os
import user
import shutil
import copy

from qtcm import Qtcm
from utilities import read_benchmark
from utilities import read_output


#- To delete old proc and log files (@@@may delete):

qi_file = os.path.join( os.getcwd(), 'rundir', 'full_365_aqua' \
                      , 'qi_full_365_aqua.nc' )
qm_file = os.path.join( os.getcwd(), 'rundir', 'full_365_aqua' \
                      , 'qm_full_365_aqua.nc' )
log_file = os.path.join( os.getcwd(), 'stdout.log' )
if os.path.exists(qi_file):  os.remove(qi_file)
if os.path.exists(qm_file):  os.remove(qm_file)
if os.path.exists(log_file):  os.remove(log_file)


#- First Run:  aqua:

inputs = {}
inputs['dt'] = 1200.
inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
                               , 'SST_Reynolds' )
inputs['outdir'] = os.path.join( os.getcwd(), 'rundir', 'full_365_aqua' )
inputs['runname'] = 'full_365_aqua'
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

model = Qtcm(**inputs)
start_time_f2py = time.time()
for i in xrange(1):
    model.run_session()
print time.time() - start_time_f2py

if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')



#- Delete old proc and log files and move old to aqua2:

#if os.path.exists(os.path.join(os.getcwd(),'rundir','full_365_aqua2')):
#    shutil.rmtree(os.path.join(os.getcwd(),'rundir','full_365_aqua2'))
#shutil.copytree(os.path.join(os.getcwd(),'rundir','full_365_aqua'), \
#                os.path.join(os.getcwd(),'rundir','full_365_aqua2'))
#os.rename( \
#   os.path.join(os.getcwd(),'rundir','full_365_aqua2','qi_full_365_aqua.nc'),
#   os.path.join(os.getcwd(),'rundir','full_365_aqua2','qi_full_365_aqua2.nc'))
#os.rename( \
#   os.path.join(os.getcwd(),'rundir','full_365_aqua2','qm_full_365_aqua.nc'),
#   os.path.join(os.getcwd(),'rundir','full_365_aqua2','qm_full_365_aqua2.nc'))
#
#
#qi_file = os.path.join( os.getcwd(), 'rundir', 'full_365_aqua' \
#                      , 'qi_full_365_aqua.nc' )
#qm_file = os.path.join( os.getcwd(), 'rundir', 'full_365_aqua' \
#                      , 'qm_full_365_aqua.nc' )
#log_file = os.path.join( os.getcwd(), 'stdout.log' )
#if os.path.exists(qi_file):  os.remove(qi_file)
#if os.path.exists(qm_file):  os.remove(qm_file)
#if os.path.exists(log_file):  os.remove(log_file)


#- Second Run:  aqua:

#model2 = Qtcm(**inputs)
#model2.run_session()
#if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


#- Third Run:  landon:

#inputs['dt'] = 1200.
#inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
#inputs['bnddir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' )
#inputs['SSTdir'] = os.path.join( os.getcwd(), 'bnddir', 'r64x42' \
#                               , 'SST_Reynolds' )
#inputs['outdir'] = os.path.join( os.getcwd(), 'rundir', 'full_365_landon' )
#inputs['runname'] = 'full_365_landon'
#inputs['landon'] = 1
#inputs['year0'] = 1
#inputs['month0'] = 11
#inputs['day0'] = 1
#inputs['lastday'] = 13
#inputs['ntout'] = 1
#inputs['ntouti'] = 1
#inputs['noout'] = 0
#inputs['mrestart'] = 0

#print '@@@@@@@@@@@@@@@@', model.get_qtcm_item('outdir') #@@@
# the problem with this is for another run, even though you've changed
#  model, it doesn't use the new values for SSTdir etc.?

#model3 = Qtcm(**inputs)
#model3.run_session()
#if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')


ivar = 'u1'
data_aq = read_benchmark( ivar, 'aquaplanet', array_type='numpy' )
data_lo = read_benchmark( ivar, 'landon', array_type='numpy' )
data1 = read_output( ivar, 'full_365_aqua', array_type='numpy' )
#data2 = read_output( ivar, 'full_365_aqua2', array_type='numpy' )
#data3 = read_output( ivar, 'full_365_landon', array_type='numpy' )
print 'aqua run meets benchmark', N.allclose(data1[0], data_aq[0])
#print 'land run meets benchmark', N.allclose(data3[0], data_lo[0])
#print 'The two aqua runs should test allclose True: ', \
#      N.allclose(data1[0], data2[0])
#print 'The aqua runs and landon should both test allclose False: ', \
#      N.allclose(data1[0], data3[0]), N.allclose(data2[0], data3[0])
print model._Qtcm__qtcm.__file__
#print model2._Qtcm__qtcm.__file__
#print model3._Qtcm__qtcm.__file__

# ====== end file ======
