#- Imports:

import numpy as N
import sys
import os
import user

from qtcm import Qtcm
from utilities import prepare_outdir


#- Make test run directory if it doesn't exist.  Delete old proc
#  and log files as needed:

rundirname = 'full_365_aqua'
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
inputs['compiled_form'] = 'full'

model = Qtcm(**inputs)
model.run_session()
#model.plotm('Qc',lat=[-40,40], time=4, lon=[0.,])
#model.plotm('Qc',lat=[-40,40], time=4)
model.plotm('Qc',lon=[20,300], time=4, approx_nlev=4, tmppreview=True)
#model.plotm('Qc',time=4, tmppreview=True)
#model.plotm('Qc',lon=[20,100], time=[0,4])  #@@@should give error
#model.plotm('Qc',lat=5.625, lon=[20,100], time=[0,4])
#model.plotm('Qc',lat=[-5.625, 50], lon=247.5, time=[0,4])
if os.path.exists('qtcm_00011115.restart'):  os.remove('qtcm_00011115.restart')




# ====== end file ======
