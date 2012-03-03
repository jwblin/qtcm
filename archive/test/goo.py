import os
import sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user
del os


#- Import other modules at module level:

import __main__
import copy
import unittest
import numpy as N
from utilities import read_benchmark
from utilities import read_output




def test_all_fields_benchmarks_are_all_different():
    """Test all benchmark outputs are different for sub-var_list.

    Only select variables are compared.
    """
    #@@@e
    sub_var_list = [ 'WD', 'GMs1', 'u1', 'cdn', 'S0', 'Ts', 'Prec' \
                    , 'top', 'advT1', 'dfsq1', 'GMq1', 'div0', 'q1' \
                    , 'Evap', 'lon', 'u0', 'OLR', 'div1', 'wet' \
                    , 'advq1', 'FLWds', 'FLWus', 'Runs', 'dfsT1', 'T1' \
                    , 'v0', 'v1', 'vort0', 'Runf', 'lat', 'FTs', 'Evapi' \
                    , 'FSWut', 'FSWus', 'us', 'vs', 'psi0', 'stype' \
                    , 'time', 'taux', 'tauy', 'cl1', 'FSWds' ]

    #sub_var_list = [ 'us', 'Ts']# 'GMs1', 'u1', 'cdn', 'S0', 'Ts', 'Prec' \
    #                , 'top', 'advT1', 'dfsq1', 'GMq1', 'div0', 'q1' \
    #                , 'Evap', 'lon', 'u0', 'OLR', 'div1', 'wet' \
    #                , 'advq1', 'FLWds', 'FLWus', 'Runs', 'dfsT1', 'T1' \
    #                , 'v0', 'v1', 'vort0', 'Runf', 'lat', 'FTs', 'Evapi' \
    #                , 'FSWut', 'FSWus', 'us', 'vs', 'psi0', 'stype' \
    #                , 'time', 'taux', 'tauy', 'cl1', 'FSWds' ]

    sub_var_list = ['u1',]
    for ivar in sub_var_list:
        data1 = read_benchmark( ivar, 'aquaplanet' \
                              , array_type='numpy' )
        data2 = read_benchmark( ivar, 'landon' \
                              , array_type='numpy' )
        data3 = read_benchmark( ivar, 'aquaplanet_params2' \
                              , array_type='numpy' )
        myoutput3 = read_output( ivar, 'full_365_aqua_multi3' \
                               , array_type='numpy' )

        #print N.allclose(data1[0], data2[0])
        #print N.allclose(data1[1], data2[1])
        #print ivar, N.allclose(data3[0], myoutput3[0])
        #print ivar, N.allclose(data3[1], myoutput3[1])
        print N.max( N.abs((data3[0]- myoutput3[0])/data3[0]) )
        print N.max( N.abs((data3[1]- myoutput3[0])/data3[1]) )
        #print N.allclose(data2[0], data3[0])
        #print N.allclose(data2[1], data3[1])
        #execfile('/Users/jlin/.stop.py')




if __name__ == "__main__":
    test_all_fields_benchmarks_are_all_different()



# ===== end file =====
