#-------------------------------------------------------------------------
# Quick view of QTCM netCDF output data.
#-------------------------------------------------------------------------


import copy
import pylab as p
from matplotlib.toolkits.basemap import Basemap
import os
import Scientific.IO.NetCDF as S


#--- Get data (assuming lon goes from 0 to 360):  Make arrays cyclic
#    in longitude:

time_idx = 5   #@@@ USER ADJUSTABLE

f = S.NetCDFFile('qm_seasonal_1yr.nc', mode='r')
lat = f.variables['lat'].getValue()
lon = f.variables['lon'].getValue()
u1_all = f.variables['u1'].getValue()
if p.allclose(lon[0], 0):
    u1 = p.zeros( (p.size(lat), p.size(lon)+1) )
    u1[:,0:-1] = u1_all[time_idx,:,:]
    u1[:,-1] = u1_all[time_idx,:,0]

tmp = copy.copy(lon)
lon = p.zeros((p.size(tmp)+1,))
lon[0:-1] = tmp[:]
lon[-1] = tmp[0]+360
del tmp

f.close()


#--- Mapping information:

map = Basemap( projection='cyl', resolution='l'
             , llcrnrlon=0, urcrnrlon=360
             , llcrnrlat=-76.875, urcrnrlat=76.875
             , lon_0=180, lat_0=0
             )
map.drawmeridians(p.arange(0,361,45), labels=[0,0,0,1])
map.drawparallels(p.arange(-90,90,30), labels=[1,0,0,1])
map.drawcoastlines()


#--- Write out contour map and view using preview:

x, y = p.meshgrid(lon,lat)
CS = map.contourf(x, y, u1, cmap=p.cm.gray)
p.text( 0.5, -0.15, 'Longitude [deg]'
      , horizontalalignment='center'
      , verticalalignment='center'
      , transform = p.gca().transAxes )
p.text( -0.11, 0.5, 'Latitude [deg]'
      , rotation='vertical'
      , horizontalalignment='center'
      , verticalalignment='center'
      , transform = p.gca().transAxes )
p.title('QTCM u1')
cbar = p.colorbar(CS, orientation='horizontal', shrink=0.6)
cbar.ax.set_xlabel('See QTCM Documentation')
p.savefig('foo.png')
os.system('preview foo.png')
