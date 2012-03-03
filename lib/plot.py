#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Functions to plot fields from Qtcm model objects and run.

Module Function:
* nice_levels:  Compute a vector of "levels" at "nice" increments.

* plot_ncdf_output:  Plot data from QTCM1 netCDF output files.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: plot.py 4 2008-06-25 01:03:28Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, Physics Department, North
#   Park University.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.4.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python plot.py".
# - See import statements throughout for non-"built-in" packages and
#   modules required.
#
# Copyright (c) 2008 by Johnny Lin.  For licensing, distribution 
# conditions, contact information, and additional documentation see
# the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
#=======================================================================




#---------------- Module General Import and Declarations ---------------

#- Import os and sys.  Also, if you're importing this module in
#  testing mode, or you're running pydoc on this module via the 
#  command line, import user-specific settings to make sure any 
#  non-standard libraries are found:

import os, sys
if (__name__ == "__main__") or \
   ("pydoc" in os.path.basename(sys.argv[0])):
    import user


#- Import package version and set module version to package version:

import package_version as _package_version
__version__ = _package_version.version
__author__  = _package_version.author
__date__    = _package_version.date
__credits__ = _package_version.credits


#- Import numpy/Numeric/numarray as appropriate:

import num_settings as num
from num_settings import N


#- Import matplotlib (must be in this order):

from matplotlib import rc as matplotlibrc
#matplotlibrc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
matplotlibrc('text', usetex=True)
import pylab
import matplotlib
from matplotlib.toolkits.basemap import Basemap


#- Other imports that then are used as module variables:

import copy
import Scientific.IO.NetCDF as S
import shutil
import tempfile
from where_close import where_close




#-------------------- Module Function:  nice_levels --------------------

def nice_levels(data, approx_nlev=10, max_nlev=28):
    """Compute a vector of "levels" at "nice" increments.

    Returns a 1-D array of "levels" (e.g., contour levels) calculated
    to give an aesthetically pleasing and human-readable interval,
    if possible.  If not, returns levels for approx_nlev levels
    between the maximum and minimum of data.  In any event, the
    function will return no more than max_nlev levels.

    Keyword Input Parameter:
    * data:  Array of values to calculate levels for.  Can be of any 
      size and shape.

    Keyword Input Parameter:
    * approx_nlev:  Integer referring to approximately how many
      levels to return.  This is the way of adjusting how "coarse"
      or "fine" to make the vector of levels.

    * max_nlev:  The maximum number of levels the function will
      permit to be returned.  The interval of levels will be adjusted
      to keep the number of levels returned under this value.  If
      approx_nlev is chosen to be greater than or equal to max_nlev,
      an exception is raised.

    Output:
    * This function returns a 1-D array of contour levels.

    Function is adaptation of parts of IDL routine contour_plot.pro
    by Johnny Lin.  This is why the capitalization conventions of
    Python are not strictly followed in this function.

    Examples:
    >>> z = N.array([-24.5, 50.3, 183.1, 20.])
    >>> out = nice_levels(z)
    >>> ['%g' % out[i] for i in range(len(out))]
    ['-30', '0', '30', '60', '90', '120', '150', '180', '210']

    >>> z = N.array([-24.5, 50.3, 183.1, 20.])
    >>> out = nice_levels(z, approx_nlev=5)
    >>> ['%g' % out[i] for i in range(len(out))]
    ['-50', '0', '50', '100', '150', '200']

    >>> z = N.array([-24.5, 50.3, 183.1, 20.])
    >>> out = nice_levels(z, approx_nlev=10)
    >>> ['%g' % out[i] for i in range(len(out))]
    ['-30', '0', '30', '60', '90', '120', '150', '180', '210']
    """
    #- Default settings and error check:

    if approx_nlev >= max_nlev:
        raise ValueError, 'max_nlev is too small'

    MAX_zd_ok = N.max(data)
    MIN_zd_ok = N.min(data)

    nlevels = N.min([approx_nlev, max_nlev])
    tmpcmax = MAX_zd_ok
    tmpcmin = MIN_zd_ok
    tmpcint = N.abs( (tmpcmax-tmpcmin)/float(nlevels) )


    #- See if the cint can be "even".  If not, return alternative
    #  contour levels vector:

    #+ Guess a possible cint.  Look for an "even" value that is
    #  closest to that:

    guesscint = N.abs( (MAX_zd_ok-MIN_zd_ok)/float(nlevels) )

    if (guesscint > 1e-10) and (guesscint < 1e+10):
        possiblecint = [1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5,
                            0.0001,  0.0002,            0.0005,
                            0.001,   0.002,             0.005,
                            0.01,    0.02,              0.05,
                            0.1,     0.2,               0.5,
                            1.,      2.,                5.,
                           10.,     20.,     30., 45., 50.,
                          100.,    200.,              500.,
                         1000.,   2000.,             5000.,
                        10000.,  20000.,            50000.,
                        1e+5, 1e+6, 1e+7, 1e+8, 1e+9, 1e+10]

        diffcint = N.abs(possiblecint-guesscint)
        tempcint = N.compress( diffcint == N.min(diffcint), possiblecint )[0]
        tcidx = N.compress( where_close( possiblecint, tempcint ),
                            N.arange(N.size(possiblecint)) )[0]


        #+ Look around at the "even" values nearby the possible option
        #  for cint.  Calculate how many contours each of those cints
	#  would give.  Dictionary ncon_count is the number of
	#  contours for a given test cint.  test_tcidxs are the indices
	#  in possiblecint to examine in detail; these index values
	#  will be the keys in ncon_count.

        if tcidx == 0:  tcidx = 1
        if tcidx == N.size(possiblecint)-1:  tcidx = N.size(possiblecint)-2

        ncon_count = {}
        test_tcidxs = [tcidx-1, tcidx, tcidx+1]
        for i in test_tcidxs:
            itcval = possiblecint[i]
            positivecon = N.arange(max_nlev+2, dtype=float)*itcval
            negativecon = (N.arange(max_nlev+2, dtype=float)+1.0)*itcval*(-1.0)
            if (MAX_zd_ok + itcval >= 0) and (MIN_zd_ok - itcval >= 0):
                ncon_count[i] = N.sum( \
                    N.logical_and(positivecon <= MAX_zd_ok + itcval,
                                  positivecon >= MIN_zd_ok - itcval) )
            elif (MAX_zd_ok + itcval < 0) and (MIN_zd_ok - itcval < 0):
                ncon_count[i] = N.sum( \
                    N.logical_and(negativecon <= MAX_zd_ok + itcval,
                                  negativecon >= MIN_zd_ok - itcval) )
            else:
                ncon_count[i] = N.sum(positivecon <= MAX_zd_ok + itcval) \
                              + N.sum(negativecon >= MIN_zd_ok - itcval)


	#+ Select the cint that has the fewest levels if it has at
	#  least nlevels-1.  Otherwise, try to find the next cint with
	#  the fewest levels that is below max_nlev.  tempcint is what
	#  you get (changed, if warranted) leaving this section:

        min_ncon_count = N.min(ncon_count.values())
        current_best_count = max_nlev
        for i in test_tcidxs:
            if (ncon_count[i] == min_ncon_count) and \
               (ncon_count[i] >= nlevels-1):
                tempcint = possiblecint[i]
                current_best_count = ncon_count[i]
                break
            elif (ncon_count[i] == min_ncon_count) and \
               (ncon_count[i] < nlevels-1):
                continue
            elif ncon_count[i] > max_nlev:
                continue
            else:
                if N.abs(ncon_count[i]-nlevels) < \
                   N.abs(current_best_count-nlevels):
                    tempcint = possiblecint[i]
                    current_best_count = ncon_count[i]
                continue


	#+ Create levels for case with neg. and pos. contours.  There
	#  is the case of all pos., all neg., and mixed pos. and neg.
	#  contours:

        positivecon = N.arange(max_nlev+2, dtype=float)*tempcint
        negativecon = (N.arange(max_nlev+2, dtype=float)+1.0)*tempcint*(-1.0)

        if (MAX_zd_ok + tempcint >= 0) and (MIN_zd_ok - tempcint >= 0):
            tmpclevels = N.compress( \
                N.logical_and(positivecon <= MAX_zd_ok + tempcint, 
                              positivecon >= MIN_zd_ok - tempcint),
                                     positivecon )
        elif (MAX_zd_ok + tempcint < 0) and (MIN_zd_ok - tempcint < 0):
            tmpclevels = N.compress( \
                N.logical_and(negativecon <= MAX_zd_ok + tempcint, 
                              negativecon >= MIN_zd_ok - tempcint),
                                     negativecon )
        else:
            uppercon = N.compress( positivecon <= MAX_zd_ok + tempcint, 
                                   positivecon )
            lowercon = N.compress( negativecon >= MIN_zd_ok - tempcint, 
                                   negativecon )
            tmpclevels = N.concatenate([lowercon, uppercon])


	#+ Sort clevels, reset number of levels, maximum, minimum,
	#  and interval of contour based on the automatic setting:

        tmpclevels = N.sort(tmpclevels)
        if (N.size(tmpclevels) <= max_nlev ) and (N.size(tmpclevels) > 0):
            nlevels = N.size(tmpclevels)
            tmpcmax = tmpclevels[-1]
            tmpcmin = tmpclevels[0]
            tmpcint = tempcint

    else:
        pass


    #- Return output:

    return N.arange(tmpcmin, tmpcmax+tmpcint, tmpcint)      




#----------------- Module Function:  mpl_latex_script1 -----------------

def mpl_latex_script1(instring):
    """Return instring to handle super/subscripts of one character.

    The string returned expresses instring so that an exponent or
    subscript with a single character after it (e.g., ^2, _s) can
    be processed as a LaTeX exponent or superscript in Matplotlib.
    Any number of these super/subscripts can be present in instring.
    See http://www.scipy.org/Cookbook/Matplotlib/UsingTex by Darrin
    Dale for some code snippets incorporated in this function.

    Positional Input Parameter:
    * instring:  String to be formatted.

    Output:
    * A string, processable by Matplotlib's LaTeX emulation module.

    Examples:
    >>> mpl_latex_script1('Precipitation [W/m^2]')
    'Precipitation [W/m$^2$]'
    >>> mpl_latex_script1('u_0^2 and u_1^2 [m^2/s^2]')
    'u$_0$$^2$ and u$_1$$^2$ [m$^2$/s$^2$]'
    """
    #- Make input string be a list of individual characters:

    instrlist = list(instring)


    #- Cycle through each ^ and _, change accordingly, and return:

    for ichar in ['^', '_']:
        num_char = instrlist.count(ichar)
        for inum in xrange(num_char):
            cidx = instrlist.index(ichar)
            instrlist[cidx] = '$' + ichar
            instrlist[cidx+1] = instrlist[cidx+1] +'$'

    return ''.join(instrlist)




#----------------- Module Function:  plot_ncdf_output ------------------

def plot_ncdf_output(id, datafn, **kwds):
    """Plot model field id from the data in netCDF file datafn.

    Positional Input Parameter:
    * id:  Name of the id of the field to plot.  String.

    * datafn:  Filename containing the output data to plot.  String.

    Input keyword parameter descriptions are found in the docstring
    for Qtcm methods ploti, plotm, and other methods that call this
    private method.  In general, those methods take the keyword
    parameters they receive and pass it along unchanged as keyword
    parameters to this function.  In that sense, this function is
    seldom used as a stand-alone function, but rather is usually
    used coupled with a Qtcm instance.

    The data fields read in from the netCDF output file are dimensioned
    (time, lat, lon).  This is different than how the data is stored
    in the compiled QTCM model fields (lon, lat, time), and at the
    Python level (lon, lat).  The reason this is the case is that
    f2py automatically makes the arrays passed between the Python
    and Fortran levels match.

    For a lat vs. lon plot, the contour plot is superimposed onto
    a cylindrical projection map of the Earth with continents drawn
    and labeled meridians and parallels.  The title also includes
    the model time, and x- and y-axis labels are not drawn.

    All numerical data used for plotting come from the netCDF output
    file for consistency (e.g., the dimensions of u1).  Currently
    this method only works for 3-D data arrays (two in space, one
    in time).
    """
    #- Accomodate other ids.  The id that this routine will use
    #  (i.e., iduse corresponds to the name in the netCDF output file)
    #  is called iduse.  Set this to id, except for the case where
    #  some aliases of ids are entered in which iduse is the alias:

    if id == 'Qc': iduse = 'Prec'
    elif id == 'FLWut': iduse = 'OLR'
    elif id == 'STYPE': iduse = 'stype'
    else: iduse = id


    #- Set defined keyword defaults.  All are set to None except for
    #  nlatlon which gets an integer:

    plotkwds_ids = ['lat', 'lon', 'time', 'fn', 'levels', 'title',
                    'xlabel', 'ylabel', 
                    'filled', 'nlatlon', 'tmppreview']

    plotkwds = {}
    for ikey in plotkwds_ids:
        if kwds.has_key(ikey):
            plotkwds[ikey] = copy.copy(kwds[ikey])
        else:
            plotkwds[ikey] = None

    if not kwds.has_key('nlatlon'):
        plotkwds['nlatlon'] = 8


    #- Get data and dimensions of iduse to plot:

    fileobj = S.NetCDFFile(datafn, mode='r')
    data = N.array(fileobj.variables[iduse].getValue())
    data_name = fileobj.variables[iduse].long_name
    data_units = fileobj.variables[iduse].units

    dim = {}
    dimname = {}
    dimunits = {}

    dim['lat'] = N.array(fileobj.variables['lat'].getValue())
    dimname['lat'] = fileobj.variables['lat'].long_name
    dimunits['lat'] = fileobj.variables['lat'].units

    dim['lon'] = N.array(fileobj.variables['lon'].getValue())
    dimname['lon'] = fileobj.variables['lon'].long_name
    dimunits['lon'] = fileobj.variables['lon'].units
    
    dim['time'] = N.array(fileobj.variables['time'].getValue())
    dimname['time'] = fileobj.variables['time'].long_name
    dimunits['time'] = fileobj.variables['time'].units

    fileobj.close()


    #- Alter data long name to remove any units.  The definition
    #  of units as the substring within the [] is the same as used in
    #  defVar in output.F90 of the compiled QTCM model.  Remove 
    #  underscores and extra whitespace in data_name and data_units, 
    #  replacing with a single whitespace character between words:

    idx1 = data_name.find('[')
    idx2 = data_name.find(']')
    if idx1 != -1 and idx2 != -1:
        data_name = data_name[:idx1] + data_name[idx2+1:]
    data_name = data_name.strip()

    data_name  = ' '.join(data_name.replace('_',' ').split())
    data_units = ' '.join(data_units.replace('_',' ').split())


    #- Alter dimension long name to remove any units.  The definition
    #  of units as the substring within the [] is the same as used in
    #  defVar in output.F90 of the compiled QTCM model.  Remove 
    #  underscores and extra whitespace in name and units, replacing 
    #  with a single whitespace character between words, and 
    #  capitalizing like a title:

    for idimkey in dim.keys():
        idimname = dimname[idimkey]
        idx1 = idimname.find('[')
        idx2 = idimname.find(']')
        if idx1 != -1 and idx2 != -1:
            idimname = idimname[:idx1] + idimname[idx2+1:]
        dimname[idimkey] = idimname.strip()

        dimname[idimkey]  = \
            ' '.join(dimname[idimkey].replace('_',' ').split()).title()
        dimunits[idimkey] = \
            ' '.join(dimunits[idimkey].replace('_',' ').split()).title()


    #- Some data checks:

    if N.rank(data) != 3:
        raise ValueError, '_plot: can only plot lat, lon, time fields'
    if not N.allclose(dim['time'], N.sort(dim['time'])):
        raise ValueError, '_plot: time not monotonically ascending'
    if not N.allclose(dim['lat'], N.sort(dim['lat'])):
        raise ValueError, '_plot: lat not monotonically ascending'
    if not N.allclose(dim['lon'], N.sort(dim['lon'])):
        raise ValueError, '_plot: lon not monotonically ascending'
    if N.shape(data)[0] != N.size(dim['time']):
        raise ValueError, '_plot: data time dim mismatch'
    if N.shape(data)[1] != N.size(dim['lat']):
        raise ValueError, '_plot: data lat dim mismatch'
    if N.shape(data)[2] != N.size(dim['lon']):
        raise ValueError, '_plot: data lon dim mismatch'


    #- Choose and describe ranges for lat, lon, and time.  The
    #  section cycles through the dictionary of dimensions.  idim is
    #  the 1-D array of the values of that dimension.  rngs is a
    #  dictionary where each entry corresponds to a dimension, and the
    #  value of the entry is the values of that dimension that are to 
    #  be plotted.  rngs_idxs are the indices in the original 
    #  dimensions array corresponding to the values in rngs.  
    #  keys_rngs_sizes_gt_1 is a list of the keys of ranges that have 
    #  sizes greater than 1:

    rngs = {}
    rngs_idxs = {}
    keys_rngs_sizes_gt_1 = []
    for idimkey in dim.keys():
        idim = dim[idimkey]

        if plotkwds[idimkey] == None:
            dim_mask = N.ones( N.size(idim), dtype=int )

        elif N.isscalar(plotkwds[idimkey]):
            dim_mask = where_close( idim, plotkwds[idimkey] )
            if N.sum(dim_mask) != 1:
                raise ValueError, 'no point chosen'

        elif (not N.isscalar(plotkwds[idimkey])) and \
             N.size(plotkwds[idimkey]) == 1:
            dim_mask = where_close( idim, plotkwds[idimkey][0] )
            if N.sum(dim_mask) != 1:
                raise ValueError, 'no point chosen'

        elif N.size(plotkwds[idimkey]) == 2:
            dim_mask = N.logical_and( idim >= plotkwds[idimkey][0],
                                      idim <= plotkwds[idimkey][-1] )

        else:
            raise ValueError, 'bad dimension range keyword entry'

        rngs[idimkey]      = N.compress( dim_mask, idim )
        rngs_idxs[idimkey] = N.compress( dim_mask, N.arange(N.size(idim)) )
        if N.size(rngs[idimkey]) > 1:
            keys_rngs_sizes_gt_1.append(idimkey)


    #- Set plot types (line or contour):

    if len(keys_rngs_sizes_gt_1) == 0:
        raise ValueError, 'cannot plot without any fixed dimension'
    elif len(keys_rngs_sizes_gt_1) == 1:
        plottype = 'line'
    elif len(keys_rngs_sizes_gt_1) == 2:
        plottype = 'contour'
    else:
        raise ValueError, 'cannot plot with > 2 varying dimensions'


    #- Set plot axis fields and axis names, depending on what sort
    #  of dimensions will be plotted.  If lon is to be plotted, it is
    #  always the x-axis.  If lat is to be plotted, it is always the
    #  y-axis.  In this section and later on in a few places, I
    #  rely on the count method for a list as a Boolean test:  If it
    #  returns 0, consider that False; > 0 is True.  The text for the
    #  title and axis labels to be passed to the xlabel, etc. methods,
    #  are called titlename, xname, and yname:

    if plottype == 'line':                       #+ Choose x-axis vector and
        x = rngs[keys_rngs_sizes_gt_1[0]]        #  x/y names for line plot
        xname = dimname[keys_rngs_sizes_gt_1[0]] + ' [' \
              + dimunits[keys_rngs_sizes_gt_1[0]] + ']'
        yname = data_name + ' [' + data_units + ']'

    elif plottype == 'contour':                  #+ Choose axis vectors and
        if keys_rngs_sizes_gt_1.count('lon'):    #  names for contour plot
            x = rngs['lon']
            xname = dimname['lon'] + ' [' + dimunits['lon'] + ']'
        if keys_rngs_sizes_gt_1.count('lat'):
            y = rngs['lat']
            yname = dimname['lat'] + ' [' + dimunits['lat'] + ']'
        if keys_rngs_sizes_gt_1.count('time'):
            if keys_rngs_sizes_gt_1.count('lon'):
                y = rngs['time']
                yname = dimname['time'] + ' [' + dimunits['time'] + ']'
            elif keys_rngs_sizes_gt_1.count('lat'):
                x = rngs['time']
                xname = dimname['time'] + ' [' + dimunits['time'] + ']'
            else:
                raise ValueError, 'bad treatment of time'

    else:
        raise ValueError, 'unrecognized plottype'


    #- Override xname, yname, and titlename with keywords, if they
    #  are not None.  titlename receives data_name and data_units
    #  by default:

    if plotkwds['xlabel'] != None:
        xname = plotkwds['xlabel']
    if plotkwds['ylabel'] != None:
        yname = plotkwds['ylabel']

    if plotkwds['title'] != None:
        titlename = plotkwds['title']
    else:
        titlename = data_name + ' [' + data_units + ']'


    #- Pick data to be plotted and plot:

    pylab.clf()                       #+ Clear any previous figures
    pylab.figure(1)                   #+ Open a pylab figure

    if plottype == 'line':            #+ Select data for a line plot
        y = data[rngs_idxs['time'],   #  and plot
                 rngs_idxs['lat'], 
                 rngs_idxs['lon']]
        pylab.plot(x, y)

    elif plottype == 'contour':       #+ Select data for a contour
        ritim = rngs_idxs['time']     #  plot and plot
        rilat = rngs_idxs['lat']
        rilon = rngs_idxs['lon']


        #* Extract subarrays depending on which two dimensions are 
        #  chosen:

        if N.size(rngs_idxs['time']) == 1:
            zgrid = num.MLab.squeeze(data[ ritim[0],
                                           rilat[0]:rilat[-1]+1, 
                                           rilon[0]:rilon[-1]+1 ])
        elif N.size(rngs_idxs['lat']) == 1:
            zgrid = num.MLab.squeeze(data[ ritim[0]:ritim[-1]+1,
                                           rilat[0],
                                           rilon[0]:rilon[-1]+1 ])
        elif N.size(rngs_idxs['lon']) == 1:
            zgrid = num.MLab.squeeze(data[ ritim[0]:ritim[-1]+1,
                                           rilat[0]:rilat[-1]+1, 
                                           rilon[0] ])
        else:
            raise ValueError, 'unrecognized configuration'


        #* Change zgrid for special case of a lat. vs. time contour 
        #  plot.  Calculate xgrid and ygrid:

        if keys_rngs_sizes_gt_1.count('time') and \
           keys_rngs_sizes_gt_1.count('lat'):
           zgrid = N.transpose(zgrid)

        xgrid, ygrid = pylab.meshgrid(x, y)
        

        #* Set contour levels:

        if plotkwds['levels'] == None:
            levels = nice_levels(zgrid)
        else:
            levels = plotkwds['levels']


        #- Plot (creating continents first if is a lat vs. lon plot)
        #  and write contour levels/color bar as appropriate:

        if keys_rngs_sizes_gt_1.count('lon') and \
           keys_rngs_sizes_gt_1.count('lat'):
            mapplot = Basemap(projection='cyl', resolution='l',
                              llcrnrlon=N.min(xgrid), llcrnrlat=N.min(ygrid),
                              urcrnrlon=N.max(xgrid), urcrnrlat=N.max(ygrid))
            mapplot.drawcoastlines()
            mapplot.drawmeridians(nice_levels(rngs['lon'], 
                                  approx_nlev=plotkwds['nlatlon']),
                                  labels=[1,0,0,1])
            mapplot.drawparallels(nice_levels(rngs['lat'],
                                  approx_nlev=plotkwds['nlatlon']),
                                  labels=[1,0,0,1])
            if plotkwds['filled']:
                plot = mapplot.contourf(xgrid, ygrid, zgrid, levels)
                pylab.colorbar(plot, orientation='horizontal', format='%g')
            else:
                plot = mapplot.contour(xgrid, ygrid, zgrid, levels)
                pylab.clabel(plot, inline=1, fontsize=10, fmt='%g')
        else:
            if plotkwds['filled']:
                plot = pylab.contourf(xgrid, ygrid, zgrid, levels)
                pylab.colorbar(plot, orientation='horizontal', format='%g')
            else:
                plot = pylab.contour(xgrid, ygrid, zgrid, levels)
                pylab.clabel(plot, inline=1, fontsize=10, fmt='%g')

    else:
        raise ValueError, 'unrecognized plottype'


    #- Add titling.  Lat vs. lon plots do not have axis labels because
    #  the map labels already make it clear, and for those plots the
    #  title also includes the time value:

    if keys_rngs_sizes_gt_1.count('lon') and \
       keys_rngs_sizes_gt_1.count('lat'):
        titlename = titlename + ' at ' \
                  + dimname['time'] + ' ' \
                  + str(rngs['time'][0]) + ' ' \
                  + dimunits['time']
        titlename = mpl_latex_script1(titlename)
        pylab.title(titlename)
    else: 
        titlename = mpl_latex_script1(titlename)
        xname = mpl_latex_script1(xname)
        yname = mpl_latex_script1(yname)
        pylab.xlabel(xname)
        pylab.ylabel(yname)
        pylab.title(titlename)


    #- Output plot to PNG file or screen.  The show command seems to
    #  have a problem on my Mac OS X, so save to a temporary file
    #  and use preview to view for fn == None and tmppreview set to
    #  True.  Note that the temporary file is not deleted by this 
    #  method:

    if plotkwds['fn'] == None:                       #+ Screen display
        if plotkwds['tmppreview'] and sys.platform == 'darwin':
            outputfn = tempfile.mkstemp('.png','qtcm_')
            pylab.savefig(outputfn[-1])
            os.system('open -a /Applications/Preview.app '+outputfn[-1])
        else:
            pylab.show()

    elif type(plotkwds['fn']) == type('a'):          #+ Write to file
        pylab.savefig(plotkwds['fn'])
        pylab.close(1)

    else:
        raise ValueError, 'cannot write to this type of file'




#--------- Addition doctest Examples as Private Module Variable --------

__test__ = { 'All positive values data':
    """
>>> z = N.array([[24.5, 50.3],[283.1, 20.]])
>>> out = nice_levels(z)
>>> ['%g' % out[i] for i in range(len(out))]
['0', '30', '60', '90', '120', '150', '180', '210', '240', '270', '300']
    """,

'Exception for max_nlev':
    """
>>> z = N.array([-24.5, 50.3, 183.1, 20.])
>>> out = nice_levels(z, approx_nlev=45)
Traceback (most recent call last):
    ...
ValueError: max_nlev is too small
    """,

'More mpl_latex_script1 examples':
    """
>>> mpl_latex_script1('u_1 variance [m^2/s^2]')
'u$_1$ variance [m$^2$/s$^2$]'
    """,

'Additional Examples of nice_levels':
    """
>>> z = N.array([-24.5, 50.3, 183.1, 20.])
>>> out = nice_levels(z, approx_nlev=45, max_nlev=46)
>>> ['%g' % out[i] for i in range(len(out))]
['-25', '-20', '-15', '-10', '-5', '0', '5', '10', '15', '20', '25', '30', '35', '40', '45', '50', '55', '60', '65', '70', '75', '80', '85', '90', '95', '100', '105', '110', '115', '120', '125', '130', '135', '140', '145', '150', '155', '160', '165', '170', '175', '180', '185']

>>> z = N.array([124.5, 150.3, 183.1, 220.])
>>> out = nice_levels(z, approx_nlev=20)
>>> ['%g' % out[i] for i in range(len(out))]
['120', '125', '130', '135', '140', '145']
        """ }




#-------------------------- Main:  Test Module -------------------------

#- Execute doctest if module is run from command line:

if __name__ == "__main__":
    """Test the module.

    Note:  To help ensure that module testing of this file works, the 
    parent directory to the current directory is added to sys.path.
    """
    import doctest, sys, os
    sys.path.append(os.pardir)
    doctest.testmod(sys.modules[__name__])




# ===== end file =====
