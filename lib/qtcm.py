#!/usr/bin/python -tt
#=======================================================================
#                        General Documentation

"""Definition of the Qtcm class.
"""

#-----------------------------------------------------------------------
#                       Additional Documentation
#
# RCS Revision Code:
#   $Id: qtcm.py 57 2008-09-12 17:27:04Z jlin $
#
# Modification History:
# - 29 May 2008:  Original by Johnny Lin, Physics Department, North
#   Park University.  Passed passably reasonable tests.
#
# Notes:
# - Written for Python 2.4.
# - Module docstrings can be tested using the doctest module.  To
#   test, execute "python qtcm.py".
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


#- Other imports that then are used as module variables:

import copy
import defaults
from field import Field
from plot import plot_ncdf_output
import Scientific.IO.NetCDF as S
import shutil
import tempfile




#--------------------------- Module Variables --------------------------

#- Field object instance used usually in type tests:

_test_field = Field('dt')


#- Dictionary of all prognostic variables and right-hand sides that
#  could be initialized and their default initial values.  List of
#  those variable keys, plus dateofmodel and title, which are all the
#  variables usually written out into a restart file:

_init_prog_dict = defaults.init_prognostic_dict
_init_vars_keys = _init_prog_dict.keys() + ['dateofmodel', 'title']




#----------------------- Custom Exception Classes ----------------------
#

class FieldNotReadableFromCompiledModel(UnboundLocalError):
    "Field variable unable to be read from QTCM compiled model."




#----------------------------- Class:  Qtcm ----------------------------
#
# Note:  With exception of the __init__ object, all the public methods
# come first followed by the private methods.  Within each category,
# the methods are listed in alphabetical order.

class Qtcm(object):
    """Qtcm model object.

    Public Instance Attributes:
    * compiled_form:  What form the compiled Fortran version of the
      QTCM model has.  This attribute is a string and can have the 
      following values:

      - 'full':  The compiled portion of the model encompasses the
	entire QTCM model.  Thus, the only compiled QTCM model
	modules or subroutines that Python should interact with is
	__qtcm.driver (which executes the entire model) and
	__qtcm.setbypy (which enables communication between the
	compiled model and the Python-level of model fields.

      - 'parts':  The compiled portion of the model encompasses
	parts of the model as separate units all the way down to
	an atmosphere timestep.  Thus, compiled QTCM model
	module/subroutines include those that are executed within
	an atmosphere timestep.

      This attribute must be set on instantiation via a keyword
      input parameter.

    * coupling_day:  Current value of atmosphere-ocean coupling
      day.  When compiled_form of type 'parts' is run, is started
      at 1 and then is updated every coupling day.  coupling_day
      not changed from instantiation value for compiled_form 'full'.
      Is a Field object with integer scalar value, but is set to
      None on instantiation (unless its value is overridden through
      keyword input parameters).  Is set to a Field object to enable
      run lists that specify it to be able to be updated.

    * init_with_instance_state:  This keyword only has an effect if 
      compiled_form = 'parts'.  Boolean scalar.

      If True, uses, in general, the values of the instance attributes
      (for prognostic variables and right-hand sides, and start
      date) are right before the run_session method call as the
      initial values for the run_session.  If False, run_session
      initialization, in general, follows the rules set by all the
      standard QTCM model input parameters (e.g., the sign of day0,
      month0, and year0, the value of the mrestart flag, etc.).
      See the docstring for the varinit method for details, as well
      as the manual.  Default is True, unless overridden, on
      instantiation.

    * QTCM fields:  This is not the name of an attribute, but a
      category of attributes.  These attributes are QTCM fields
      whose names correspond to those given in the _qtcm_fields_ids
      attribute, and whose values are Field objects.  Initially set
      on instantiation.

      Override of default values on instantiation is accomplished
      by setting the value to a keyword input parameter with the
      Field object id as the keyword key.  The __init__ method takes
      all keyword input and sets it to instance attributes with the
      key as the attribute name and attribute values being Field
      objects.  Note that the instantiation keyword parameter does
      not have to itself be a Field object; if the keyword is set
      to a plain value, the __init__method will form a Field object
      around that value and set the Field object to the instance
      attribute.

    * runlists:  Dictionary of lists of methods, routines, and other
      run lists that can be executed by the run_list method.  Altering
      lists in this dictionary enables one to simply change order
      or execution or insert/delete methods.  Runlist's keys are
      strings, and are names that describe that the category of
      routines does.  Run list names should not be preceeded by two
      underscores (though runlist elements may be very private
      variables), nor should runlist names be the same as any
      instance attribute.  See list of run lists below for details
      about the lists.

    * sodir:  Name of temporary directory containing all the copies
      of shared object files for this instance of Qtcm.  This
      temporary directory is located in the default temporary
      directory, whose location is platform dependent, using the
      rules of the module tempfile.  The name of the temporary
      directory is generated randomly, and should not conflict with
      any existing temporary directory.  String.  Set on model
      instantiation.


    Public Instance Methods:
    * get_qtcm_item:  Get field from the compiled QTCM model (same
      as get_qtcm1_item).

    * get_qtcm1_item:  Get field from the compiled QTCM model.

    * make_snapshot:  Make copy of current state of run session
      variables that would be used for a restart.

    * more_first_method_at_atm_oc_step:  More of first method
      executed at the atmosphere-coupling timestep.  This gives an
      easy way to change the model at the atmosphere-coupling
      timestep, after the update of the interval attribute:  Just
      overload this method.

    * plotm:  Plot mean output for model field id.  Matplotlib is
      used.

    * qtcm:  Run the qtcm atmosphere over a coupling interval step.

    * run_list:  Run runlist of run lists and/or instance methods.

    * run_session:  Run model session.

    * set_qtcm_item:  Set Python-accessible compiled QTCM model 
      fields (same as set_qtcm1_item).

    * set_qtcm1_item:  Set Python-accessible compiled QTCM model
      fields.

    * sync_set_py_values_to_snapshot:  Set Python attributes to
      snapshot values from a dictionary created by make_snapshot.

    * sync_set_qtcm_items_to_all_py_values:  Synchronize so that
      any Python attribute that corresponds to a compiled model
      Python-changable variable is set, on both the Python attribute
      side and the compiled model side, to the value of the Python
      attribute.

    * sync_set_all_py_values_to_qtcm_items:  Set all Python-level
      attributes to QTCM compiled model values.  Synchoronize so
      that any Python attribute that corresponds to a compiled model
      Python-changable variable is set to the value of the QTCM
      compiled variable.

    * varinit:  Initialize model variables in a run session.


    List of Run Lists:
    * 'atm_oc_step':  Methods to run at an atmosphere-ocean coupling
      time step.  As default, set to methods that set the calendar,
      run the ocean (make SST for atmo), run the atmosphere, and
      output all output variables.  Set on model instantiation.
      See method __run_parts for details.

    * 'init_model':  Methods to initialize model variables.  Not
      all model variables are initalized in these methods.  Note
      that many of these methods use the values of QTCM parameters
      as defined as instance attributes (and the defaults of which
      are set by the _set_all_qtcm_scalar_fields_to_defaults method).
      As default, set to initialize the atmosphere, ocean, and
      output files.  Set on model instantiation.  See method
      __run_parts for details.

    * 'qtcminit':  Duplicates functionality of the Fortran subroutine
      qtcminit, for compiled_form 'parts'.


    Private Instance Attributes (partial list):
    * _monlen:  Integer array of the number of days in each month.
      Assumes a 365 day year.

    * __qtcm:  The compiled Fortran model for this instance.  Set on
      instantiation.

    * _qtcm_fields_ids:  Field ids for all default QTCM fields.
      NB: The model assumes this is a list of all the default fields.
      Set on instantiation.

    * _runlists_long_names:  Descriptions of the standard run lists.


    Private Instance Methods (partial list):
    * _bartropic_mode_at_atm_step:  Calculate the atmosphere
      barotropic mode at atmosphere timestep.

    * _first_method_at_atm_oc_step:  First method executed at the 
      atmosphere-coupling timestep.

    * __run_parts:  Run parts model starting at the atmosphere-ocean
      coupling level.  This method duplicates the functionality of
      the driver subroutine in the original compiled QTCM model.

    * _set_qtcm_array_item_in_model:  Set Python-accessible QTCM
      array settings in compiled model.

    * _set_qtcm_scalar_item_in_model:  Set Python scalar variable
      in the compiled QTCM model.


    Description of Compiled QTCM Modules Available to Class:
    * _qtcm_full_365.so:  QTCM shared object file for the full
      model, compiled assuming a year is 365 days.  This is essentially
      the default version of QTCM, with netCDF output used.

    * _qtcm_parts_365.so:  QTCM shared object file for the model
    separated so units within the atmosphere timestep are separate
    compiled units.
    """




#----------------- Class Qtcm Private Method:  __init__ ----------------

    def __init__(self, **kwds):
        """Initialize Qtcm instance.

	Items in kwds are attributes of the model instance to be
	set.  There are two types of attributes:  (1) Field objects,
	and (2) Non-Field objects.  In both cases you can pass in
	the name and value in the input kwds dictionary.  For Field
	object attributes, the kwds key is the id, and the kwds
	value can be either the value of the Field object or a Field
	object itself.  For Non-Field objects, the kwds key:value
	pair directly provides the instance attribute name and
	value.  A more detailed description of this input and other
	input is provided in the class documentation.

        Examples:
          model = Qtcm(dt=1200., compiled_form='parts')
          model = Qtcm(dt=Field('dt'), compiled_form='full')
          model = Qtcm(dt=Field('dt', 1200), compiled_form='full')
          inputs = {}
          inputs['dt'] = 1200.
          inputs['title'] ='QTCM spinup part 1 test (aquaplanet)'
          inputs['compiled_form'] = 'parts'
          model = Qtcm(**inputs)

	NB:  Not all compiled QTCM model variables are accessible
	at the Python level, and vice versa.  Only the variables
	in defaults (whose names are specified by the _qtcm_fields_ids
	attribute) can be passed between the compiled QTCM model
	and the Python model instance.

	Additionally, there is only a full-as-possible syncronization
	between compiled QTCM model fields and Python fields at the
	beginning and end of the execution of a Python method.  If
	the compiled QTCM model changes a variable as the model is
	running, the Python field that variable corresponds to will
	in general not also be changed, since you need to execute
	get_qtcm_item at the Python level to retrieve it.  This
	only occurs with Python methods (e.g., run_session, varinit).
        """
	#- Set instance private attributes:

        self._qtcm_fields_ids = defaults.qtcm_fields_ids


        #- Set compiled_form keyword, which is required input:

        if kwds.has_key('compiled_form'):
            setattr( self, 'compiled_form', kwds['compiled_form'] )
        else:
            raise ValueError, 'must pass in complied_form on instantiation'


        #- Set compiled model as private variable and set default 
        #  values for QTCM parameters:

        self._set_compiled_qtcm_attr()
        self._set_all_qtcm_scalar_fields_to_defaults()


	#- For QTCM variables entered in as keyword input parameters,
	#  overwrite the values with the ones passed in on instantiation.
	#  For all other keyword input parameters, set as attributes.
	#  (Recall the set_qtcm_item method also makes sure the Qtcm
	#  instance attribute matches the compiled model value and
	#  that set_qtcm_item takes as input Field variables or Field
	#  variable values:

        for ikey in kwds.keys():
            if ikey in self._qtcm_fields_ids:
                self.set_qtcm_item(ikey, kwds[ikey])
            else:
                setattr(self, ikey, kwds[ikey])


        #- Set other attributes:
        #  * coupling_day:  If not yet set via kwds, initialized to 
        #    None.
	#  * init_with_instance_state:  If not yet set via kwds, set
	#    to default value of True.

        if not hasattr(self, 'coupling_day'):
	    self.coupling_day = Field( 'coupling_day', None, units='dy' \
				     , long_name='current coupling day')
        else:
	    self.coupling_day = Field( 'coupling_day', self.coupling_day \
                                     , units='dy' \
				     , long_name='current coupling day')

        if not hasattr(self, 'init_with_instance_state'):
            self.init_with_instance_state = True

        if not hasattr(self, '_monlen'):
            self._monlen = N.array([31,28,31,30,31,30,31,31,30,31,30,31])


        #- Set run lists:

        self.runlists = {}

        #+ Runlist for initializng the entire model (atmos. and ocean):
        self.runlists['init_model'] = \
            [ 'qtcminit',
              '__qtcm.wrapcall.woceaninit',
              '__qtcm.wrapcall.woutpinit' ]

        #+ Runlist for initializing the atmos. part of the model:
        self.runlists['qtcminit'] = \
            [ '__qtcm.wrapcall.wparinit',
              '__qtcm.wrapcall.wbndinit',
              'varinit',
              {'__qtcm.wrapcall.wtimemanager':[1,]},
              'atm_physics1' ]

        #+ Runlist for calculating atmospheric physics at one instant:
        self.runlists['atm_physics1'] = \
            [ '__qtcm.wrapcall.wmconvct',
              '__qtcm.wrapcall.wcloud',
              '__qtcm.wrapcall.wradsw',
              '__qtcm.wrapcall.wradlw', 
              '__qtcm.wrapcall.wsflux' ]

        #+ Runlist for coupling atmos./ocean at a coupling timestep:
        self.runlists['atm_oc_step'] = \
            [ '_first_method_at_atm_oc_step',
              {'__qtcm.wrapcall.wtimemanager':[self.coupling_day,]},
              {'__qtcm.wrapcall.wocean': [self.interval, self.coupling_day]},
              'qtcm',
              '__qtcm.wrapcall.woutpall' ]

        #+ Runlist for calculating atmos. at one atmos. timestep:
        self.runlists['atm_step'] = \
            [ 'atm_physics1',             #* physics assoc. w/ temp. mode T1
              '__qtcm.wrapcall.wsland1',  #* land module
              '__qtcm.wrapcall.wadvctuv', #* advection of mom'm u1,v1,u0,v0 
              '__qtcm.wrapcall.wadvcttq', #* advection of T1,q1
              '__qtcm.wrapcall.wdffus',   #* diffusion of u1,v1,u0,v0,T1,q1
              '__qtcm.wrapcall.wbarcl',   #* baroclinic mode
              '_bartropic_mode_at_atm_step', #* barotropic mode
              '__qtcm.wrapcall.wvarmean' ]   #* time mean (sum up variables, 
                                             #  output as spec'd. in output.F)

        #+ Runlist for calculate the atmos. barotropic mode at the
        #  barotropic timestep:
        self.runlists['atm_bartr_mode'] = \
            [ '__qtcm.wrapcall.wsavebartr',  #+ save u0,v0 for surf. geopot.
              '__qtcm.wrapcall.wbartr',      #+ barotropic mode
              '__qtcm.wrapcall.wgradphis']   #+ surf. geopot. gradient


        #- Set run lists descriptions:

        self._runlists_long_names = {}
        self._runlists_long_names['init_model'] = \
            'initialize the entire model, i.e., the atmosphere and ocean ' + \
            'components and output'

        self._runlists_long_names['qtcminit'] = \
            'initialize the atmosphere portion of the entire model'

        self._runlists_long_names['atm_physics1'] = \
            'calculate atmospheric physics at one instant'

        self._runlists_long_names['atm_oc_step'] = \
            'calculate the atmosphere and ocean models at a coupling timestep'

        self._runlists_long_names['atm_step'] = \
            'calculate the entire atmosphere at one atmosphere timestep'

        self._runlists_long_names['atm_bartr_mode'] = \
            'calculate the atmospheric barotropic mode at the barotropic ' + \
            'timestep'


        #- Check run list keys for descriptions and the actual lists
        #  match, at least on instantiation:

        if len(self.runlists) != len(self._runlists_long_names):
            raise ValueError, 'number of run lists incorrect'
        for ikey in self.runlists.keys():
            if not self._runlists_long_names.has_key(ikey):
                raise ValueError, 'description and runlist key mismatch'




#----------------- Class Qtcm Private Method:  __del__ -----------------

    def __del__(self):
	"""Method to execute before garbage collection.

	Delete the sodir temporary directory.
        """
        shutil.rmtree(self.sodir)




#-------------- Class Qtcm Public Method:  get_qtcm_item ---------------

    def get_qtcm_item(self, key):
        """Get field from the compiled QTCM model.

	Gets the value of the variable key from the compiled QTCM
	model.  (Note that not all compiled QTCM model variables
	are accessible to Python.)  Method does nothing else (such
	as saving that value as an object attribute).

	For scalar variables, in order to know which Fortran routine
	to use to access the item, this method needs to know the
	type of the variable named key.  For arrays, it needs to
	know the rank of the array.  It assumes that type is the
	same as the corresponding Field object default value, and
	rank is the same as the corresponding Field object default
	value.

	If the compiled QTCM model variable is unreadable, the
	custom exception FieldNotReadableFromCompiledModel is thrown.

        Positional Input Parameter:
        * key:  Name of QTCM variable.  String.  Scalar.

	Output:
	* Returns the value of the variable named key.  Scalar
	  string or numeric value, or numeric array.  The value
          returned is a copy of the values, not a reference to the
          array in memory.
        """
	#- Get information about the field defined by key, using
	#  the default values:

        tmpfield = Field(key)
        value_type = tmpfield.typecode()
        value_rank = tmpfield.rank()


	#- Execute compiled model subroutine to extract the variable
	#  from the compiled model.  Which compiled model subroutine
	#  is used for the extraction depends on value_type and
	#  value_rank.  For arrays, the __qtcm.setbypy module
	#  attribute corresponding to an array (e.g., real_rank1_array)
	#  is set to None, thus deallocating the attribute in
	#  __qtcm.setbypy, at the end of extracting the value from the
	#  compiled model.  Note that if the variable referenced by key
        #  is not readable in the compiled module (which is given via a
        #  module logical variable is_readable), a custom exception is
        #  raised.  There is a single return line for this method.  The
        #  return is set to the local variable tmp, which is what the
        #  value of the variable obtained from the compiled model is
        #  set to:

        if value_type in N.typecodes['Float']:
            if value_rank == 0:
                tmp = copy.copy(self.__qtcm.setbypy.getitem_real(key))
                if not self.__qtcm.setbypy.is_readable:
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
            else:
                self.__qtcm.setbypy.getitem_real_array(key)
                if not self.__qtcm.setbypy.is_readable:
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
                else:
                    if value_rank == 1:
                        tmp = copy.copy(self.__qtcm.setbypy.real_rank1_array)
                        self.__qtcm.setbypy.real_rank1_array = None
                    elif value_rank == 2:
                        tmp = copy.copy(self.__qtcm.setbypy.real_rank2_array)
                        self.__qtcm.setbypy.real_rank2_array = None
                    elif value_rank == 3:
                        tmp = copy.copy(self.__qtcm.setbypy.real_rank3_array)
                        self.__qtcm.setbypy.real_rank3_array = None
                    else:
                        raise ValueError, 'bad rank value'

        elif value_type in N.typecodes['Integer']:
            if value_rank == 0:
                tmp = copy.copy(self.__qtcm.setbypy.getitem_int(key))
                if not self.__qtcm.setbypy.is_readable:
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
            else:
                raise ValueError, 'variable option does not exist'

        elif value_type in num.typecodes['Character']:
            if value_rank == 0:
                tmp = self.__qtcm.setbypy.getitem_str(key)
                if not self.__qtcm.setbypy.is_readable:
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
                tmp = tmp.strip()
            else:
                raise ValueError, 'variable option does not exist'

        else:
            raise TypeError, 'Value type does not exist'


	#- Return value obtained from compiled model (should already
	#  be a copy):

        return tmp




#-------------- Class Qtcm Public Method:  get_qtcm1_item --------------

    def get_qtcm1_item(self, *args, **kwds):
        """Get field from the compiled QTCM model.

	Identical to get_qtcm_item.  See get_qtcm_item description.
        """
        return self.get_qtcm_item(*args, **kwds)




#--------------- Class Qtcm Public Method:  make_snapshot --------------

    def make_snapshot(self):
        """Make copy of current state of run session variables.

	The snapshot is a copy of the variables of the current state
	of the run session that would be used for a restart.  The
	variables are saved in a dictionary, and the dictionary is
	saved in the attribute snapshot.  As most of the variables
	are Field objects, most of the values in the dictionary are
	Field objects.  Any variables that don't exist are just
	left out of the dictionary.
        """
        self.snapshot = {}
        for ikey in _init_vars_keys:
            if hasattr(self, ikey):
                self.snapshot[ikey] = copy.deepcopy(getattr(self, ikey))




#------------------ Class Qtcm Private Method:  _plot ------------------

    def _plot(self, id, datafn, **kwds):
        """Plot model field id from the data in netCDF file datafn.

        See docstrings for submodule plot in package qtcm, module
        function plot_ncdf_output.
        """
        plot_ncdf_output(id, datafn, **kwds)




#------------------ Class Qtcm Public Method:  plotm -------------------

    def plotm(self, id, **kwds):
        """Plot mean output for model field id.

	Plotting is done using Matplotlib, and can accomodate either
	line plots or contour plots.  This method is designed for
	a quick look at the output, and thus only a few plotting
	parameters can be adjusted.

        Positional Input Parameter:
        * id:  Name of the id of the field to plot.  String.

        Input Keyword Parameters:
	* lat:  Latitude range to plot.  List of the min and max
	  value, or a scalar, in degrees.  If the keyword is not 
          passed in, the entire range is plotted.

	* lon:  Longitude range to plot.  List of the min and max
  	  value, or a scalar, in degrees.  If the keyword is not
          passed in, the entire range is plotted.

	* time:  Time range to plot.  List of the min and max value,
	  or a scalar, in days.  If the keyword is not passed in,
          the entire range is plotted.

        * fn:  Filename to write figure out to.  String.  If set to
          None, plot is output to screen.

	* levels:  If a contour plot, the min, max, and interval
	  of the contour levels to be plotted.  If set to None, 
          contour levels are calculated automatically.

	* filled:  If is True, and the plot is a contour plot, a
	  filled contour plot is drawn along with a horizontal
	  colorbar.  If not True, a line contour plot is drawn.  If
	  plot is a line plot, this keyword is ignored.

        * title:  Title of plot.  String.

        * xlabel:  Label for the x-axis.

        * ylabel:  Label for the y-axis.

	* nlatlon:  For lat vs. lon contour plots, this sets the
	  approximate number of meridions and parallels to annotate.
	  Ignored for all other types of plots.  Default is 8.

	* tmppreview:  If True, and your platform is 'darwin',
	  instead of screen display the method will write the plot
	  to a temp file and display that file using Preview.  You
	  must delete the temp file yourself (it's located probably
	  in /tmp, and is created using the tempfile module default
	  settings).  The variable is boolean.

	The default of all input keywords is None, unless otherwise
	noted.  Depending on which series of keywords is chosen,
	the range is chosen accordingly.  Title, x-axis and y-axis
	labels are automatically added, unless overridden by the
	title, xlabel, and ylabel keywords.

	For a lat vs. lon plot, the contour plot is superimposed
	onto a cylindrical projection map of the Earth with continents
	drawn and labeled meridians and parallels.  The title also
	includes the model time, and x- and y-axis labels are not
	drawn.

        Examples for model an instance of Qtcm:
	* model.plotm('Qc', lat=1.875):  A time vs. longitude contour
	  plot is made for the full range of time and longitude,
	  at the latitude 1.875 deg N, for mean precipitation.  The
	  period over which the mean is taken is self.ntout.

	* model.plotm('Evap', lat=1.875, lon=[100,200]):  A time
	  vs. longitude contour plot of evaporation is made for the
	  longitude points between 100 and 200 degrees east, at the
	  latitude 1.875 deg N.  The period over which the mean is
	  taken is self.ntout.

	* model.plotm('cl1', lat=1.875, lon=[100,200], time=20):
	  A deep cloud amount vs. longitude line plot is made for
	  the longitude points between 100 and 200 degrees east,
	  at the latitude 1.875 deg N, at day 20 of the model run.
	  The period over which the mean is taken is self.ntout.
        """ 
        datafn = os.path.join( self.outdir.value, 
                               "qm_" + self.runname.value + ".nc" )
        self._plot(id, datafn, **kwds) 




#------------------- Class Qtcm Public Method:  qtcm -------------------

    def qtcm(self):
        """Run the qtcm atmosphere over a coupling interval step.

	This method duplicates the functionality of the qtcm
	subroutine in the original compiled QTCM model.  It assumes
	all fields that can be synchronized between the Python and
	compiled QTCM model levels have been so before the calling
	of this method.  The coupling interval step is given by the
	interval attribute.  It is meant to be called once in a run
	session, and assumes that the compiled_form is 'parts'.
        """
        #- Get time changng boundary data, SST, albedo etc.  Check that
        #  the number of atmos. time steps within the atmos.-ocean
        #  coupling interval is a divisor of 1 day (nastep is the 
        #  number of atmospheric time steps within one air-sea coupling 
        #  interval and set the value of nastep in the compiled QTCM
        #  model:

        self.run_list(['__qtcm.wrapcall.wgetbnd',])
        nastep = self.interval.value * (86400./self.dt.value)
        if N.mod(nastep, self.mt0.value) != 0:
           raise ValueError, 'Fatal error: mt0*dt not a divisor of 1 day'
        self.set_qtcm_item('nastep', int(nastep))


        #- The land-atmosphere coupling loop:

        for self.it.value in xrange(1, self.nastep.value+1):
            self.run_list(['atm_step',])




#----------------- Class Qtcm Public Method:  run_list -----------------

    def run_list(self, runlist):
        """Run runlist of run lists and/or instance methods.

	Run through list of elements that specify other run lists
	or instance method names to execute through in runlist
	order.  Methods with private attribute names are automatically
	mangled as needed to become executable by the method.  Note
	that if an item in runlist is an instance method, it should
	be the entire name (after "self.") of the callable method,
	separated by periods as appropriate.

        Input Via Arguments:

	* runlist:  List whose elements are 1-element dictionaries
	  or strings.  The list can contain a mix of 1-element
	  dictionaries and strings, or just one of those types:

	  + If 1-element dictionaries:  The key of the dictionary
	    is a string, and is the name of the method to execute.
	    The value of the entry is a list, which gives the
	    positional input parameters, or is None, if there are
	    no input parameters.

	  + If strings:  Each string is the name of a run list or
	    the name of the method to execute.  All methods are
	    assumed to not require any positional input parameters.

	The methods in runlist are called in the order given in
	runlist.  For each element, we first check if the key name
	corresponds to the key of an entry in self.runlists.  If
	so, run_list is executed using that run list (recursive
	call).  If the key name is not a run list, we check if it
	is a method of the instance, and if so the method is called.
	Any other value throws an exception.

	If input parameters for a method are of class Field, we
	first try to pass the parameters into the method as is
	(i.e., as Field object(s)).  If that fails, we pass its
	parameters in as the value of the Field object.

	Examples:
	  a = Qtcm(...) 
          a.run_list( ['qtcminit' \
                    , {'sync_set_all_py_values_to_qtcm_items':None} \
		    , '__qtcm.driver' \ 
                    , {'set_qtcm_item': ['outdir', '/home/jlin']}] )
          a.run_list( [{'sync_set_all_py_values_to_qtcm_items':None} \
		    , {'__qtcm.driver':None} \ 
                    , {'set_qtcm_item': ['outdir', '/home/jlin']} )
          a.run_list(['sync_set_all_py_values_to_qtcm_items',])
	"""
        #- Comparison cases for string and dict types:

        str_type = type('a')
        dict_type = type({'a':None})


        #- Check runlist is a list.  This prevents undetected errors
        #  if another iterable type (e.g., a string) is passed in as 
        #  runlist:

        if type(runlist) != type([]):
            raise TypeError, 'runlist must be a list'


        #- Loop through all methods in runlist:

        for imethod in runlist:


	    #+ Case for runlist item being a string:  Access name of
	    #  run list or method, and mangle private name as needed
	    #  if is a method.  If item is itself the name for a 
            #  run list, recursively call run_list.  If is a method, 
            #  recursively find the callable method, and call the 
            #  method:

            if type(imethod) == str_type:
                imethodname = imethod
                if self.runlists.has_key(imethodname):
                    if hasattr(self, imethodname):
                        raise AttributeError, \
                            'run list cannot be same name as instance attrib.'
                    self.run_list(self.runlists[imethodname])
                else:
                    if imethodname[0:2] == '__':
                        imethodname = '_' + self.__class__.__name__ \
                                    + imethodname
                    imethodname_seplist = imethodname.split('.')
                    f = self
                    for subname in imethodname_seplist:
                        f = getattr(f, subname)
                    f()


            #+ Case for runlist item being a 1-entry dictionary:

            elif type(imethod) == dict_type:

                #* Access name of method, and mangle private name if needed:

                imethodname = imethod.keys()[0]
                if imethodname[0:2] == '__':
                    imethodname = '_' + self.__class__.__name__ + imethodname
                imethodname_seplist = imethodname.split('.')


                #* Recursively find the callable method and extract the
                #  list of input parameters for the method:

                f = self
                for subname in imethodname_seplist:
                    f = getattr(f, subname)
                inputparams = imethod.values()[0]


		#* Call the method of the runlist item:  For the case of
		#  no input parameters, call without any parameters.  If
		#  the input parameters cause a TypeError while being
		#  passed into the method, assume it's because one
		#  or more of the parameters is a Field object.  Replace
		#  the Field input parameters with their values, and
		#  retry calling the method:

                if inputparams == None:
                    f()
                else:
                    try:
                        f(*inputparams)
                    except TypeError:
                        inputparams_vals = []
                        for iparam in inputparams:
                            if type(iparam) == type(_test_field):
                                inputparams_vals.append(iparam.value)
                        try:
                            f(*inputparams_vals)
                        except:
                            print "Unexpected error (A):", sys.exc_info()[0]
                            raise
                    except:
                        print "Unexpected error (B):", sys.exc_info()[0]
                        raise


            #+ Case for runlist item of bad input type:

            else:
                raise TypeError, 'bad run_list element input type'




#---------------- Class Qtcm Public Method:  run_session ---------------

    def run_session(self, **kwds):
        """Run model session.

        Run the QTCM model in a model session.  A model session is
        defined as a "complete" model run, at the end of which
        restart files are written and the Python Qtcm instance is
        synchronized to the Fortran model.
        
        The following tasks are done by this method, in this order:
        * The attribute snapshot, if present, is deleted.
        * The compiled model is synchronized with the Python model 
          instance so that any Python attribute that corresponds to a 
          compiled model Python-changable variable is set to the value 
          of the Python attribute.  If the compiled model variable is
          not accessible, nothing is done, and the Python attribute is
          left unchanged.
        * The model is run until lastday.  Restart and output files 
          are written accordingly.
        * The Python-changable model attributes are set to the compiled
          model's values.  If the compiled model variable is not 
          accessible, an exception is raised.
        * A snapshot of the Python attributes that store the variables
          that would be used in a restart run is taken and stored as
          the attribute snapshot.

        Input Via Keywords Arguments:

	* cont:  If set to False, the run session is not a continuation
	  of the previous run, but a new run session.  If set to
	  True, the run session is a continuation of the previous
	  run, and whatever the field values are in the Python
	  instance are used in the model (if init_with_instance_state
	  is True).  If set to an integer greater than zero, the run
	  session is a continuation just like cont=True, but the
	  value set to cont is used for lastday and replaces
	  lastday.value.  This keyword has no effect if
	  compiled_form='full'.  Default is False.  (Note whatever
	  cont is set to in this method call is stored as attribute
	  _cont, in case you really need to access it elsewhere.)
        """
        #- Set default values for arguments and keywords:

        if kwds.has_key('cont'):
            cont = kwds['cont']
        else:
            cont = False
        self._cont = cont

        if type(cont) == type(1):
            self.lastday.value = cont
        elif type(cont) == type(False):
            pass
        else:
            raise TypeError, 'cont keyword must be integer or boolean'


# @@@this option doesn't work right.  only python stdout is
#    written while compiled code generate stdout seems to go to
#    /dev/null.  and the file isn't closed at the end of this
#    method executing.  stdout doesn't come back to screen.
#
#        * stdout:  String.  If set to:
#             'screen':  Writes to screen only.
#             'logfile':  Writes to logfile only; file opened in mode "a+".
#        * logfile:  String.  Name of logfile.  Set to os.devnull to send it
#          to /dev/null (or equivalent in your OS).
#        if kwds.has_key('stdout'):
#            stdout = kwds['stdout']
#        else:
#            stdout = 'screen'
#
#        if kwds.has_key('logfile'):
#            logfile = kwds['logfile']
#        else:
#            logfile = 'stdout.log'
#
#
#        #- Open stdout redirect, if needed:
#
#        if stdout=='screen':
#            pass
#        elif stdout=='logfile':
#            all_stdout_log_file = open(logfile, 'a+')
#            os.dup2(all_stdout_log_file.fileno(), 1)


	#- Delete snapshot attribute, if present.  Set the values
	#  of Python-setable variables in the compiled code to the
	#  value at the Python-level:

        if hasattr(self, 'snapshot'):  del self.snapshot
        self.sync_set_qtcm_items_to_all_py_values()


        #- Run model:

        if self.compiled_form == 'full':
            self.run_list(['__qtcm.driver',])
        elif self.compiled_form == 'parts':
            self.run_list(['__run_parts',])
        else:
            raise ValueError, 'Compiled form not recognized'


        #- Set the values of Python-setable variables at the Python-level
        #  to the values in the compiled code.  Take snapshot of model
        #  variables that would be used in another session:

        self.sync_set_all_py_values_to_qtcm_items()
        self.make_snapshot()


#@@@this doesn't seem to work right, so comment out.
#        #- Close stdout redirect, if needed:
#
#        if stdout=='screen':
#            pass
#        elif stdout=='logfile':
#            all_stdout_log_file.close()




#-------------- Class Qtcm Public Method:  set_qtcm_item ---------------

    def set_qtcm_item(self, *args):
        """Set Python-accessible compiled QTCM model fields.

	Sets the value in the compiled QTCM model as well as at the
	Python level.  If no value given, the default value is used.
	When the compiled model variable is set, a copy of the
	Python value input via the parameter list is passed to the
	Fortran model, not a reference.

	If this method is called with a single positional input
	argument, and that argument is a string, the compiled QTCM
	model variable is set to the default value of the Python
	counterpart of that same name.  If the single positional
	argument is a Field variable, the compiled QTCM model
	variable is set to that Field variable (the id and value
	attributes of the Field variable thus corresponding to the
	name and value of the compiled QTCM variable, respectively).

        Positional Input Parameters (for 1 argument):
	* Name of QTCM variable (scalar string) or a Field variable.

        Positional Input Parameters (for 2 arguments):
        * key:  Name of QTCM variable.  String.  Scalar.
	* value_in:  Value of variable or a Field object.  Scalar
	  or numeric array.  If you want to set an array to a single
	  quantity, value can be a scalar.  Type of value must be
	  the same as the type of Field(key), which is the same as
	  in the compiled QTCM model; the routine doesn't check for
	  this, however, as the compiled will return a fatal error
	  if this mismatch happens.

	Some compiled QTCM model variables are not ready to be set.
	An example is a compiled QTCM model pointer variable prior
	to the pointer being associated with a target (this would
	result in a bus error).  In such cases, this method will
	throw a FieldNotReadableFromCompiledModel exception, nothing
	will be set in the compiled QTCM model, and the Python
	counterpart variable (if it previously existed) would be
	left unchanged.  Otherwise, both the compiled QTCM model
	variable and its Python attribute counterpart are set by
	the method to the same value, overwriting any previous
	values held by either.
        """
	#- Set key and value based on key and value_in.  After this
	#  section, key will be a string and value_in either a scalar
	#  or a numeric array:

        if len(args) == 1:
            if type(args[0]) == type(_test_field):
                key = args[0].id
                value_in = args[0].value
            else:
                key = args[0]
                value_in = Field(key).value   #- This sets value to the default
        elif len(args) == 2:
            key = args[0]
            value_in = args[1]
            if type(value_in) == type(_test_field):   #- value_in a Field obj.
                if key != value_in.id: 
                    raise ValueError, 'inconsistent input args.'
                value_in = args[1].value
        else:
            raise ValueError, 'set_qtcm_item uses only 1 or 2 args'


	#- Make value a copy of the input value, not a reference 
        #  assignment.  Find type and rank of value.  Use value 
        #  instead of value_in for rest of method:

        value = copy.copy(value_in)
        value_rank = N.rank(value)
        if value_rank == 0:
            value_dtype = type(value)
        else:
            value_dtype = num.typecode(value)

        # (Here I provide code to check type of value is good, if you
        # want to do this at the Python level.  I commented it out,
        # however, for speed gains, since the compiled QTCM model
        # already checks this, returning a fatal error if there is a
        # mismatch):
        #
        # field_default_dtype = Field(key).typecode()
        # if ( (value_dtype in N.typecodes['Integer']) and \
        #      (field_default_dtype not in N.typecodes['Integer']) ) or \
        #    ( (value_dtype in N.typecodes['Float']) and \
        #      (field_default_dtype not in N.typecodes['Float']) ) or \
        #    ( (value_dtype in num.typecodes['Character']) and \
        #      (field_default_dtype not in num.typecodes['Character']) ):
        #     raise TypeError, 'Type of value does not match default'


	#- Set compiled QTCM variable using methods chosen based
	#  upon upon field rank (we use field rank to enable you to
	#  set all values in an array to the same value).  If
	#  FieldNotReadableFromCompiledModel is returned, nothing is
	#  set in the compiled QTCM model and the exception is sent 
        #  to continue upwards.  Set Python attribute accordingly:

        field_rank = Field(key).rank()


        #+ For rank of the field 0 (i.e. a scalar).  The first line
        #  after the "try" is to just see if the variable can be set:

        if field_rank == 0:
            try:
                tmp = self.get_qtcm_item(key)
                self._set_qtcm_scalar_item_in_model(key, value)
                setattr( self, key, Field(key, value) )
            except FieldNotReadableFromCompiledModel:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not writeable'
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise


        #+ For rank of the field non-zero (i.e., an array).  Note that
        #  in this section, we ensure that the Field value is also an 
        #  array as in the compiled model:

        else:
            try:
                ashape = self._set_qtcm_array_item_in_model(key, value)
                if value_rank == 0:         
                    valarr = N.empty(ashape, dtype=value_dtype)
                    valarr.fill(value)         
                else:
                    valarr = value
                setattr( self, key, Field(key, valarr) )
            except FieldNotReadableFromCompiledModel:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not writeable'
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise




#-------------- Class Qtcm Public Method:  set_qtcm1_item --------------

    def set_qtcm1_item(self, *args, **kwds):
        """Set Python-accessible compiled QTCM model fields.

	Identical to set_qtcm_item.  See set_qtcm_item description.
        """
        return self.set_qtcm_item(*args, **kwds)




#------ Class Qtcm Public Method:  sync_set_py_values_to_snapshot ------

    def sync_set_py_values_to_snapshot(self, snapshot=None):
        """Set Python attributes to snapshot values.

	The snapshot is copies of the variables of the current state
	of the run session that would be used for a restart.  This
	method sets the Python attributes corresponding to the
	snapshot variables to the snapshot values.  It does not set
	anything on the compiled QTCM model side.

        Keyword Input Parameter:
	* snapshot:  The snapshot (a dictionary, following the rules
	  of method make_snapshot) that is the source of the value
	  for the syncronization.  If set to None, the instance
	  attribute snapshot is used.  Default is None.
        """
        if snapshot == None:
            for ikey in self.snapshot.keys():
                setattr(self, ikey, self.snapshot[ikey])
        else:
            for ikey in snapshot.keys():
                setattr(self, ikey, snapshot[ikey])




#--- Class Qtcm Public Method:  sync_set_qtcm_items_to_all_py_values ---

    def sync_set_qtcm_items_to_all_py_values(self):
	"""Set QTCM items in compiled model to all Python-level values.

	Synchronize so that any Python attribute that corresponds
	to a compiled model Python-changable variable is set, on
	both the Python attribute side and the compiled model side,
	to the value of the Python attribute.  Note this method
	only sets attributes that are already defined in the object;
	it does not create new attributes.  If a compiled QTCM model
	variable is not ready to be set (e.g., a pointer variable
	is not yet associated, that variable is not set in the
	compiled QTCM model, and its Python counterpart is left
	unchanged.
        """
        for ikey in self._qtcm_fields_ids:
            if hasattr(self, ikey):
                try:
                    self.set_qtcm_item(getattr(self, ikey))
                except FieldNotReadableFromCompiledModel:
                    pass
                except:
                    print "Unexpected error:", sys.exc_info()[0]
                    raise




#--- Class Qtcm Public Method:  sync_set_all_py_values_to_qtcm_items ---

    def sync_set_all_py_values_to_qtcm_items(self):
	"""Set all Python-level attributes to QTCM compiled model values.

	Synchoronize so that any Python attribute that corresponds
	to a compiled model Python-changable variable is set to the
	value of the QTCM compiled variable.  Note this method goes
        through all items listed in self._qtcm_fields_ids and sets 
        those values as object attributes.

	If a compiled QTCM model variable is not ready to be read
	(e.g., a pointer variable is not yet associated, a
	FieldNotReadableFromCompiledModel exception is raised,
	because for the situations where this method is called
	(usually after the model has run for some time), that
	situation should not occur.
        """
        for ikey in self._qtcm_fields_ids:
            try:
                setattr( self, ikey, Field(ikey, self.get_qtcm_item(ikey)) )
            except FieldNotReadableFromCompiledModel:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not readable'
            except:
                print "Unexpected error:", sys.exc_info()[0]
                raise




#------ Part of Run Private Method:  _first_method_at_atm_oc_step ------

    def _first_method_at_atm_oc_step(self):
        """First method executed at the atmosphere-coupling timestep.

	This is a private method.  If you wish to add more computations
	at the beginning of the atmosphere-ocean coupling timestep,
	overload more_first_method_at_atm_oc_step.
        """
        self.interval.value = self.get_qtcm_item('interval')
        self.more_first_method_at_atm_oc_step()




#----- Part of Run Public Method:  more_first_method_at_atm_oc_step ----

    def more_first_method_at_atm_oc_step(self):
        """More of first method executed at the atmo.-coupling timestep.

	This gives an easy way to change the model at the
	atmosphere-coupling timestep, after the update of the interval
        attribute:  Just overload this method.
        """
        pass




#------- Part of Run Private Method:  _bartropic_mode_at_atm_step ------

    def _bartropic_mode_at_atm_step(self):
        """Calculate the atmos. barotropic mode at atmos. timestep.

	The calculation is made at self.it timestep for the atmosphere,
	which is the time of day in terms of atmospheric timesteps.
	It assumes that the compiled_form is 'parts'.
	""" 
        if N.mod(self.it.value, self.mt0.value) == 0:
	    self.run_list(['atm_bartr_mode',])




#----------------- Part of Run Public Method:  varinit -----------------

    def varinit(self):
	"""Initialize model variables in a run session.

	Method duplicates the functionality of the Fortran QTCM1
	subroutine varinit, but with changes that will enable us
	to handle restarts at the Python-level in a dynamic way.
        This method only works with compiled_form set to 'parts'.

	If init_with_instance_state is False, this method just
	executes the compiled Fortran QTCM model varinit subroutine,
	and thus run_session initialization follows the rules set
	by all the standard QTCM model input parameters (e.g., the
	sign of day0, month0, and year0, the value of the mrestart
	flag, etc.).  If init_with_instance_state is True, variable
	initialization uses the algorithm described below.

	First, all prognostic pointer variables are associated.
	Next, initialization of prognostic variables and right-hand
	sides to default values (which for most of the variables
	is 0) occur for the cases where the corresponding Python
	attribute is not defined.  If the corresponding Python
	attribute is defined, the compiled QTCM model variable is
	set to the value that the attribute already has (the mrestart
	flag, given as the instance attribute mrestart, is ignored).
	The exception is day0, month0, and year0, which are overwritten
	with values derived from dateofmodel to set the run to start
	the day after dateofmodel.  If dateofmodel is less than or
	equal to 0, day0, month0, and year0 are set to their
	respective instance values, if valid (for invalid values,
	day0, month0, and year0 are all set to 1).

	Note for init_with_instance_state True or False, at the end
	of this method, day0, month0, and year0, may all be changed,
	and dateofmodel may be inconsistent with will be updated
	to match the values of day0, month0, and year0.
        """
	#- Dictionary of all prognostic variables and right-hand
	#  sides that could be initialized and their default initial
	#  values (note that these default values should be the same
	#  as given in the first part of the Fortran QTCM1 subroutine 
        #  varinit, prior to the mrestart test, but there is no check 
        #  in this method for this:

        init_dict = copy.copy(_init_prog_dict)


        #- If do not initialize QTCM variables with the instance
        #  state:  Run __qtcm.wrapcall.wvarinit, set all Python values 
        #  of these initialized variables to the compiled QTCM model 
        #  value, recalculate dateofmodel and put updated value into
        #  both the Python and compiled model sides:

        if not self.init_with_instance_state:
            self.run_list(['__qtcm.wrapcall.wvarinit',])

            update_list = init_dict.keys() + ['day0', 'month0', 'year0']
            for ikey in update_list:
                setattr(self, ikey, Field(ikey, self.get_qtcm_item(ikey)) )

	    self.dateofmodel.value = self.year0.value*10000 \
			           + self.month0.value*100 \
                                   + self.day0.value
            self.set_qtcm_item(self.dateofmodel)


        #- If do initialize QTCM variables with the instance state:

        else:
            #+ Associate pointer variables:

            if not self._cont:
                self.run_list(['__qtcm.varptrinit',])


	    #+ Calculate WD at 70% saturation and make that initial
	    #  value as part of init_dict.  STYPE equals 0 for ocean,
	    #  1 for forest, 2 for grass, and 3 for desert.  WD0 is
	    #  indexed in Fortran from 0 to NSTYPE-1, which matches
	    #  Python indexing so there doesn't need to be any index
	    #  shifting below:

	    tmpstype = self.get_qtcm_item('STYPE') 
            tmpwd0 = self.get_qtcm_item('WD0') 
            tmpwd = N.choose(tmpstype.astype(int), tmpwd0) * 0.7 
            init_dict['WD'] = tmpwd


            #+ Initialize all prognostic variables and right hand 
            #  sides if the attributes do not currently exist:

            for ikey in init_dict.keys():
                if not hasattr(self, ikey):
                    self.set_qtcm_item(ikey, init_dict[ikey])
                else:
                    self.set_qtcm_item(getattr(self, ikey))


	    #+ Calculate the day0, month0, year0 values for the day
	    #  after dateofmodel.  Checks that the day after
            #  dateofmodel is correct, i.e., if the next day is in
            #  the next month, the month, day, and year are changed
            #  as needed:

            yearr = self.dateofmodel.value / 10000
            monthr = N.mod(self.dateofmodel.value, 10000) / 100
            dayr = N.mod(self.dateofmodel.value, 100) + 1

            if self._monlen[monthr-1] < dayr:
                dayr = 1
                monthr = monthr + 1
                if monthr > 12:
                    monthr = 1
                    yearr = yearr + 1


	    #+ Set day0, month0, year0.  If dateofmodel is > zero,
	    #  change day0, etc. so the model will run the day after
	    #  dateofmodel.  If dateofmodel is <= 0, keep day0, etc.
	    #  unchanged if day0, etc. are valid (for invalid cases,
	    #  set day0, etc. all to 1).

            if self.dateofmodel.value > 0:
                self.year0.value = yearr
                self.month0.value = monthr
                self.day0.value = dayr
            else:
                if (self.day0.value < 1) or (self.day0.value > 31):
                    self.day0.value = 1
                if (self.month0.value < 1) or (self.month0.value > 12):
                    self.month0.value = 1
                if self.year0.value < 1:
                    self.year0.value = 1

	    #+ Update dateofmodel with the changed year0, month0,
	    #  day0 values:

	    self.dateofmodel.value = self.year0.value*10000 \
			           + self.month0.value*100 \
                                   + self.day0.value


            #+ Update the QTCM side with the changed Python date
            #  attributes:

            update_list = ['day0', 'month0', 'year0', 'dateofmodel']
            for ikey in update_list:
                    self.set_qtcm_item(getattr(self, ikey))




#--------------- Class Qtcm Private Method:  __run_parts ---------------

    def __run_parts(self):
        """Run parts model starting at the atmos.-oc. coupling level.

	This method duplicates the functionality of the driver
	subroutine in the original compiled QTCM model.  It assumes
	all fields that can be synchronized between the Python and
	compiled QTCM model levels have been so before the calling
	of this method.  It is meant to be called once in a run
	session, and assumes that the compiled_form is 'parts'.
        """
	#- Initialize variables.  If the run session is not a
	#  continuation of the previous run, initialize the model
	#  using the init_model run list.  Otherwise, do not do
        #  that initialization but just varinit which mainly sets
        #  the date-related fields for continuation run sessions:

        interval = self.interval.value
        lastday = self.lastday.value
        if not self._cont:
            self.run_list(self.runlists['init_model'])
        else:
            self.run_list(['varinit',])


	#- Main loop of ocean-atmosphere coupling.  Note that with
	#  the current way of handling netCDF output of the dayofmodel
	#  in the compiled QTCM model, startday has to be 1 or the
	#  timemanager function will feed the wrong value to the output
	#  routine, resulting in totally messed up output:

        startday = 1
        endday = lastday+interval

        for self.coupling_day.value in xrange(startday, endday, interval):
            self.run_list(self.runlists['atm_oc_step'])
            print 'Driver: Running for %i days at model date %i ' \
                  % (self.coupling_day.value, self.get_qtcm_item('dateofmodel'))


        #- Write restart file:

        self.run_list(['__qtcm.wrapcall.woutrestart',])


	#- Because startday has to be at 1, as noted in the comment
	#  earlier where startday is set, if this is a continuation
	#  run, we need to rewrite the values of the time output file
	#  variable in both the instantaneous and mean output files
	#  to make it correct.  Note that the compiled QTCM model
        #  output algorithm will output a first time value of 0; this
        #  method checks that that is the case.  This section also
        #  deals with the peculiar case where oldvalue[0] will not be
        #  a scalar, and when oldvalue is a short integer and Python
        #  will not automatically make a long to short int conversion:

        if self._cont:
            testoutdir = self.get_qtcm_item('outdir')
            if self.outdir.value != testoutdir:
                raise ValueError, 'bad outdir value for contiguous run'

            path = self.outdir.value
            suffix = self.runname.value
            fn_inst = os.path.join( path, "qi_" + suffix + ".nc" )
            fn_mean = os.path.join( path, "qm_" + suffix + ".nc" )

            for ifn in [fn_inst, fn_mean]:
                fileobj = S.NetCDFFile(ifn, mode='r+')
                oldvalue = fileobj.variables['time'].getValue()
                varshape = N.shape(oldvalue)
                if len(varshape) != 1:  raise ValueError, 'bad time shape'
                if oldvalue[0] != 0:
                    if N.squeeze(oldvalue[0]) != 0:
                        raise ValueError, 'bad time first value'
                newvalue = N.arange(varshape[0], dtype=num.typecode(oldvalue))
                fileobj.variables['time'].assignValue(newvalue)
                fileobj.close()




#- Class Qtcm Private Method:  _set_all_qtcm_scalar_fields_to_defaults -

    def _set_all_qtcm_scalar_fields_to_defaults(self):
        """Set all scalar QTCM fields to their default values.

	Sets the values in the compiled QTCM model as well as at
	the Python level to their default values.  Scalar values
	are defined as those with rank 0.
        """
        for ikey in self._qtcm_fields_ids:
            if N.rank(Field(ikey).value) == 0:
                self.set_qtcm_item(ikey)




#--------- Class Qtcm Private Method:  _set_compiled_qtcm_attr ---------

    def _set_compiled_qtcm_attr(self):
	"""Set compiled QTCM attribute in Qtcm instance.

	This method makes a copy of the needed compiled shared
	object file in a unique hidden subdirectory of the current
	working directory, imports that .so file, and sets that
	imported shared object library to the compiled QTCM attribute
	self.__qtcm.
        """
	#- Create temporary shared object directory for instance,
	#  self.sodir.  Extract the random part of the sodir
	#  directory name, put a '_' in front of it, and set to local
	#  variable.  Note how the sodir prefix is hard-wired in:

        self.sodir = tempfile.mkdtemp(prefix='qtcm_sodir_')


	#- Selected name of the shared object library of interest
	#  based on the compiled_form attribute:

        if self.compiled_form == 'full':
            soname = '_qtcm_full_365'
        elif self.compiled_form == 'parts':
            soname = '_qtcm_parts_365'
        else:
            raise ValueError, 'Compiled form not recognized'


	#- Assuming that the directory path of the original shared
	#  object is the same as that for the package_version module,
	#  the original shared object is copied to the temporary
	#  directory, and that copy (which will now be unique for
	#  this instance, if the temporary directory is put at the
	#  beginning of sys.path) is loaded and set to __qtcm.  The
	#  shared object name is removed from sys.modules, to prevent
	#  future instances of this class setting __qtcm to the
	#  modules in other instance sodirs.  To make sure nothing
	#  bad happens, at the end I check that the path for the
	#  __qtcm attribute is sodir:

        origsofile = \
            os.path.join( os.path.split(_package_version.__file__)[0] \
                        , soname + '.so' )
        shutil.copy2(origsofile, self.sodir)
        sys.path.insert(0, self.sodir)
        self.__qtcm = __import__(soname)
        sys.path.remove(self.sodir)
        del sys.modules[soname]

        if os.path.split(self.__qtcm.__file__)[0] != self.sodir:
            raise ValueError, 'Incorrect import of .so library'




#------ Class Qtcm Private Method:  _set_qtcm_array_item_in_model ------

    def _set_qtcm_array_item_in_model(self, key, value):
        """Set Python-accessible QTCM array settings in compiled model.

	Sets the value of arrays in the compiled QTCM model only.
	Custom exception FieldNotReadableFromCompiledModel is raised
	if the compiled QTCM model variable is not readable/writable.

        Positional Input Parameters:
        * key:  Name of QTCM variable.  String.  Scalar.

	* value:  Value of variable.  Must be numeric scalar or
	  real numeric array.  If a scalar, all values in the array
	  named key are set to that scalar value.  If value cannot
          be mapped onto the compiled QTCM model array successfully,
          you'll receive a bus error or other unexpected error.

        Output:
	* Sets value into the compiled QTCM model array variable,
	  and also returns the shape of the array that was set.
        """
        #- Set type and rank information:

        default_type = Field(key).typecode()
        value_type = num.typecode(value)
        field_rank = Field(key).rank()


	#- Test the value_type and default type are the same for
	#  floating point arrays:

        if (value_type in N.typecodes['Float']) and \
           (default_type not in N.typecodes['Float']):
            raise TypeError, 'value type different from default for key'


	#- Use the Fortran routine getitem_real_array to create the
        #  memory needed and to see if the compiled QTCM model variable 
        #  is readable.  Raise an exception if it is not readable.
        #  Set the Fortran module variable real_rank?_array to value,
        #  write that to the compiled QTCM model variable, and deallo-
        #  cate module variable real_rank?_array:

        if value_type in N.typecodes['Float']:
            if field_rank == 1:
                self.__qtcm.setbypy.getitem_real_array(key)
                if not self.__qtcm.setbypy.is_readable:
                    self.__qtcm.setbypy.real_rank1_array = None
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
                else:
                    if N.rank(value) == 0:
                        ashape = N.shape(self.__qtcm.setbypy.real_rank1_array)
                        tmpa = N.zeros(ashape, float)
                        tmpa.fill(value)
                        self.__qtcm.setbypy.real_rank1_array = tmpa
                    else:
                        ashape = N.shape(value)
                        self.__qtcm.setbypy.real_rank1_array = value
                    self.__qtcm.setbypy.setitem_real_array(key)
                    self.__qtcm.setbypy.real_rank1_array = None

            elif field_rank == 2:
                self.__qtcm.setbypy.getitem_real_array(key)
                if not self.__qtcm.setbypy.is_readable:
                    self.__qtcm.setbypy.real_rank2_array = None
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
                else:
                    if N.rank(value) == 0:
                        ashape = N.shape(self.__qtcm.setbypy.real_rank2_array)
                        tmpa = N.zeros(ashape, float)
                        tmpa.fill(value)
                        self.__qtcm.setbypy.real_rank2_array = tmpa
                    else:
                        ashape = N.shape(value)
                        self.__qtcm.setbypy.real_rank2_array = value
                    self.__qtcm.setbypy.setitem_real_array(key)
                    self.__qtcm.setbypy.real_rank2_array = None

            elif field_rank == 3:
                self.__qtcm.setbypy.getitem_real_array(key)
                if not self.__qtcm.setbypy.is_readable:
                    self.__qtcm.setbypy.real_rank3_array = None
                    raise FieldNotReadableFromCompiledModel, \
                          'Compiled model variable not readable'
                else:
                    if N.rank(value) == 0:
                        ashape = N.shape(self.__qtcm.setbypy.real_rank3_array)
                        tmpa = N.zeros(ashape, float)
                        tmpa.fill(value)
                        self.__qtcm.setbypy.real_rank3_array = tmpa
                    else:
                        ashape = N.shape(value)
                        self.__qtcm.setbypy.real_rank3_array = value
                    self.__qtcm.setbypy.setitem_real_array(key)
                    self.__qtcm.setbypy.real_rank3_array = None

            else:
                raise ValueError, 'bad rank value'

        elif value_type in N.typecodes['Integer']:
            raise TypeError, 'array type not yet supported'

        elif value_type in num.typecodes['Character']:
            raise TypeError, 'array type not yet supported'

        else:
            raise TypeError, 'value is of unsupported type'


        #- Return the shape of the array that was set:

        return ashape




#----- Class Qtcm Private Method:  _set_qtcm_scalar_item_in_model ------

    def _set_qtcm_scalar_item_in_model(self, key, value):
        """Set Python scalar variable in the compiled QTCM model.

	Sets value of Python scalar variable key in the compiled
	QTCM model.  Nothing else is done (e.g., nothing is set on
	the Python side).  Exception FieldNotReadableFromCompiledModel
	is raised if the compiled QTCM model variable is not
	readable/writable.

        Positional Input Parameters:
        * key:  Name of QTCM variable.  String.  Scalar.
        * value:  Value of variable.  String or numeric value.  Must be 
          a scalar.
        """
        default_type = Field(key).typecode()
        value_type = num.typecode(value)

        if   value_type in N.typecodes['Float']:
            if default_type not in N.typecodes['Float']:
                raise TypeError, 'value type different from default for key'
            tmp = copy.copy(self.__qtcm.setbypy.getitem_real(key))
            if not self.__qtcm.setbypy.is_readable:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not readable'
            self.__qtcm.setbypy.setitem_real(key, value)

        elif value_type in N.typecodes['Integer']:
            if default_type not in N.typecodes['Integer']:
                raise TypeError, 'value type different from default for key'
            tmp = copy.copy(self.__qtcm.setbypy.getitem_int(key))
            if not self.__qtcm.setbypy.is_readable:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not readable'
            self.__qtcm.setbypy.setitem_int(key, value)

        elif value_type in num.typecodes['Character']:
            if default_type not in num.typecodes['Character']:
                raise TypeError, 'value type different from default for key'
            tmp = self.__qtcm.setbypy.getitem_str(key)
            if not self.__qtcm.setbypy.is_readable:
                raise FieldNotReadableFromCompiledModel, \
                      'Compiled model variable not readable'
            self.__qtcm.setbypy.setitem_str(key, value)

        else:
            raise TypeError, 'value is of unsupported type'




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
