#!/bin/csh
# QTCM model benchmark run:  aquaplanet, 15 days, daily mean and
#    instantaneous output, netCDF output files
# Produces a restart file starting in November (MONTH0 = 11)
#
# ============================================================

# Set Macro switches for compilation (see makefile)
# -------------------------------------------------
#
set DEFS = " -DNETCDFOUT"

# String to document the macros in effect
set defs_str = `echo "$DEFS" | sed -e 's/ -D/_/g' `
# adjust RUNNAME 
set RUNNAME = `basename $0 .sh`
set MONTH0 = 11

# Set directories 
# ---------------
#
set SCRIPTDIR = `pwd`                    # directory of the script
set QTCMROOT  = $SCRIPTDIR/..            # root directory
set INIDIR    = $QTCMROOT/inidata        # directory of restart file
set SRCDIR    = $QTCMROOT/src            # Set the source directory 
set BNDDIR    = $QTCMROOT/bnddata/r64x42 # Boundary data directory
set OUTDIR    = $QTCMROOT/proc/${RUNNAME} #${defs_str}  # Output dirctory
set LOGFILE = ${RUNNAME}.log             # model output log file
set SCRIPT  = `basename $0`              # script name 

# Let the user know the settings:
#
echo Script Settings:
echo Output directory: $OUTDIR
echo RUNNAME : "$RUNNAME"
echo Preprocessor macro definitions: $DEFS
echo -----------------------
echo ' '


# Compile qtcm in source dir
# --------------------------
#
cd $SRCDIR
# remove old object files if preprocessor flags are set 
if ( "$DEFS" != '' )  make clean 
make DEFS="${DEFS}" qtcm || exit


if ( ! -d $OUTDIR )  mkdir $OUTDIR

# Store script and source code in an archive in the output dir
#
tar -uf $OUTDIR/${RUNNAME}_src.tar *.f90 makefile
cd $SCRIPTDIR
tar -uf $OUTDIR/${RUNNAME}_src.tar $SCRIPT
gzip -f $OUTDIR/${RUNNAME}_src.tar

cd $OUTDIR

# Get restart file and executable
# ----------------------------------------------
# 
#cp $INIDIR/qtcm*restart .
cp -fp $SRCDIR/qtcm .


# Copy the script into the output log file
# ----------------------------------------
#
cp $SCRIPTDIR/$SCRIPT $LOGFILE
chmod -x $OUTDIR/$LOGFILE


# Generate user input file
# ------------------------
#
sed 's/\!.*$//' > driver.in << EOF
 &driverdata
 title='QTCM spinup part 1 (aquaplanet)'
 bnddir='$BNDDIR' 
 SSTdir='$BNDDIR/SST_Reynolds'
 outdir='$OUTDIR'
 runname='$RUNNAME'
 landon=0
 year0=1
 month0 = $MONTH0  
 day0=1
 lastday = 15
 ntout=1
 ntouti=1
 noout=0
  mrestart=0
  viscT=12.0e4
  viscQ=12.0e4
  viscU=7.0e4
  ziml=100.
 !
 !----------------------------------------------------------------
 ! default values for input specified in driver.f90:
 !  title='QTCM default title'          ! A decrciptive title 
 !  bnddir='../bnddata'                 ! bnd. data other than SST
 !  SSTdir='../bnddata/SST_Reynolds'    ! where SST files are
 !  outdir='../proc/qtcm_output'        ! where output go
 !  runname='runname'                   ! string used to create an output filename
 !  landon=1                            ! if not 1: land is like ocean with fake 'SST'
 !  SSTmode='seasonal'                  ! decide what kind of SST to use
 !  year0 = -1                          ! starting year if <0 use year in restart file
 !  month0 = -1                         ! starting month if <0 use month in restart file
 !  day0 = -1                           ! -1 and mrestart=1: Use date saved in the restart file
 !  lastday = daysperyear               ! last day of integration
 !  interval = 1                        ! coupling interval
 !  noout=0                             ! no output for the first noout days
 !  ntout=-30                           ! monthly mean output
 !  ntouti=0                            ! monthly instantaneous data output
 !  ntoutr=-30                          ! monthly restart file
 !  mrestart=1                          ! =1: restart using qtcm_????????.restart
 !  dt=1200.                            ! time step [seconds]
 !  mt0=1                               ! barotropic timestep every mt0 timesteps
 !  viscT=12.0e5                        ! temperature diffusion parameter [m^s/s]
 !  viscQ=12.0e5                        ! humidity diffusion parameter [m^s/s]
 !  viscU=7.0e5                         ! viscocity parameter [m^s/s]
 !  ziml=400                            ! Atmos. Mixed layer depth [m]
 !  weml=0.02                           ! Mixed layer entrainment velocity [m]
 !  Wsmin=4.5                          ! Minumum wind speed for flux calculation [m/s]
 !  V1b=-0.2204077                      ! V1 projection coeffcicient at top of mixed layer
 !  arr1name='?'                        ! Auxiliary output array names 1...8
 !  arr2name='?'                        !='?' Array not included in output.
 !  arr3name='?'
 !  arr4name='?'
 !  arr5name='?'
 !  arr6name='?'
 !  arr7name='?'
 !  arr8name='?'
 ! Example:
 ! arr1name='dps1dx Mode 1 contribution to dpsdx [m/s^2]'
 !----------------------------------------------------------------
 &end
EOF
# Run qtcm using nice
# -------------------
# Timing the run. Output to screen and logfile
#
/usr/bin/time nice ./qtcm |& tee -a $LOGFILE


#
# ---------- end of c-shell script -------------- 

