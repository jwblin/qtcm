!******************************************************************************
! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/driver.f90,v $
! $Author: munnich $
! $Revision: 1.8 $
! $Date: 2002/07/25 00:02:00 $
!
! File: driver.f90
!
! Drives qtcm and ocean, structured like a coupled ocean-atmosphere
!   model
!
! Code history:
!   Original Version:                 Bill Weibel    Nov 1996
!   Partial revision for SST input :  Ning Zeng      Jan 1997
!   QTCM1 release (QTCM1 Version2.0)                 June 1998
!******************************************************************************

Program Driver

  Use Calendar, only : dateofmodel
  Use Input, only : lastday, interval
  Use Prognostic, only : u1
  Implicit none
  Integer   :: day      

  !
  ! Initializing 
  call driveInit               !setup user-control parameters
  call qtcmInit                !atmo
  call oceanInit               !ocean
  call outpInit                !initialize output (open global files)      
  print'(a,i8.8)', 'Driver: Starting at model date ', dateofmodel
  !
  ! Main loop of ocean-atmosphere coupling
  do day = 1, lastday, interval
     if (day .eq. 11) then
        u1(53,13) = u1(53,13) * 2.0
        u1(53,14) = u1(53,14) * 2.0
        u1(53,15) = u1(53,15) * 2.0
        u1(54,13) = u1(54,13) * 2.0
        u1(54,14) = u1(54,14) * 2.0
        u1(54,15) = u1(54,15) * 2.0
        u1(55,13) = u1(55,13) * 2.0
        u1(55,14) = u1(55,14) * 2.0
        u1(55,15) = u1(55,15) * 2.0
     endif
     call timeManager(day)     ! set the calendar
     call ocean(interval,day)  ! ocean; generate SST for atmo
     call qtcm(interval)       ! atmo
     call outpAll              ! output all
     print '(a,i6,a,i8.8)', 'Driver: Running for',day, &
         &             ' days at model date ', dateofmodel

     ! call flush(6) ! flush not known on some machines
  enddo
  !
  ! Write Restart file
  call out_restart

  print*,'QTCM finished normally'
end program Driver
!
! -----------------------------------------------------------------------------
!
Subroutine DriveInit 
  !  Read input parameters from driver.in
  !
  !  Use Dimensions
  Use Input
  Use Calendar
  Use AuxVars

  Implicit None
  Real :: viscT, viscQ, viscU

  Namelist /driverdata/ title,bnddir,SSTdir,outdir ,runname         &
       &  ,landon,SSTmode,year0,month0,day0,lastday,interval        &
       &  ,noout,ntout,ntouti,mrestart,dt,mt0,ntoutr                &
       &  ,viscT, viscQ, viscU                                      &
       &  ,ziml     &  ! mixed layer depth [m] ~ cloud base
       &  ,weml     &  ! mixed layer entrainment velocoty [m/s]
       &  ,VVsmin  &  ! minimum speed for surface winds [m/s]
       &  ,arr1name,arr2name,arr3name,arr4name  &
       &  ,arr5name,arr6name,arr7name,arr8name  &
       &;

  !
  ! default values for input:
  title='QTCM default title'      ! A decrciptive title 
  bnddir='../bnddata'             ! bnd. data other than SST
  SSTdir='../bnddata/SST_Reynolds'! where SST files are
  outdir='../proc/qtcm_output'    ! where output go
  runname='runname'               ! string for an output filename
  landon=1                        ! if not 1: land = ocean w. fake 'SST'
  SSTmode='seasonal'              ! decide what kind of SST to use
  year0 = 0                       ! starting year if <0 use year in restart
  month0 = -1                     ! starting month if <0 use mo in restart
  day0 = -1                       ! starting day if <0 use day in restart
  lastday = daysperyear           ! last day of integration
  interval = 1                    ! coupling interval
  noout=0                         ! no output for the first noout days
  nooutr=0                        ! no restart file for the first nooutr days
  ntout=-30                       ! monthly mean output
  ntouti=0                        ! monthly instantaneous data output
  ntoutr=0                        ! restart file only at end of model run
  mrestart=1                      ! =1: restart using qtcm.restart
  dt=1200.                        ! time step [seconds]
  mt0=1                           ! barotropic timestep every mt0 timesteps
  viscT=12.0e5                    ! temperature diffusion parameter [m^s/s]
  viscQ=12.0e5                    ! humidity diffusion parameter [m^s/s]
  viscU=7.0e5                     ! viscocity parameter [m^s/s]
  ziml=500                        ! Atmos. Mixed layer depth [m]
  weml=0.01                       ! Mixed layer entrainment velocity [m]
  VVsmin=4.5                      ! Minumum wind speed for fluxes [m/s]
  arr1name='?'                    ! Auxiliary output array names 1...8
  arr2name='?'                    !='?' Array not included in output.
  arr3name='?'
  arr4name='?'
  arr5name='?'
  arr6name='?'
  arr7name='?'
  arr8name='?'
  ! Example:
  ! arr1name='dps1dx Mode 1 contribution to dpsdx [m/s^2]'
  !
  ! read user input
  Open(13, file='driver.in', status='old')
  Read (13,driverdata)
  Close(13)
  Write (*,driverdata)

  viscxu0=viscU
  viscyu0=viscU
  visc4x=viscU
  visc4y=viscU
  viscxu1=viscU
  viscyu1=viscU
  viscxT=viscT
  viscyT=viscT
  viscxq=viscQ
  viscyq=viscQ
  Return
End Subroutine DriveInit
!
! =============================================================================

