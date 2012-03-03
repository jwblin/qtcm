! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: physics1.F90
!
! Code history:
!   Original (QTCM1 Version 1a):  Ning Zeng, David Neelin, et al. Jan. 1997
!   Beta version (QTCM1 V1b):     June 1997
!   QTCM1 release (QTCMes1 V2.0):   June 1998
!   QTCM1 V2.1 :                  Jan  1999
!   QTCM1 V2.2f90 :              Nov 2000 
!
! Some routines possibly called but not in this file:
!   getSST          in getbnd.F90
!   varmean         in output.F90
!   sland           in sland.F90
!   radsw \
!   radlw  >        in clrad.F90
!   cloud /
!   sflux/mlsflux   in sflux.F90/mlsflux.F90
! and if chosen by a preprocessor option      
!   cplmean         in cplmean.F90
!   set_qclock      in driver_util.F90
!
! File I/O:  There are files remain open throughout (global), and others
!   are closed immediately in the corresponding opening routine before
!   exit (local).
!   Convention: units 20-29 are reserved for global files (always output);
!   units 10-19 for local files (usually input except for qtcm.restart.out);
!   no other units are used.
!
! ---------------------------------------------------------------------
!
Subroutine physics1
  !
  ! Model physics package associated with temperature mode 1:
  !   moist convection, cloud, longwave and shortwave radiation,
  !   surface fluxes of heat,moisture and momentum
  ! 
  Implicit None

  ! print*,'TOP physics1'
  Call mconvct         !moist convection; precipitation
  ! print*,'mconvct finished'
  Call cloud           !cloud prediction
  Call radsw           !shortwave scheme
  Call radlw           !longwave radiation scheme
  Call Sflux           !surface fluxes
  ! call other_flux      !e.g., stochastic forcing; add to rhs

  Return
End Subroutine physics1


!======================================================================
!            End of physics1.F90 file
!======================================================================
