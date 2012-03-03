! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: qtcm.F90
!
! This file contains the subroutine qtcm.
!
! ---------------------------------------------------------------------
!
Subroutine qtcm(interval)
  !
  ! Main routine for the Quasi-Equilibrium Tropical Circulation Model.
  !
  Use Input, Only : dt, mt0, nastep, it  , noout
  Use Calendar, Only : dayofmodel  ! for outpInst
  Use QtcmUtilMod, Only : barcl, bartr, advctTq, dffus, &
      xfilter, gradphis, SaveBartr, advctuv

  Implicit None

  ! --Arguments:
  Integer interval         ! interval of ocean-atmos.
  ! coupling, unaltered on exit

  ! --Local Variables
  !  Integer it
  !----------------------------------

  Call getbnd          !get time evolving bndry data, SST, albedo etc.

  !- number of atmospheric time steps within one air-sea coupling interval
  nastep = interval* (86400./dt)

  If(Mod(nastep,mt0).Ne.0) Then
     Print*,'qtcm: Fatal error: mt0*dt not a divisor of 1 day'
     Stop
  Endif

  !- land-atmosphere coupling loop
  !    land runs at the same timestep as atmo., easy modification otherwise,
  !    but if land runs less often than [sflux] in [physics1], Evap over land 
  !    needs to be fixed in [sflux]
  Do it=1,nastep
!!$#if def DIURNAL
!!$     Call set_qclock(it)   !advance the diurnal clock
!!$#endif
     Call physics1         !physics associated with temperature mode T1
     Call sland1           !land

     Call advctuv       !advection of momentum u1,v1,u0,v0
     Call advctTq       !advection of T1, q1
     Call dffus         !diffusion of u1,v1,u0,v0,T1,q1

     Call barcl         !the baroclinic mode
     If(Mod(it,mt0).Eq.0) then
#ifdef NO_ABL
        Call bartr      !barotropic mode
#else
        Call SaveBartr !save u0,v0 for surface geopotential
        Call bartr      !barotropic mode
        Call gradphis  !surface geopotential gradient
#endif
     End If
     !If(Mod(it,3)==0) Call varmean          !time mean (sum up variables, output in output.F)
     Call varmean          !time mean (sum up variables, output in output.F)
#ifdef MXL_OCEAN
     ! ocean mixed layer runs need mean sfc fluxes
#define CPLMEAN
#endif

#ifdef BLEND_SST
#define CPLMEAN
#endif

#ifdef CPLMEAN
     Call cplmean(it)      !optional mean sfc. flux for air-sea coupling
#endif

     !
     ! output within the air-sea coupling interval, useful for model diagnosis.
     !   you may want to set the user-control parameter ntouti=1 (day)
     ! if(dayofmodel< 1) call outpInst

  End Do
  Return
End Subroutine qtcm
!
!======================================================================
