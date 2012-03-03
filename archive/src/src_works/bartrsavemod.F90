! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: bartrsavemod.F90
!======================================================================
!
!
Module BartrSave
!
! Saved barotropic velocity 
! Needed for geopotential gradient calculation
!
  Use Dimensions
  Implicit None

  Real, Dimension(nx,ny) :: u0sav, v0sav
End Module BartrSave


