! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: qtcminit.F90
!
! ---------------------------------------------------------------------
!
Subroutine qtcminit()
  !
  ! Initialize qtcm atmospheric model variables, datasets, and
  ! diagnostic fields.
  !
  Implicit None

  !(begin comment out)
  !   In Version 2.3 - current-release (August 2002):
  !   This is the original order of the routine.  However, the problem
  !   with this order is that STYPE is not read in for landon until
  !   bndinit, but in varinit you use STYPE to calculate the original
  !   values of WD for the non-restart case.  The manual is also gives
  !   two conflicting (p. 29 and 32) orders.
  !Call parinit            !Initialize model parameters
  !Call varinit            !Initialize variables
  !Call TimeManager(1)     !mm set model time 
  !Call bndinit            !input boundary datasets
  !Call physics1           !diagnostic fields for initial condition
  !(end comment out)

  ! JWL:  I've altered the above order to this order in order for
  ! varinit to work correctly.  bndinit doesn't seem to have any
  ! dependencies such that it has to come after TimeManager or
  ! varinit.
  Call parinit            !Initialize model parameters
  Call bndinit            !input boundary datasets
  Call varinit            !Initialize variables
  Call TimeManager(1)     !mm set model time 
  Call physics1           !diagnostic fields for initial condition

  Return
End Subroutine qtcminit
!
!======================================================================
