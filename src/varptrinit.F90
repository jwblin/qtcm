! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: varptrinit.F90
!
! This file contains the variable initialization subroutine.
!
! Code history:
!   Original (QTCM1 Version 1a):  Ning Zeng, David Neelin, et al. Jan. 1997
!   Beta version (QTCM1 V1b):     June 1997
!   QTCM1 release (QTCMes1 V2.0):   June 1998
!   QTCM1 V2.1 :                  Jan  1999
!   QTCM1 V2.2f90 :              Nov 2000, Matthias Munnich
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
Subroutine varptrinit
  !
  ! Associate model pointer variables.
  !

  Use Prognostic
  Implicit None

  ! --Local Variables & Functions:  None
  !
  ! --File I/O:  None
  !
  ! --Altered Module Variables:
  ! all the prognostic pointer variables
  !

  u1  =>u1workA
  v1  =>v1workA
  T1  =>T1workA
  q1  =>q1workA

  Return
End Subroutine varptrinit


!======================================================================
!            End of varptrinit.F90 file
!======================================================================
