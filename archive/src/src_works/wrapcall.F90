!===========================================================================
! Name:
!   WrapCall - Module of wraps of subroutine calls for Python access
!
! Description:
!   This module contains subroutines whose sole purpose is to call other
!   single subroutines in the QTCM model.  These wrapper routines are
!   needed because f2py, for some reason I can't figure out, will not
!   properly wrap Fortran routines (that are then callable at the Python
!   level) that create local arrays using parameters obtained through a
!   use statment.  Thus, a routine foo with the following definition
!
!       subroutine foo
!         use dimensions
!         real a(nx,ny)
!         [...]
!       end subroutine foo
!
!   where nx and ny are defined in dimensions, will return an error,
!   either resulting in the extension module being unable to be created,
!   or in results that are different from running the pure-Fortran
!   version of the model.
!
! Modification History:
! * 28 May 2008:  Johnny Lin, Physics Department, North Park University.
!
! Copyright (c) 2008 by Johnny Lin.  For licensing, distribution condi-
! tions, contact information, and additional documentation see the URL
! http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
!===========================================================================

Module WrapCall
  Implicit None

Contains
  Subroutine wadvcttq
    Call advcttq
  End Subroutine wadvcttq

  Subroutine wadvctuv
    Call advctuv
  End Subroutine wadvctuv

  Subroutine wbarcl
    Call barcl
  End Subroutine wbarcl
            
  Subroutine wbartr
    Call bartr
  End Subroutine wbartr

  Subroutine wbndinit
    Call bndinit
  End Subroutine wbndinit

  Subroutine wcloud
    Call cloud
  End Subroutine wcloud

  Subroutine wdffus
    Call dffus
  End Subroutine wdffus

  Subroutine wgradphis
    Call gradphis
  End Subroutine wgradphis

  Subroutine wmconvct
    Call mconvct
  End Subroutine wmconvct

  Subroutine wocean(ndays,it)
    implicit none
    integer, intent(in) :: ndays,it
    Call ocean(ndays,it)
  End Subroutine wocean

  Subroutine woceaninit
    Call oceaninit
  End Subroutine woceaninit

  Subroutine wradlw
    Call radlw
  End Subroutine wradlw

  Subroutine wradsw
    Call radsw
  End Subroutine wradsw

  Subroutine wsavebartr
    Call savebartr
  End Subroutine wsavebartr

  Subroutine wsflux
    Call sflux
  End Subroutine wsflux

End Module WrapCall
