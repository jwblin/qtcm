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
!   By wrapping these calls into this file, I also avoid having to
!   separate out the QTCM1 subroutines in order to access them at
!   the Python level.  That separation is needed if I'm to use
!   f2py to create wrappers for those subroutines.
!
!   All subroutines in this module begin with "w", with the rest of the
!   name being the Fortran QTCM1 subroutine name.  The calling interface
!   for the "w" version is the same as the Fortran QTCM1 original version.
!   There are no subroutines in this module that do not have an exact
!   counterpart in the Fortran QTCM1 code.
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


  Subroutine wgetbnd
    Call getbnd
  End Subroutine wgetbnd


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


  Subroutine woutpall
    Call outpall
  End Subroutine woutpall


  Subroutine woutpinit
    Call outpinit
  End Subroutine woutpinit


  Subroutine woutrestart
    Call outrestart
  End Subroutine woutrestart


  Subroutine wparinit
    Call parinit
  End Subroutine wparinit


  Subroutine wqtcminit
    Call qtcminit
  End Subroutine wqtcminit


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


  Subroutine wsland1
    Call sland1
  End Subroutine wsland1


  Subroutine wtimemanager(iday)
    implicit none
    integer, intent(in) :: iday
    Call timemanager(iday)
  End Subroutine wtimemanager


  Subroutine wvarinit
    Call varinit
  End Subroutine wvarinit


  Subroutine wvarmean
    Call varmean
  End Subroutine wvarmean

End Module WrapCall
