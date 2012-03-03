! File:  outpinit.F90
!
! ---------------------------------------------------------------------
!
Subroutine outpInit 
  !
  ! NetCDF version of outpInit
  !
  ! For netCDF output the files are open at the first output time.
  ! With isInit='.true.' All that is done by defOutVars is
  ! setting up the (linked) list of output arrays to store the
  ! accumulated values of the variables for mean value output.
  ! 
  Use ncFileData, Only: isinit
  !
  isinit=.True.
  ! defOutVars calls defVar which for isinit=.true.
  ! only sets the variable var0, which indicates the fields to average
  Call defOutVars 
  isinit=.False.
  Return
End Subroutine outpInit
!
! ===================================================================
