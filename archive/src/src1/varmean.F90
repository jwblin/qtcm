! This file the varmean output subroutine for qtcm model
!
! ---------------------------------------------------------------------
!
Subroutine varMean
  !
  ! accumulate (sum up) all variables in the linked list
  ! of output variables
  !
  Use OutputVars
  Use Input,    Only: noout, ntout, it, mt0,nastep
  Use Calendar, Only: dayofmodel

  Implicit None
  Type(outvar), Pointer :: var

  If(dayofmodel <= noout .Or. ntout==0 ) Return
  var=>first
  ! print*,'varMean: first var name = ',var%name
  Do  
  ! loop over list of output variables
     Call oacc(var)
     If(Associated(var%var,last%var)) Exit ! var was last in linked list
     var=> var%next
  End Do
  Return
End Subroutine varMean
!
! ===================================================================

