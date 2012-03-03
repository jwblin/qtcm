! File:  outpall.F90
!
! ---------------------------------------------------------------------
!
Subroutine outpAll
  !
  ! Output means and instantaneous values;
  ! called by the driver.
  !
  Use Input, Only : ntouti, ntout, ntoutr, noout, nooutr
  Use Calendar

  Implicit None

  Logical outTime

  If (outTime(ntouti,noout)) Call outpInst ! output instanteneous values 
  If (outTime(ntout ,noout)) Call outpMean ! output mean values
  If (outTime(ntoutr,nooutr)) Call outrestart  ! output restart file
!  If (dayofyear==DAYSPERYEAR) Call outrestart ! End of year restart



  Return
End Subroutine outpAll
!
! ===================================================================
