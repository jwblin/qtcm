! File:  timemanager.F90
! ---------------------------------------------------------------------
!
Subroutine TimeManager(iday)
  !     
  !  timemanager  -   set the model calendar
  !     
  !  A "date" is coded in an integer as yyyymmdd.
  !  inputs:
  !  it  -  integer equal to the current day of the  model integration 
  !
  !  Using a 360 or 365 day calendar 
  !
  Use Input, Only: SStmode
  Use Calendar

  Implicit None

  Integer, Intent(In) :: iday
  Integer             :: day

  dayofmodel = iday
  if(SSTmode=='perpetual') then
     yearofmodel=year0
     dayofyear=Mod((day0-1)+cummonth(month0),DAYSPERYEAR)+1
     monthofyear=month0
     dayofmonth=day0
  else
     day = dayofmodel -1+ cummonth(month0) + day0-1 !day since 00:00ZJan.1 of year0
     yearofmodel = day/DAYSPERYEAR + year0
     dayofyear=Mod((day0-1)+cummonth(month0)+dayofmodel-1,DAYSPERYEAR)+1
     monthofyear=12
     Do While(cummonth(monthofyear)>=dayofyear)
        monthofyear = monthofyear-1
     Enddo
     dayofmonth = dayofyear - cummonth(monthofyear)
  endif

  dateofmodel = yearofmodel*10000 + monthofyear*100 + dayofmonth

  Return
End Subroutine TimeManager
!
! ===================================================================== 
