! calendar and time variables
Module Calendar
#ifdef YEAR360
  !
  !  360 day calendar no leap year
  Integer, Dimension(12), Parameter ::        &
       &       monlen(12)   = (/30,30,30,30,30,30   & ! length of month
       &                       ,30,30,30,30,30,30/),           & ! length of month
       &       cummonth(12) = (/0,30,60,90,120,150,& ! cum. days of prior months
       &                        180,210,240,270,300,330/)
  Integer, Parameter :: daysperyear=360
  Real,    Parameter :: rdaysperyear=360.
#else
  !
  ! Default: 365 day calendar no leap year
  Integer, Dimension(12), Parameter ::         &
       &       monlen(12)   = (/31,28,31,30,31,30   & ! length of month
       &                       ,31,31,30,31,30,31/) &
       &     , cummonth(12) = (/0,31,59,90,120,151  & ! cum. days of prior months
       &                       ,181,212,243,273,304,334/)

  Integer, Dimension(0:13), Parameter ::         &
       &    midmonth=(/-16, 15,  45, 74, 105, 135, 166, & ! integer for TimeInterp
       &                          196, 227, 258,288, 319, 349, 380/)
  !  real, dimension(12), parameter :: &   ! julian day of monthly means
       !       &    midmonth=/( 15.5,  45, 74.5, 105, 135.5, 166, &
  !       &                          196.5, 227.5, 258, 288.5, 319, 349.5 /), 
  Integer, Parameter :: daysperyear=365    
  Real,    Parameter :: rdaysperyear=365.
#endif

  ! #if def DIURNAL
  !         integer*8 itsteps
  !         integer timeofday
  ! #endif

  Integer dayofmodel, dateofmodel, yearofmodel
  Integer dayofmonth, dayofyear, monthofyear
  Integer year0, month0, day0
End Module Calendar
!
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
! ---------------------------------------------------------------------
!
! #if def DIURNAL
! Subroutine set_qclock(it)
!
! Advance the internal clock of QTCM, using the model step as a reference.
! The clock is intended for components of the QTCM which have diurnal cycles.
!
! timeofday is set to the current time of day, in seconds.
!
! Note: it is assumed that it is 00Z when this is first called, otherwise
!   an offset should be  used here.  NZ
! 
!   Use Input,    Only : dt
!   Use Calendar, Only : timeofday
! 
!   Implicit None
!   Integer it
! 
!   timeofday = Mod(Int(dt)*(it-1),86400)
!   Return
! End Subroutine set_qclock
! #endif
!
! ---------------------------------------------------------------------
!
Integer Function julian(date)
  !
  !  return the "date" as the number of days since the model reference year.
  !
  Use Calendar
  Implicit None
  Integer date, day, month, year
  year   = date/10000 - year0
  month  = Mod(date,10000)/100 
  day    = Mod(date,100)-day0+1
  julian = day+cummonth(month)-cummonth(month0)+year*DAYSPERYEAR
  Return
End Function julian
!
! ===================================================================== 
