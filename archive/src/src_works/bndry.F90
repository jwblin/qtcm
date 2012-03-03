! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/bndry.f90,v $
! $Author: munnich $
! $Revision: 1.5 $Date: 2002/07/18 20:34:35 $
!
! Declarations for bndry.f
!
!-- Model variables
Module bndry
  Use dimensions
  Integer, Parameter :: NVAR = 7
  Real var_now(nx,ny),var_prior(nx,ny,nvar),var_next(nx,ny,nvar)
  Integer time1, time2
End Module bndry
!
!*********************************************************************
! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/bndry.f90,v $
! $Author: munnich $
! $Revision: 1.5 $
! $Date: 2002/07/18 20:34:35 $
!
! File: bndry.f
!
! Generic routines for reading multi-variable bndary data files 
!   (currently read monthly climatology only); 
! Two different types of bndary data are handled currently:
!  bndry1: type 1 data, assumed centered on the 15th of each month, 
!     linear interpolation in between; example: monthly SST
!  bndry2: type 2 data from 15th to the 15th (e.g.,00000101.var covers 
!     from dec15 to jan15th); no interpolation. example: SST tendency
!
! The date information is converted into file names:
!   In this version, bc file names include the date and variable name
!   in this way:
!               file name = yyyymmdd.var
!   where yyyy = year , mm = month, dd=dayofmonth, var = a variable name.
!   if seasonal, yyyy=0000
!   dd=15: type 1 data; dd=01: type 2 data 
!
! Code history:
!   Original version:  Ning Zeng/Bill Weibel  May 1998, adopted from ocean.f
!   Modified to interpolate Aug 1998 by K. Hales
!
! Routines called but not in this file
!   julian             in driver_util.F
!   TimeInterp         in utilities.F
!*********************************************************************
!
!
Subroutine bndry1(var, vardir, vardummy)
  !*******************************************************************
  ! Time evolving boundary conditions other than real-time SST:
  !   prescribed fields such as observed soil moisture, surface
  !   radiative fluxes, cloud cover, etc.
  !******************************************************************
  !
  Use dimensions
  Use bndry
  Use Calendar

  Implicit None

  Real, Dimension(nx,ny),Intent(out) :: vardummy(nx,ny)
  Character(len=3), Intent(in)       :: var
  Character(len=*), Intent(in)       ::  vardir
  !
  !- Altered
  !    vardummy
  !
  !- Local Variables and Functions
  Integer julian
  Integer orig_time1,orig_time2,prior_dofm
  Integer month
  Character(len=3), Dimension(nvar) :: var0
  Integer i,j,var_index,m,num_var
  Integer day1,day2
  Real thistime,ndays
  Logical firstcall(NVAR)
  !
  Data var0/NVAR*'   '/        !array of var names
  Data firstcall/NVAR*.True./
  Data prior_dofm/0/
  Data num_var/0/              !variable number corres. to
  !position in var0
  Data var_index/1/
  !
  Save prior_dofm,orig_time1,orig_time2
  Save firstcall
  Save var_index,num_var
  !
  !
  !
  ! *************************************************************
  ! Assign an index for each variable and initialze the variables
  ! *************************************************************
  !
  ! Assign a value of 'var_index' corresponding to each var called 
  ! (like num_var above)
  ! 'var_index' is what we called 'k'
  Do m=1,NVAR
     If (var.Eq.var0(m)) Then
        var_index=m
        firstcall(var_index)=.False.
     End If
  End Do
  !
  ! this loads up the array var0 with variable names, only first
  ! time per variable; 
  If(firstcall(var_index)) Then
     num_var=num_var+1
     If (num_var.Gt.NVAR) Then
        Write(*,*) 'bndry: Error. capacity of array "var0" exceeded.'
        Write(*,*) 'Remedy: increase size of NVAR in bndry'
        Stop
     End If
     var_index=num_var
     var0(var_index)=var
     Call bndry1_init(var,vardir,var_index) !initialize the var
  End If
  !
  !
  ! **************************************************************
  ! Get values for var_prior and var_next
  ! **************************************************************
  !     Basically the structure is the same as 
  !     used for reading and interpolating sst
  !
  ! Reset time1 and time2 if bndry is called more than once 
  ! in a month (ie used to get more than one variable)
  If (dateofmodel.Eq.prior_dofm) Then
     time1=orig_time1
     time2=orig_time2
  End If
  !
  ! Read in the variable of the next month 
  ! This is altered from subroutine ocean (time1,time2 equivalent to
  ! date1,date2):
  If (dateofmodel .Ge. time2) Then
     Do j=1,ny
        Do i=1,nx
           var_prior(i,j,var_index) = var_next(i,j,var_index)
        End Do
     End Do
     orig_time1=time1        !keep the unadvanced time1 and
     orig_time2=time2        !time2 for the next var called
     time1 = time2
     time2 = time1 + 100                     !advance time2 by 1 month
     month=Mod(time1-15,10000)/100           !if month=12
     If(month.Eq.12)time2= time1+10000-1100   !  set to Jan of next year
     Call readvar(time2,var,vardir, var_next(1,1,var_index))
     prior_dofm=dateofmodel         !set prior_dateofmodel
  Endif
  !
  !
  !
  ! *************************************************************
  ! Interpolate Variable
  ! *************************************************************
  ! do interpolation
  ! below is modified from subroutine Ocean
  ! ndays (=interval) needs to be passed!
  !
  !     timeinterp requires julian dates or some other monotonic time.
  day1 = julian(time1)
  day2 = julian(time2)
  !     !!!!! ndays needs to be passed somehow
  ndays=1.
  thistime = dayofmodel + ndays/2.       ! compute timemean of sst
  Call TimeInterp(var_prior(1,1,var_index),var_next(1,1,var_index)  &
       &     ,vardummy(1,1) ,day1,day2,thistime,nx,ny)
  !
  !
  !
  Return
End Subroutine bndry1
!
!
!
! ***********************************************************
! Subroutine readvar: reads the variable from monthly
! climatological data files
! ****************************************************************
!
Subroutine readvar(date,var_name,vardir,vardummy)
  !
  Use dimensions
  Use bndry
  Use Calendar

  Implicit None

  Real, Dimension(nx,ny), Intent(Out) ::  vardummy(nx,ny)
  Character(len=3), Intent(In)        :: var_name
  Character(len=*), Intent(In)        :: vardir
  Integer ios
  !
  !- Local Variables and Functions
  Integer            :: lnblnk, month, date
  Character(len=2)   :: monthc
  Character(len=12)  :: varfile
  Character(len=130) :: fname
  !
  !
  ! assign the file name 
  month=Mod(date,10000)/100
  Write(monthc,'(i2.2)') month
  varfile='0000'//monthc//'15.'//var_name    !file format: 0000mm15.var
  !
  ! read the var
  fname=vardir(1:lnblnk(vardir))//'/'//varfile
  Open(13,file=fname,form='formatted',status='old',iostat=ios)
  If(ios.Ne.0) Then
     Write(6,*) 'readvar: ERROR opening data file: ', Trim(fname)
     Write(6,*) 'Readvar inputs: '
     Write(6,*) '   date    = ',date
     Write(6,*) '   varname = ',var_name
     Write(6,*) '   vardir  = ',vardir
     Stop
  Endif
  !
  Read(13,*) vardummy
  Close(13)
  !
  !
  Write (*,*) 'bndry: Read ',var_name,' at model date', date        &
       &   ,' from file'
  Write(*,*) '         "',vardir(1:lnblnk(vardir))//'/'//varfile,'"'
  !
  !
  Return
End Subroutine readvar
!
!*********************************************************************
!
!*********************************************************************
!
!
Subroutine bndry1_init(var,vardir,var_index)
  !
  !
  !  Get the initial boundary data for the dates on either side of an
  !  arbitrary starting date
  !
  ! this is from subroutine oceaninit in mxlayer.F modified for bndry
  !
  Use dimensions
  Use bndry
  Use Calendar

  Implicit None

  Integer, Intent(In)           ::  var_index
  Character(Len=*), Intent(In)  :: vardir
  Character(Len=3), Intent(In) :: var

  Integer :: i,j
  !
  Write (*,*) "BndryInit: get the initial boundary data"
  !
  !-  get the dates(centered on the 15th) on either side of the starting date.
  ! before the 15th
  If(day0.Ge.1.And.day0.Le.15) Then 
     time1 = year0*10000 + (month0-1)*100 + 15
     If(month0.Eq.1)time1= (year0-1)*10000 + 12*100 + 15
     time2= year0*10000 + month0*100 + 15
  End If
  ! after the 15th
  If(day0.Ge.15.And.day0.Le.31) Then
     time1= year0*10000 + month0*100 + 15
     time2 = year0*10000 + (month0+1)*100 + 15
     If(month0.Eq.12)time2= (year0+1)*10000 + 1*100 + 15
  End If
  !
  !
  Call readvar(time1,var,vardir,var_prior(1,1,var_index))
  !
  Do j=1,ny
     Do i=1,nx
        var_now(i,j) = var_prior(i,j,var_index)      !Fudge the initial 
     End Do
  End Do
  !
  Call readvar(time2,var,vardir,var_next(1,1,var_index))
  !
  Return
End Subroutine bndry1_init
!
!*********************************************************************
!
Subroutine bndry2(var, vardir, vardummy)
  !************
  ! read in seasonal boundary data centered on the 1st of a month.
  ! NO time interpolation needed 
  ! for q-flux (temperature tendency term) calculation  
  !*********************************
  !
  Use Dimensions
  Use Calendar

  Implicit None

  !
  Real, Dimension(nx,ny) , Intent(Out) :: vardummy
  Character(Len=3), Intent(In) :: var
  Character(Len=*), Intent(In) :: vardir
  !
  !- Altered
  !    vardummy
  !
  !- Local Variables and Functions
  Integer lnblnk, month, day, month_prior
  Integer month_read
  Character(Len=2)   :: monthc
  Character(Len=12)  :: varfile
  Character(Len=130) :: fname
  !
  Data month_prior/0/
  Save month_prior
  !
  month=monthofyear  ! mod(dateofmodel,10000)/100
  day=dayofmonth     ! mod(dateofmodel,10000)-month*100
  month_read=month
  If(day.Ge.15) month_read=month+1  
  If(month_read.Eq.13) month_read=1
  If(month_read.Eq.month_prior) Return
  !
  Write(monthc,'(i2.2)') month_read              !convert into a string
  varfile='0000'//monthc//'01'//'.'//var    !file format: 0000mm01.var
  !
  ! read the var
  fname=vardir(1:lnblnk(vardir))//'/'//varfile
  Open(13,file=fname,form='formatted',status='old')
  Read(13,*) vardummy
  Close(13)
  !
  Write (*,*) 'bndry: Read ',var,' at model date', dateofmodel      &
       &   ,' from file'
  Write(*,*) '         "',Trim(fname),'"'
  !
  month_prior=month_read
  !
  Return
End Subroutine bndry2
!*********************************************************************
! ********************************************************************
! end bndry.F
! ********************************************************************
