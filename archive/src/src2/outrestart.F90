! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: outrestart.F90
!
! ---------------------------------------------------------------------
!
Subroutine outrestart
  !
  ! Save model variables for restart. 
  ! Be sure to save ALL the prognostic variables for atmo and land
  !
  Use Prognostic
  Use Barotropic, Only : vort0,psi0,rhsvort0,u0bar,rhsu0bar
  Use Surface,    Only : Ts
  Use Land,       Only : WD
  Use Input,      Only : title
  Use Calendar,   Only : dayofmodel, dateofmodel

  Implicit None

  ! --File I/O:  Subroutine opens and closes following files:
  !     unit 12 = qtcm.restart.out
  !   units closed prior to entry, and on exit of the routine
  !
  ! --Local Variables & Functions:
  Character(len=30)  :: fname
  Integer,Save       :: icount=1
  !
  write(fname,'(a,i8.8,a)') 'qtcm_',dateofmodel,'.restart'
#ifdef   ENSEMB_INI
  write(fname,"('qtcm.ini_',i2.2)") icount
#endif
  Open(12,file=fname,form='unformatted')
  Write(12)dateofmodel,title
  Write(12)u1(:,1:ny),v1,T1(:,1:ny),q1(:,1:ny)    & ! baroclinic mode
       &      ,u0(:,1:ny),v0,vort0,psi0,rhsvort0  &
       &      ,u0bar,rhsu0bar                     & ! barotropic mode
       !    2      ,u2,v2                           ! ekman
       &      ,Ts,WD                                     ! land
  Close(12)
  Write(6,'(a,i8.8)') ' Restart file written at end of ', dateofmodel
  icount=icount+1

  Return
End Subroutine outrestart
!
!======================================================================
