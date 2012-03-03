! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: varinit.F90
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
Subroutine varinit
  !
  ! Initialize model variables.  Will read data for restart run
  ! if mrestart=1
  !

  Use Prognostic
  Use Barotropic, Only : vort0,psi0,rhsvort0,u0bar,rhsu0bar
  Use Land,       Only : WD, WD0
  Use Surface,    Only : Ts, STYPE
  Use Input,      Only : mrestart
  Use Calendar

  Implicit None

  ! --Local Variables & Functions:
  !
  !
  ! --File I/O:  Subroutine opens and closes following files:
  !     unit 12 = qtcm.restart
  !   units closed prior to entry, and on exit of the routine
  !
  ! --Altered Module Variables:
  ! all the prognostic variables
  !
  ! --Local Variables & Functions:
  Integer :: idateofmodel = 0101
  Integer i,j
  Character(len=130) retitle
  Integer dayr, monthr, yearr

  Call varptrinit       ! associate prognostic pointer variables

  ! initialize (array syntax)
  u1=0.0
  v1=0.0
  T1=-100.0
  q1=-50.0
  u0=0.0
  v0=0.0
  vort0=0.0
  rhsvort0=0.0
  u0bar=0.0
  rhsu0bar=0.0
  Ts=295.0
  WD=0.0
  psi0=0.0    !- JWL add:  seems like it was inadvertently left out

  If(mrestart.Eq.1) Then
     ! Read in initial condition (all prognostic variables/rhs and a few other 
     ! variables for initializing physics) if restart

     Open(12,file='qtcm.restart'                         &
          &         ,form='unformatted',status='old')
     Read(12)idateofmodel,retitle
     Read(12)u1(:,1:ny),v1,T1(:,1:ny),q1(:,1:ny)    & ! baroclinic mode
          &      ,u0(:,1:ny),v0,vort0,psi0,rhsvort0 &
          &      ,u0bar,rhsu0bar                    & !barotropic mode
          !    2      ,u2,v2                                      !ekman
     &      ,Ts,WD                                    !land
     Close(12)
     Write(*,*)'varinit: restart the previous run ',retitle,         &
          &       'at the end of',idateofmodel

     If( (day0<1) .Or. (day0>31) .Or. (month0<1) .Or. (month0>12)     &
          &   .Or. (year0<0) ) Then
        !  calculate start date from restart file
        yearr=idateofmodel/10000
        monthr=Mod(idateofmodel,10000)/100
        dayr=Mod(idateofmodel,100)+1
        If (monlen(monthr) < dayr) Then
           dayr=1
           monthr=monthr+1
           If(monthr.Gt.12) Then
              monthr=1
              yearr=yearr+1
           Endif
        Endif
        ! set the start date
        If (year0<0 )  year0=yearr
        If ( (month0<1) .Or. (month0>12) ) month0=monthr
        If ( (day0<1) .Or. (day0>31) ) day0=dayr
        Write(6,'(a,i4.4,i2.2,i2.2)')                                 &
             &              'Using starting date from restart file:', &
             &               year0,month0,day0
     Endif
  Else
     !- Initialize land points for non-restart run with land on
     Do j=1,ny
        Do i=1,nx
           Ts(i,j)=295.
           WD(i,j)=WD0(Int(STYPE(i,j)))*.7      !70% saturation to start with
           !         everything else = zero
        End Do
     End Do
     day0=1
     If(month0<1 .Or. month0>12) month0=1
     If(year0<1) year0=1
  End If  ! mrestart.eq.1

  Return
End Subroutine varinit


!======================================================================
!            End of varinit.F90 file
!======================================================================
