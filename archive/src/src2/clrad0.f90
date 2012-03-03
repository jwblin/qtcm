!*********************************************************************
! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/clrad0.f90,v $
! $Author: munnich $
! $Revision: 1.6 $
! $Date: 2002/07/18 20:34:36 $
!
! File:  clrad0.f
!
! This contains a simple cloud and radiation package. Retaining much
!   accuracy, its "analytical" ability makes it easy to understand what's
!   going on.
!
! 1. cloud:     Interactive (High+Mid) clouds depend on precipitation through
!               empirical relationship; low+Nimbo clouds fixed to a reference
!               value; radiation are treated correspondingly
!
! 2. solartop:  Diurnally averaged ; easy modification for orbital
!               parameters
!
! 3. Shortwave: Single cloud-atmosphere layer, two-stream, analytical
!
! 4. Longwave:  Linear longwave radiation scheme
!
! Code history:
!   Original version (QTCM1 V1a):  Zeng/Chou/Neelin  Jan 1997
!   Beta version (QTCM1 V1b):                        Jun 1997
!   QTCM1 release (QTCM1 V2.0):                      Jun 1998
!   QTCM1 V2.1: new b1 profile/reclassification
!     of interactive/noninteractive cloud types      Dec 1998
!*********************************************************************
Module Clrad
  Use Dimensions
  Real, Dimension(nx,ny) :: cl1
  ! Simple radiation
  Real :: cllowref=0.36,cl1ref=0.21
End Module Clrad

Subroutine cloud
  !************
  ! cloud prediction; deep convective cloud cover proportional to convective
  !   precipitation 
  ! 
  ! We lump clouds into two categories: interactive and noninteractive, with
  ! the former proportional to rainfall (predicted by QTCM), the latter
  ! fixed at its tropical mean.
  ! In the following, we lump high clouds (types 5,6,7 of ISCCP) and AsAc (3) as
  !  a single interactive type (called High+Middle), while 1,2,4 are lumped
  !  as one noninteractive cloud type (loosely called low clouds);
  !  the reason we put Nimbostratus as noninteractive is because it has
  !  little correlation with rainfall.
  Use Dimensions
  Use Clrad, Only : cl1
  Use SurfaceFluxes, Only: Qc

  Implicit None

  ! --Altered Common Variables:
  ! cl1          !cloud cover          fraction

  ! --Local Variables & Functions:
  Integer i,j
  Real :: cl1P=0.0018

  !*****************************************************************

  Do j=1,ny
     Do i=1,nx
        !proportional to deep convective prec.
        cl1(i,j)=Min(cl1P*Qc(i,j),1.0)   !cloud cover not allowed to exceed 100%
        !        cl1(i,j)=cl1ref                !@ mean cloud
     End Do
  End Do

  Return
End Subroutine cloud


Subroutine solartop
  !************
  ! Solar radiation at TOA, diurnally averaged
  !
  ! Astronomical formulae used:
  ! 1.  solar declination angle as a function of dayofyear     Zeng(1994)
  !     derived by assuming a circular orbit (tiny error)
  ! 2.  correction(<3.4%) to solar constant due to the 
  !     ellipticity of the earth orbit, function of dayofyear  Dickinson(1983)
  ! 3.  Diurnally averaged solar zenith angle             Sellers(1965),Zeng
  !
  ! Length of year =DAYSPERYEAR to accomodate qtcm coupling
  !
  Use Dimensions
  Use Fluxes, Only: S0
  Use Grid, Only: pi, dyi
  Use PhysicalConstants, Only: Rearth
  Use Calendar, Only: dayofyear, rdaysperyear
  Implicit None

  ! --Altered Common Variables:
  ! S0

  ! --Local Variables & Functions:
  Real sint,cost,dtheta
  Real cosZbar
  Real cosHor,sinHor,Hor
  Real sint0,cost0,theta0,theta0max
  Real SolarC
  Real dayspring,dayperihelion,ecce
  Integer i,j
  !*****************************************************************

  !- orbital parameters
  theta0max=23.447                 !degree; max solar dec. angle
  dayspring = 81.                  !Julian day of Vernal Equinox
  dayperihelion = 3.               !day of year closest to the sun
  ecce = 0.034                     !related to eccentricity of earth orbit

  sint0  = Sin(theta0max/180.*pi) *       &          !Zeng (1994)
       &         Sin(2.*pi*(dayofyear-dayspring)/rdaysperyear)
  theta0 = Asin(sint0)                    !solar declination angle
  cost0  = Cos(theta0)
  SolarC = 1370.*(1.+ecce* &
       &                Cos(2.*pi*(dayofyear-dayperihelion)/rdaysperyear))
  dtheta=Rearth*dyi

  Do j=1,ny
     sint=Sin((j-(ny+1)*0.5)*dtheta)               ! sint=sin(lat)
     cost=Sqrt(1.0-sint**2)
     cosHor=-sint*sint0/(cost*cost0)
     ! winter and summer polar regions; not needed for a tropical model
     cosHor=Max(cosHor,-1.0)       !polar; no night; halfday length=pi
     cosHor=Min(cosHor,1.0)        !polar; no day; halfday length=0
     Hor=Acos(cosHor)
     sinHor=Sqrt(1.0-cosHor**2)
     cosZbar=(sint*sint0*Hor+cost*cost0*sinHor)/pi   !Sellers(1966)
     S0(1,j)=SolarC*cosZbar
     Do i=2,nx                             !no dependence on i
        S0(i,j)=S0(1,j)
     End Do
  End Do

  Return
End Subroutine solartop


Subroutine radsw
  !************
  ! A rudimentary shortwave radiation scheme; Zeng/Neelin (J.Climate 1998)
  ! Analytically solved surface and top fluxes for a single layer of 
  !   cloud-atmo with reflectivity 'albdc' (including cloud
  !   albedo and atmospheric backscattering) and absorptivity 'absp', over a 
  !   surface with albedo 'ALBDs'. The cloud albedo dependence on cloud cover
  !   is computed with solar zenith angle effects absorbed in.
  !
  Use Dimensions
  Use Fluxes, Only: FSWds, FSWus, FSWut, FSW, S0
  Use Surface, Only: ALBDs
  Use Clrad
  Implicit None

  ! --Altered Common Variables:
  ! FSWs         !net (absorbed by ground) surface solar         W/m**2
  ! FSW          !total SW absorbed by column   W/m**2

  ! --Local Variables & Functions:
  Integer i,j
  Real absp0,abspc1,abspc3,absp,albdc,albdc0,albdc1,albdc3
  Real ALBDst,frac
  !  Real cllowref,cl1ref
  !  Common/clcom/ cllowref,cl1ref       !@can we use total cloud for LW ?
  !
  Data absp0/0.19/     !atmo. absorptivity;   dimensionless
  Data abspc1/0.00/    !low cloud absorptivity;   dimensionless
  Data abspc3/0.08/    !H+M cloud absorptivity; <- 4% absorption at cl=.5
  Data albdc0/0.00/    !clear sky back scattering
  Data albdc1/0.30/    !converting low cloud cover into albedo  1/100%
  Data albdc3/0.40/    !converting H+M cloud cover into albedo  1/100%
  Data frac/0.5/       !fraction absorbed before reflection occurs at top
  !*****************************************************************

  Call solartop        !get incoming solar at TOA

  Do j=1,ny
     Do i=1,nx

        albdc=albdc0+albdc1*cllowref+albdc3*cl1(i,j)     !cloud/atmo albedo

        absp=absp0+abspc1*cllowref+abspc3*cl1(i,j)       !total absorptivity

        ALBDst=ALBDs(i,j)

        FSWds(i,j)=S0(i,j)*(1.0-albdc)*(1.0-absp)      !downward surface
        FSWus(i,j)=FSWds(i,j)*ALBDst                   !upward surface
        FSWut(i,j)=S0(i,j)*( (1.0-frac*absp)*albdc  &   !upward top
             &       +((1.0-albdc)*(1.0-absp))**2*ALBDst )

        FSW(i,j)=S0(i,j)+FSWus(i,j)-FSWds(i,j)-FSWut(i,j)

     End Do
  End Do

  Return
End Subroutine radsw



Subroutine radlw
  !************
  ! Longwave radiation; linear scheme (LINRAD) of Chou and Neelin (1996)
  ! linearization coefficients calculated by NZ using a tropical mean profile
  ! 
  Use Dimensions
  Use Clrad
  Use Fluxes, Only: FLWds, FLWus, FLWut, FLW
  Use Surface, Only: Ts
  Use Qtcmpar, Only: Tsref
  Use Prognostic, Only: T1, q1
  Implicit None

  ! --Altered Common Variables:
  ! FLWds         !downward longwave flux at surface         W/m**2
  ! FLWus         !upward longwave flux at surface           W/m**2
  ! FLWut         !upward longwave flux at top, i.e., OLR    W/m**2
  ! FLW           !FLW=FLWus-FLWds-FLWut; total LW absorbed by column   W/m**2

  ! --Local Variables & Functions:
  Integer i,j
  Real FLWdsref,FLWusref,FLWutref &
       & ,eps_rT1t,eps_rT1s,eps_rq1t,eps_rq1s,eps_rclt,eps_rcls &
       & ,eps_rTst,eps_rTss

  !
  !  Real cllowref,cl1ref
  !  Common/clcom/ cllowref,cl1ref       !@can we use total cloud for LW ?

  ! Values derived from CC(a1,b1 to 200mb, Chia's weakly nonlinear coefficients
  !   are used by weighted averaging high+middle clouds)
  Data FLWdsref/431./       !reference state, LW downward surface
  Data FLWusref/475./       !reference state, upward surface
  Data FLWutref/253./       !reference state, upward top  (OLR)
  Data eps_rT1t/1.23/       !top upward IR due to T1, W/m**2/K  -> -42 day
  Data eps_rT1s/1.29/       !surface downward IR due to T1, W/m**2/K
  Data eps_rq1t/-0.41/      !top upward IR due to q1(in Kelvin) ->75 day
  Data eps_rq1s/1.23/       !surface downward IR due to q1(in Kelvin)
  Data eps_rclt/-81.1/      !top upward IR due to cloud,W/m2/100%
  Data eps_rcls/15.9/       !surface downward IR due to cloud,W/m2/100%
  Data eps_rTst/0.533/      !top upward due to Ts    W/m^2/K
  Data eps_rTss/6.283/      !surface upward due to Ts W/m^2/K -> 17 day 
  Save FLWdsref,FLWusref,FLWutref &
       & ,eps_rT1t,eps_rT1s,eps_rq1t,eps_rq1s,eps_rclt,eps_rcls &
       & ,eps_rTst,eps_rTss
  !*****************************************************************

  Do j=1,ny
     Do i=1,nx

        FLWds(i,j)=FLWdsref+eps_rT1s*T1(i,j)+eps_rq1s*q1(i,j) &
             &       +eps_rcls*(cl1(i,j)-cl1ref)
        FLWus(i,j)=FLWusref+eps_rTss*(Ts(i,j)-Tsref)
        FLWut(i,j)=FLWutref+eps_rT1t*T1(i,j)+eps_rq1t*q1(i,j) &
             &       +eps_rclt*(cl1(i,j)-cl1ref)&
             &       +eps_rTst*(Ts(i,j)-Tsref)

        FLW(i,j)=FLWus(i,j)-FLWds(i,j)-FLWut(i,j)

     End Do
  End Do

  Return
End Subroutine radlw
