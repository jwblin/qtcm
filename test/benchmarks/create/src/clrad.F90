!
! This contains the cloud and radiation package of Chou and Neelin (Chou, Ph.D
!   Dissertation, UCLA 1997)
!
! Code history:
!   version 2.0 (QTCM1):  Chou/Neelin  June 1997
!   version 2.3: cloud parameterization modified 
!                by J. D. Neelin, H. Su,  July 2002  
!
!**************************************************************************
Module Clrad
  ! for Chia Chou cloud-radiation package
  Use Dimensions

  Integer, Parameter              :: clt=4, cltm1=clt-1 !cloud types
  Real, Dimension(clt), Parameter :: cldref = &
       &    (/0.1051000,      0.1047000,     0.1096000,     0.2234000/)

  Real                            :: cosZ !solar zenith angle
  Real, Dimension(0:clt,nx,ny)    :: cld  !radiatively active cloud cover
  Real, Dimension(nx,ny)          :: cl1  !cloud 1 output array
#ifdef DELTA_CO2
  Real, Dimension(nx,ny)          :: co2  !carbon dioxide
#endif
End Module Clrad

Subroutine cloud
  !
  ! cloud prediction scheme
  !
  Use Fluxes
  Use Clrad
  Use Input, Only : bnddir
  Use AuxVars
  Implicit None
  !
  !
  !
  ! --Input Variables:
  !   Qc, cldref
  !
  ! --output Variables:
  !   cldref: reference (mean) cloud cover
  !   cld: radiatively active cloud cover after overlap calculation
  !   cldtot: total cloud fraction of each type (i.e. before overlap calculation)
  !                   type 0: clear sky
  !                   type 1: Deep cloud +CsCc
  !                   type 2: Cirrus
  !                   type 3: Stratus
  !                   type 4: AsAc+CuSc (approx constant in space and time)
  !   Note: Cloud cover of cloud type 4 varies only due to cloud overlap
  !
  ! --Local Variables & Functions:
  Integer i,j
  Character(Len=120) :: clddir
  Real  :: cl1P=7.76275869E-4              ! coeff. for cloud prediction scheme
  Real  :: cl2fac=1.5                      ! ratio of cloud type 2 to cloud type 1, 
  Real, Dimension(clt,nx,ny), Save :: cldtot ! Total cloud fraction for each type 

#ifdef COUNTCAP          !count times when cld(1) cap is reached 
  Integer, Dimension(nx,ny)  :: icount = 0
#endif
  !
  !*****************************************************************
  !
#ifdef OBSCLD
  !     read ISCCP observed cloudiness
  clddir=Trim(bnddir)//'/CLOUD_ISCCP' 
  ! note: cldtot(1), cldtot(2) are parameterized below, 
  ! cldtot(4) not from data so currently only  cldtot(3) used 
  ! NB: reads only once per month
  Call readobscloud(cldtot,clddir)
  ! assign preliminary values of cloud 4, to be modified below
  ! by overlapping cloud 1 and 2
  Do j=1,ny
     Do i=1,nx
        cldtot(4,i,j)=min(cldref(4)/(1-cldref(1)-cldref(2)),1-cldtot(3,i,j))
     End Do
  End Do
  
#else
  ! assign preliminary values of cloud 3 and 4, to be modified below
  ! by overlapping cloud 1 and 2 
  Do j=1,ny
     Do i=1,nx
        cldtot(3,i,j)=cldref(3)/(1-cldref(1)-cldref(2))
        cldtot(4,i,j)=min(cldref(4)/(1-cldref(1)-cldref(2)),1-cldtot(3,i,j))
     End Do
  End Do
#endif
  ! 
  !Parameterization of cloud fractions:
  !cld(1,i,j) deep+CsCc: depends on deep convective precip
  !approx proportional for typical values, capped at 1.
  !cld(2,i,j) cirrus: For typical values proportional to cld(1) 
  !with factor (cl2fac) based on ISCCP data. Modified by
  !cloud overlap condition that ensures total cloud fraction .lt. 1 
  !Overlap region of cld(1) with cld(2) uses cld(1) properties, 
  !hence (1-cld(1)) factor (equiv to random overlap assumption)
  !For cld(3) and cld(4) overlap with either cld(1) or cld(2)
  !uses cld(1)/cld(2) properties hence (1-cld(1)-cld(2)).
  !If using spatially constant cld(3)=cldref(3)/(1-cldref(1)-cldref(2)),
  !helps keep spatial mean similar to cldref(3); likewise for cld(4).
  !
  !For this paramzn sum_k cld(k).le.1 if cld(1).le.1, cldtot(2).le.1
  !and cldtot(3)+cldtot(4).le.1. For latter condition,
  !reduce cldtot(4) if cldtot(3)+cldtot(4).gt.1
  !
  !cld(0) clear sky: as residual. 
  !
  !
  Do j=1,NY
     Do i=1,NX
        cldtot(1,i,j)=min(cl1P*Qc(i,j),1.)   !depends on precip
        cld(1,i,j)=cldtot(1,i,j)  !when cld(1) overlaps it dominates
#ifdef COUNTCAP
        icount(i,j)=0
        if(cl1P*Qc(i,j).gt.1.) then
           icount(i,j)=1
           arr1(i,j)=icount(i,j)
        endif
#endif
        cldtot(2,i,j)=min(cld(1,i,j)*cl2fac,1.)   !before overlap calculation
        cld(2,i,j)=cldtot(2,i,j)*(1.-cld(1,i,j))  !non-overlapped by cld(1)
        cld(3,i,j)=cldtot(3,i,j)*(1.-cld(1,i,j)-cld(2,i,j)) 
        cld(4,i,j)=cldtot(4,i,j)*(1.-cld(1,i,j)-cld(2,i,j)) 
        cld(0,i,j)=1.-cld(1,i,j)-cld(2,i,j)-cld(3,i,j)-cld(4,i,j) 
        cl1(i,j)=cld(1,i,j)    !for output
!!$        arr2(i,j)=cld(2,i,j)
        arr3(i,j)=cld(3,i,j)
!!$        arr4(i,j)=cld(4,i,j)
     End Do
  End Do
  Return
End Subroutine cloud
!
!
!*********************************************************************
!
Subroutine readobscloud(cld, dir)
  !
  ! Read CLOUD COVER data files from ISCCP C2 (monthly)
  !
  Use Dimensions
  Use Clrad,    Only : clt
  Use Calendar
  Use AuxVars
  Implicit None
  Integer, Parameter :: nyOBS=32
  Integer, Parameter :: jskip=(ny-nyOBS)/2

  Real, Dimension(clt,nx,ny), Intent(Out) :: cld
  Character(Len=*),           Intent(In)  :: dir
  !
  !
  !- Local Variables and Functions
  Integer lnblnk, month_prior, n, i, j
  Character(Len=2)   :: monthc
  Character(Len=12)  :: file
  Character(len=305) :: fname
  Real :: cldmn
  !
  Data month_prior/0/
  Save month_prior
  !
  If(monthofyear.Eq.month_prior) Return
  !  jskip=(ny-nystd)/2
  !
  Write(monthc,'(i2.2)') monthofyear       !convert into a string
  file='0000'//monthc//'15'//'.cld'     !file format: 0000mm15.cld
  !
  ! read the var
  fname=dir(1:lnblnk(dir))//'/'//file
  Open(13,file=fname,form='formatted',status='old')
  Do n=1,3
     Do j=1+jskip,ny-jskip
        Do i=1,nx
           Read(13,*) cld(n,i,j)
        End Do
     End Do
     ! Compute zonal mean for fill in
     cldmn=0.
     Do i=1,nx
        cldmn=cldmn+cld(n,i,jskip+1)
     Enddo
     cldmn=cldmn/nx
     ! Fill in southern region
     Do j=1,jskip
        Do i=1,nx
           cld(n,i,j)=cldmn
        End Do
     End Do
     ! Compute zonal mean for fill in
     cldmn=0.
     Do i=1,nx
        cldmn=cldmn+cld(n,i,ny-jskip)
     Enddo
     cldmn=cldmn/nx
     Do j=ny-jskip+1,ny
        Do i=1,nx
           cld(n,i,j)=cldmn
        End Do
     End Do

  End Do
  Close(13)

  Write (*,*) 'readobscld: Read cloud fraction at model date', &
       &  dateofmodel,' from file:'
  Write (*,*) dir(1:lnblnk(dir))//'/'//file

  month_prior=monthofyear

!!$  do j=1,ny
!!$     do i=1,nx
!!$        arr5(i,j)=cld(1,i,j)
!!$        arr6(i,j)=cld(2,i,j)
!!$        arr7(i,j)=cld(3,i,j)
!!$     enddo
!!$  enddo

  Return
End Subroutine readobscloud


Subroutine radsw
  !************
  ! a simplified solar radiation scheme derived from Fu-Liou scheme (1993).
  !             FSW=S0*Ca(1)(1+bs(1))
  !             FSWds=S0*Ca(2)(1+bs(2))
  !      where
  !             Ca(i=1,2)=a(i=1,2)*cosZ+b(i=1,2
  !             bs(i=1,2)=c(i=1,2)*ALBDs
  !
  !     cosZ    : cosine solar zenith angle
  !     FSWds   : downward solar radiation at the surface
  !     FSWus   : solar reflection at the surface
  !     FSWut   : upward solar radiation at the top of the atmosphere
  !     SolarC  : solar constant
  !     S0      : S0 = SolarC * cosZ
  !     ALBDs   : surface albedo
  !     Ca(1)   : atmoshperic column absorption including clouds
  !     Ca(2)   : ratio of surface solar irradiance and S0
  !     bs(1)   : first order contribution due to surface reflection
  !     bs(2)   : first order contribution due to surface reflection
  !
  !************
  ! Solar radiation at TOA (S0), without diurnal cycle
  ! Astronomical formulae used:
  ! 1.  solar declination angle as a function of dayofyear
  ! 2.  correction(<3.4%) to solar constant due to the
  !     ellipticity of the earth orbit, function of dayofyear  Dickinson(1983)
  !
  ! --Local Variables & Functions:
  !    SolarC : solar constant
  !    cosZ   :  zenith angle in cosine, cosZ=0 at night
  !    cosZ1  : day average of cosZ
  !    cosZ2  : day average of cosZ**2
  !    lamda  : latitude
  !    delta  : solar inclination
  !    H      : hour angle at sunset and sunrise
  !
  ! --Altered Common Variables:
  ! FSW          !total SW absorbed by column   W/m**2
  !
  Use Fluxes
  Use Clrad
  Use Surface,           Only : ALBDs
  Use Grid,              Only : dy, pi
  Use PhysicalConstants, Only : Rearth
  Use Calendar
  Implicit None

  Integer i,j,k
  Real cosZ1,cosZ2,SolarC
  Real lamda,sinlamda,sindelta,coslamda                             &
       &,      cosdelta,sinH,sin2H,cosH                                   &
       &,      dtheta,delta,H
  Real dayspring,dayperigee,ecce,theta0max
  Real Ca(0:clt,2),bs(0:clt,2),a(0:clt,2),b(0:clt,2),c(0:clt,2)
  ! #if def SOLAR_CORR
  ! ! scorra=1-scorrs Correction factors for solar absorbtion
  !   Real, parameter :: scorrs=0.9, scorra=0. !, scorra=1.-scorrs
  !   logical :: firstcall = .true.
  ! #endif
  !
  !     The coefficients are calculated from Fu-Liou (1993) scheme with
  !     a prescribed ocean aerosol distribution.
  !
  Data (a(i,1),i=0,clt)/                                            &
       &              -0.7223241E-01,  0.2799124E-01, -0.3336635E-01,     &
       &              -0.3025872E-01, -0.5196485E-01/
  Data (b(i,1),i=0,clt)/                                            &
       &               0.2511654E+00,  0.1754200E+00,  0.2181227E+00,     &
       &               0.2593093E+00,  0.2454770E+00/
  Data (a(i,2),i=0,clt)/                                            &
       &               0.1867886E+00,  0.1466294E+00,  0.3125147E+00,     &
       &               0.2819863E+00,  0.3791754E+00/
  Data (b(i,2),i=0,clt)/                                            &
       &               0.5848927E+00,  0.8151811E-01,  0.3584743E+00,     &
       &               0.1479825E+00,  0.2931454E+00/
  Data (c(i,1),i=0,clt)/                                            &
       &               0.1662798E+00,  0.6341881E-01,  0.1743809E+00,     &
       &               0.6970352E-01,  0.6061328E+00/
  Data (c(i,2),i=0,clt)/                                            &
       &               0.1160093E+00,  0.7485661E+00,  0.2741095E+00,     &
       &               0.5867739E+00,  0.1347113E+01/
  !*****************************************************************
  ! #if def SOLAR_CORR
  !   if(firstcall) then
  !      firstcall=.false.
  !      print*,' Using solar correction: factor: scorrs=',scorrs
  !   endif
  ! #endif
  !
  !
  !- orbital parameters
  theta0max=23.447                 !degree; max solar dec. angle
  dayspring = 81.                  !Julian day of Vernal Equinox
  dayperigee = 3.                  !day of year closest to the sun
  ecce = 0.034                     !related to eccentricity of earth orbit
  !
  sindelta=Sin(theta0max/180.*pi) *                 &!Zeng (1994)
       &       Sin(2.*pi*(dayofyear-dayspring)/rdaysperyear)
  delta=Asin(sindelta)
  cosdelta=Cos(delta)
  SolarC =1370.*(1.+ecce*                 &! Dickinson(1983)
       &              Cos(2.*pi*(dayofyear-dayperigee)/rdaysperyear))
  dtheta=dy/Rearth
  Do j=1,ny
     lamda=(j-(ny+1)/2)*dtheta
     sinlamda=Sin(lamda)
     coslamda=Sqrt(1.-sinlamda**2)
     cosH=-sinlamda*sindelta/(coslamda*cosdelta)
     cosH=Max(cosH,-1.)       !polar; no night; halfday length=pi
     cosH=Min(cosH,1.)        !polar; no day; halfday length=0
     H=Acos(cosH)             
     sinH=Sqrt(1.-cosH**2)
     sin2H=2*sinH*cosH
     cosZ1=(sinlamda*sindelta*H                                      &
          &             +coslamda*cosdelta*sinH)/pi
     cosZ2=((sinlamda*sindelta)**2*H                                 &
          &           +2*sinlamda*coslamda*sindelta*cosdelta*sinH    &
          &           +(coslamda*cosdelta)**2*(H/2.-sin2H/4.))/pi
     Do k=0,clt
        Ca(k,1)=a(k,1)*cosZ2+b(k,1)*cosZ1
        Ca(k,2)=a(k,2)*cosZ2+b(k,2)*cosZ1
     Enddo
     Do i=1,nx
        ! empirical formula from Fu-Liou scheme
        Do k=0,clt
           bs(k,1)=c(k,1)*ALBDs(i,j)
           bs(k,2)=c(k,2)*ALBDs(i,j)
        Enddo
        FSW(i,j)=0.
        FSWds(i,j)=0.
        Do k=0,clt
           FSW(i,j)=FSW(i,j)+cld(k,i,j)*Ca(k,1)*(1.+bs(k,1))
           FSWds(i,j)=FSWds(i,j)+cld(k,i,j)*Ca(k,2)*(1.+bs(k,2))
        Enddo
        S0(i,j)=SolarC*cosZ1
        ! #if def SOLAR_CORR
        ! scorra=1-scorrs Correction factors for solar absorbtion
        !      FSW(i,j)=SolarC*(FSW(i,j)+FSWds(i,j)*(1.-ALBDs(i,j))*scorra)
        ! atmospheric absorption
        !      FSWds(i,j)=SolarC*(FSWds(i,j)*scorrs ) ! surface solar irradiance
        ! #else
        FSW(i,j)=SolarC*FSW(i,j)                    !atmospheric absorption
        FSWds(i,j)=SolarC*FSWds(i,j)                ! surface solar irradiance
        ! #endif
        FSWus(i,j)=FSWds(i,j)*ALBDs(i,j)
        FSWut(i,j)=S0(i,j)-FSWds(i,j)-FSW(i,j)+FSWus(i,j)
     End Do
  End Do

  Return
End Subroutine radsw


Subroutine radlw
  !
  ! Longwave radiation; weakly nonlinear scheme (LINRAD) of Chou & Neelin (1996)
  ! linearization coefficients calculated by Chou using a tropical mean profile
  !
  Use Prognostic, Only : T1,q1
  Use Surface,    Only : Ts,ALBDs
  Use Grid,       Only : dy
  Use QtcmPar,    Only : Tsref
  Use Fluxes
  Use Clrad
  Implicit None

  ! --Altered Common Variables:
  ! FLWds         !downward longwave flux at surface         W/m**2
  ! FLWus         !upward longwave flux at surface           W/m**2
  ! FLWut         !upward longwave flux at top, i.e., OLR    W/m**2
  ! FLW           !FLW=FLWus-FLWds-FLWut; total LW absorbed by column   W/m**2
  !
  ! --Local Variables & Functions:
  ! eps_rxxt    coeficients for calculating total column LW contributed
  !        by   x = T     temperature variation
  !               = q     moisture variation
  !               = c     cloudiness variation
  !               = co2   carbon dioxide
  !             x = s     surface temperature variation
  !        at   x = t     top of the atmosphere
  !               = s     the surface
  ! eps_rxx     coeficients for calculating column LW of each cloud type
  !
  !
  Integer i,j,n
  Real eps_rTst,eps_rqst,eps_rTtt               &!LW radiation,coeff.
       &,      eps_rqtt,eps_rstt                                          &
       &,      eps_rcs(clt),eps_rct(clt)                                  &
       &,      eps_rTt(0:clt),eps_rqt(0:clt)                              &
       &,      eps_rst(0:clt),eps_rTs(0:clt)                              &
       &,      eps_rqs(0:clt),eps_rss                                     &
#ifdef DELTA_CO2
       &,      FLWdsref,FLWusref,FLWutref                                 &
       &,      eps_rco2s(0:clt),eps_rco2t(0:clt),eps_rco2st,eps_rco2tt    &
       &,      co2m                                                       
#else
       &,      FLWdsref,FLWusref,FLWutref                                 
#endif
  !
  !*****************************************************************
  Data (eps_rct(n),n=1,clt) /                                       &
       &     -0.100751E+03, -0.616890E+02, -0.122923E+02,  -0.276311E+02/
  Data (eps_rcs(n),n=1,clt) /                                       &
       &      0.230918E+02,  0.843160E+01,  0.325455E+02,   0.205806E+02/
  Data (eps_rTt(n),n=0,clt) /                                       &
       &      0.133549E+01,  0.152938E+01,  0.145419E+01,  0.141225E+01,  &
       &      0.140892E+01/
  Data (eps_rqt(n),n=0,clt) /                                       &
       &     -0.806434E+00, -0.644083E-01, -0.352098E+00, -0.471530E+00,  &
       &     -0.421301E+00/
  Data (eps_rst(n),n=0,clt) /                                       &
       &      0.535400E+00,  0.109650E-02,  0.208200E+00,  0.900000E-03,  &
       &      0.115160E+00/
  Data (eps_rTs(n),n=0,clt) /                                       &
       &      0.129186E+01,  0.174937E+01,  0.151779E+01,  0.183933E+01,  &
       &      0.168531E+01/
  Data (eps_rqs(n),n=0,clt) /                                       &
       &      0.259074E+01,  0.656097E+00,  0.173063E+01,  0.118402E+00,  &
       &      0.862426E+00/
  Data eps_rss   /                                                  &
       &      0.628300E+01/
  Data FLWdsref /443.4288/,    FLWusref /475.3227/,                 &
       &     FLWutref /240.3406/
#ifdef DELTA_CO2
  Data (eps_rco2t(n),n=0,clt) /                                     &
       &     -0.564960E+01, -0.326653E+01, -0.419050E+01, -0.559220E+01,  &
       &     -0.522538E+01/
  Data (eps_rco2s(n),n=0,clt) /                                     &
       &      0.188300E+00,  0.131395E+00,  0.165100E+00,  0.112100E+00,  &
       &      0.139727E+00/
  Data co2m /330.00/
#endif
  !*****************************************************************
  !
  Do j=1,ny
     Do i=1,nx
#ifdef DELTA_CO2
        co2(i,j)=330.
#endif
        eps_rTst=0.
        eps_rqst=0.
        eps_rTtt=0.
        eps_rqtt=0.
        eps_rstt=0.
#ifdef DELTA_CO2
        eps_rco2st=0.
        eps_rco2tt=0.
#endif

        Do n=0,clt                              !calculating eps_rxxt
           eps_rTst=eps_rTst+eps_rTs(n)*cld(n,i,j)  
           eps_rqst=eps_rqst+eps_rqs(n)*cld(n,i,j)
           eps_rTtt=eps_rTtt+eps_rTt(n)*cld(n,i,j)
           eps_rqtt=eps_rqtt+eps_rqt(n)*cld(n,i,j)
           eps_rstt=eps_rstt+eps_rst(n)*cld(n,i,j)
#ifdef DELTA_CO2
           eps_rco2tt=eps_rco2tt+eps_rco2t(n)*cld(n,i,j)
           eps_rco2st=eps_rco2st+eps_rco2s(n)*cld(n,i,j)
#endif
        Enddo
        !*********   LW contributed by T, q, Ts and cloudiness
#ifdef DELTA_CO2
        FLWds(i,j)=FLWdsref+eps_rTst*T1(i,j)+eps_rqst*q1(i,j)          &
             &              +eps_rco2st*(co2(i,j)-co2m)/330.           
#else
        FLWds(i,j)=FLWdsref+eps_rTst*T1(i,j)+eps_rqst*q1(i,j)          
#endif
        FLWut(i,j)=FLWutref+eps_rTtt*T1(i,j)+eps_rqtt*q1(i,j)          &
#ifdef DELTA_CO2
             &              +eps_rstt*(Ts(i,j)-Tsref)                  &
             &              +eps_rco2tt*(co2(i,j)-co2m)/330.           
#else
             &              +eps_rstt*(Ts(i,j)-Tsref)                  
#endif
        Do n=1,clt                         !LW contributed by cloudiness
           FLWds(i,j)=FLWds(i,j)+eps_rcs(n)*(cld(n,i,j)-cldref(n))
           FLWut(i,j)=FLWut(i,j)+eps_rct(n)*(cld(n,i,j)-cldref(n))
        Enddo
        !if(Ts(i,j).lt.226.6) then
        !   print *, 'Ts < Tsmin, causing FLWus < 0'
        !   stop
        !endif
        FLWus(i,j)=max(FLWusref+eps_rss*(Ts(i,j)-Tsref),0.)  !LW contributed Ts
        FLW(i,j)=FLWus(i,j)-FLWds(i,j)-FLWut(i,j)          !sum of total LW
     End Do
  End Do

  Return
End Subroutine radlw
