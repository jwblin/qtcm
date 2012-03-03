!
! File  clrad1d.f90
!
! This contains the cloud and radiation package of Chou and Neelin (Chou, Ph.D
!   Dissertation, UCLA 1997)
!
! Code history:
!   beta version(QTCM2.0):  Chou/Neelin  June 1997
!
!**************************************************************************
Subroutine cloud
  !************
  ! cloud prediction scheme
  !
  Implicit None
  Include 'hgrid.h'
  Include 'qtcm.h'
  Include 'cld.h'
  Include 'pmean.h'


  Integer itdriver,it
  Common/modelday/itdriver,it

  ! --Input Variables:
  !   Qc, cldref

  ! --output Variables:
  !   cldref: reference (mean) cloud cover
  !   cld: cloud cover
  !   cld: cloud cover
  !                   type 0: clear sky
  !                   type 1: Deep cloud +CsCc
  !                   type 2: Cirrus
  !                   type 3: Stratus

  ! --Local Variables & Functions:
  Integer i,j,n
  Real*8 cl1P     !coefficients for cloud prediction scheme
  Data cl1P/7.76275869D-4/
  Data (cldref(n),n=1,clt)/
1      0 .1051000,      0.1047000,     0.1096000/

  !*****************************************************************

  !      call readobscloud('cld','../bnddata/CLOUD_ISCCP')
  !                                             !read ISCCP observed cloudiness
  Do j=1,NY
     Do i=1,NX
        cld(i,j,1)=cl1P*Qc(i,j)  !proportional to deep convective precipitation
        !        cld(i,j,2)=cldref(2)
        !        cld(i,j,3)=cldref(3)
        cld(i,j,0)=1.-cld(i,j,1)-cld(i,j,2)-cld(i,j,3)
        cl3(i,j)=cld(i,j,1)
     End Do
  End Do


  Return
End Subroutine cloud

!
!*********************************************************************
!
! Read CLOUD COVER data files from ISCCP C2 (monthly)
!
!*********************************************************************
! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/clrad1d.f90,v $
! $Author: munnich $
! $Revision: 1.5 $
! $Date: 2002/07/16 22:30:12 $
!
! File: ocean.f
!
! Read bndary data files (monthly climatology)
!
! The date information is converted into file names:
!   In this version, bc file names include the date and variable name
!   in this way:
!
!               file name = yyyymmdd.var
!
!   where yyyy = year , mm = month, dd=dayofmonth, var = a variable name.
!
!   if SSTmode is seasonal, yyyy=0000
!
! Code history:
!   Original version:  Ning Zeng/Bill Weibel  May 1998, adopted from ocean.f
!
! Routines called but not in this file
!   date_of_model      in driver.f
!*********************************************************************

Subroutine readobscloud(var, vardir)
  !************
  ! Time evolving boundary conditions other than SST:
  !   prescribed fields such as observed soil moisture, surface
  !   radiative fluxes, cloud cover, etc.
  ! It can be called from the ocean-atmo coupling or land-atmo coupling
  !   loops depending on the input data's time resolution, or other
  !   coupling intervals if desired.
  !*********************************

  Implicit None
  Include 'hgrid.h'
  Include 'calendar.h'
  Include 'cld.h'

  Real*8 vardummy(NX,NY,3)
  Character*(3) var
  Character*(*) vardir

  !- Altered
  !    vardummy

  !- Local Variables and Functions
  Integer lnblnk, date_of_model, month, month_prior, n, i, j
  Character*(2) monthc
  Character*(12) varfile

  Data month_prior/0/
  Save month_prior
  month=Mod(date_of_model(),10000)/100
  If(month.Eq.month_prior) Return

  Write(monthc,'(i2.2)') month              !convert into a string
  varfile='0000'//monthc//'15'//'.'//var    !file format: 0000mm15.var

  ! read the var
  Open(13,file=vardir(1:lnblnk(vardir))//'/'//varfile
  &       ,form='formatted',status='old')
  Read(13,*) vardummy
  Close(13)
  Do n=1,3
     Do j=1,ny
        Do i=1,nx
           cld(i,j,n)=vardummy(i,j,n)
        Enddo
     Enddo
  Enddo

  Write (*,*) 'bndry: Read ',var,' at model date', date_of_model(),
  &' from file',vardir(1:lnblnk(vardir))//'/'//varfile

  month_prior=month

  Return
End Subroutine readobscloud


Subroutine radsw
  !************
  ! a simplified solar radiation scheme derived from Fu-Liou scheme (1993).
  !             FSW=S0*Ca(1)(1+bs(1))
  !             FSWds=S0*Ca(2)(1+bs(2))
  !      where
  !             Ca(i=1,2)=a(i=1,2)*cosZ+b(i=1,2)
  !             bs(i=1,2)=c(i=1,2)*ALBDs
  !
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
  !    cosZ   :  zenith angle in cosine, cosZ=0 at night
  !    lamda  : latitude
  !    delta  : solar inclination
  !    hh     : hour angle, h=0 at noon
  !
  ! --Altered Common Variables:
  ! FSW          !total SW absorbed by column   W/m**2
  !
  Implicit None
  Include 'hgrid.h'
  Include 'qtcm.h'
  Include 'calendar.h'
  Include 'cld.h'
  Include 'pmean.h'
  Include 'driver.h'

  !
  ! --Altered Common Variables:
  ! FSW          !total SW absorbed by column   W/m**2

  ! --Local Variables & Functions:
  Integer i,j,k,m,maxb
  Parameter(maxb=72)
  Real*8 lamda,delta,hh
  Real*8 sinlamda(ny),sindelta,coslamda(ny)
2 ,      cosdelta,coshh(nx,maxb),SolarC
  Real*8 dayspring,dayperigee,ecce,theta0max,dtheta
  Real*8 Ca(0:clt,2),bs(0:clt,2),a(0:clt,2),b(0:clt,2),c(0:clt,2)

  !      data (a(i,1),i=0,clt)
  !     1            / -6.5183550D-02,  3.5192758D-02, -2.4936408D-02,
  !     2              -1.1718750D-02/
  !      data (b(i,1),i=0,clt)
  !     1            /  0.2267621    ,  0.1460207    ,  0.1901388    ,
  !     2               0.2174571    /
  !      data (a(i,2),i=0,clt)
  !     1            /  0.1291060    ,  0.1611494    ,  0.2892808    ,
  !     2               0.3057677    /
  !      data (b(i,2),i=0,clt)
  !     1            /  0.6696874    ,  9.0886146D-02,  0.4165963    ,
  !     2               0.1688761    /
  !      data (c(i,1),i=0,clt)
  !     1            /  0.1694734    ,  7.8532398D-02,  0.1926342    ,
  !     2               9.1768831D-02/
  !      data (c(i,2),i=0,clt)
  !     1            /  7.5819805D-02,  0.7279596    ,  0.2416583    ,
  !     2               0.5395555    /
  Data (a(i,1),i=0,clt)/
  *              -0.6991568D-01,  0.2862155D-01, -0.3275377D-01,
  *              -0.3075491D-01/
  Data (b(i,1),i=0,clt)/
  *               0.2489715D+00,  0.1758697D+00,  0.2177709D+00,
  *               0.2603912D+00/
  Data (a(i,2),i=0,clt)/
  *               0.1344528D+00,  0.1476796D+00,  0.2950021D+00,
  *               0.2853160D+00/
  Data (b(i,2),i=0,clt)/
  *               0.6461399D+00,  0.8235896D-01,  0.3858907D+00,
  *               0.1503189D+00/
  Data (c(i,1),i=0,clt)/
  *               0.1659058D+00,  0.6383526D-01,  0.1763771D+00,
  *               0.6906892D-01/
  Data (c(i,2),i=0,clt)/
  *               0.7897843D-01,  0.7440951D+00,  0.2550164D+00,
  *               0.5797699D+00/
  !*****************************************************************


  !- orbital parameters
  theta0max=23.447                 !degree; max solar dec. angle
  dayspring = 81.                  !Julian day of Vernal Equinox
  dayperigee = 3.                  !day of year closest to the sun
  ecce = 0.034                     !related to eccentricity of earth orbit

  If(timeofday .Eq. 0.)Then
     sindelta=dsin(theta0max/180.*pi) *                 !Zeng (1994)
     &         dsin(2.*pi*(dayofyear-dayspring)/360.)    !1yr=360day for QTCM
     delta=dasin(sindelta)
     cosdelta=dcos(delta)
     SolarC =1370.*(one+ecce*                 ! Dickinson(1983)
     &                dcos(2.*pi*(dayofyear-dayperigee)/365.))
     dtheta=dy/Rearth
     sindelta=dsin(delta)
     cosdelta=dsqrt(one-sindelta**2)
     If(dayofyear .Eq. 1)Then
        If((86400/dt) .Gt. maxb)
2       Write(*,*)'hour angle array is too small, change maxb'
        Do j=1,ny
           lamda=(j-(NY+1)/2)*dtheta
           sinlamda(j)=dsin(lamda)
           coslamda(j)=dsqrt(one-sinlamda(j)**2)
        Enddo
        m=86400./dt
        Do k=1,m
           Do i=1,NX
              hh=2*pi*k/m-pi+(i-1/2.)*dx/Rearth
              coshh(i,k)=dcos(hh)
           Enddo
        Enddo
     Endif
  Endif
  m=timeofday/dt+1
  Do j=1,NY
     Do i=1,NX
        cosZ=sinlamda(j)*sindelta+coslamda(j)*cosdelta*coshh(i,m)
        S0(i,j)=0.
        FSWds(i,j)=0.                              ! at night
        FSW(i,j)=0.
        FSWus(i,j)=0.
        FSWut(i,j)=0.
        If(cosZ .Gt. 1.D-8)Then               ! at day
           !******************************* empirical formula from Fu-Liou scheme
           Do k=0,clt
              Ca(k,1)=a(k,1)*cosZ+b(k,1)
              Ca(k,2)=a(k,2)*cosZ+b(k,2)
           Enddo
           Do k=0,clt
              bs(k,1)=c(k,1)*ALBDs(i,j)
              bs(k,2)=c(k,2)*ALBDs(i,j)
           Enddo
           !c
           FSW(i,j)=0
           FSWds(i,j)=0
           Do k=0,clt
              FSW(i,j)=FSW(i,j)+cld(i,j,k)*Ca(k,1)*(1.+bs(k,1))
              FSWds(i,j)=FSWds(i,j)+cld(i,j,k)*Ca(k,2)*(1.+bs(k,2))
           Enddo
           S0(i,j)=SolarC*cosZ
           FSW(i,j)=S0(i,j)*FSW(i,j)               !atmospheric absorption
           FSWds(i,j)=S0(i,j)*FSWds(i,j)           ! surface solar irradiance
           FSWus(i,j)=FSWds(i,j)*ALBDs(i,j)
           FSWut(i,j)=S0(i,j)-FSWds(i,j)-FSW(i,j)+FSWus(i,j)
           ! re-arrange the portion of SW absorbed by atmosphere and surface
           !        FSW(i,j)=FSW(i,j)+50.
           !        FSWds(i,j)=FSWds(i,j)-50./(1-ALBDs(i,j))
           !        FSWus(i,j)=FSWds(i,j)*ALBDs(i,j)
           !        FSWut(i,j)=S0(i,j)-FSWds(i,j)-FSW(i,j)+FSWus(i,j)
        Endif
     End Do
  End Do

  Return
End Subroutine radsw


Subroutine radlw
  !***
  ! Longwave radiation; weakly nonlinear scheme (LINRAD) of Chou & Neelin (1996)
  ! linearization coefficients calculated by Chou using a tropical mean profile
  !
  Implicit None
  Include 'hgrid.h'
  Include 'qtcm.h'
  Include 'cld.h'
  Include 'pmean.h'


  ! --Altered Common Variables:
  ! FLWds         !downward longwave flux at surface         W/m**2
  ! FLWus         !upward longwave flux at surface           W/m**2
  ! FLWut         !upward longwave flux at top, i.e., OLR    W/m**2
  ! FLW           !FLW=FLWus-FLWds-FLWut; total LW absorbed by column   W/m**2

  ! --Local Variables & Functions:
  ! eps_rxxt    coeficients for calculating total column LW contributed
  !        by   x = T     temperature variation
  !               = q     moisture variation
  !               = c     cloudiness variation
  !             x = s     surface temperature variation
  !        at   x = t     top of the atmosphere
  !               = s     the surface
  ! eps_rxx     coeficients for calculating column LW of each cloud type
  !


  Integer i,j,n
  Real*8 eps_rTst,eps_rqst,eps_rTtt               !LW radiation,coeff.
1 ,      eps_rqtt,eps_rstt
2 ,      eps_rcs(clt),eps_rct(clt)
3 ,      eps_rTt(0:clt),eps_rqt(0:clt)
4 ,      eps_rst(0:clt),eps_rTs(0:clt)
5 ,      eps_rqs(0:clt),eps_rss
7 ,      FLWdsref,FLWusref,FLWutref

  !      data (eps_rct(n),n=1,clt) /
  !     1    -100.8353,   -61.6350,   -12.2580/
  !      data (eps_rcs(n),n=1,clt) /
  !     1      24.8918,     8.3830,    36.1260/
  !      data (eps_rTt(n),n=0,clt) /
  !     1      0.926000,      0.553503,      0.698000,      0.998000/
  !      data (eps_rqt(n),n=0,clt) /
  !     1     -1.477000,     -0.000502,     -0.575000,     -1.079000/
  !      data (eps_rst(n),n=0,clt) /
  !     1      0.533000,      0.000000,      0.207000,      0.001000/
  !      data (eps_rTs(n),n=0,clt) /
  !     1      1.014000,      1.397605,      1.200000,      1.4760000/
  !      data (eps_rqs(n),n=0,clt) /
  !     1      2.345000,      0.507853,      1.508000,      0.0390000/
  !      data eps_rss   /6.283000/
  !ccccccc
  Data (eps_rct(n),n=1,clt) /
  *     -0.100388D+03, -0.616890D+02, -0.118121D+02/
  Data (eps_rcs(n),n=1,clt) /
  *      0.248086D+02,  0.843160D+01,  0.348483D+02/
  Data (eps_rTt(n),n=0,clt) /
  *      0.925083D+00,  0.555392D+00,  0.697919D+00,  0.994707D+00/
  Data (eps_rqt(n),n=0,clt) /
  *     -0.147942D+01, -0.836689D-02, -0.575450D+00, -0.109708D+01/
  Data (eps_rst(n),n=0,clt) /
  *      0.535400D+00,  0.303675D-02,  0.208200D+00,  0.218000D-01/
  Data (eps_rTs(n),n=0,clt) /
  *      0.101384D+01,  0.139465D+01,  0.120000D+01,  0.145793D+01/
  Data (eps_rqs(n),n=0,clt) /
  *      0.234464D+01,  0.520346D+00,  0.150826D+01,  0.128629D+00/
  Data eps_rss   /
  *      0.628300D+01/
  !ccccccc   from nonlinear scheme
  Data FLWdsref /441.7868/,    FLWusref /478.1020/,
2 FLWutref /245.9178/
  !*****************************************************************

  Do j=1,NY
     Do i=1,NX
        eps_rTst=0  
        eps_rqst=0
        eps_rTtt=0
        eps_rqtt=0
        eps_rstt=0
        Do n=0,clt                              !calculating eps_rxxt
           eps_rTst=eps_rTst+eps_rTs(n)*cld(i,j,n)  
           eps_rqst=eps_rqst+eps_rqs(n)*cld(i,j,n)
           eps_rTtt=eps_rTtt+eps_rTt(n)*cld(i,j,n)
           eps_rqtt=eps_rqtt+eps_rqt(n)*cld(i,j,n)
           eps_rstt=eps_rstt+eps_rst(n)*cld(i,j,n)
        Enddo
        !*********   LW contributed by T, q, Ts and cloudiness
        FLWds(i,j)=FLWdsref+eps_rTst*T1(i,j)+eps_rqst*q1(i,j)
        FLWut(i,j)=FLWutref+eps_rTtt*T1(i,j)+eps_rqtt*q1(i,j)
        *                 +eps_rstt*(Ts(i,j)-Tsref)
        Do n=1,clt                         !LW contributed by cloudiness
           FLWds(i,j)=FLWds(i,j)+eps_rcs(n)*(cld(i,j,n)-cldref(n))
           FLWut(i,j)=FLWut(i,j)+eps_rct(n)*(cld(i,j,n)-cldref(n))
        Enddo
        FLWus(i,j)=FLWusref+eps_rss*(Ts(i,j)-Tsref)        !LW contributed Ts
        FLW(i,j)=FLWus(i,j)-FLWds(i,j)-FLWut(i,j)          !sum of total LW
     End Do
  End Do

  Return
End Subroutine radlw

