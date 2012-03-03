! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: qtcm.F90
!
! This file contains most of the subroutines and functions called by
!   subroutine qtcm, which forms the heart of the Quasi-Equilibrium
!   Tropical Circulation Model (QTCM). The cloud and radiation subroutines
!   are in a separate file and so is the land surface model
!
! Code history:
!   Original (QTCM1 Version 1a):  Ning Zeng, David Neelin, et al. Jan. 1997
!   Beta version (QTCM1 V1b):     June 1997
!   QTCM1 release (QTCMes1 V2.0):   June 1998
!   QTCM1 V2.1 :                  Jan  1999
!   QTCM1 V2.2f90 :              Nov 2000 
!
! Routines called but not in this file:
!   getSST          in getbnd.F90
!   varmean         in output.F90
!   sland           in sland.F90
!   radsw \
!   radlw  >        in clrad.F90
!   cloud /
!   sflux/mlsflux   in sflux.F90/mlsflux.F90
! and if chosen by a preprocessor option      
!   cplmean         in cplmean.F90
!   set_qclock      in driver_util.F90
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
Subroutine qtcminit()
  !
  ! Initialize qtcm atmospheric model variables, datasets, and
  ! diagnostic fields.
  !
  Implicit None

  !(begin comment out)
  !   In Version 2.3 - current-release (August 2002):
  !   This is the original order of the routine.  However, the problem
  !   with this order is that STYPE is not read in for landon until
  !   bndinit, but in varinit you use STYPE to calculate the original
  !   values of WD for the non-restart case.  The manual is also gives
  !   two conflicting (p. 29 and 32) orders.
  !Call parinit            !Initialize model parameters
  !Call varinit            !Initialize variables
  !Call TimeManager(1)     !mm set model time 
  !Call bndinit            !input boundary datasets
  !Call physics1           !diagnostic fields for initial condition
  !(end comment out)

  ! JWL:  I've altered the above order to this order in order for 
  ! varinit to work correctly.  bndinit doesn't seem to have any
  ! dependencies such that it has to come after TimeManager or
  ! varinit.
  Call parinit            !Initialize model parameters
  Call bndinit            !input boundary datasets
  Call varinit            !Initialize variables
  Call TimeManager(1)     !mm set model time
  Call physics1           !diagnostic fields for initial condition

  Return
End Subroutine qtcminit
!
! ---------------------------------------------------------------------
!
Subroutine qtcm(interval)
  !
  ! Main routine for the Quasi-Equilibrium Tropical Circulation Model.
  !
  Use Input, Only : dt, mt0, nastep, it  , noout
  Use Calendar, Only : dayofmodel  ! for outpInst

  Implicit None

  ! --Arguments:
  Integer interval         ! interval of ocean-atmos.
  ! coupling, unaltered on exit

  ! --Local Variables
  !  Integer it
  !----------------------------------

  Call getbnd          !get time evolving bndry data, SST, albedo etc.

  !- number of atmospheric time steps within one air-sea coupling interval
  nastep = interval* (86400./dt)

  If(Mod(nastep,mt0).Ne.0) Then
     Print*,'qtcm: Fatal error: mt0*dt not a divisor of 1 day'
     Stop
  Endif

  !- land-atmosphere coupling loop
  !    land runs at the same timestep as atmo., easy modification otherwise,
  !    but if land runs less often than [sflux] in [physics1], Evap over land 
  !    needs to be fixed in [sflux]
  Do it=1,nastep
!!$#if def DIURNAL
!!$     Call set_qclock(it)   !advance the diurnal clock
!!$#endif
     Call physics1         !physics associated with temperature mode T1
     Call sland1           !land

     Call advctuv       !advection of momentum u1,v1,u0,v0
     Call advctTq       !advection of T1, q1
     Call dffus         !diffusion of u1,v1,u0,v0,T1,q1

     Call barcl         !the baroclinic mode
     If(Mod(it,mt0).Eq.0) then
#ifdef NO_ABL
        Call bartr      !barotropic mode
#else
        Call Savebartr !save u0,v0 for surface geopotential
        Call bartr      !barotropic mode  @@@JWL: remove underscore f/ name
        Call gradphis  !surface geopot. gradient  @@@JWL: rm underscore f/ name
#endif
     End If
     !If(Mod(it,3)==0) Call varmean          !time mean (sum up variables, output in output.F)
     Call varmean          !time mean (sum up variables, output in output.F)
#ifdef MXL_OCEAN
     ! ocean mixed layer runs need mean sfc fluxes
#define CPLMEAN
#endif

#ifdef BLEND_SST
#define CPLMEAN
#endif

#ifdef CPLMEAN
     Call cplmean(it)      !optional mean sfc. flux for air-sea coupling
#endif

     !
     ! output within the air-sea coupling interval, useful for model diagnosis.
     !   you may want to set the user-control parameter ntouti=1 (day)
     ! if(dayofmodel< 1) call outpInst

  End Do
  Return
End Subroutine qtcm
!
! ---------------------------------------------------------------------
!
Subroutine parinit
  !
  ! Initialize physical constants and model parameters
  !
  Use PhysicalConstants
  Use Qtcmpar
  Use Grid
#ifdef SPONGES
  Use Sponges
#endif
  Use Input
  Use Calendar

  Implicit None

  ! --File I/O:  Subroutine opens and closes following files:
  !     unit 12 = qtcmpar.in
  !   units closed prior to entry, and on exit of the routine
  !
  ! --Local Variables:
  Real thetau,thetav
#ifdef SPONGES
  Real sig,y0,yj,y0dfs
#endif
  Integer i,j,k

  pi=Asin(1.)*2.
  !
  ! set up saturation humidity look-up table
  Call humtable
#ifdef LINEAR_T1C
  !
#else
  ! set up nonlinear T1c look-up table
  Call t1ctable
#endif
  !
  !- set up model domain, grid size, geometric factor, Coriolis etc.
  dx=2.*pi*Rearth/nx           !whole globe, periodic in x
  dxi=1./dx
  dy=YB/90.*pi*Rearth/ny              
  dyi=1./dy

  Do j=0,ny      !for staggered variables along y; pole points not needed
     thetav=dy/Rearth*(j+.5-(ny+1.)/2.)
     fv(j)=2.*Omega*Sin(thetav)  !Coriolis parameter at half grid (j+1/2)
     cosv(j)=Cos(thetav)         !cos(lat) at half (j+1/2) (staggered)
     !         cosv(j)=1.                  !turns off spherical geometry
     cosvi(j)=1./cosv(j)
     dxvi(j)=cosvi(j)*dxi 
     dyvi(j)=cosvi(j)*dyi 
  End Do

  Do j=1,ny        !for centered variables
     thetau=dy/Rearth*(j-(ny+1.)/2.)
     fu(j)=2.*Omega*Sin(thetau)  !Coriolis parameter at full grid points
     cosu(j)=Cos(thetau)         !cos(lat) at full grid   j
     !         cosu(j)=1.                  !turns off spherical geometry
     cosui(j)=1./cosu(j)                  
     dxui(j)=cosui(j)*dxi 
     dyui(j)=cosui(j)*dyi 
  End Do
!- consider spherical coordinates; coef. for diffusion terms  
!- all weights must satisfy sum_j(cos(j)sum_k(weight4(j,k)T(j+k-3))) = 0.  
!- for interior points and if cos(j)=1, weights are (1,-4,6,-4,1) as in
!- Cartesian Coordinates
  Do j=2,ny-1
     weight2u(j,1)=cosv(j)*cosui(j)
     weight2u(j,2)=-(cosv(j)+cosv(j-1))*cosui(j)
     weight2u(j,3)=cosv(j-1)*cosui(j)
     weight4u(j,1)=cosu(j+1)*cosui(j)
     weight4u(j,2)=-2.*(cosu(j+1)+cosu(j))*cosui(j)
     weight4u(j,3)=(2.*(cosu(j+1)+cosu(j)*2.+cosu(j-1))              &
    &                  -cosu(j+1)-cosu(j-1))*cosui(j)
     weight4u(j,4)=-2.*(cosu(j)+cosu(j-1))*cosui(j)
     weight4u(j,5)=cosu(j-1)*cosui(j)
  End Do
  Do j=1,ny-1
     weight2v(j,1)=cosu(j+1)*cosvi(j)
     weight2v(j,2)=-(cosu(j+1)+cosu(j))*cosvi(j)
     weight2v(j,3)=cosu(j)*cosvi(j)
     weight4v(j,1)=cosv(j+1)*cosvi(j)
     weight4v(j,2)=-2.*(cosv(j+1)+cosv(j))*cosvi(j)
     weight4v(j,3)=(2.*(cosv(j+1)+cosv(j)*2.+cosv(j-1))              &
    &                  -cosv(j+1)-cosv(j-1))*cosvi(j)
     weight4v(j,4)=-2.*(cosv(j)+cosv(j-1))*cosvi(j)
     weight4v(j,5)=cosv(j-1)*cosvi(j)
  End Do
  Do j=1,ny
     latt(j) = (YB/ny)*(2*j-ny-1)
  End Do
  Do i=1,nx
     lont(i)=(360./Real(nx))*(i-1)
  End Do

#ifdef SPONGES
  Do j=1,ny
     !
     !c-sponge boundary outside y0 N/S with a linear dropoff; 
     !with a lower cutoff 
     ! @Only spngh1 and spngh3 are used currently to reduce excessive 
     !bnd. precip. 
     !sig=10.0               !in degree latitude
     !y0=45.0
     !     y0dfs=31.0
     !yj=thetau*180.0/pi
     !spngh0(j)=0.               !mask sponge boundary (for Ts land)
     !spngh1(j)=1.                !reduce heat/moisture flux
     !spngh2(j)=1.                !increase diffusion for T,U,V,Q 
     !spngh3(j)=1.                !increase momentum damping
     !spngh4(j)=1.                !reduce precipitation Qc only
     !
     !c-model results are somewhat sensitive to the use of spngh
     !c-ZNC rad1 used spngh0 and spngh1; ZNC rad0 used spngh0,spngh1 and spngh3 
     !c-Su et al. used spngh0 and spngh2  
     !if(abs(yj).gt.y0) then
     !   spngh0(j)=min(1.,(abs(yj)-y0)/sig)
     !   spngh0(j)=min(10,(abs(yj)-y0))
     !   spngh0(j)=min(15,(abs(yj)-y0))
     !   spngh1(j)=max(.5,1.-(abs(yj)-y0)/sig)
     !   spngh2(j)=1./max(.25,1.-(abs(yj)-y0)/sig)
     !   spngh3(j)=1./max(.1,1.-(abs(yj)-y0)/sig)
     !   spngh4(j)=max(.10,1.-(abs(yj)-y0)/sig)
     !end if
     !          if(abs(yj).gt.y0dfs) then
     !            spngh2(j)=1./max(.250,1.-(abs(yj)-y0dfs)/sig)
     !          endif
  End Do
  Write(*,*)'parinit: spngh0,1,2,3,4='                              &
       & ,spngh0(1),spngh1(1),spngh2(1),spngh3(1),spngh4(1)
#endif  
  !
  !c- To handle the periodic boundary condition in finite differencing along x,
  ! we use a 'double referencing' approach;  the following is the indexing:
  Do i=1,nx
     im1(i)=i-1
     ip1(i)=i+1
     im2(i)=i-2
     ip2(i)=i+2
  End Do
  im1(1)=nx
  ip1(nx)=1
  im2(1)=nx-1
  im2(2)=nx
  ip2(nx)=2
  ip2(nx-1)=1


  Return
End Subroutine parinit
!
! ---------------------------------------------------------------------
!
Subroutine bndinit
  !
  ! Time-independent boundary conditions; Initializing time-
  ! dependent boundaries
  !

  Use Surface, Only : nx,ny, ALBDs, STYPE, TOP, CDN
  Use Input,   Only : bnddir, landon
  Use Land,    Only : Z0

  Implicit None

  ! --File I/O:  Subroutine opens and closes following files:
  !     unit 13 = SST01
  !     unit 14 = STYPE 
  !     unit 15 = TOP
  !   units closed prior to entry, and on exit of the routine
  !
  ! --Altered Module Variables:
  ! Ts,CDN,STYPE,ALBDs
  !
  ! --Local Variables & Functions:
  Real VONKAR,hPBL
  Integer, Dimension (nx,ny) :: iSTYPE
  Integer i,j
  Character(len=305) :: fname
  ! Functions
  Integer lnblnk
  !
  !c input  !return the directory name with the trailing blanks removed
  fname=bnddir(1:lnblnk(bnddir))//'/ALBD_Darnell/00001315.alb'
  Open(13,file=fname,status='old')
  Read(13,*) ALBDs  ! annual mean surface albedo
  Close(13)

  fname=bnddir(1:lnblnk(bnddir))//'/STYPE'
  Open(13,file=fname, status='old')
  Read(13,*) STYPE   ! surface type 
  Close(13)
  ! print*,'Surface type:'
  ! print '(64i1)',((ifix(STYPE(i,j)),i=1,nx),j=ny,1,-1)

#ifdef TOPO
  fname=bnddir(1:lnblnk(bnddir))//'/TOP'
  Open(13,file=fname,status='old')
  Read(13,*) TOP   ! relative topography, height/10km
  Close(13)

  Where (top<0.1) top=0.0
  ! Print*,'Topography [500m]:'
  ! Print '(64i1)',((ifix(top(i,j)*20),i=1,nx),j=ny,1,-1)
#endif


  !  Treat land surface as ocean if land is turned off
  If(landon.Ne.1) Then
     STYPE=0.  !array
  Else

  End If
  !
  != Compute CDN and ALBDs using the assigned surface type
  VONKAR=0.4      !Von Karman constant
  Do j=1,ny
     Do i=1,nx
        !
        !- neutral drag
        !  Deardorff (1972) (CCM1 formulation, Williamson et al. 1987; 
        !NCAR TN285)
        hPBL=2000.      !2km PBL depth
        CDN(i,j)= 1./                                                 &
             &      (Log(.025*hPBL/Z0(Int(STYPE(i,j))))/VONKAR+8.4)**2
        !
        !  BATS/CCM2 formulation; CDN larger than CCM1 especially over forest
        !         hPBL=75.        !PBL depth
        !         CDN(i,j)=(VONKAR/DLOG(hPBL/Z0(int(STYPE(i,j)))))**2
        !
        !
        !- surface albedo
        !         ALBDs(i,j)=albdveg(int(STYPE(i,j)))   !get albedo from sfc type
        !

        !        CDN(i,j)= 1./                                                 &
        !             &      (Log(100./Z0(Int(STYPE(i,j))))/VONKAR)**2

     End Do
  End Do



  Return
End Subroutine bndinit
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
  Character(len=305) retitle
  Integer dayr, monthr, yearr

  u1  =>u1workA
  v1  =>v1workA
  T1  =>T1workA
  q1  =>q1workA

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
!
! ---------------------------------------------------------------------
!
!@@@JWL:  Rename outrestart from out_restart
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
! ---------------------------------------------------------------------
!
Subroutine physics1
  !
  ! Model physics package associated with temperature mode 1:
  !   moist convection, cloud, longwave and shortwave radiation,
  !   surface fluxes of heat,moisture and momentum
  ! 
  Implicit None

  ! print*,'TOP physics1'
  Call mconvct         !moist convection; precipitation
  ! print*,'mconvct finished'
  Call cloud           !cloud prediction
  Call radsw           !shortwave scheme
  Call radlw           !longwave radiation scheme
  Call Sflux           !surface fluxes
  ! call other_flux      !e.g., stochastic forcing; add to rhs

  Return
End Subroutine physics1
!
! ---------------------------------------------------------------------
!
Subroutine mconvct
  !
  ! Moist convection; Betts-Miller scheme projected on temperature mode 1
  !
  Use Prognostic,  Only : T1,q1,nx,ny
  Use Surface,     Only : Qc
#ifdef SPONGES
  Use Sponges,     Only : spngh4
#endif
  Use Qtcmpar
  Use AuxVars
  Use Grid, Only: latt

  Implicit None

  ! --Altered Module Variables:
  ! Qc           !convective heating, also precipitation     W/m**2
  !
  ! --Local Variables & Functions:
  Real CAPE1
  Real dTrefhat
  Real dqrefhat
  Real T1c(nx,ny)
  Integer i,j

  Real x,nlt1c

  dTrefhat=Tcrefhat-Trefhat        ! difference in Tcrhat-Trhat 
  dqrefhat=qcrefhat-qrefhat        ! difference in qcrhat-qrhat

  Do j=1,ny
     Do i=1,nx
#ifdef LINEAR_T1C
        ! linear closure of version V2.1
        !
        !- calculate projected CAPE using qc closure 1
        T1c(i,j)=(a1hat*T1(i,j)+b1hat*q1(i,j)                  &
             &            -dTrefhat-dqrefhat) / (a1hat+bb1hat)
#else
        ! default in V2.2 nonlinear closure
        !
        !- calculate projected CAPE using nonlinear T1c table
        x=Trefhat-Tcrefhat+a1hat*T1(i,j)+qrefhat+b1hat*q1(i,j)
        ! print*,'mconvct: i,j = ',i,j,', x = ',x
        ! call flush(6)
        T1c(i,j)=nlt1c(i,j,x)
#endif

        CAPE1=a1hat*(T1c(i,j)-T1(i,j))+dTrefhat
        CAPE1=Max(CAPE1,0.)            ! reset to zero if negative
        Qc(i,j)=eps_c*CAPE1          & ! convective heating               
#ifdef SPONGES
             &            *Cpg       & ! convert to flux W/m^2
             &            *spngh4(j)   ! sponge bnd.
#else
             &            *Cpg         ! convert to flux W/m^2
#endif
     End Do
  End Do
  Call xfilter(Qc)
  Do j=1,ny
     Do i=1,nx
       Qc(i,j)=Max(Qc(i,j),0.)
     End Do
  End Do
  Return
End Subroutine mconvct
!
! ---------------------------------------------------------------------
!
Subroutine barcl
  !
  ! Baroclinic mode V1; explicit (Euler-backward in a sense)
  !
  Use Prognostic
  Use Baroclinic
  Use Barotropic,        Only : div0
  Use Grid
  Use Qtcmpar
  Use Fluxes,            Only : FSW, FLW, Evap, Qc, FTs, taux, tauy
  Use PhysicalConstants, Only : Rair, gpTi, Hlatent, Cp
  Use Input,             Only : dt
  Use Calendar,          Only : dayofmodel  ! for debugging
  Use Surface,           Only : STYPE
#ifdef SPONGES
  Use Sponges,           Only : spngh1,spngh3
#endif
  Use AuxVars
   
  Implicit None

  ! --Altered module Variables:
  ! u1,v1,T1,q1
  !
  ! --Local Variables & Functions:
  Real Cpgi,a1hati,b1hati 
  Real GMs, t1max,t1avg ! ,GMs2,GMq2
  Real uatv,vatu,rhs
  Real epstau,rdyi
  Real, Dimension(ny) :: rdxui
  !GMq now as GMq1 in defined in qtcm.h and part of barclcom
  Integer i,j
  !     integer iop
  Real, Dimension(:,:), Pointer :: psave
  Real GMsq,qsm,q1m

#ifdef CHIWIND
  Real iru1(nx,ny),irv1(nx,0:ny)
#endif

  Do j=1,ny
     rdxui(j)=Rair*dxui(j)
  Enddo
  rdyi=Rair*dyi
  Cpgi=1./Cpg
  a1hati=1./a1hat
  b1hati=1./b1hat
  epstau=gpTi*V1s/V1sqhat

  !
  != x-momentum equation
  Do j=1,ny
     Do i=1,nx
        vatu=.250*(v1(i,j)+v1(ip1(i),j)+v1(i,j-1)+v1(ip1(i),j-1))
        rhs = -epstau*taux(i,j)                    & !damping
             &   -eps_i1*u1(i,j)                   & ! vertical diffusion 
#ifdef SPONGES
        &    *spngh3(j)                       & !sponge
#endif
        &    +fu(j)*vatu                      & !Coriolis
             &    -rdxui(j)*(T1(ip1(i),j)-T1(i,j)) & !Pressure gradient
             &    +advu1(i,j)                      & !advection of u
             &    +dfsu1(i,j)                           !diffusion of u

        u1(i,j)=u1(i,j)+dt*rhs
     End Do
  End Do


  Call xfilter(u1(:,1:ny))          !spatial filter at high latitudes
  !
  != y-momentum equation
  Do j=1,ny-1
     Do i=1,nx
        uatv=.25*(u1(i,j)+u1(im1(i),j)+u1(i,j+1)+u1(im1(i),j+1))
        rhs=-epstau*tauy(i,j)                     & !surface stress
             &   -eps_i1*v1(i,j)                  & !vertical diffusion 
#ifdef SPONGES
        &   *spngh3(j)                       & !sponge
#endif 
        &    -fv(j)*uatv                     & !Coriolis                 
             &    -rdyi*(T1(i,j+1)-T1(i,j))  & !pressure gradient force  
             &    +advv1(i,j)                & !advection of v           
             &    +dfsv1(i,j)                  !diffusion of v

        v1(i,j)=v1(i,j)+dt*rhs
     End Do
  End Do
  !
  ! boundary conditions in y
  ! pole values are interpolated from points one grid away; these points
  !   are in a sense ghost points; note the sign 
  !     do i=1,nx/2
  !         iop=i+nx/2           !locate the opposite side of the poles
  !
  !         v1(i,0)=half*(v1(i,1)-v1(iop,1))            !South Pole
  !         v1(iop,0)=-v1(i,0)    !point on the other side has opposite sign
  !
  !         v1(i,ny)=.5*(v1(i,ny-1)-v1(iop,ny-1))     !North Pole
  !         v1(iop,ny)=-v1(i,ny)
  !     end do
  !
  ! This is not needed if the curvature term is kept.

  Do i=1,nx
     v1(i,0)=0.
     v1(i,ny)=0.
  Enddo

  Call xfilter(v1(:,1:ny))
  !
  != divergence; diagnostic; curvature term neglected; using u,v at t+1
  Do j=1,ny
     Do i=1,nx
        div1(i,j)=(u1(i,j)-u1(im1(i),j))*dxui(j)     &
       &    +(v1(i,j)*cosv(j)-v1(i,j-1)*cosv(j-1))*dyui(j)
        !not considering the spherical coordinates  
        !     &              +(v1(i,j)-v1(i,j-1))*dyi
        ! mode 2, not implemented 
        !   div2(i,j)=(u2(i,j)-u2(im1(i),j))*dxui(j)+(v2(i,j)*cosv(j)    &
        ! &    -v2(i,j-1)*cosv(j-1))*dyui(j)
        !not considering the spherical coordinates
        !   div2(i,j)=(u2(i,j)-u2(im1(i),j))*dxui(j)+(v2(i,j)-v2(i,j-1))*dyi
     End Do
  End Do
#ifdef CHIWIND
  != solve for the velocity potential for experiments using irrotational winds;
  !
       call fatdfe_neu(dx,dy,cosu,cosv,div1(1,1),chi1)      ! Poisson solver
  Do i=1,nx
     irv1(i,0)=0.
     irv1(i,NY)=0.
     iru1(i,NY)=(chi1(ip1(i),NY)-chi1(i,NY))*dxui(NY)     
     arr2(i,NY)=iru1(i,NY)
     arr3(i,NY)=irv1(i,NY)
  End Do

  Do j=1,ny-1
     Do i=1,nx
        iru1(i,j)=(chi1(ip1(i),j)-chi1(i,j))*dxui(j)     
        irv1(i,j)=(chi1(i,j+1)-chi1(i,j))*dyi
        arr2(i,j)=iru1(i,j)
        arr3(i,j)=irv1(i,j)
     End Do
  End Do

  Do j=1,ny-1
     Do i=1,nx
        arr4(i,j)=(iru1(i,j)-iru1(im1(i),j))*dxui(j)     &
       &    +(irv1(i,j)*cosv(j)-irv1(i,j-1)*cosv(j-1))*dyui(j)
        arr5(i,j)=(irv1(ip1(i),j)-irv1(i,j))*dxvi(j)     &
       &    -(iru1(i,j+1)*cosu(j+1)-iru1(i,j)*cosu(j))*dyvi(j)
     End Do
  End Do
  Do i=1,nx
     arr4(i,ny)=(iru1(i,ny)-iru1(im1(i),ny))*dxui(ny)     &
   &    +(irv1(i,ny)*cosv(ny)-irv1(i,ny-1)*cosv(ny-1))*dyui(ny)
     arr5(i,ny)=0.
  End Do
#endif


  !
  != Temperature equation
  !GMs=GMsr

  GMsq=GMqp
  qsm=19.0*0.001*Hlatent/Cp   ! qsm-the threshold value of surface moisture
  ! (19 g/kg) for cloud-top effect on Gms
  q1m=(qsm-qrefs)/b1s   

  Do j=1,ny
     Do i=1,nx
        !GMs=GMsr+GMsp*max(T1(i,j),-10.0)  !nonlinear correction to GMs
        !GMs=GMsr+GMqp*max(q1(i,j),-6.0)  !nonlinear correction to GMs
        GMs=GMsr+GMsq*Max(q1(i,j),q1m)  !nonlinear correction to GMs
        GMs1(i,j)=GMs
        !         GMs=GMsr+GMsp*T1(i,j)  !nonlinear correction to GMs; @not used
        !        GMs=GMsr
        !          GMs2=GMs2r
        rhs=                             &
             &       -GMs*div1(i,j)      &!adiabatic warming
#ifdef TOPO 
             & !      -GMs0r*div0(i,j)    &!topographic lifting
             & !    2       - GMs2*div2(i,j)
#endif
             &       +(                  &
             &         +Qc(i,j)           &!convective heating
             &         +FSW(i,j)          &!shortwave heating
             &         +FLW(i,j)          &!longwave heating
             &         +FTs(i,j)          &!sensible heat
#ifdef SPONGES
!            &         +other_flux(i,j)   &!other heat fluxes
             &        )*Cpgi             &!flux -> heating rate;
             &         *spngh1(j)         !sponge boundary
#else
!            &         +other_flux(i,j)   &!other heat fluxes
             &        )*Cpgi              !flux -> heating rate;
#endif
        !
        ! Divide by a1hat to account for the vertical averaging/projection
        !   note: this effect is already absorbed in adv (through DTijk) and dfs.
        !        arr3(i,j)=rhs*a1hati
        rhs=rhs*a1hati+advT1(i,j)+dfsT1(i,j)

        T1(i,j)=T1(i,j)+dt*rhs
     End Do
  End Do

  Call xfilter(T1(:,1:ny))
  ! Moisture equation
  Do j=1,ny
     Do i=1,nx
        GMq1(i,j)=GMqr+GMqp*q1(i,j)  !nonlinear correction to GMq
        !         GMq=GMqr
        !          GMq2=GMq2r
        rhs=                                 &
             &       +GMq1(i,j)*div1(i,j)    & !moisture convergence
#ifdef TOPO
        &    !   -GMq0r*Div0(i,j)      & !topographic lifting
        &    !    &       -GMq2*div2(i,j)  &
#endif
             &       +(               &
             &        -Qc(i,j)        &!convective moistening (drying, of course)
             &        +Evap(i,j)      &!moistening due to evaporation
#ifdef SPONGES
!            &        +other_mflux(i,j)  &!other moisture fluxes
             &        )*Cpgi             & !flux -> heating rate
             &          *spngh1(j)         ! sponge boundary
#else
!            &        +other_mflux(i,j)  &!other moisture fluxes
        &        )*Cpgi               !flux -> heating rate
#endif
        !        arr1(i,j)=rhs*b1hati
        rhs=rhs*b1hati+advq1(i,j)+dfsq1(i,j)  !adv,dfs

        q1(i,j)=q1(i,j)+dt*rhs
     End Do
  End Do

  Call xfilter(q1(:,1:ny))          !spatial filter at high latitudes

  Return 
End Subroutine barcl
!
! ---------------------------------------------------------------------
!
Subroutine bartr
  !
  ! Barotropic mode V0
  ! Non-divergent barotropic vorticity equation; vort0 is the only prognostic
  !   variable; streamfunction psi0 is solved with a Poisson solver;
  !   then u0, v0 are diagnosed from psi0; time step can be larger than
  !   baroclinic mode due to the lack of gravity waves, but be careful with
  !   advection
  ! Note: the dynamic effect of topography is implemented and is found to
  !   improve northern hemisphere rainfall and flow pattern significantly.
  !   Commented out here. To use it, uncomment fv*div0 in vort0 equation;
  !   you can neglect the small contribution to wind field (chi0).
  !
  Use Prognostic, Only : u0,v0,u1,v1
  Use Barotropic
  Use Qtcmpar
  Use Grid
  Use Surface,    Only : CV, taux, tauy , STYPE , TOP 
  Use PhysicalConstants, Only : gpTi
  Use Input,      Only : dt,mt0

  Use AuxVars
#ifdef SPONGES
  Use Sponges,     Only : spngh3
#endif


  Implicit None

  ! --Altered Module Variables:
  ! vort0,rhsvort0,u0bar,psi0,u0,v0
  !
  ! --Local Variables & Functions:
  Real zbar,hbar                 ! function
  Real advu0bar,dfsu0bar
  Real V1st,we
  Real adams1,adams2,adams3
  Real tauxbar 
  Integer i,j,k,k_1,k_2,ktemp,kw
  !     REAL res,res2
  !
  Data k/1/,k_1/2/,k_2/3/
  Save k,k_1,k_2

  adams1=23.0/12.
  adams2=-16.0/12.
  adams3=5.0/12.

  !  barotropic vorticity equation
  !  rhs for bartr
  Do j=1,ny-1
     Do i=1,nx
#ifdef TOPO
        !@ use winds at the terrain-following surface (interpolated using V1z)
        we=TOP(i,j)*NZ+1.
        kw=we
        we=we-kw
        V1st=(1.-we)*V1z(kw)+we*V1z(kw+1)
        !
        ! divergence generated by topographic lifting; @not used in standard 
        ! QTCM1 V2.0. div0 centered at (i+1/2,j+1/2), so an average is done; 
        ! however, this makes it smoother (can potentially make significant 
        ! difference, at large topo. gradient and low resolution).

        div0(i,j)=                       &!centered at (i+1/2,j+1/2)
             &   ((u0(i,j)+V1st*u1(i,j))*(TOP(ip1(i),j)-TOP(i,j))        &
             & +(u0(i,j+1)+V1st*u1(i,j+1))*(TOP(ip1(i),j+1)-TOP(i,j+1))) &
             &    *.5*dxvi(j)                                            &
             &  +((v0(i,j)+V1st*v1(i,j))*(TOP(i,j+1)-TOP(i,j))           &
             &    +(v0(ip1(i),j)+V1st*v1(ip1(i),j))*(TOP(ip1(i),j+1)     &
             &      -TOP(ip1(i),j)))*.5*dyi
#endif
        rhsvort0(i,j,k)=                                                 &
#ifdef TOPO           
        ! diverg. due to TOPO
        &   -fv(j)*div0(i,j)                                             &
#endif
        &   -gpTi*(  (tauy(ip1(i),j)-tauy(i,j))*dxvi(j)                  &
        &           -(taux(i,j+1)*cosu(j+1)-taux(i,j)*cosu(j))*dyvi(j) ) &
        &   -(fu(j+1)-fu(j))*dyi*(v0(i,j)+v0(ip1(i),j))*.5               &
        &   +( (advv0(ip1(i),j)-advv0(i,j))*dxvi(j)                      &
        &   -(advu0(i,j+1)*cosu(j+1)-advu0(i,j)*cosu(j))*dyvi(j))        &
        &   +( (dfsv0(ip1(i),j)-dfsv0(i,j))*dxvi(j)                      &
        &   -(dfsu0(i,j+1)*cosu(j+1)-dfsu0(i,j)*cosu(j))*dyvi(j))
        !not considering spherical geometry
        !     &   -gpTi*(  (tauy(ip1(i),j)-tauy(i,j))*dxvi(j)                  &
        !     &     -(taux(i,j+1)   -taux(i,j))*dyi     )                      &
        !     &   -(fu(j+1)-fu(j))*dyi*(v0(i,j)+v0(ip1(i),j))*.5               &
        !     &   +( (advv0(ip1(i),j)-advv0(i,j))*dxvi(j)                      &
        !     &     -(advu0(i,j+1)-advu0(i,j))*dyi )                           &
        !     &   +( (dfsv0(ip1(i),j)-dfsv0(i,j))*dxvi(j)                      &
        !     &     -(dfsu0(i,j+1)-dfsu0(i,j))*dyi )                   

        vort0(i,j)=vort0(i,j)+dt*mt0*(adams1*rhsvort0(i,j,k)              &
             &      +adams2*rhsvort0(i,j,k_1)+adams3*rhsvort0(i,j,k_2))

     End Do
  End Do
  Call xfilter(vort0)

  !
  != advance u0bar(u0 domain average); needed as a bnd. condition for psi0
  !use area-weighted domain average 
  tauxbar=hbar(taux,cosu,nx,ny)

  !no-area weighting considered
  !tauxbar=zbar(taux,nxy)
  !      u2CVbar=zbar2(u2(1,1),CV(1,1),nxy)

  !area-weighted domain average
  advu0bar=hbar(advu0,cosu,nx,ny)
  dfsu0bar=hbar(dfsu0,cosu,nx,ny)

  !no-area weighting considered
  !advu0bar=zbar(advu0(1,1),nxy)
  !dfsu0bar=zbar(dfsu0(1,1),nxy)

  rhsu0bar(k)=                        &
       &    -gpTi*tauxbar             &
       &    +advu0bar                 &
       &    +dfsu0bar
  !    &    -eps_20*u2CVbar           &
  u0bar=u0bar+dt*mt0*(adams1*rhsu0bar(k)                &
       &    +adams2*rhsu0bar(k_1)+adams3*rhsu0bar(k_2))
  !
  != solve for streamfunction psi0 by inverting Laplacian
  ! y boundary: v=0 => psi0=const.    (no x dependence)
  !                        =-u0bar*Dy at north (set psi0(south)=0)
  Do i=1,nx
     psi0(i,ny)=-ny*dy*u0bar     !not ny-1 due to the staggered grid
     psi0(i,0)=0.
  End Do

  Call fatdfe_di(dx,dy,cosu,cosv,vort0(1,1),psi0(1,0))  ! poisson solver

  ! check residual; Note: curvature term has been neglected (Holton 1992).
!!$     res=0.
!!$     res2=0.
!!$     do j=1,ny-1
!!$     do i=1,nx
!!$       res=res+abs(vort0(i,j)                                  &
!!$    &    - (psi0(ip1(i),j)-2.0*psi0(i,j)+psi0(im1(i),j))       &
!!$    &      *dxvi(j)**2                                         &
!!$    &    - (psi0(i,j+1)-2.0*psi0(i,j)+psi0(i,j-1))             &
!!$    &      *dyi**2)                     !curvature term neglected
!!$       res2=res2+abs(vort0(i,j))
!!$     end do
!!$     end do
!!$     write(*,*)'res,res/rhs',res,res/res2
  !====================================================================
  !
  !
  != calculate u0,v0
  Do j=1,ny
     Do i=1,nx
        u0(i,j)=(psi0(i,j-1)-psi0(i,j))*dyi
        !    1           +(chi0(ip1(i),j)-chi0(i,j))*dxui(j)  
        v0(i,j)=(psi0(i,j)-psi0(im1(i),j))*dxvi(j)
        !    1           +(chi0(i,j+1)-chi0(i,j))*dyi
     End Do
  End Do
  !
  ! force the boundary condition
  Do i=1,nx
     v0(i,0)=0.
     v0(i,ny)=0.
  End Do


  ! shift the pointers for rhs by 1. timestep
  ktemp=k_2
  k_2=k_1
  k_1=k
  k=ktemp
  Return
End Subroutine bartr
!
! ---------------------------------------------------------------------
!
Subroutine advctuv
  !
  ! Advection for momentum u,v,u0,v0, ...
  ! Centered second order differencing except duy at N/S boundaries
  !
  ! Array adv* is used to pass variables.  This is useful for
  !   the barotropic mode where the curl of advection is needed later
  !   (because it is solved in streamfunction-vorticity formulation).
  ! Otherwise the procedure can be applied pointwise without saving
  !   in adv*.  Since memory is not a major concern here we do
  !   this for diffusion as well as advection/diffusion for other
  !   modes.  This always looks neater.
  ! Note the three-level looping due to the nonlinear advection and
  !   projection.
  ! Moderately optimized for speed
  !

  Use Prognostic, Only : u0,v0,u1,v1
  Use Baroclinic, Only : div1,advu1,advv1,advwu1 !,advwv1 
  Use Barotropic, Only : div0,advu0,advv0,advwu0 !,advwv0 
  Use Grid
  Use Qtcmpar
  Use AuxVars

  Implicit None

  !
  !
  ! --Altered Module Variables:
  ! advu1,advv1,advu0,advv0
  !
  ! --Local Variables & Functions:
  Real u1atu,v1atu,u1atv,v1atv
  Real u0atu,v0atu,u0atv,v0atv
  Real du1xdxi,du1ydyi,dv1xdxi,dv1ydyi
  Real du0xdxi,du0ydyi,dv0xdxi,dv0ydyi                   
  Real ajp1(ny),ajm1(ny),aj(ny)
  Real vdv1,vdv2,vdv3,vdv4!,Vijkt(8)
  Real divv1,divv2,divv3,divv4!,Vwijkt(8)
  Real dxih,dyih       ! mm speed
  Integer i,j

  dyih=.5*dyi 
  !
  !c- handle 1st-order differencing for duy at y boundaries
  Do j=2,ny-1            !centered differencing in interior
     ajp1(j)=.5       
     aj(j)=0.
     ajm1(j)=-.5
  End Do
  ajp1(1)=1.             !1st order differencing at bnd.
  aj(1)=-1.
  ajm1(1)=0.
  ajp1(ny)=0.
  aj(ny)=1.
  ajm1(ny)=-1.
  !
  !=  x momentum equation; advection for u1, u0
  !
  Do j=1,ny
     dxih=0.5*dxui(j)
     Do i=1,nx
        ! mode V1
        du1xdxi=(u1(ip1(i),j)-u1(im1(i),j))*dxih       !centered differencing
        du1ydyi=(ajp1(j)*u1(i,j+1)+aj(j)*u1(i,j)      &
             &         +ajm1(j)*u1(i,j-1))*dyi
        u1atu=u1(i,j)
        v1atu=.250*(v1(i,j)+v1(ip1(i),j)+v1(i,j-1)    &
             &                     +v1(ip1(i),j-1))
        ! mode V0
        du0xdxi=(u0(ip1(i),j)-u0(im1(i),j))*dxih       !centered differencing
        du0ydyi=(ajp1(j)*u0(i,j+1)+aj(j)*u0(i,j)      &
             &         +ajm1(j)*u0(i,j-1))*dyi
        u0atu=u0(i,j)
        v0atu=.250*(v0(i,j)+v0(ip1(i),j)+v0(i,j-1)    &
             &                     +v0(ip1(i),j-1))
        !
        ! Loop for different modes; 
        ! explicit summation makes the routine [advctuv]
        !   more than twice as fast. However, looping is easier to implement
        !   if you need to add more modes 
        ! (simply change the parameter value NVMOD).
        !       do k=0,NVMOD-1
        !         advu(k)=0.
        !         do n=0,NVMOD-1
        !         do m=0,NVMOD-1
        !           advu(k)=advu(k)-Vijk(m,n,k)
        !    1       *(uatu(m)*dux(n)*dxi+vatu(m)*duy(n)*dyi)
        !    2       -Vwijk(m,n,k)*(div(m)*uatu(n)) !symbolic form w-advection 
        !         end do
        !         end do
        !       end do
        !
        vdv1=u0atu*du0xdxi+v0atu*du0ydyi
        vdv2=u1atu*du0xdxi+v1atu*du0ydyi
        vdv3=u0atu*du1xdxi+v0atu*du1ydyi
        vdv4=u1atu*du1xdxi+v1atu*du1ydyi
        !
        != divergence; curvature term neglected; using u,v at time t
        div1(i,j)=(u1(i,j)-u1(im1(i),j))*dxui(j)     &
             &   +(v1(i,j)*cosv(j)-v1(i,j-1)*cosv(j-1))*dyui(j)
        !not considering spherical coordinates
        !     &              +(v1(i,j)-v1(i,j-1))*dyi
        divv1=0.
        divv3=0.
#ifdef TOPO
        divv1=div0(i,j)*u0atu
        divv3=div0(i,j)*u1atu
#endif
        divv2=div1(i,j)*u0atu      
        divv4=div1(i,j)*u1atu 

!!$        advu0(i,j)=                        &
!!$             &      -Vijkt(1)*vdv1         &
!!$             &      -Vijkt(2)*vdv2         &
!!$             &      -Vijkt(3)*vdv3         &
!!$             &      -Vijkt(4)*vdv4         &
        !             &      -Vijkt(2)*vdv2         &  ! Vijkt(2)=0
        !             &      -Vijkt(3)*vdv3         &  ! Vijkr(3)=0
        !             &      -Vwijkt(1)*divv1       &
        !             &      -Vwijkt(2)*divv2       &
        advu0(i,j)=                        &
             &      -vdv1                   &  ! Vijkt(1)=Vijk(0,0,0)=1.0
             &      -Vijkt(4)*vdv4         &
#ifdef V1V0ADVH
        &       *0.5                  &
#endif
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(3)*divv3       &
             &      -Vwijkt(4)*divv4       &
#ifdef V1V0ADVH
        &       *0.5                  &
#endif
        &;
#endif
        advwu0(i,j)=                       &
             &      -Vwijkt(3)*divv3       &
             &      -Vwijkt(4)*divv4
        !             &      -Vwijkt(1)*divv1       & ! Vwijkt(1)=0
        !             &      -Vwijkt(2)*divv2       & ! Vwijkt(2)=0

        !    1      -Vijk(0,0,0)*vdv1    !this is about 30% slower than above
        !    2      -Vijk(1,0,0)*vdv2
        !    3      -Vijk(0,1,0)*vdv3
        !    4      -Vijk(1,1,0)*vdv4
        !    1      -Vwijk(0,0,0)*divv1
        !    2      -Vwijk(1,0,0)*divv2
        !    3      -Vwijk(0,1,0)*divv3
        !    4      -Vwijk(1,1,0)*divv4
        !
        !             &      -Vijkt(5)*vdv1         & !Vijkt(5)=0
        !             &      -Vijkt(6)*vdv2         & !Vijkt(6)=1
        !             &      -Vijkt(7)*vdv3         &  !Vijkt(7)=1
        !             &      -Vwijkt(5)*divv1       & !Vijkt(5)=0
        !             &      -Vwijkt(6)*divv2       & !Vijkt(6)=0

        advu1(i,j)=                        &
             &      -( vdv2                 &  !Vijkt(6)=1
             &         +vdv3                &  !Vijkt(7)=1
             &         +Vijkt(8)*vdv4 )     &
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4   
#endif

        advwu1(i,j)=                       &
             &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4   
     End Do
  End Do

  !
  !=  y momentum equation; advection for v0, v1
  Do j=1,ny-1
     dxih=0.5*dxvi(j)
     Do i=1,nx
        ! mode V1
        dv1xdxi=(v1(ip1(i),j)-v1(im1(i),j))*dxih
        dv1ydyi=(v1(i,j+1)-v1(i,j-1))*dyih
        u1atv=.250*(u1(i,j)+u1(im1(i),j)+u1(i,j+1)    &
             &                     +u1(im1(i),j+1))
        v1atv=v1(i,j)
        ! mode V0
        dv0xdxi=(v0(ip1(i),j)-v0(im1(i),j))*dxih
        dv0ydyi=(v0(i,j+1)-v0(i,j-1))*dyih
        u0atv=.250*(u0(i,j)+u0(im1(i),j)+u0(i,j+1)    &
             &                     +u0(im1(i),j+1))
        v0atv=v0(i,j)
        !
        vdv1=u0atv*dv0xdxi+v0atv*dv0ydyi
        vdv2=u1atv*dv0xdxi+v1atv*dv0ydyi
        vdv3=u0atv*dv1xdxi+v0atv*dv1ydyi
        vdv4=u1atv*dv1xdxi+v1atv*dv1ydyi
        !
        divv1=0.
        divv3=0.
#ifdef TOPO
        divv1=div0(i,j)*v0atv               !div0 is 0
        divv3=div0(i,j)*v1atv               !div0 is 0
#endif
        divv2=div1(i,j)*v0atv
        divv4=div1(i,j)*v1atv

        !             &      -Vijkt(1)*vdv1          & !   =-vdv1
        !             &      -Vijkt(2)*vdv2          & !   =0
        !             &      -Vijkt(3)*vdv3          & !   =0
        !             &      -Vijkt(3)*vdv3          & !   =0
        !             &      -Vijkt(5)*vdv1         &  ! = 0
        !             &      -Vijkt(6)*vdv2         &  ! = vdv2
        !             &      -Vijkt(7)*vdv3         &  ! = vdv3
        !             &      -Vwijkt(5)*divv1       &  ! = 0
        !             &      -Vwijkt(6)*divv2       &  ! = 0
        advv0(i,j)=                         &
             &      -vdv1          &
             &      -Vijkt(4)*vdv4          &
#ifdef V1V0ADVVH
        &       *0.5                  &
#endif
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(3)*divv3        &
             &      -Vwijkt(4)*divv4        &
#endif
#ifdef V1V0ADVVH
        &       *0.5                  &
#endif
        &;


        advv1(i,j)=                        &
             &      - (  vdv2              &
             &         +vdv3               &
             &         +Vijkt(8)*vdv4 )    & 
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4
#endif

     End Do
  End Do
  Return
End Subroutine advctuv
!
! ---------------------------------------------------------------------
!
Subroutine advctTq
  !
  ! Advection of T and q; structured to allow more modes in the future
  ! Moderately optimized for speed
  !

  Use Prognostic
  Use Baroclinic, Only : advT1,advq1
  Use Grid
  Use Qtcmpar
  Implicit None

  !
  ! --Altered Module Variables:
  ! adv
  !
  ! --Local Variables & Functions:
  Real u0atC,v0atC,u1atC,v1atC
  Real dT1x,dT1y,dq1x,dq1y
  Real ajp1(ny),ajm1(ny),aj(ny) !,dxi
  Integer i,j

  !  dxi=1./dx
  !
  !c- handle 1st-order differencing for dTy,dqy at y boundaries
  Do j=2,ny-1            !centered differencing in interior
     ajp1(j)=.5       
     aj(j)=0.
     ajm1(j)=-.5
  End Do
  ajp1(1)=1.             !1st order differencing at bnd.
  aj(1)=-1.
  ajm1(1)=0.
  ajp1(ny)=0.
  aj(ny)=1.
  ajm1(ny)=-1.
  !
  Do j=1,ny
     Do i=1,nx
        ! advection due to mode V1
        u1atC=(u1(i,j)+u1(im1(i),j))*.5
        v1atC=(v1(i,j)+v1(i,j-1))*.5
        ! mode V0
        u0atC=(u0(i,j)+u0(im1(i),j))*.5
        v0atC=(v0(i,j)+v0(i,j-1))*.5
        !
        ! dTx,dTy; centered differencing
        dT1x=(T1(ip1(i),j)-T1(im1(i),j))*.5
        dT1y=ajp1(j)*T1(i,j+1)+aj(j)*T1(i,j)                          &
             &         +ajm1(j)*T1(i,j-1)
        !
        ! dqx,dqy
        dq1x=(q1(ip1(i),j)-q1(im1(i),j))*.5
        dq1y=ajp1(j)*q1(i,j+1)+aj(j)*q1(i,j)                          &
             &         +ajm1(j)*q1(i,j-1)
        !
        ! sum up contributions from different modes 
        !       do k=1,NTMOD
        !         advT(k)=0.
        !         advq(k)=0.
        !         end do
        !         do n=1,NTMOD
        !         do m=0,NVMOD-1
        !             advT(k)=advT(k)-DTijk(m,n,k)
        !    1               *(uatC(m)*dTx(n)*dxi+vatC(m)*dTy(n)*dyi)
        !             advq(k)=advq(k)-Dqijk(m,n,k)
        !    1               *(uatC(m)*dqx(n)*dxi+vatC(m)*dqy(n)*dyi)
        !         end do
        !         end do
        !       end do
        !
        advT1(i,j)=                                                   &
             &    -DTijkt(1)*(u0atC*dT1x*dxui(j)+v0atC*dT1y*dyi)      &
             &    -DTijkt(2)*(u1atC*dT1x*dxui(j)+v1atC*dT1y*dyi)
        !    1    -DTijk(0,1,1)*(u0atC*dT1x*dxui(j)+v0atC*dT1y*dyi)
        !    2    -DTijk(1,1,1)*(u1atC*dT1x*dxui(j)+v1atC*dT1y*dyi)
        advq1(i,j)=                                                   &
             &    -Dqijkt(1)*(u0atC*dq1x*dxui(j)+v0atC*dq1y*dyi)      &
             &    -Dqijkt(2)*(u1atC*dq1x*dxui(j)+v1atC*dq1y*dyi)
     End Do
  End Do

  Return
End Subroutine advctTq
!
! ---------------------------------------------------------------------
!
Subroutine dffus
  !
  ! Diffusion for u,v,u0,v0,T1,q1
  ! Arrays dfs* are used to pass variables. 
  ! N/S boundary: flux zero (handled by dj(),...)
  !
  Use Prognostic
  Use Baroclinic, Only : dfsu1,dfsv1,dfsT1,dfsq1
  Use Barotropic, Only : dfsu0,dfsv0
  Use Grid
  Use Input, Only : viscxu0,viscyu0,viscxq,viscyq,visc4x,visc4y &
       &            ,viscxu1,viscyu1,viscxT, viscyT
#ifdef SPONGES
  Use Sponges,     Only : spngh2
#endif
  
  Implicit None
#ifdef SPONGES
  Integer :: i,j
#endif

  ! --Altered Module Variables:
  ! dfsu1,dfsv1,dfsu0,dfsv0,dfsT1,dfsq1
  !

  !
  != u1, v1
#ifdef VI2U1
  Call nabla2(viscxU1,viscyU1,u1(:,1:ny),dxui,1,weight2u,dfsu1)
  Call nabla2(viscxU1,viscyU1,v1,dxvi,0,weight2v,dfsv1)
#else
  Call nabla4mm5(viscxU1,viscyU1,u1(:,1:ny),dxui,1,weight4u,dfsu1)
  Call nabla4mm5(viscxU1,viscyU1,v1,dxvi,0,weight4v,dfsv1)
#endif

  !
  != u0, v0
#ifdef VI2U0
  Call nabla2(viscxU0,viscyU0,u0(:,1:ny),dxui,1,weight2u,dfsu0)
  Call nabla2(viscxU0,viscyU0,v0,dxvi,0,weight2v,dfsv0)
#else
  Call nabla4mm5(visc4x,visc4y,u0(:,1:ny),dxui,1,weight4u,dfsu0)
  Call nabla4mm5(visc4x,visc4y,v0,dxvi,0,weight4v,dfsv0)
#endif

  !
  != T1
  ! call nabla4mm5(viscxT,viscyT,T1,dxui,1,weight4u,dfsT1)
  Call nabla2(viscxT,viscyT,T1(:,1:ny),dxui,1,weight2u,dfsT1)

  !
  != q1
  ! call nabla4mm5(viscxq,viscyq,q1,dxui,1,weight4u,dfsq1)
  Call nabla2(viscxq,viscyq,q1(:,1:ny),dxui,1,weight2u,dfsq1)

#ifdef SPONGES
  Do j=1,ny
     Do i=1,nx
        dfsu1(i,j)=dfsu1(i,j)*spngh2(j)
        dfsv1(i,j)=dfsv1(i,j)*spngh2(j)
        dfsu0(i,j)=dfsu0(i,j)*spngh2(j)
        dfsv0(i,j)=dfsv0(i,j)*spngh2(j)
        dfsT1(i,j)=dfsT1(i,j)*spngh2(j)
        dfsq1(i,j)=dfsq1(i,j)*spngh2(j)
     End Do
  End Do
#endif

  Return
End Subroutine dffus
!
! ---------------------------------------------------------------------
!
Subroutine nabla2(viscx,viscy,fld,dxi,js,weight2,d2fld)
  Use Grid, Only: dyi, im1,im2,ip1,ip2, nx,ny
  Implicit None
  ! Inputs
  Integer, Intent(in)                   :: js
  Real, Intent(in)                      :: viscx, viscy
  Real, Dimension(js:ny) , Intent(in)   :: dxi
  Real, Dimension(nx,js:ny), Intent(in) ::  fld
  ! Local
  !  Real, Dimension(nx,js:ny)             :: nabla2
  Real, Dimension(nx,js:ny)             :: d2fld
  Real                                  :: vydy2i,vxdx2i
  Integer                               :: i,j,jsp1
  Real weight2(ny,3)

  vydy2i=viscy*dyi**2

 Do j=js+1,ny-1
     vxdx2i=viscx*dxi(j)**2
     Do i=1,nx
        d2fld(i,j)= vxdx2i*( fld(ip1(i),j)-2.*fld(i,j)+fld(im1(i),j) )   &
#ifdef NO_SPHERDFS
            &     + vydy2i*( fld(i,j+1)   -2.*fld(i,j)+fld(i,j-1)    )
#else
             &    + vydy2i*( fld(i,j+1)*weight2(j,1)+                    &
             &               fld(i,j)*weight2(j,2)+                      &
             &               fld(i,j-1)*weight2(j,3) )
#endif
     Enddo
  Enddo

  ! First and last latitude bands
  jsp1=js+1
  !c- j=js  
  vxdx2i=viscx*dxi(js)**2
  Do i=1,nx
     d2fld(i,js)= vxdx2i*( fld(ip1(i),js)-2.*fld(i,js)+fld(im1(i),js) )  &  
          &     + vydy2i*( fld(i,jsp1)-fld(i,js) )
  Enddo
  !c- j=ny
  vxdx2i=viscx*dxi(ny)**2
  Do i=1,nx
     d2fld(i,ny)= vxdx2i*( fld(ip1(i),ny)-2.*fld(i,ny)+fld(im1(i),ny) )  &  
          &      +vydy2i*( -fld(i,ny)+fld(i,ny-1) )
  Enddo

  Return
End Subroutine nabla2
!
! ---------------------------------------------------------------------
!
Subroutine nabla4(visc,visc2,fld,dxi,js,d4fld)
  Use Grid, Only: dyi, dx, im1,im2,ip1,ip2, nx,ny
  Implicit None
  Integer, Intent(in)                   :: js    ! start index of fld
  Real, Dimension(js:ny) , Intent(in)   :: dxi   ! lon grid distances
  Real, Dimension(nx,js:ny), Intent(in) ::  fld  ! field
  Real,                      Intent(in) ::  visc ! viscocity
  Real,                      Intent(in) ::  visc2 ! viscocity (not used)

  Real, Dimension(nx,js:ny)             :: d4fld 
  Real  :: dy2i,dy4i,dx2i,dx4i, vidx2
  Integer :: i,j,jsp1,jsp2,nym1,nym2


  dy2i=dyi**2
  dy4i=dy2i**2
  Do j=js+2,ny-2
     dx2i=dxi(j)**2
     vidx2=visc/dy2i  ! convert input viscosity it in m^2/s to nabla4 viscosity
     dx4i=dx2i**2
     Do i=1,nx
        d4fld(i,j)= -vidx2                                                 &
             &     *(   dx4i   *( fld(ip2(i),j)+fld(im2(i),j) )            &
             &       +  dy4i   *( fld(i,j+2)   +fld(i,j-2)    )            &
             &       + 2.*(dy2i*dx2i)*( fld(ip1(i),j+1)+fld(ip1(i),j-1)    &
             &                         +fld(im1(i),j+1)+fld(im1(i),j-1) )  &
             &       - 4.*(dy4i+dy2i*dx2i)*( fld(ip1(i),j)+fld(im1(i),j) ) &
             &       - 4.*(dx4i+dy2i*dx2i)*( fld(i,j+1)   +fld(i,j-1) )    &
             &       + ( 6.*(dy4i+dx4i)+8.*dy2i*dx2i) *fld(i,j)  )
     Enddo
  Enddo

  jsp1=js+1
  jsp2=js+2
  nym1=ny-1
  nym2=ny-2
  Do i=1,nx
     d4fld(i,js)=d4fld(i,jsp2)
     d4fld(i,jsp1)=d4fld(i,jsp2)
     d4fld(i,ny)=d4fld(i,nym2)
     d4fld(i,nym1)=d4fld(i,nym2)
  Enddo

  Return
End Subroutine nabla4
!
! ---------------------------------------------------------------------
!
Subroutine nabla4mm5(viscx,viscy,fld,dxi,js,weight4,d4fld)
  Use Grid, Only: dyi, im1,im2,ip1,ip2, nx,ny
  Implicit None
  ! Inputs
  Integer, Intent(in)                   :: js
  Real, Intent(in)                      :: viscx, viscy
  Real, Dimension(js:ny) , Intent(in)   :: dxi
  Real, Dimension(nx,js:ny), Intent(in) ::  fld
  ! Local
  Real, Dimension(nx,js:ny)             :: d4fld
  Real                                  :: vydy2i,vxdx2i
  Integer                               :: i,j,jsp1,jsp2
  Real weight4(ny,5)

  vydy2i=viscy*dyi**2

  Do j=js+2,ny-2
     vxdx2i=viscx*dxi(j)**2
     Do i=1,nx
        d4fld(i,j)= -vxdx2i*( fld(ip2(i),j)+fld(im2(i),j)             &
             &                   -4.*(fld(ip1(i),j)+fld(im1(i),j))    &
             &                  +6.*fld(i,j)  )                       &
#ifdef NO_SPHERDFS 
             &          -vydy2i*( fld(i,j+2)+fld(i,j-2)               &
             &                   -4.*(fld(i,j+1)+fld(i,j-1))          &
             &                   +6.*fld(i,j) )
#else
             &      -vydy2i*( fld(i,j+2)*weight4(j,1)+                &
             &                fld(i,j+1)*weight4(j,2)+                &
             &                fld(i,j)*weight4(j,3)+                  &
             &                fld(i,j-1)*weight4(j,4)+                &
             &                fld(i,j-2)*weight4(j,5))               
#endif
     Enddo
  Enddo

  ! First and last 2 latitude bands
  !- j=js
  vxdx2i=viscx*dxi(js)**2
  jsp1=js+1
  jsp2=js+2
  Do i=1,nx
     d4fld(i,js)= vxdx2i*( fld(ip1(i),js)-2.*fld(i,js)+fld(im1(i),js) )  &
          &     +vydy2i*( fld(i,jsp1)-fld(i,js) )
  Enddo
  !- j=js+1
  vxdx2i=viscx*dxi(jsp1)**2
  Do i=1,nx
     d4fld(i,jsp1)= vxdx2i*( fld(ip1(i),jsp1)-2.*fld(i,jsp1)+fld(im1(i),jsp1) ) &
          &     +vydy2i*( fld(i,jsp2)+fld(i,js)-2.*fld(i,jsp1) )
  Enddo
  !- j=ny-1
  vxdx2i=viscx*dxi(ny-1)**2
  Do i=1,nx
     d4fld(i,ny-1)= vxdx2i*( fld(ip1(i),ny-1)-2.*fld(i,ny-1)+fld(im1(i),ny-1) ) &
          &        +vydy2i*( fld(i,ny)+fld(i,ny-2)-2.*fld(i,ny-1) )
  Enddo
  !- j=ny
  vxdx2i=viscx*dxi(ny)**2
  Do i=1,nx
     d4fld(i,ny)= vxdx2i*( fld(ip1(i),ny)-2.*fld(i,ny)+fld(im1(i),ny) ) &
          &      +vydy2i*( -fld(i,ny)+fld(i,ny-1) )
  Enddo

  Return
End Subroutine nabla4mm5
!
! ---------------------------------------------------------------------
!
Subroutine fatdfe_di(dx,dy,cosu,cosv,rhs,u)
  !
  ! Front-end for the FATD solver
  ! 
  !! Solves 2-dimensional, second-order, periodic (in X) elliptic
  !! problem with coefficients given as functions of Y (constant in X).
  !! with periodicity in X s.t.  u(0,j) <- u(nx,j)  and  u(nx+1,j) <- u(1,j)
  !  Spherical geometry and Dirichlet boundary conditions in Y.
  !
  ! Arguments:
  !  integer         dx,dy        - grid spacing (unchanged)
  !  real(nx,ny)     rhs        - vorticity (unchanged); undefined at rhs(i,ny)
  !  real(nx,0:ny) u      - stream function (returned)
  !                         u(i,0) is used as the southern B.C. (unchanged)
  !                         u(i,ny) is used as the northern B.C. (unchanged)
  !
  ! written: Mar 24, 1997  A.J.Adcroft, UCLA (adcroft@atmos.ucla.edu)
  ! modified: to accomodate the spherical geometry    4/98  N.Zeng
  !
  Use Dimensions
  Implicit None

  ! Arguments
  Real dx,dy
  Real rhs(nx,ny)
  Real u(nx,0:ny)
  Real cosu(ny),cosv(0:ny)
  ! Local
  Integer i,j,nym1
  Real dx2,dy2i,dxc2
  Real a(ny-1),b(ny-1),c(ny-1),d(ny-1)
  Real bet(nx,ny-1),gam(nx,ny-1)
  Real wsave(2*(nx)+15,ny-1)
  Real r(nx,ny-1),uu(nx,ny-1)
  Logical firstcall
  Data firstcall/.True./
  Save nym1,dx2,dy2i,a,bet,gam,wsave,firstcall
  ! -----------------------------------------------------------------------------
  !
  If (firstcall) Then
     !      print *,'Doing FATD initialisation'
     nym1=ny-1
     dx2=dx**2
     dy2i=1.0/dy**2
     Do j=1,nym1
        a(j)=dy2i*cosv(j)**2*dx2
        b(j)=-2.0-2.0*dy2i*cosv(j)**2*dx2
        c(j)=dy2i*cosv(j)**2*dx2
        d(j)=1.0
     Enddo
     !
     ! Dirichlet boundary conditions: u=const.
     a(1)=0.
     c(nym1)=0.
     !
     Call inifatd(nx,nym1,a,b,c,d,bet,gam,wsave)
     !      print *,'fatdpkg: Finished FATD initialisation'
     firstcall=.False.
  Endif
  !
  Do j=1,nym1
     dxc2=dx2*cosv(j)**2
     Do i=1,nx
        r(i,j)=rhs(i,j)*dxc2     !multiply rhs through by (dx*cos)**2
     Enddo
  Enddo
  !
  Do i=1,nx            !mv bnd values to rhs
     r(i,1)=r(i,1)-cosv(1)**2*dx2*dy2i*u(i,0)
     r(i,nym1)=r(i,nym1)-cosv(nym1)**2*dx2*dy2i*u(i,ny)
  Enddo
  !
  Call fatd(nx,nym1,a,bet,gam,wsave,r,uu)
  !     print *,'Returned from FATD'
  !
  Do j=1,nym1
     Do i=1,nx
        u(i,j)=uu(i,j)
     Enddo
  Enddo
  !
  Return
End Subroutine fatdfe_di
!
! ---------------------------------------------------------------------
!
Subroutine fatdfe_neu(dx,dy,cosu,cosv,rhs,u)
! Solve for velocity potential chi, Neuman boundary in y (correspondes to v=0)
!
! Front-end for the FATD solver
!
! This interfaces  2^N+1 sized arrays with the
! 2^N sized arrays used in the FATD solver.
!
! Arguments:
!  integer      dx,dy   - grid spacing (unchanged)
!  real(NX,NY)  rhs     - divergence (unchanged)
!  real(NX,NY)  u       - velocity potential (returned)
!                         u(i,NY) is used as the northern B.C. (unchanged)
!
! written: Mar 24, 1997  A.J.Adcroft, UCLA (adcroft@atmos.ucla.edu)
! Revised: Sam Burns, 2002; Adopted by H. Su, August 2002  
!
      Use Dimensions
      Implicit None
! Arguments
      REAL dx,dy
      REAL cosu(ny),cosv(0:ny)
      REAL rhs(nx,ny)
      REAL u(nx,ny)
! Local
      integer i,j
      REAL rdx2,rdy2
      REAL a(ny),b(ny),c(ny),d(ny)
      REAL bet(nx,ny),gam(nx,ny)
      REAL wsave(2*(nx)+15)
      REAL r(nx,ny),uu(nx,ny)
      logical firstcall
      data firstcall/.true./
      save rdy2,a,bet,gam,wsave,firstcall
! ---------------------------------------------------------------------------
      if (firstcall) then
!       print *,'Doing FATD initialisation'
        rdx2=1.0/(dx**2)
        rdy2=1.0/(dy**2)
        do j=2,ny-1
          a(j)=rdy2*(cosu(j)*dx)**2
          b(j)=-2.0-2.0*rdy2*(cosu(j)*dx)**2
          c(j)=rdy2*(cosu(j)*dx)**2
          d(j)=1.0
        enddo

! Southern boundary condition (Neuman)
        c(1)=rdy2*(cosu(1)*dx)**2
        b(1)=-2.0-rdy2*(cosu(1)*dx)**2
        a(1)=0.0
        d(1)=1.0

! Northern boundary condition (Neuman)
        c(ny)=0.0
        b(ny)=-2.0-rdy2*(cosu(ny)*dx)**2
        a(ny)=rdy2*(cosu(ny)*dx)**2
        d(ny)=1.0

        call inifatd(nx,ny,a,b,c,d,bet,gam,wsave)
!       print *,'fatdpkg: Finished FATD initialisation'
        firstcall=.false.
      endif

      do j=1,ny
        do i=1,nx
          r(i,j)=rhs(i,j)*(cosu(j)*dx)**2  !multply rhs by (cosu(j)*dx)**2
        enddo
      enddo

!     print *,'Calling FATD'
      call fatd(nx,ny,a,bet,gam,wsave,r,uu)
!     print *,'Returned from FATD'
      do j=1,ny
       do i=1,nx
        u(i,j)=uu(i,j)
       enddo
      enddo

      Return
End Subroutine fatdfe_neu 
!
! ---------------------------------------------------------------------
!
Subroutine xfilter(x)
  !
  ! spatial filter of high wave numbers at high latitudes; 
  !   using Arakawa/Lamb filter and FFTPACK (shared with the poisson solver
  !   of Adcroft) 
  !
  Use Dimensions
  Use Grid
  Implicit None
  !
  ! --Local Variables & Functions:
  Real wsave(2*nx+15),x(nx,ny),ara(nx,ny)
  Integer i,j,m,js,jn
  Logical firstcall
  Data firstcall/.True./
  Save ara,firstcall,js,jn,wsave

  If(firstcall) Then
     Call rffti(nx,wsave)       !initialize FFTPACK
     !
     Do j=1,ny                  !low-pass filter of Arakawa/Lamb(1977);
        Do i=2,nx                  !  no filtering on wave-number zero(i=1);
           m=i/2                    !m=wave number
           !ara(i,j)=cosu(j)/Sin(m*pi/nx)*.9     !a reduction of .9
           !ara(i,j)=Cos(latt(j)*pi/180.)/Sin(m*pi/nx)*.9     !a reduction of .9
!less filter at midlats1
           ara(i,j)=0.9*sqrt(cos(latt(j)*pi/180.))/Sin(m*pi/NX)
           ara(i,j)=Min(ara(i,j),1.)
        End Do
     End Do
     !
     !js=(1.-45./YB)*ny/2           !apply the filter outside 45 N/S
     js=(1.-60./YB)*ny/2           !apply the filter outside 60 N/S
     jn=ny-js+1
     !
     firstcall=.False.
  Endif
  !
  !
  !Do j=1,ny
  Do j=1,js
     Call rfftf(nx,x(1,j),wsave)  !Fourier transform; rfftf returns in order
     ! of increasing wave number: 0,1,1,2,2,...
     Do i=2,nx
        x(i,j)=x(i,j)*ara(i,j)   !filtering; not on wavenumber zero
     End Do
     !
     Call rfftb(nx,x(1,j),wsave)  !transform back
  End Do
  Do j=jn,ny
     Call rfftf(nx,x(1,j),wsave)  !Fourier transform; rfftf returns in order
     ! of increasing wave number: 0,1,1,2,2,...
     Do i=2,nx
        x(i,j)=x(i,j)*ara(i,j)   !filtering; not on wavenumber zero
     End Do
     !
     Call rfftb(nx,x(1,j),wsave)  !transform back
  End Do
  !
  Return
End Subroutine xfilter
!
! ---------------------------------------------------------------------
!
Module BartrSave
!
! Saved barotropic velocity 
! Needed for geopotential gradient calculation
!
  Use Dimensions
  Implicit None

  Real, Dimension(nx,ny) :: u0sav, v0sav
End Module BartrSave
!
! ---------------------------------------------------------------------
!
!@@@JWL:  remove underscore from name
Subroutine gradphis
! Compute geopotential gradients using NZ (4.6), (4.10), (4.11) and A1. 
! See NZ appendix A
!
! if macro switch DPHIS1 is set the baroclinic component of the surface
! geopotential is computed separately and stored in arr3.
!
  Use Dimensions
  Use Prognostic,        Only : u0,v0,u1,v1,T1
  Use Barotropic,        Only : advu0,advv0,dfsu0,dfsv0
  Use PhysicalConstants, Only : Rair, gpTi
  !  Use QtcmPar,           Only : a1phat
  Use Grid,              Only : dy,fu,fv,dx,cosu   !, dxui
  Use Input,             Only : dt,mt0
  Use Surface,           Only : taux,tauy, dphisdx,dphisdy,ps    

  Use QtcmPar
  Use AuxVars

  Use BartrSave

  Implicit None

  Real, Dimension(ny),Save  :: radxi
  Real,Save                 :: dti, radyi, ra1phat

  Integer              :: i,j
  Logical,Save              :: firstcall=.True.

  If(firstcall) Then
     firstcall=.False.
     dti=1/(dt*mt0)
     ra1phat=Rair*a1phat
!     ra1phat=Rair*( a1phat-a1phatb )&  ! Mean dphidx,dphidy over subcloud layer, NZ (A1)
     Do j=1,ny
        radxi(j)=ra1phat/dx/cosu(j)    ! *dxui(j)
     Enddo
     radyi=ra1phat/dy
     Return  ! as long as u0sav is not in the restart file.
  Endif
  !
  ! zonal gradient of geopotential  see NZ appendix 
  Do j=1,ny  
     Do i=1,nx-1
        dphisdx(i,j)=-(u0(i,j)-u0sav(i,j))*dti & ! u0 tendency  
             &       +advu0(i,j)               & ! advection
             &       +dfsu0(i,j)               & ! diffusion
             &       +fu(j)*.250*(v0(i,j)+v0(i+1,j)+&
             &                    v0(i,j-1)+v0(i+1,j-1)) & ! Coriolis
             &       -gpTi*taux(i,j)           & ! surface stress  mode 0 NZ(4.6)
             &       -(T1(i+1,j)-T1(i,j))*radxi(j) ! Mode 1 contribution NZ(4.11)
     Enddo
     ! last lon (i=nx)  
     dphisdx(nx,j)=-(u0(nx,j)-u0sav(nx,j))*dti & ! u0 tendency  
          &        +advu0(nx,j)                & ! Advection
          &        +dfsu0(nx,j)                & ! Diffusion
          &        +fu(j)*.250*(v0(nx,j)+v0(1,j)+&
          &                     v0(nx,j-1)+v0(1,j-1)) & ! Coriolis
          &        -gpTi*taux(nx,j)            & ! surface stress mode 0 NZ(4.6)
          &        -(T1(1,j)-T1(nx,j))*radxi(j)  ! mode 1 contribution NZ (4.11)
  Enddo
  !
  ! meridional gradient of geopotential
  !
  ! first lat
  Do i=1,nx-1
     dphisdy(i,1)=-(v0(i,1)-v0sav(i,1))*dti   & ! v0 tendency  
          &       +advv0(i,1)                 & ! advection
          &       +dfsv0(i,1)                 & ! diffusion
          &       -fv(1)*.50*(u0(i,1)+u0(i+1,1)) & ! Coriolis
          &       -gpTi*tauy(i,1)             & ! surface stress mode 0 NZ(4.6)
          &       -(T1(i,1+1)-T1(i,1))*radyi    ! Mode 1 contribution NZ (4.11)
  Enddo
  !last lon (i=nx)
  dphisdy(nx,1)=-(v0(nx,1)-v0sav(nx,1))*dti   & ! v0 tendency  
       &        +advv0(nx,1)                  & ! advection
       &        +dfsv0(nx,1)                  & ! diffusion
       &        -fv(1)*.50*(u0(nx,1)+u0(1,1)) & ! Coriolis
       &        -gpTi*tauy(nx,1)              & ! surface stress  mode 0 NZ (4.6)
       &        -(T1(nx,1+1)-T1(nx,1))*radyi    ! Mode 1 contribution NZ (4.11)
  Do j=2,ny-1 
     Do i=1,nx-1
        dphisdy(i,j)=-(v0(i,j)-v0sav(i,j))*dti  & ! v0 tendency  
             &       +advv0(i,j)                & ! advection
             &       +dfsv0(i,j)                & ! diffusion
             &       -fv(j)*.250*(u0(i,j)+u0(i+1,j)+&
             &           u0(i,j-1)+u0(i+1,j-1)) & ! Coriolis
             &       -gpTi*tauy(i,j)          & ! surface stress mode 0 NZ (4.6)
             &       -(T1(i,j+1)-T1(i,j))*radyi ! Mode 1 contribution NZ (4.11)
     Enddo
     !last lon (i=nx)
     dphisdy(nx,j)=-(v0(nx,j)-v0sav(nx,j))*dti & ! v0 tendency  
          &        +advv0(nx,j)                & ! advection
          &        +dfsv0(nx,j)                & ! diffusion
          &        -fv(j)*.250*(u0(nx,j)+u0(1,j)+&
          &              u0(nx,j-1)+u0(1,j-1)) & ! Coriolis
          &        -gpTi*tauy(nx,j)            & ! surface stress mode 0 NZ (4.6)
          &        -(T1(nx,j+1)-T1(nx,j))*radyi  ! Mode 1 contribution NZ (4.11)
  Enddo
  ! extrapolate to boundary rows
  Do i=1,nx
     dphisdy(i,0)=2*dphisdy(i,1)-dphisdy(i,2)
     dphisdy(i,ny)=2*dphisdy(i,ny-1)-dphisdy(i,ny-2)
  Enddo

  Call dphiint(dphisdx,dphisdy,ps)

End Subroutine gradphis
!
!------------------------------------------------------------------
!
! save barotropic velocities for tendency estimate
!@@@JWL: remove underscore from name
Subroutine Savebartr
  Use Prognostic, Only : u0,v0
  Use BartrSave
  Integer   :: i,j

  Do j=1,ny
     Do i=1,nx
        u0sav(i,j)=u0(i,j)
        v0sav(i,j)=v0(i,j)
     Enddo
  Enddo
  Return
End Subroutine Savebartr
!
!------------------------------------------------------------------
!
! Integration of Grad phi to recover pressure (trapezoidal rule) 
!
Subroutine dphiint(phix,phiy,ps)

  Use Dimensions
  Use Grid
  !  Use BartrSave
  !  Use Surface, only : ps
  Implicit None

  Real, Dimension(nx,ny), Intent(In) :: phix,phiy! , p2
  Real, Dimension(nx,ny),Intent(Out) :: ps
  Real, Dimension(ny)    :: dxuh
  Real                   :: dyh
  Real, Parameter        :: pref=101325  ! [Pa]
  Real, Parameter        :: rhoair=1.2
  Real                   :: avg1, avg2 !

  Integer :: i,j, nyh=ny/2

  Do j=1,ny
     dxuh(j)=0.5*cosu(j)*dx
  Enddo
  dyh=0.5*dy



  ! method (a): first along equator than meridional 
  ps(1,nyh)=0.  ! [Pa]

  Do i=2,nx
     ps(i,nyh)= ps(i-1,nyh)+(phix(i,nyh)+phix(i-1,nyh))*dxuh(nyh)
  Enddo

  Do j=nyh+1,ny
     Do i=1,nx
        ps(i,j)=ps(i,j-1)+(phiy(i,j)+phiy(i,j-1))*dyh
     Enddo
  Enddo
  Do j=nyh-1,1,-1
     Do i=1,nx
        ps(i,j)=ps(i,j+1)-(phiy(i,j)+phiy(i,j+1))*dyh
     Enddo
  Enddo

  ! Normalize:
  avg1=0.
  avg2=0.
  Do j=1,ny
     Do i=1,nx
        avg1=avg1+ps(i,j)
     Enddo
  Enddo

  avg1=avg1/nxy
  avg2=avg2/nxy

  Do j=1,ny
     Do i=1,nx
        ps(i,j)=rhoair*(ps(i,j) -avg1) +pref
     Enddo
  Enddo

  Return
End Subroutine dphiint
!            End of qtcm.F file
!======================================================================
