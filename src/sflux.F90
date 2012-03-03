Subroutine Sflux
  ! Surface fluxes of sensible heat, moisture, and momentum;
  ! Computes:
  ! CV=1/ra, ra=1./(CDN(i,j)*VVsE(i,j))        aerodynamic resistance, for land
  ! FTs       !sensible heat                     W/m**2
  ! Evap      !evaporation                       W/m**2
  ! taux      !surface stress along x            N/m**2 (=10 dyn/cm**2)
  ! tauy      !surface stress along y            N/m**2 (=10 dyn/cm**2)
  ! 
  Use Surface,           Only : CV,FTs,Evap,taux,tauy,us,vs &
       &                      ,Ts,STYPE,CDN ! not changed
  Use Prognostic,        Only : T1,q1,nx,ny
  Use PhysicalConstants, Only : Cp, Rhoair, Hlatent, gpTi
  Use Qtcmpar,           Only : a1s,b1s,qrefs,Trefs

  Implicit None

  Integer  :: i,j
  Real     :: HlCpi,rhoCp ! derived constants for speedup
  Real     :: hsat        ! hsat - function
#ifdef NO_ABL
  Real, Dimension(nx,ny) :: VVs & ! effective wind speed for surface stress
       & , VVsE  ! NO_ABL effective wind speed for Evap and FTs
#else
  Real, Dimension(nx,ny), Save :: VVs  ! effective wind spd for sfc fluxes
#endif

#ifdef NO_ABL
  Call SfcWindNoABL(us,vs,VVs,VVsE)
#else
  ! Momentum ABL
  Call SfcWindABL(us,vs,VVs) 
#endif

  Do j=1,ny
     Do i=1,nx
        CV(i,j)=CDN(i,j)*VVs(i,j)
     Enddo
  Enddo
  !
  ! surface stress
  ! taux on u-grid, tauy on v-grid
  Do j=1,ny
     Do i=1,nx-1
        taux(i,j)=rhoair*(CV(i,j)+CV(i+1,j))*0.5*us(i,j)
     End Do
     taux(nx,j)=rhoair*(CV(nx,j)+CV(1,j))*0.5*us(i,j)
  End Do
  Do j=1,ny-1 ! tauy(:,ny) is never used.
     Do i=1,nx
        tauy(i,j)=rhoair*(CV(i,j)+CV(i,j+1))*0.5*vs(i,j)
     End Do
  End Do

#ifdef NO_ABL
  !  use a different the effective wind for Evap and FTs
  Do j=1,ny
     Do i=1,nx
        CV(i,j)=CDN(i,j)*VVsE(i,j)
     Enddo
  Enddo
#endif

  !
  !- evaporation; pseudo-potential evaporation over land;    W/m^2
  rhoCp=rhoair*Cp
  HlCpi=Hlatent/Cp
  Do j=1,ny
     Do i=1,nx
        Evap(i,j)=rhoCp*CV(i,j)           & !Cp: q is in Kelvin
             &     *( hsat(Ts(i,j))*HlCpi & !saturation humidity at surface
             &        -qrefs-b1s*q1(i,j) )  ! - surface humidity 
     End Do
  End Do
  !
  !- sensible heat flux;  W/m^2
  Do j=1,ny
     Do i=1,nx
        FTs(i,j)=rhoCp*CV(i,j)*(Ts(i,j)-Trefs-a1s*T1(i,j))
        ! Limit the heat flux into the ocean to 5 W/m^2.
        if  (STYPE(i,j)==0.) FTs(i,j)=max(FTs(i,j),-5.)
     End Do
  End Do

  Return
End Subroutine Sflux
!
!-----------------------------------------------------------
!
Subroutine SfcWindNoABL(us,vs,VVs,VVsE)
  ! Note: the wind speed for evaporation/sensible heat is allowed to be
  !   computed differently from that for momentum; without a
  !   PBL parameterization q,T at surface are not simulated accurately enough.
  !   This leads to somewhat too sensitive dependence of evaporation on wind
  !   speed.  The evaporation-wind feedback, therefore MJO is also quite
  !   sensitive to this parameterization.
  ! 
  Use Prognostic,  Only : nx,ny,u0,v0,u1,v1
  Use Qtcmpar,     Only : V1s

  Implicit None
  Real, Dimension(nx,ny), Intent(out) :: us,vs,VVs,VVsE
  Real, Parameter ::  &
       eta=0.6        &
  & , VVsEmin=5.0     & ! minimum surface wind speed for Evap and FTs
  & , VVsmin=4.0      & ! minimum surface wind speed for momentum 
  & , V1sE=-0.17        ! winds proj. coeff. for Evap (QTCM1V2.2 value)
  Integer         :: i,j
  Real            :: VVsEsq,u1_T,u0_T,v1_T,v0_T

  Do j=1,ny
     ! i=1
     u1_T=(u1(1,j)+u1(nx,j))*0.5
     u0_T=(u0(1,j)+u0(nx,j))*0.5
     v1_T=(v1(1,j)+v1(1,j-1))*0.5
     v0_T=(v0(1,j)+v0(1,j-1))*0.5
     !
     ! zonal and meridional surface winds on T-grid
     us(1,j)=u1_T*V1s+u0_T
     vs(1,j)=v1_T*V1s+v0_T
     !
     ! effective surface wind speed for Evaporation and Sensible heat
     VVsEsq=(u1_T*V1sE+u0_T)**2+(v1_T*V1sE+v0_T)**2
     VVsE(1,j)=Sqrt(VVsEmin**2+eta**2*VVsEsq)
     !
     ! effective surface wind speed for Stress 
     VVs(1,j)=Sqrt(us(1,j)**2+vs(1,j)**2+VVsmin**2)
     Do i=2,nx
        u1_T=(u1(i,j)+u1(i-1,j))*0.5
        u0_T=(u0(i,j)+u0(i-1,j))*0.5
        v1_T=(v1(i,j)+v1(i,j-1))*0.5
        v0_T=(v0(i,j)+v0(i,j-1))*0.5
        !
        ! zonal and meridional surface winds on T-grid
        us(i,j)=u1_T*V1s+u0_T
        vs(i,j)=v1_T*V1s+v0_T
        !
        ! effective surface wind speed for Evaporation and Sensible heat
        VVsEsq=(u1_T*V1sE+u0_T)**2+(v1_T*V1sE+v0_T)**2
        VVsE(i,j)=Sqrt(VVsEmin**2+eta**2*VVsEsq)
        !
        ! effective surface wind speed for Stress 
        VVs(i,j)=Sqrt(us(i,j)**2+vs(i,j)**2+VVsmin**2)
     End Do
  End Do

  Return
End Subroutine SfcWindNoABL
!
!===========================================================
