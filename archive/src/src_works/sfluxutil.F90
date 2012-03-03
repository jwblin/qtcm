!=========================================================================
! sflux utilities
!=========================================================================


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


!====== end file ======
