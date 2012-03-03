! Single soil layer (thin (.1m) for Ts; large field capacity for water);
!   a major difference between SLAND1 and the bucket model is that SLAND1 
!   has an evapotranspiration limited by a bulk surface resistance.
!   Other features include: 
!   interception loss is computed using an 'interception function' which
!   incorporates the stochastic property of rainfall in time/space;
!   soft runoff (fast and slow components);

Module Land

  Use Dimensions

  Implicit None
  Real, Dimension(nx,ny), Target  :: Evapi, WD, Runf, wet, Runs
  Integer, Dimension(nxy)     :: iland, jland
  !
  !            surface type; ocean or vegetation type over land
  Integer, Parameter          :: NSTYPE=4   

  Real, Dimension(0:NSTYPE-1) ::         &
       &  rsmin    =  (/0., 150., 200.,200./) &!min bulk surface resistance (BATS)
       & , albdveg =  (/.07, .12, .19,.30 /)  &!veg. albedo (BATS)
       & , Z0      =  (/.0024,2.,.1,.05/)     &!roughness length (BATS)
       & , xla     =  (/0.,6.,3.,1./)         &!Leaf Area Index (BATS)
       & , WD0     =  (/0.,500.,400.,300./)    !field capacity SIB2/CSU (approximately)
  !defined as the part above wilting point
End Module Land
!
!------------------------------------------------------------------------
!
Subroutine bucket
  !************
  ! An implementation of the bucket model (not identical to Manabe et al. 1965).
  ! Land-surface; Bucket model hydrology except for a 'softer'
  ! runoff scheme.  May need to be solved implicitly when soil
  ! depth is thin
  !
  Use Surface
  Use Land
  Use Input, Only : dt
  Use PhysicalConstants, Only : Hlatent

  Implicit None
  !
  ! --Altered Common Variables:
  !
  ! --Local Variables & Functions:
  Real WD00,soilC,beta
  Integer i,j
  !*****************************************************************
  !
  !- Soil heat Capacity * Depth; assume .1m water-like soil
  soilC=4.18e3*1.e3*.1       !1cal/g/K * 1g/cm3 * .1m;  J/KM2
  !
  Do j=1,ny
     Do i=1,nx
        If(Int(STYPE(i,j)).Ne.0) Then
           !c- hydrology  !W/m**2->mm(depth of water;=kg/m2)/s
           WD00=150.              !field capacity=150mm everywhere
           wet(i,j)=WD(i,j)/WD00
           Runf(i,j)=wet(i,j)**4*Qc(i,j)      !depending on P only
           beta=wet(i,j)
           Evap(i,j)=beta*Evap(i,j)       !actual ET = beta * PotentialE
           WD(i,j)=WD(i,j)+dt*(Qc(i,j)-Evap(i,j)-Runf(i,j))/Hlatent
           !c- ground temperature; note conversion factor Cpg
           Ts(i,j)=Ts(i,j)+ dt*(FSWds(i,j)-FSWus(i,j)+FLWds(i,j)           &
                &            -FLWus(i,j)-Evap(i,j)-FTs(i,j))/soilC
        End If
     End Do
  End Do

  Return
End Subroutine bucket
