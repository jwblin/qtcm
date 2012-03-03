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

!************************************************
!
Subroutine sland1
  !************
  ! Simple-LAND(SLAND) version 1B           Zeng and Neelin (1998)
  ! In the following implementation, sensible heat flux and a 'swamp' evaporation
  !   (evaporation like over ocean except using the actual land roughness)
  !   are assumed given (calculated in QTCM1 routine [sflux]).
  !
  Use Surface
  Use PhysicalConstants, Only : Hlatent
#ifdef SPONGES
  Use Sponges,           Only : spngh0
#endif
  Use AuxVars
  Use Land
  Use Input,             Only : dt

  Implicit None

  !
  ! --Altered Common Variables:
  ! Ts,qs,WD,Runf                   over land only
  !
  ! --Local Variables & Functions:
  Real Rsnet,Fsnet,FWnet
  Real Evapi0,tau_0,Fitc
  Real eta,wra,soilC,Rung,tiny , wet4
  Real Wimax,rints,tau_r,Rung0
#ifdef SPONGES
  Real eps_Ts
#endif
  ! Integer Clapp       !@integer for speed; modify if the value is real
  Integer i,j,iS
  !
  Data tiny/0.0000000001/   !tiny to avoid division by zero
  !
  !- non-surface type dependent parameters
  Data Wimax/0.1/     !max intercepted water per leaf area; kg/m2 (mm)
  Data rints/1.06e-3/ !storm intensity, 3.8mm/hr(from ARME)->mm/s  
  Data tau_r/4320./   !storm duration,  1.2hrs ->s
  Data Rung0/4.e-4/   !subsurface runoff at saturation; kg/m2(mm)/s
  ! Data Clapp/4/       !Clapp-Hornberger B
  Save Wimax,rints,tau_r,Rung0 ! ,Clapp
  !
  !*****************************************************************
  !  time scale to relax Ts to Ts0 at sponge boundary (2 day here)
#ifdef SPONGES
  eps_Ts=1./(86400.*.2)
#endif
  !
  !  Soil heat Capacity * density * Depth; assume .1m water-like soil
  soilC=4.18e3*1.e3*.1       !1cal/g/K * 1g/cm3 * .1m;  J/KM2
  !
  Do j=1,ny
     Do i=1,nx
        iS=STYPE(i,j)
        If(iS.Ne.0.) Then
           !
           != Diagnostics;  surface hydrology
           !
           !- Interception loss; stochastic rainfall effects included
           Rsnet= FSWds(i,j)-FSWus(i,j)+FLWds(i,j)-FLWus(i,j)  !net sfc. radiation
           !
           Evapi0=Rsnet-FTs(i,j)        !simply set Evapi0 to available energy
           !
           Evapi0=Max(tiny,Evapi0)    !minimum of tiny if too small
           !
           tau_0=Wimax*xla(iS)/Evapi0   &!time to evaporate a saturated  can
                &        *Hlatent                !  in seconds (Evapi0 in W/m2)
           !
           Fitc=(tau_r+tau_0*.8)        &!interception function F2 of Zeng (&
                &    *Qc(i,j)/(Hlatent*rints*tau_r)        !  Qc in W/m2
           !
           Evapi(i,j)=Evapi0*Fitc            !interception loss

           Evapi(i,j)=min(Evapi(i,j),0.5*Qc(i,j)) ! cap Evapi


           !           arr1(i,j)=Evapi
           !
           !---------------------------------------------
           !
           !- evaporation; aerodynamic resistance ra from [sflux]
           !  stomatal resistance rs=rsmin/wet**alpha with: 
           !   alpha=1 - linear dependence 
           !   alpha<1 - mimic deep roots water usage such as during dry season  
           !   alpha=0 - not limited by water, but by a mimimum stomatal resistance
           !
           wet(i,j)=WD(i,j)/WD0(iS)           !wetness

           !        arr3(i,j)=wet

           !       wra=one                       !alpha=0 - no wetness dependence
           !       wra=wet                       !alpha=1
           !       wra=sqrt(wet)                !alpha=1/2
           wra=Sqrt(Sqrt(wet(i,j)))         !alpha=1/4
           !
           wra=wra/CV(i,j)               ! CV=CDN*wind=1/ra
           eta=wra/(rsmin(iS)+wra)       !eta=ra/(rs+ra)

           Evap(i,j)=eta*Evap(i,j)       !actual ET = eta * E(swamp)
           !E(swamp) was calculated in [sflux]
           !
           Evap(i,j)=Evap(i,j)+Evapi(i,j)     !total Evap = ET + Interception loss
           !---------------------------------------------
           !
           !- runoff; note: interception is not available for runoff
           !
           wet4=wet(i,j)**4
           Runs(i,j)=(Qc(i,j)-Evapi(i,j))*wet4   !surface runoff, BATS formulation
           !
           !       Rung=Hlatent*Rung0*wet**(2*Clapp+3)   !subsurface runoff; in W/m2
           !       Rung=Hlatent*Rung0*(wet**2)**Clapp*wet**3   !optimization for speed
           Rung=Hlatent*Rung0*(wet4**2)*wet(i,j)**3   !optimization for speed
           !
           Runf(i,j)=Runs(i,j)+Rung           !total runoff; in W/m2
           !        arr4(i,j)=Runs
           !=========================================================
           !
           != Prognostic equations
           !
           !c- soil moisture equation
           !
           FWnet=(Qc(i,j)-Evap(i,j)-Runf(i,j))/Hlatent  !net water flux; mm/s
           !
           WD(i,j)=WD(i,j)+dt*FWnet
           !
           WD(i,j)=Max(WD(i,j),0.)      !zero if negative (due to numerics)
           !
           !- An initially saturated surface will be quickly drained out to a subsaturated
           !    wetness (typically .8 for tropical rainforest) by the power law
           !    subsurface drainage term; So supersaturation rarely occurs
           !       if(WD(i,j)/WD0(iS).gt.1.00) write(*,*) 'warning! wet='
           !    &   ,WD(i,j)/WD0(iS),i,j
           !---------------------------------------------
           !
           !- ground temperature; 10cm water-like soil gives a 20min/K damping rate
           !    with typical flux; this mimics the surface soil layer response. 
           Fsnet= Rsnet-Evap(i,j)-FTs(i,j)      !net energy absorbed by surface
           !
           !
#ifdef SPONGES
           Ts(i,j)=Ts(i,j)+dt*Fsnet/soilC                   &
                &  +dt*eps_Ts*(Ts0(i,j)-Ts(i,j))*spngh0(j)   !relax to Ts0
#else
           Ts(i,j)=Ts(i,j)+dt*Fsnet/soilC                   
#endif
        End If
     End Do
  End Do
  !
  Return
End Subroutine sland1
!
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
