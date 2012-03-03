! 
! Modules for QTCM1 
!
! ---------------------------------------------------------------------
!
! horizontal grid dimensions 
Module Dimensions
  Integer, Parameter :: nx=64      & ! long. grid size
       &              , ny=42      & ! latit. grid size
       &              , nxy=nx*ny    ! array size
End Module Dimensions
!
! ---------------------------------------------------------------------
!
Module Prognostic
  Use Dimensions
  ! 1st baroclinic:
  !  Storage arrays for 2 time slices for the baroclinic prognostics
  ! One could use storage allocation, but this was easier
  Real, Dimension(nx,0:ny+1),Target   ::  u1workA,u1workB  & ! [m/s]
       &                               ,  T1workA,T1workB  & ! [K]
       &                               ,  q1workA,q1workB    ! [K]
  Real, Dimension(nx,0:ny),Target :: v1workA,v1workB      ! [m/s]
  ! u1,T1,q,v1 point to field at current time step
  Real,Dimension(:,:),  Pointer  :: u1,    & ! [m/s]
       &                            T1,    & ! [K]
       &                            q1       ! [K]
  Real, Dimension(:,:), Pointer  :: v1       ! [m/s]
  ! barotropic:  (actual prognostics are vort0, u0bar
  Real, Dimension(nx,0:ny+1)   :: u0      ! [m/s]
  Real, Dimension(nx,0:ny), Target :: v0      ! [m/s]



  ! Note: surface prognostics 
  !       Ts, WD are in the Surface and Land modules, respectively

End Module Prognostic
!
! ---------------------------------------------------------------------
!
! baroclinic mode
Module Baroclinic
  Use  Dimensions
  Real, Dimension(nx,ny)   :: div1   &
       &                    , chi1   &
       &                    , advu1  &
       &                    , advT1  &
       &                    , advq1  &
       &                    , dfsu1  &
       &                    , dfsT1  &
       &                    , dfsq1  &
       &                    , GMq1   &
       &                    , advwu1 &
       &                    , advwv1 &
       &                    , GMs1
  Real, Dimension(nx,0:ny) :: advv1  &
       &                    , dfsv1
End Module Baroclinic
!
! ---------------------------------------------------------------------
!
! barotropic mode V0; note only vort0 and u0bar are prognostic
Module Barotropic
  Use  Dimensions
  Real, Dimension(nx,0:ny), Target :: psi0    
  Real, Dimension(nx,ny)      :: vort0  &
       &                       , div0       &
       &                       , chi0   &
       &                       , advu0  &
       &                       , dfsu0  &
       &                       , advwu0 &
       &                       , advwv0
  Real, Dimension(nx,0:ny)    :: advv0 &
       &                       , dfsv0 
  Real u0bar,rhsu0bar(3),rhsvort0(nx,ny,3)  
End Module Barotropic
!
! ---------------------------------------------------------------------
!
! surface fluxes
Module SurfaceFluxes
  Use  Dimensions
  Real, Dimension(nx,ny) :: Qc          & !prec.
       &                  , FLWds,FLWus & !down/up LW fluxes at surface
       &                  , FSWds,FSWus & !down/up SW fluxes at surface
       &                  , Evap,FTs    & !lat./sens heat flux
       &                  , taux,tauy     !surface stress 
End Module SurfaceFluxes
!
! ---------------------------------------------------------------------
!
! heat and moisture fluxes
Module Fluxes
  Use  SurfaceFluxes
  Real, Dimension(nx,ny) :: FLWut, FLW &!LW at top; LW column total
       &                  , S0         &!incoming solar at top     
       &                  , FSWut,FSW   !SW at top; SW absorbed by 
End Module Fluxes
!
! ---------------------------------------------------------------------
!
!switch #if def MODE2
!switch Module Mode2
!switch   ! ekman   @not used in this version
!switch   !
!switch   Use  Dimensions
!switch   Real u2(nx,ny),v2(nx,ny),div2(nx,ny),rhsu2(nx,ny,3),rhsv2(nx,ny,3)
!switch End Module Mode2
!switch #endif
!
! ---------------------------------------------------------------------
!
! auxiliary arrays for output
Module AuxVars
  Use Dimensions
  Real, Dimension(nx,ny) ::arr1,arr2,arr3,arr4,arr5,arr6,arr7,arr8
  Character(len=130) :: arr1name='?', arr2name='?'  &
       &, arr3name='?', arr4name='?'  &
       &, arr5name='?', arr6name='?'  &
       &, arr7name='?', arr8name='?'
End Module AuxVars
!
! ---------------------------------------------------------------------
!
! surface variables. Ts is prognostic over land.
Module Surface
  Use  SurfaceFluxes
  Real, Dimension(nx,ny) ::     & ! All on T-grid
       &                    Ts  & ! surface temperature [K]
       &                  , Ts0 & ! surface temperature [K]
       &                  , stype & ! STYPE real for output [-]
       &                  , CDN & ! neutral drag coefficient [-]
       &                  , CV  & ! CDN * [effective wind velocity] [m/s]
       &                  , ALBDs & ! surface albedo [-]
       &                  , TOP & ! topography [1000m]
       &                  , ps  & ! surface pressure [Pa]
       &                  , us  & ! surface zonal velocity [m/s]
       &                  , vs  & ! surface merid. velocity [m/s]
       &                  , ub  & ! top of mixed layer zonal velocity [m/s]
       &                  , vb  & ! top of mixed layer merid. velocity [m/s]
       &                  , dphisdx  ! surface geopotential zonal  grad [m/s^2]
  Real, Dimension(nx,0:ny) :: dphisdy ! surface geopotential merid. grad [m/s^2]
End Module Surface
!
! ---------------------------------------------------------------------
!
! constants
Module PhysicalConstants
  Real, Parameter ::          &
       &    Rearth=6.37d6     & ! earth radius           [m]
       &  , rhoair=1.20       & ! air density            [kg/m^3]
       &  , Rair=287.040      & ! gas constant air       [J/kg/K]
       &  , Omega=7.292d-5    & ! earth angular velocity [1/s]
       &  , Hlatent=2.43d6    & ! latent heat            [J/kg]
       &  , Cp=1.004d3        & ! heat capacity          [J/kg/K] 
       &  , gravity=9.8       & ! gravitational acceler. [m/s^2]
       &  , delp=80000        & ! total pressure diff    [Pa]
       &  , gpTi = gravity/delp !  NZ $g/p_T$            [m/s^2/Pa]
End Module PhysicalConstants
!
! ---------------------------------------------------------------------
!
! grid layout
Module Grid
  Use  Dimensions
  Real pi                     ! =2*acos(0)
  Real :: dx,dy, dxi,dyi      ! grid sizes and inverse at equator [m],[1/m]
  Real, Dimension(0:ny) ::  & ! at v-points:
       &   fv               & ! Coriolis parameter v-points [1/s]
       & , cosv             & ! cos latitude v-points
       & , cosvi            & ! 1/cosv
       !       & , dxv              & ! zonal grid distance [m]
  & , dxvi,dyvi               ! 1/(cosv*dxv),1/(cosv*dyv)   [1/m]
  Real, Dimension(ny)   ::  & ! at u-points
       &   fu               & ! Coriolis parameter [1/s]
       & , cosu             & ! cos latitude u-points
       & , cosui            & ! 1/cosu
       !       & , dxu              & ! zonal grid distance [m]
  & , dxui,dyui               ! 1/(cosu*dxu),1/(cosu*dyu)   [1/m]
  Real, Dimension(ny):: latt ! latitude of T-points
  Real, Dimension(nx):: lont ! longitude of T-points
  Real, Dimension(ny,3):: weight2u,weight2v ! 2nd-order diffusion weight factor
  ! for spherical geometry
  Real, Dimension(ny,5):: weight4u,weight4v ! 4th-order diffusion weight factor
  ! for spherical geometry
  Real, Parameter :: YB = 78.75 !domain spans YB S - YB N [degree]
  !
  ! indexing for finite differencing at the periodic x boundary
  Integer im1(nx),ip1(nx),im2(nx),ip2(nx)
End Module Grid
!
! ---------------------------------------------------------------------
!
#ifdef SPONGES
Module Sponges
  ! sponge boundary for heat flux, diffusion, momentum 
  !
  Use  Dimensions
  Real spngh0(ny),spngh1(ny),spngh2(ny),spngh3(ny) !sponge boundary
  Real spngh4(ny)
End Module Sponges
#endif
!
! ---------------------------------------------------------------------
!
!    user input via driver.in
Module Input
  Real dt,viscxu0,viscyu0,viscxu1,viscyu1   !atmo timestep; del2 viscosity
  Real visc4x,visc4y                  !del4 viscosity
  Real viscxT,viscyT,viscxq,viscyq    ! T,q diffusivities 
  Real V1b                            ! top of ABL V1 projection [-]
  Real weml, ziml, VVsmin         ! for mixed layer ABL (mlabl.F90)
  Integer it ! time of day in time steps
  Integer noout,nooutr,ntout,ntouti,ntoutr,mrestart,landon,nastep,mt0
  Integer lastday, interval
  Character(len=130) bnddir,outdir
  Character(len=130) title,runname
  !
  !- for ocean (mixed layer, data ocean)
  Character(len=130) SSTdir
  Character(len=20) SSTmode
End Module Input
!
! =======================================================================
