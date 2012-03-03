! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: qtcmutil.F90
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
! Some routines called but not in this file:
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
Subroutine SaveBartr
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
End Subroutine SaveBartr
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
!
!------------------------------------------------------------------
!
