! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File:  parinit.F90
!
! Code history:
!   Original (QTCM1 Version 1a):  Ning Zeng, David Neelin, et al. Jan. 1997
!   Beta version (QTCM1 V1b):     June 1997
!   QTCM1 release (QTCMes1 V2.0):   June 1998
!   QTCM1 V2.1 :                  Jan  1999
!   QTCM1 V2.2f90 :              Nov 2000, Matthias Munnich
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


!======================================================================
!            End of parinit.F90 file
!======================================================================
