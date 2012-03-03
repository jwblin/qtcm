! *****************************************************************************
!  Code history:
!     Hui Su oct. 1998
!     Matthias Munnich may 1999 adjusted call to bndry1 to use bnddir
!
! *****************************************************************************
Subroutine getbnd
  !  Get boundary data such as SST, albedo, surface fluxes for atmosphere model 
  !
  Use Surface, Only : Ts, Ts0, STYPE, ALBDs, nx, ny
  Use Input,   Only : bnddir, landon

  Implicit None
  Integer       :: i,j,lnblnk
  Integer, Save :: ncall=0
  Integer       :: nwet

  If(ncall==0) Then
     ncall=1 ! initially use Ts from restart file   
  Else
     Call getSST(Ts0)         ! get data SST
     Do j=1,ny
        Do i=1,nx
           If(STYPE(i,j)==0.) Ts(i,j)=Ts0(i,j) !update Ts over ocean
           ! faked SST works only for 
           ! interpolated ocean SST over land, e.g., Reynolds SST 
           If((landon==0).And.(stype(i,j)>0.)) Ts(i,j)=Ts0(i,j)
        Enddo
     Enddo
  Endif
  !
  !- get the climatological sfc albedo; if not called, will use albedo 
  !  according to surface type (see [bndinit])
  !  Note: other prescribed clim. fields can be done similarly
  Call bndry1('alb',bnddir(1:lnblnk(bnddir))//'/ALBD_Darnell',ALBDs)
  Return
End Subroutine getbnd
!
! =============================================================================

