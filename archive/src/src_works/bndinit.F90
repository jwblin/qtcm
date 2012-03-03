! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: bndinit.F90
!
! Code history:
!   Original (QTCM1 Version 1a):  Ning Zeng, David Neelin, et al. Jan. 1997
!   Beta version (QTCM1 V1b):     June 1997
!   QTCM1 release (QTCMes1 V2.0):   June 1998
!   QTCM1 V2.1 :                  Jan  1999
!   QTCM1 V2.2f90 :              Nov 2000, Matthias Munnich
!
! Some routines possibly called but not in this file:
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
  Character(len=130) :: fname
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


!======================================================================
!            End of bndinit.F90 file
!======================================================================
