! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/qtcm.f90,v $
! $Author: munnich $
! $Revision: 1.7 $
! $Date: 2002/07/19 01:33:35 $
!
! File: mconvct.F90
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


!======================================================================
!            End of mconvct.F90 file
!======================================================================
