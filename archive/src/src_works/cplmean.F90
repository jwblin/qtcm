!*********************************************************************
! $Source: /long/munnich/CVSROOT/UCLA/qtcmf90/src/cplmean.f90,v $
! $Author: munnich $
! $Revision: 1.5 $
! $Date: 2002/07/18 20:34:36 $
!
! File: cplmean.f90
!
! Cplmean accumulates coupling variables (fluxes) from the qtcm over
! time.  The qtcm time step is supplied as the argument "it".  If
! "it" equals 1, cplmean clears the flux buffers before accumulating.
! The fluxes are returned as time means by the routines 
! GetStress and GetSfcHeat.
!
! Code History:
!   Original version:  Bill Weibel  Spring 1997
!
!*********************************************************************
!
Subroutine CplMean(it)
  !
  ! Accumulate coupling variables over the length of the coupling time step
  ! 
  Use Surface

  Implicit None

  !
  ! Arguments:
  Integer it  ! qtcm time step
  !
  ! Commons
  Real samples, FLWdscum(Nx,Ny),FLWuscum(Nx,Ny)
  Real FSWdscum(Nx,Ny),FSWuscum(Nx,Ny)
  Real Evapcum(Nx,Ny),FTscum(Nx,Ny)
  Common /cplcount/ samples
  Common /SfcHeat/ FSWdscum,FSWuscum,FLWdscum,FLWuscum,             &
       &                 Evapcum,FTscum
  Real tauxcum(nx,ny),tauycum(nx,ny)
  Real Qccum(nx,ny)
  Common /SfcStress/  tauxcum,tauycum
  Common /WaterFlux/  Qccum
  Save /cplcount/, /SfcHeat/, /SfcStress/
  !
  ! Local variables
  Integer i,j
  !
  ! Sum the surface fluxes
  !   
  If (it .Eq. 1) Then
     Call as(FSWdscum,Nx*Ny,0.)
     Call as(FSWuscum,Nx*Ny,0.)
     Call as(FLWdscum,Nx*Ny,0.)
     Call as(FLWuscum,Nx*Ny,0.)
     Call as(Evapcum,Nx*Ny,0.)
     Call as(FTscum,Nx*Ny,0.)
     Call as(tauxcum,Nx*Ny,0.)
     Call as(tauycum,Nx*Ny,0.)
     Call as(Qccum,Nx*Ny,0.)
     samples = 0.0
  End If
  !
  samples = samples + 1.0
  Do j = 1,Ny
     Do i=1,Nx
        FSWdscum(i,j) = FSWdscum(i,j) + FSWds(i,j)
        FSWuscum(i,j) = FSWuscum(i,j) + FSWus(i,j)
        FLWdscum(i,j) = FLWdscum(i,j) + FLWds(i,j)
        FLWuscum(i,j) = FLWuscum(i,j) + FLWus(i,j)
        Evapcum(i,j) = Evapcum(i,j) + Evap(i,j)
        FTscum(i,j) = FTscum(i,j) + FTs(i,j)
        tauxcum(i,j) = tauxcum(i,j) + taux(i,j)
        tauycum(i,j) = tauycum(i,j) + tauy(i,j)
        Qccum(i,j) = Qccum(i,j) + Qc(i,j)
     End Do
  End Do
  !
  Return
End Subroutine CplMean
!
!
Subroutine getStress(taux_out,tauy_out)
  !
  ! Coupling interface routine.
  ! Return time-averaged wind stress up through argument list
  !
  Use dimensions
  Implicit None
  !
  ! --Local Variables:
  Real taux_out(nx,ny),tauy_out(nx,ny)
  Integer i,j
  ! -- Globals
  Real divisor, samples
  Common /cplcount/ samples
  Real tauxcum(nx,ny),tauycum(nx,ny)
  Common /SfcStress/  tauxcum,tauycum
  !*****************************************************************
  !
  divisor = 1.0/samples
  Do j=1,ny
     Do i=1,nx
        taux_out(i,j) = tauxcum(i,j) * divisor
        tauy_out(i,j) = tauycum(i,j) * divisor
     End Do
  End Do
  !
  Return
End Subroutine getStress
!
Subroutine getPrec(Qc)
  !
  ! Coupling interface routine.
  ! Return time-averaged precipitation up through argument list
  ! Evaporation can be obtained via getSfcHeat
  ! NOTE: Unit are energy unit
  !
  Use Dimensions

  Implicit None

  !
  ! --Local Variables:
  Real Qc(nx,ny)
  Integer i,j
  ! -- Globals
  Real divisor, samples
  Common /cplcount/ samples
  Real Qccum(nx,ny)
  Common /WaterFlux/  Qccum
  !*****************************************************************
  !
  divisor = 1.0/samples
  Do j=1,ny
     Do i=1,nx
        Qc(i,j) = Qccum(i,j) * divisor
     End Do
  End Do
  !
  Return
End Subroutine getPrec
!
Subroutine getSfcHeat(FSWds,FSWus,FLWds,FLWus,Evap,FTs)
  !
  ! GetSfcHeat is an interface for Ocean-Atmosphere coupling.
  ! There are no inputs.  The fluxes, appearing as arguments, are synonymous
  ! with variables defined in qtcm.h.  But, in this case, they are averages
  ! over some number of time steps.  The number of steps is given by 'samples'
  ! in the common block 'SfcHeat'
  !
  ! Bugs:
  !     The driver calls Ocean before QTCM.  So this routine 
  !     may be called before the flux variables are defined.
  !
  Use Dimensions

  Implicit None

  Integer i,j
  Real    :: divisor, samples=1.0 ! avoid divide by 0.:
  Real FLWdscum(Nx,Ny),FLWuscum(Nx,Ny)
  Real FSWdscum(Nx,Ny),FSWuscum(Nx,Ny)
  Real Evapcum(Nx,Ny),FTscum(Nx,Ny)
  Common /cplcount/ samples
  Common /SfcHeat/ FSWdscum,FSWuscum,FLWdscum,FLWuscum,             &
       &                 Evapcum,FTscum
  Real FSWds(Nx,Ny),FSWus(Nx,Ny)
  Real FLWds(Nx,Ny),FLWus(Nx,Ny)
  Real Evap(Nx,Ny),FTs(Nx,Ny)
  !
  divisor = 1.0/samples
  Do j=1,Ny
     Do i=1,Nx
        FSWds(i,j) = FSWdscum(i,j) * divisor
        FSWus(i,j) = FSWuscum(i,j) * divisor
        FLWds(i,j) = FLWdscum(i,j) * divisor
        FLWus(i,j) = FLWuscum(i,j) * divisor
        Evap(i,j) = Evapcum(i,j) * divisor
        FTs(i,j) = FTscum(i,j) * divisor
     End Do
  End Do
  Return
End Subroutine getSfcHeat
