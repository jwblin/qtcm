!
! Module definition for module clrad for Chia Chou cloud-radiation package
!
! Code history:
!   version 2.0 (QTCM1):  Chou/Neelin  June 1997
!   version 2.3: cloud parameterization modified 
!                by J. D. Neelin, H. Su,  July 2002  
!
!**************************************************************************

Module Clrad
  ! for Chia Chou cloud-radiation package
  Use Dimensions

  Integer, Parameter              :: clt=4, cltm1=clt-1 !cloud types
  Real, Dimension(clt), Parameter :: cldref = &
       &    (/0.1051000,      0.1047000,     0.1096000,     0.2234000/)

  Real                            :: cosZ !solar zenith angle
  Real, Dimension(0:clt,nx,ny)    :: cld  !radiatively active cloud cover
  Real, Dimension(nx,ny)          :: cl1  !cloud 1 output array
#ifdef DELTA_CO2
  Real, Dimension(nx,ny)          :: co2  !carbon dioxide
#endif
End Module Clrad


!====== end file ======
