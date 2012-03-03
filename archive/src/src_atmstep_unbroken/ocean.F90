! File: ocean.F90
!  
! This file is just the ocean subroutine that then calls other ways of  
!   obtaining SST:  by reading in observed SST, or using a mixed-layer
!   ocean model to compute SST, or combining observed SST in masked regions
!   and mixed-layer SST in non-masked regions       
!
! Code history:
!   For readsst: Original version:  Bill Weibel  Nov 1996
!              Partial revision:  Ning Zeng    Jan 1997
!
!   For mxlocean: Original version:  H. Su/N. Zeng Nov. 1998
!              
!   For blendsst: Original version:  H. Su July 2002  
!
! ---------------------------------------------------------------------
!
      Subroutine Ocean(ndays,it)
!
  Use Dimensions
  Use SSTout 

      implicit none

      integer, intent(in) :: ndays,it

#if defined BLEND_SST
! BLEND_SST and MXL_OCEAN are exclusive macros
#undef MXL_OCEAN

      Real, Dimension(NX,NY) :: Tsobs,Tsmxl 
      
      call mxlocean(ndays,it)
      Tsmxl = Tnow
      call readsst(ndays,it)
      Tsobs = Tnow 
!
! combines observed SST in masked regions and mixed-layer ocean in
! non-masked regions
!
      call blendsst(Tsobs,Tsmxl)

#elif defined MXL_OCEAN
      call mxlocean(ndays,it)      

#else 
!
! By default, read in observed SST
!
      call readsst(ndays,it)
#endif      

      Return
End subroutine Ocean
!
!=====================================================================
