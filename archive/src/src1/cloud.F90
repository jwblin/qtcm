!
! Subroutine cloud from the cloud and radiation package of Chou and 
! Neelin (Chou, Ph.D Dissertation, UCLA 1997)
!
! Code history:
!   version 2.0 (QTCM1):  Chou/Neelin  June 1997
!   version 2.3: cloud parameterization modified 
!                by J. D. Neelin, H. Su,  July 2002  
!
!**************************************************************************

Subroutine cloud
  !
  ! cloud prediction scheme
  !
  Use Fluxes
  Use Clrad
  Use Input, Only : bnddir
  Use AuxVars
  Implicit None
  !
  !
  !
  ! --Input Variables:
  !   Qc, cldref
  !
  ! --output Variables:
  !   cldref: reference (mean) cloud cover
  !   cld: radiatively active cloud cover after overlap calculation
  !   cldtot: total cloud fraction of each type (i.e. before overlap calculation)
  !                   type 0: clear sky
  !                   type 1: Deep cloud +CsCc
  !                   type 2: Cirrus
  !                   type 3: Stratus
  !                   type 4: AsAc+CuSc (approx constant in space and time)
  !   Note: Cloud cover of cloud type 4 varies only due to cloud overlap
  !
  ! --Local Variables & Functions:
  Integer i,j
  Character(Len=120) :: clddir
  Real  :: cl1P=7.76275869E-4              ! coeff. for cloud prediction scheme
  Real  :: cl2fac=1.5                      ! ratio of cloud type 2 to cloud type 1, 
  Real, Dimension(clt,nx,ny), Save :: cldtot ! Total cloud fraction for each type 

#ifdef COUNTCAP          
  !count times when cld(1) cap is reached 
  Integer, Dimension(nx,ny)  :: icount = 0
#endif
  !
  !*****************************************************************
  !
#ifdef OBSCLD
  !     read ISCCP observed cloudiness
  clddir=Trim(bnddir)//'/CLOUD_ISCCP' 
  ! note: cldtot(1), cldtot(2) are parameterized below, 
  ! cldtot(4) not from data so currently only  cldtot(3) used 
  ! NB: reads only once per month
  Call readobscloud(cldtot,clddir)
  ! assign preliminary values of cloud 4, to be modified below
  ! by overlapping cloud 1 and 2
  Do j=1,ny
     Do i=1,nx
        cldtot(4,i,j)=min(cldref(4)/(1-cldref(1)-cldref(2)),1-cldtot(3,i,j))
     End Do
  End Do
  
#else
  ! assign preliminary values of cloud 3 and 4, to be modified below
  ! by overlapping cloud 1 and 2 
  Do j=1,ny
     Do i=1,nx
        cldtot(3,i,j)=cldref(3)/(1-cldref(1)-cldref(2))
        cldtot(4,i,j)=min(cldref(4)/(1-cldref(1)-cldref(2)),1-cldtot(3,i,j))
     End Do
  End Do
#endif
  ! 
  !Parameterization of cloud fractions:
  !cld(1,i,j) deep+CsCc: depends on deep convective precip
  !approx proportional for typical values, capped at 1.
  !cld(2,i,j) cirrus: For typical values proportional to cld(1) 
  !with factor (cl2fac) based on ISCCP data. Modified by
  !cloud overlap condition that ensures total cloud fraction .lt. 1 
  !Overlap region of cld(1) with cld(2) uses cld(1) properties, 
  !hence (1-cld(1)) factor (equiv to random overlap assumption)
  !For cld(3) and cld(4) overlap with either cld(1) or cld(2)
  !uses cld(1)/cld(2) properties hence (1-cld(1)-cld(2)).
  !If using spatially constant cld(3)=cldref(3)/(1-cldref(1)-cldref(2)),
  !helps keep spatial mean similar to cldref(3); likewise for cld(4).
  !
  !For this paramzn sum_k cld(k).le.1 if cld(1).le.1, cldtot(2).le.1
  !and cldtot(3)+cldtot(4).le.1. For latter condition,
  !reduce cldtot(4) if cldtot(3)+cldtot(4).gt.1
  !
  !cld(0) clear sky: as residual. 
  !
  !
  Do j=1,NY
     Do i=1,NX
        cldtot(1,i,j)=min(cl1P*Qc(i,j),1.)   !depends on precip
        cld(1,i,j)=cldtot(1,i,j)  !when cld(1) overlaps it dominates
#ifdef COUNTCAP
        icount(i,j)=0
        if(cl1P*Qc(i,j).gt.1.) then
           icount(i,j)=1
           arr1(i,j)=icount(i,j)
        endif
#endif
        cldtot(2,i,j)=min(cld(1,i,j)*cl2fac,1.)   !before overlap calculation
        cld(2,i,j)=cldtot(2,i,j)*(1.-cld(1,i,j))  !non-overlapped by cld(1)
        cld(3,i,j)=cldtot(3,i,j)*(1.-cld(1,i,j)-cld(2,i,j)) 
        cld(4,i,j)=cldtot(4,i,j)*(1.-cld(1,i,j)-cld(2,i,j)) 
        cld(0,i,j)=1.-cld(1,i,j)-cld(2,i,j)-cld(3,i,j)-cld(4,i,j) 
        cl1(i,j)=cld(1,i,j)    !for output
!!$        arr2(i,j)=cld(2,i,j)
        arr3(i,j)=cld(3,i,j)
!!$        arr4(i,j)=cld(4,i,j)
     End Do
  End Do
  Return
End Subroutine cloud


!====== end file ======
