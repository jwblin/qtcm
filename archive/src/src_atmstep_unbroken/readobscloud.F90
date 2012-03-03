!
! From the cloud and radiation package of Chou and Neelin (Chou, Ph.D
!   Dissertation, UCLA 1997)
!
! Code history:
!   version 2.0 (QTCM1):  Chou/Neelin  June 1997
!   version 2.3: cloud parameterization modified 
!                by J. D. Neelin, H. Su,  July 2002  
!
!**************************************************************************

Subroutine readobscloud(cld, dir)
  !
  ! Read CLOUD COVER data files from ISCCP C2 (monthly)
  !
  Use Dimensions
  Use Clrad,    Only : clt
  Use Calendar
  Use AuxVars
  Implicit None
  Integer, Parameter :: nyOBS=32
  Integer, Parameter :: jskip=(ny-nyOBS)/2

  Real, Dimension(clt,nx,ny), Intent(Out) :: cld
  Character(Len=*),           Intent(In)  :: dir
  !
  !
  !- Local Variables and Functions
  Integer lnblnk, month_prior, n, i, j
  Character(Len=2)   :: monthc
  Character(Len=12)  :: file
  Character(len=130) :: fname
  Real :: cldmn
  !
  Data month_prior/0/
  Save month_prior
  !
  If(monthofyear.Eq.month_prior) Return
  !  jskip=(ny-nystd)/2
  !
  Write(monthc,'(i2.2)') monthofyear       !convert into a string
  file='0000'//monthc//'15'//'.cld'     !file format: 0000mm15.cld
  !
  ! read the var
  fname=dir(1:lnblnk(dir))//'/'//file
  Open(13,file=fname,form='formatted',status='old')
  Do n=1,3
     Do j=1+jskip,ny-jskip
        Do i=1,nx
           Read(13,*) cld(n,i,j)
        End Do
     End Do
     ! Compute zonal mean for fill in
     cldmn=0.
     Do i=1,nx
        cldmn=cldmn+cld(n,i,jskip+1)
     Enddo
     cldmn=cldmn/nx
     ! Fill in southern region
     Do j=1,jskip
        Do i=1,nx
           cld(n,i,j)=cldmn
        End Do
     End Do
     ! Compute zonal mean for fill in
     cldmn=0.
     Do i=1,nx
        cldmn=cldmn+cld(n,i,ny-jskip)
     Enddo
     cldmn=cldmn/nx
     Do j=ny-jskip+1,ny
        Do i=1,nx
           cld(n,i,j)=cldmn
        End Do
     End Do

  End Do
  Close(13)

  Write (*,*) 'readobscld: Read cloud fraction at model date', &
       &  dateofmodel,' from file:'
  Write (*,*) dir(1:lnblnk(dir))//'/'//file

  month_prior=monthofyear

!!$  do j=1,ny
!!$     do i=1,nx
!!$        arr5(i,j)=cld(1,i,j)
!!$        arr6(i,j)=cld(2,i,j)
!!$        arr7(i,j)=cld(3,i,j)
!!$     enddo
!!$  enddo

  Return
End Subroutine readobscloud


!====== end file ======
