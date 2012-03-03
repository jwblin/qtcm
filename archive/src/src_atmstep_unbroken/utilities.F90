Subroutine TimeInterp(Prior, Upcoming, Result ,t1,t2,now, nx,ny)
  !
  !  A generic time interpolator.
  !
  Implicit None
  Integer t1,t2
  Real now
  Integer nx,ny
  Real Prior(nx,ny),Upcoming(nx,ny),Result(nx,ny)
  Real divisor,fraction
  !
  Integer i,j
  !
  divisor = t2 - t1
  If ( divisor .Gt. 0.0 ) Then
     fraction = (now - t1)/divisor
     Do j=1,ny
        Do i=1,nx
           Result(i,j)=Prior(i,j)+fraction*(Upcoming(i,j)-Prior(i,j))
        End Do
     End Do
  Else
     !        Result = Prior
     Do j=1,ny
        Do i=1,nx
           Result(i,j) = Prior(i,j)
        End Do
     End Do
  Endif
  Return
End Subroutine TimeInterp
!
!
Subroutine humtable
  !************
  ! Set up look-up table for saturation humidity every .1K
  ! between 200K and 400K
  !
  Implicit None
  !
  ! --Buried Common Variables:
  Real hTmin,hTmax,hdT,hdTi,hums
  Common/hsatcom/ hTmin,hTmax,hdT,hdTi,hums(2001)
  !
  ! --Local Variables & Functions:
  Real T
  Real hsat0         ! function
  Integer NT
  Integer k
  !*****************************************************************
  !
  hTmin=200.
  hTmax=400.
  hdT=0.1
  hdTi=1./hdT
  NT=(hTmax-hTmin)/hdT+1
  Do k=1,NT
     T=hTmin+hdT*(k-1)
     hums(k)=hsat0(T)
  End Do
  !
  Return
End Subroutine humtable
!
!
!
Real Function hsat(T)
  !************
  ! Calculates hsat using a table. Linear interpolation
  ! used at points in-between.
  !
  Implicit None
  !
  ! --Arguments:
  Real T
  !
  ! --Buried Common Variables:
  Real hTmin,hTmax,hdT,hdTi,hums
  Common/hsatcom/ hTmin,hTmax,hdT,hdTi,hums(2001)
  !
  ! --Local Variables & Functions:
  Real w
  Integer k
  !*****************************************************************
  !
  If(T.Lt.hTmin) Then
     !       write(*,*) 'terror !!!! T below 200K, hsat has been set to 200K'
     T=hTmin
  End If
  If(T.Gt.hTmax) Then
     !       write(*,*) 'terror !!!! T above 400, hsat has been set to 400K'
     T=hTmax-hdT
  End If
  w=1.0+(T-hTmin)*hdTi
  k=w
  w=w-k
  hsat=(1.0-w)*hums(k)+w*hums(k+1)
  !
  Return
End Function hsat
!
! 
Subroutine t1ctable
  !************
  ! Set up look-up table for nonlinear T1c
  ! T1c between -300 and 300 at the interval of 1.
  !
  Use T1cTableIn
  Use PhysicalConstants, only : Hlatent, Cp
  Use Qtcmpar, only : a1hat
  Implicit None
  !
  ! --Buried Common Variables:
  Real T1cs(601),rhsx(601)
  Common/T1cscom/ T1cs,rhsx
  !
  ! --Local Variables & Functions:
  !  Integer NP
  Real prs(np),alpha(np),Tcref(np),a1(np)  ! np vertical layers
  Real qcp(np),qcphat,hsat,Cpi ! , hlatent,a1hat,
  !
  Integer i,k
  !*****************************************************************
  !  Open(31,file='vtpar.in',form='formatted',status='old')
  !  NP=11
  print*, 'tqctable: prs(k),alpha(k),Tcref(k),a1(k)'
  Do k=1,NP
     !     Read(31,*) prs(k),alpha(k),Tcref(k),a1(k)
     ! table is formally a 1D array for easier initialization in qtcmpar.f90
     prs(k)=table(1+4*(k-1))
     alpha(k)=table(2+4*(k-1))
     Tcref(k)=table(3+4*(k-1))
     a1(k)=table(4+4*(k-1))
  print*,prs(k),alpha(k),Tcref(k),a1(k)
  Enddo
  !  Read(31,*) a1hat
  !  Close(31)
  !
  !  hlatent=2.43d6
  Cpi=1./Cp
  !
  Do i=1,601
     T1cs(i)=-300.+(i-1)*1.
     Do k=1,NP
! old        qcp(k)=alpha(k)*hsat(Tcref(k)+a1(k)*T1cs(i))*Hlatent*Cpi
        qcp(k)=alpha(k)*hsat(Tcref(k)+a1(k)*T1cs(i))     & 
     &         *Hlatent*Cpi*1.0e3/prs(k)

     Enddo
     qcphat=0.0
     Do k=1,NP-1
        qcphat=qcphat+(qcp(k)+qcp(k+1))*0.5*(prs(k)-prs(k+1))
     Enddo
     qcphat=qcphat/(prs(1)-prs(NP))
     rhsx(i)=a1hat*T1cs(i)+qcphat
  Enddo
  Return
End Subroutine t1ctable
!
Real Function nlt1c(i,j,x)
  !************
  ! Calculates nonlinear T1c using a table. Linear interpolation
  ! used at points in-between.
  !
  Use dimensions
  Use Grid, Only: lont,latt ! For blow-up info
  Implicit None
  !
  ! --Arguments:
  Integer i,j
  Real x 
  !
  ! --Buried Common Variables:
  Real T1cs(601),rhsx(601),alpha
  Common/T1cscom/ T1cs,rhsx
  !
  ! --Local Variables & Functions:
  Integer k,ii,jj
  Integer k0(nx,ny)   ! guess of neighboring
  Data k0 /nxy*0/
  Save k0  
  !*****************************************************************
  !
  If(x.Lt.rhsx(1)) Then
     Write(*,*) 'terror !!!! T1c < -300. extropolate T1c from -300.'
     Write(*,*) 'x=',x
     nlt1c=T1cs(1)+(T1cs(2)-T1cs(1))/(rhsx(2)-rhsx(1))*               &
          &       (x-rhsx(1)) 
     Write(*,*) 'T1cs=',nlt1c
     Call FatalInfo(i,j)
     Stop
  End If
  If(x.Gt.rhsx(601)) Then
     Write(*,*) 'terror !!!! T1c > 300. extropolate T1c from 300.'
     Write(*,*) 'x=',x
     nlt1c=T1cs(601)+(T1cs(601)-T1cs(600))/(rhsx(601)-rhsx(600))*     &
          &       (x-rhsx(601)) 
     Write(*,*) 'T1cs=',nlt1c
     Call FatalInfo(i,j)
     Stop
  End If
  !
  Call embrace(rhsx,601,x,k0(i,j))
  !
  k=k0(i,j)
  alpha=(x-rhsx(k))/(rhsx(k+1)-rhsx(k))
  nlt1c=T1cs(k+1)*alpha+T1cs(k)*(1.0-alpha)
  !
  Return
End Function nlt1c
!
Subroutine FatalInfo(i,j)
  Use dimensions
  Use Grid, Only: lont,latt ! For blow-up info
  Use Prognostic, Only: T1
  Use Surface, Only: stype
  Use Input, Only: ntouti, ntout, noout
  Implicit None
  Integer, intent(in) :: i,j
  ! --Local Variables & Functions:
  Integer ii,jj
  Write(*,*) 'Location:     (i,j) = (',i,',',j,')'
  Write(*,*) '          (lon,lat) = (',lont(i),',',latt(j),')'
  !Write(*,*) 'Surface type neighborhood:'
  !Write(*,'(3(1x,i1.1))') ((ifix(stype(ii,jj)),ii=max(i-2,1),min(i+2,nx)) &
  !                          & ,jj=min(j+1,ny),max(j-1,1),-1) 
  !Write(*,'(3(1x,f6.0))') ((T1(ii,jj),ii=i-2,i+2),ii=max(i-2,1),min(i+2,nx)) &
  !                          & ,jj=min(j+1,ny),max(j-1,1),-1) 
  ! write current stage  output to instaneous output
  ntouti=1
  ! write mean up to this point to mean output.
  if(ntout/=0) ntout=1
  noout=0
  call outpAll
  Return
End Subroutine FatalInfo
!
Subroutine embrace(xx,n,x,k0)
  !************
  !  search for embracing k0 with xx(k0)< x < xx(k0+1)
  !  k0 on input as initial guess. 
  !  Based on Num. Rec. "hunt" routine,
  !  "xx" must be increasing.
  !
  Implicit None
  ! --Arguments:
  Integer n  ! dimension of xx
  Integer k0 ! Input: Initial guess,
  ! Output lower embracing index
  Real xx(n) ! lookup table
  Real x     ! value to embrace
  !
  ! --Local Variables & Functions:
  Integer k1,k2,inc,bisect
  !*****************************************************************
  !
  If(k0.Le.0) Then ! no initial guess -> do bisection search
     k0=bisect(xx,n,x) 
     Return
  Else If (xx(k0).Le. x) Then !search upward
     inc=1
     k1=k0+inc
     Do While (xx(k1).Le. x) 
        k0=k1
        inc=2*inc
        k1=Min(k1+inc,n)
     Enddo
  Else If (xx(k0).Gt.x) Then !search downward
     inc=1
     k1=k0
     k0=k1-inc
     Do While (xx(k0).Gt.x)
        inc=2*inc
        k1=k0
        k0=Max(k0-inc,1)
     Enddo
  Endif
  Do While (k1-k0.Gt.1) !bisection search
     k2=(k1+k0)/2
     If(x.Gt.xx(k2)) Then
        k0=k2
     Else
        k1=k2
     Endif
  Enddo
  Return
End Subroutine embrace
!
!
Integer Function bisect(xx,n,x)
   !
   !
   !  bisection search for embracing 
   !  xx(bisect) < x < xx(bisect+1)
   !
   Implicit None
   Integer n, k1, k2
   Real xx(n),x
   !*****************************************************************
   bisect=0
   k1=n+1
   Do While (k1-bisect.Gt.1)
      k2=(k1+bisect)/2
      If(x.Gt.xx(k2)) Then
         bisect=k2
      Else
         k1=k2
      Endif
   Enddo
   Return
End Function bisect
!
!
!
Real Function hsat0(T1)
  !************
  ! Saturation specific humidity at 1000 mb; empirical formula from
  !  J. Shuttleworth applicable to normal range of temperature.
  !
  Implicit None
  !
  ! --Arguments:
  Real T1      !input T1 in Kelvin
  !
  ! --Local Variables & Functions:
  Real T,esat
  !*****************************************************************
  !
  T=T1-273.15                           ! convert to Celsius
  esat=.6108*Exp(17.270*T/(237.3+T))  ! unit: kpa
  hsat0=.622*esat*0.01                  ! unit= kg/kg, assuming
  !                                           ! Pressure=100kpa
  !
  Return
End Function hsat0
!
Real Function hsats(T)
  !************
  ! Slope of saturation humidity; unit of hsat/K
  Implicit None
  !
  ! --Arguments:
  Real T      !input T in Kelvin
  !
  ! --local
  Real hsat
  !*****************************************************************
  !
  hsats = 10.0*(hsat(T+.050)-hsat(T-.050))
  !
  Return
End Function hsats
!
Real Function zbar(a,n)
  !
  !************
  ! Mean of an array of size n
  !
  Implicit None
  !
  ! --Arguments:
  !     zbar          ! function:  mean of array a(n)
  Integer n     ! array dimension
  Real a(n)   ! array.  not altered on exit.
  !
  ! --Local Variables & Functions:
  Integer i
  !*****************************************************************
  !
  zbar=a(1)
  Do i=2,n
     zbar=zbar+a(i)
  End Do
  zbar=zbar/n
  !
  Return
End Function zbar
!
!
Real Function zbar2(a,b,n)
  !************
  ! Mean of arrays a*b, each of size n
  !
  Implicit None
  !
  ! --Arguments:
  Integer n
  Real a(n),b(n)
  !
  ! --Local Variables & Functions:
  Integer i
  !*****************************************************************
  !
  !
  zbar2=a(1)*b(1)
  Do i=2,n
     zbar2=zbar2+a(i)*b(i)
  End Do
  zbar2=zbar2/n
  !
  Return
End Function zbar2
!
!
!
Real Function hbar(a,wgh,nx,ny)
  !
  !************
  ! area-weighted domain mean of an array of size (nx,ny)
  ! weights are wgh(ny) 
  Implicit None
  !
  ! --Arguments:
  !     hbar              ! function:  mean of array a(n)
  Integer nx,ny           ! array dimension
  Real a(nx,ny),wgh(ny)   ! array.  not altered on exit.
  !
  ! --Local Variables & Functions:
  Integer i,j
  !*****************************************************************
  !
  hbar=0.
  Do j=1,NY
  Do i=1,NX
     hbar=hbar+a(i,j)*wgh(j)
  End Do
  End Do
  hbar=hbar/(nx*ny)
  !
  Return
End Function hbar
!
!
!
Subroutine as(a,n,val)
  !************
  ! Assign a single value val to array a of length n; no need if FORTRAN90
  !
  Implicit None
  !
  ! --Arguments:
  Integer n
  Real a(n),val
  !
  ! --Local Variables & Functions:
  Integer i
  !*****************************************************************
  !
  Do i=1,n
     a(i)=val
  End Do
  !
  Return
End Subroutine as
!
!************************************************************
!
Integer Function lnblnk(c)
  Character(Len=*), Intent(In) ::  c
  Integer i
  lnblnk=Len(c)
  Do i=lnblnk,1,-1
     If(c(i:i).Ne.' ') Then 
        lnblnk=i
        Return
     Endif
  Enddo
  lnblnk=0
End Function lnblnk
!
!************************************************************
!
Integer Function fswlen(c)
  ! Return length of first word
  Character(Len=*), Intent(in) ::  c
  Integer lenc,i
  lenc=Len(c)
  Do i=1,lenc
     If(c(i:i)==' ')  Then
         fswlen=i-1
         Return 
     endif
  Enddo
  fswlen=lenc
  Return
End Function fswlen
!
!************************************************************
!
Subroutine maxloc2d(arr,stype,nx,ny,max,mi,mj)
  !
  ! Find the maximum of a 2D array arr
  ! return
  implicit none
  Integer, intent(in)                :: nx,ny
  real, dimension(nx,ny), intent(in) :: arr,stype
  real, intent(out)                  :: max
  Integer, intent(out)               :: mi,mj
  Integer                            :: i,j
  
  

  max=arr(1,1)
  mi=1
  mj=1
  do j=1,ny
     do i=1,nx
        if((arr(i,j)>max) .and. (stype(i,j)<0.01)) then 
           max=arr(i,j)
           mi=i
           mj=j
        endif
     enddo
  enddo
  return
End Subroutine maxloc2d
!
!************************************************************
!
Subroutine avg2d(arr,stype,nx,ny,avg)
  !
  ! compute average avg of array arr
  ! return
  implicit none
  Integer, intent(in)                :: nx,ny
  real, dimension(nx,ny), intent(in) :: arr,stype
  real, intent(out)                  :: avg
  Integer                            :: i,j,navg
  
  

  avg=0.0
  navg=0
  do j=1,ny
     do i=1,nx
        if(stype(i,j) <0.01) then
           avg=avg+arr(i,j)
           navg=navg+1
        endif
     enddo
  enddo
  avg=avg/navg
  return
End Subroutine avg2d
           
        
