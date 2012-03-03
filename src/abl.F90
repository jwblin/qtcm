!
! SfcWind_ABL: Solves for steady states of the mixed layer momentum equations
!              Using 2-D Newton-Ralphson

Subroutine SfcWindABL(u,v,VVs)

  Use Prognostic, Only : u0,v0,u1,v1,nx,ny
  Use Grid,       Only : fu     ! Coriolis parameter at T points
  Use Surface,    Only : cdn  & ! drag coeff
       &  , ub,vb             & ! top of ABL wind
       &  , dphisdx,dphisdy   & ! gradient surface geopotential
       &  , STYPE               ! surface type 
  ! Use Input,      Only : V1b  & ! Mode 1 projection coeff. at top of ABL
  !     &  , weml              & ! entrainment velocity
  Use Input,      Only : weml & ! entrainment velocit
       &  , ziml              & ! mixed layer height, minimum wind speed
       &  , VVsmin  

  Implicit None

  Real, Dimension(nx,ny), Intent(inout) :: u    & ! zonal surface wind
       & , v    & ! meridional surface wind
       & , VVs   ! wind speed for surface fluxes

  Real   :: f          & ! remainder of zonal wind equation 
       &  , dfdu       & ! df/du
       &  , dfdv       & ! df/dv
       &  , g          & ! remainder of meridional wind equation 
       &  , dgdu       & ! dg/du
       &  , dgdv       & ! dg/dv
       &  , determ     & ! deterimant of jacobian
       &  , dphisdxT,dphisdyT & ! geopot gradient at T-point
       &  , cdzi       & !cd/zi
       &  , cdziwind   & ! cd/zi*VVs
       &  , cdziwindi  & ! cd/zi/VVs
       &  , usav, vsav & ! save old us, vs, VVs for cases without convergence
       &  , VVssav
  Real, Save :: wezi     & ! we/zi
       &      , VVsminsq & ! VVsmin**2
       &      , zii      & ! 1/zi
       &      , V1b        ! Mode 1 projection coeff. at top of ABL
  Real :: uij,vij,svij,fuj, V1zi, V1interpol

  Integer       :: i,j, iterate       ! Indices
  Logical, Save :: firstcall = .True. ! Initialization flag

  If(firstcall) Then
     firstcall=.False.
     ! compute mode 1 projection at top of ABL
     V1b=V1interpol(ziml)
     wezi=weml/ziml
     zii=1./ziml
     VVsminsq=VVsmin**2
     Print*,' Abl: weml = ',weml,',  ziml=',ziml, ',  VVsmin=',VVsmin &
          & ,',  V1b = ',V1b
     Do j=1,ny
        Do i=1,nx
           !
           ! Initially use top of ABL as first guess
           ! Note: This make QTCM not strictly restartable
           u(i,j)=ub(i,j)
           v(i,j)=vb(i,j)
           VVs(i,j)=Sqrt(VVsminsq+u(i,j)**2+v(i,j)**2)
           ! Double land cdn to compensate mountain drag
           If(STYPE(i,j)>0.5) Then
              cdn(i,j)=cdn(i,j)*2.0
           Else
              cdn(i,j)= 0.0011
           Endif
        End Do
     End Do
     Print*,'Land cdn doubled.'
     Print*,'Ocean cdn= ',cdn(nx/2,ny/2) ! A water point
  End If

  !
  ! Cloud base winds:
  Do j=1,ny
     ! i=1
     ub(1,j)=0.5*( (u1(1,j)+u1(nx,j))*V1b  & ! at T grid
          &        +u0(1,j)+u0(nx,j)   )        
     vb(1,j)=0.5*( (v1(1,j)+v1(1,j-1))*V1b &
          &       + v0(1,j)+v0(1,j-1)     )
     Do i=2,nx
        ub(i,j)=0.5*( (u1(i,j)+u1(i-1,j))*V1b & ! at T grid
             &        +u0(i,j)+u0(i-1,j)   )        
        vb(i,j)=0.5*( (v1(i,j)+v1(i,j-1))*V1b &
             &       + v0(i,j)+v0(i,j-1)     )
     End Do
  End Do

  Do j=1,ny
     fuj=fu(j)
     Do i=1,nx
        uij=u(i,j)
        vij=v(i,j)
        svij=VVs(i,j)
        usav=uij
        vsav=vij
        VVssav=svij
        !
        ! Interpolate geopot. gradients to T-grid
        If(i==1) Then
           dphisdxT=0.5*(dphisdx(1,j)+dphisdx(nx,j))
        Else
           dphisdxT=0.5*(dphisdx(i-1,j)+dphisdx(i,j))
        Endif
        dphisdyT=0.5*(dphisdy(i,j-1)+dphisdy(i,j))

        cdzi=cdn(i,j)*zii
        Do iterate = 1,10
           cdziwind=cdzi*svij
           cdziwindi=cdzi/svij
           !
           ! Remainder of steady state equation for zonal wind
           f=ub(i,j)*wezi+fuj*vij-dphisdxT-uij*(cdziwind+wezi)
           !
           ! Remainder of state equation for meridional wind
           g=vb(i,j)*wezi-fuj*uij-dphisdyT-vij*(cdziwind+wezi)
           If(Abs(f)+Abs(g)<1e-9) Exit ! Convergence reached
           !
           ! Derivatives
           dfdu=-wezi - cdziwind - uij**2*cdziwindi 
           dfdv=fuj - uij*vij*cdziwindi 
           dgdv=-wezi - cdziwind - vij**2*cdziwindi
           dgdu=-fuj - uij*vij*cdziwindi 

           determ = dfdu*dgdv - dfdv*dgdu
           If (determ /= 0.) Then
              !
              ! New Newton-Ralphson guess of u,v
              uij=uij-(dgdv*f-dfdv*g)/determ
              vij=vij-(dfdu*g-dgdu*f)/determ
              svij=Sqrt(VVsminsq+uij**2+vij**2)
           Else
              Print*,'Abl: WARNING: determinant zero at (i,j)=(',i,',',j,')'
           End If
        End Do
        If(iterate >9 ) Then
           Print*,'No convergence for i,j,ub,vb,f,g = ',i,j,ub(i,j),vb(i,j),f,g
           Print*,'Using prior values for us,vs = ',usav,vsav
           u(i,j)=usav
           v(i,j)=vsav
           VVs(i,j)=VVssav
        Else
           u(i,j)=uij
           v(i,j)=vij
           VVs(i,j)=svij
        Endif
     End Do
  End Do

  Return
End Subroutine SfcWindABL
!
! -----------------------------------------------------------------------------
!
Function V1interpol(zi)
  !
  ! Interpolate linearly V1z to top of ABL level
  Use Qtcmpar, Only : V1z1d,nz
  Implicit None
  Real, Dimension(nz) :: zV1,V1z

  Integer :: i
  Real, Intent(In) :: zi
  Real :: V1interpol, we, zp

  ! V1z1d is formally a 1D array 
  ! extract height zV1, and V1z from it
  Do i=1,nz
     zV1(i)=V1z1d(2*i-1)
     V1z(i)=V1z1d(2*i)
  Enddo

  ! Find first tabled value above zi
  i=1
  Do i=1,nz
     If (zV1(i)>=zi) Exit
  Enddo
  if(i==1 .or. i==nz+1) then
     Print*,' ABL height too big or small: ziml = ',zi,' m'
     Stop
  Endif

  ! Interpolate
  we=(zV1(i)-zi)/(zV1(i)-zV1(i-1))
  V1interpol=V1z(i)*(1.0-we)+V1z(i-1)*we
  Return

End Function V1interpol
!                        End of abl.f90
! =============================================================================
