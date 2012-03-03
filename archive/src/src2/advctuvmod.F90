Module advctuvmod
  Use Dimensions
  !@@@Use Prognostic, Only : u0,v0,u1,v1
  Use Prognostic
  Use Baroclinic
  Use Barotropic
  !@@@Use Baroclinic, Only : div1,advu1,advv1,advwu1 !,advwv1 
  !@@@Use Barotropic, Only : div0,advu0,advv0,advwu0 !,advwv0 
  Use Grid
  Use Qtcmpar
  Use AuxVars

Contains
!------------------------------------------------------------------
!
Subroutine advctuv
  !
  ! Advection for momentum u,v,u0,v0, ...
  ! Centered second order differencing except duy at N/S boundaries
  !
  ! Array adv* is used to pass variables.  This is useful for
  !   the barotropic mode where the curl of advection is needed later
  !   (because it is solved in streamfunction-vorticity formulation).
  ! Otherwise the procedure can be applied pointwise without saving
  !   in adv*.  Since memory is not a major concern here we do
  !   this for diffusion as well as advection/diffusion for other
  !   modes.  This always looks neater.
  ! Note the three-level looping due to the nonlinear advection and
  !   projection.
  ! Moderately optimized for speed
  !


  Implicit None

  !
  !
  ! --Altered Module Variables:
  ! advu1,advv1,advu0,advv0
  !
  ! --Local Variables & Functions:
!f2py call advctuv()

  Real u1atu,v1atu,u1atv,v1atv
  Real u0atu,v0atu,u0atv,v0atv
  Real du1xdxi,du1ydyi,dv1xdxi,dv1ydyi
  Real du0xdxi,du0ydyi,dv0xdxi,dv0ydyi                   
  Real ajp1(ny),ajm1(ny),aj(ny)
  Real vdv1,vdv2,vdv3,vdv4!,Vijkt(8)
  Real divv1,divv2,divv3,divv4!,Vwijkt(8)
  Real dxih,dyih       ! mm speed
  Integer i,j

  dyih=.5*dyi 
  !
  !c- handle 1st-order differencing for duy at y boundaries
  Do j=2,ny-1            !centered differencing in interior
     ajp1(j)=.5       
     aj(j)=0.
     ajm1(j)=-.5
  End Do
  ajp1(1)=1.             !1st order differencing at bnd.
  aj(1)=-1.
  ajm1(1)=0.
  ajp1(ny)=0.
  aj(ny)=1.
  ajm1(ny)=-1.
  !
  !=  x momentum equation; advection for u1, u0
  !
  Do j=1,ny
     dxih=0.5*dxui(j)
     Do i=1,nx
        ! mode V1
        du1xdxi=(u1(ip1(i),j)-u1(im1(i),j))*dxih       !centered differencing
        du1ydyi=(ajp1(j)*u1(i,j+1)+aj(j)*u1(i,j)      &
             &         +ajm1(j)*u1(i,j-1))*dyi
        u1atu=u1(i,j)
        v1atu=.250*(v1(i,j)+v1(ip1(i),j)+v1(i,j-1)    &
             &                     +v1(ip1(i),j-1))
        ! mode V0
        du0xdxi=(u0(ip1(i),j)-u0(im1(i),j))*dxih       !centered differencing
        du0ydyi=(ajp1(j)*u0(i,j+1)+aj(j)*u0(i,j)      &
             &         +ajm1(j)*u0(i,j-1))*dyi
        u0atu=u0(i,j)
        v0atu=.250*(v0(i,j)+v0(ip1(i),j)+v0(i,j-1)    &
             &                     +v0(ip1(i),j-1))
        !
        ! Loop for different modes; 
        ! explicit summation makes the routine [advctuv]
        !   more than twice as fast. However, looping is easier to implement
        !   if you need to add more modes 
        ! (simply change the parameter value NVMOD).
        !       do k=0,NVMOD-1
        !         advu(k)=0.
        !         do n=0,NVMOD-1
        !         do m=0,NVMOD-1
        !           advu(k)=advu(k)-Vijk(m,n,k)
        !    1       *(uatu(m)*dux(n)*dxi+vatu(m)*duy(n)*dyi)
        !    2       -Vwijk(m,n,k)*(div(m)*uatu(n)) !symbolic form w-advection 
        !         end do
        !         end do
        !       end do
        !
        vdv1=u0atu*du0xdxi+v0atu*du0ydyi
        vdv2=u1atu*du0xdxi+v1atu*du0ydyi
        vdv3=u0atu*du1xdxi+v0atu*du1ydyi
        vdv4=u1atu*du1xdxi+v1atu*du1ydyi
        !
        != divergence; curvature term neglected; using u,v at time t
        div1(i,j)=(u1(i,j)-u1(im1(i),j))*dxui(j)     &
             &   +(v1(i,j)*cosv(j)-v1(i,j-1)*cosv(j-1))*dyui(j)
        !not considering spherical coordinates
        !     &              +(v1(i,j)-v1(i,j-1))*dyi
        divv1=0.
        divv3=0.
#ifdef TOPO
        divv1=div0(i,j)*u0atu
        divv3=div0(i,j)*u1atu
#endif
        divv2=div1(i,j)*u0atu      
        divv4=div1(i,j)*u1atu 

!!$        advu0(i,j)=                        &
!!$             &      -Vijkt(1)*vdv1         &
!!$             &      -Vijkt(2)*vdv2         &
!!$             &      -Vijkt(3)*vdv3         &
!!$             &      -Vijkt(4)*vdv4         &
        !             &      -Vijkt(2)*vdv2         &  ! Vijkt(2)=0
        !             &      -Vijkt(3)*vdv3         &  ! Vijkr(3)=0
        !             &      -Vwijkt(1)*divv1       &
        !             &      -Vwijkt(2)*divv2       &
        advu0(i,j)=                        &
             &      -vdv1                   &  ! Vijkt(1)=Vijk(0,0,0)=1.0
             &      -Vijkt(4)*vdv4         &
#ifdef V1V0ADVH
        &       *0.5                  &
#endif
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(3)*divv3       &
             &      -Vwijkt(4)*divv4       &
#ifdef V1V0ADVH
        &       *0.5                  &
#endif
        &;
#endif

        advwu0(i,j)=                       &
             &      -Vwijkt(3)*divv3       &
             &      -Vwijkt(4)*divv4
        !             &      -Vwijkt(1)*divv1       & ! Vwijkt(1)=0
        !             &      -Vwijkt(2)*divv2       & ! Vwijkt(2)=0

        !    1      -Vijk(0,0,0)*vdv1    !this is about 30% slower than above
        !    2      -Vijk(1,0,0)*vdv2
        !    3      -Vijk(0,1,0)*vdv3
        !    4      -Vijk(1,1,0)*vdv4
        !    1      -Vwijk(0,0,0)*divv1
        !    2      -Vwijk(1,0,0)*divv2
        !    3      -Vwijk(0,1,0)*divv3
        !    4      -Vwijk(1,1,0)*divv4
        !
        !             &      -Vijkt(5)*vdv1         & !Vijkt(5)=0
        !             &      -Vijkt(6)*vdv2         & !Vijkt(6)=1
        !             &      -Vijkt(7)*vdv3         &  !Vijkt(7)=1
        !             &      -Vwijkt(5)*divv1       & !Vijkt(5)=0
        !             &      -Vwijkt(6)*divv2       & !Vijkt(6)=0

        advu1(i,j)=                        &
             &      -( vdv2                 &  !Vijkt(6)=1
             &         +vdv3                &  !Vijkt(7)=1
             &         +Vijkt(8)*vdv4 )     &
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4   
#endif

        advwu1(i,j)=                       &
             &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4   
     End Do
  End Do

  !
  !=  y momentum equation; advection for v0, v1
  Do j=1,ny-1
     dxih=0.5*dxvi(j)
     Do i=1,nx
        ! mode V1
        dv1xdxi=(v1(ip1(i),j)-v1(im1(i),j))*dxih
        dv1ydyi=(v1(i,j+1)-v1(i,j-1))*dyih
        u1atv=.250*(u1(i,j)+u1(im1(i),j)+u1(i,j+1)    &
             &                     +u1(im1(i),j+1))
        v1atv=v1(i,j)
        ! mode V0
        dv0xdxi=(v0(ip1(i),j)-v0(im1(i),j))*dxih
        dv0ydyi=(v0(i,j+1)-v0(i,j-1))*dyih
        u0atv=.250*(u0(i,j)+u0(im1(i),j)+u0(i,j+1)    &
             &                     +u0(im1(i),j+1))
        v0atv=v0(i,j)
        !
        vdv1=u0atv*dv0xdxi+v0atv*dv0ydyi
        vdv2=u1atv*dv0xdxi+v1atv*dv0ydyi
        vdv3=u0atv*dv1xdxi+v0atv*dv1ydyi
        vdv4=u1atv*dv1xdxi+v1atv*dv1ydyi
        !
        divv1=0.
        divv3=0.
#ifdef TOPO
        divv1=div0(i,j)*v0atv               !div0 is 0
        divv3=div0(i,j)*v1atv               !div0 is 0
#endif
        divv2=div1(i,j)*v0atv
        divv4=div1(i,j)*v1atv

        !             &      -Vijkt(1)*vdv1          & !   =-vdv1
        !             &      -Vijkt(2)*vdv2          & !   =0
        !             &      -Vijkt(3)*vdv3          & !   =0
        !             &      -Vijkt(3)*vdv3          & !   =0
        !             &      -Vijkt(5)*vdv1         &  ! = 0
        !             &      -Vijkt(6)*vdv2         &  ! = vdv2
        !             &      -Vijkt(7)*vdv3         &  ! = vdv3
        !             &      -Vwijkt(5)*divv1       &  ! = 0
        !             &      -Vwijkt(6)*divv2       &  ! = 0
        advv0(i,j)=                         &
             &      -vdv1          &
             &      -Vijkt(4)*vdv4          &
#ifdef V1V0ADVVH
        &       *0.5                  &
#endif
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(3)*divv3        &
             &      -Vwijkt(4)*divv4        &
#endif
#ifdef V1V0ADVVH
        &       *0.5                  &
#endif
        &;


        advv1(i,j)=                        &
             &      - (  vdv2              &
             &         +vdv3               &
             &         +Vijkt(8)*vdv4 )    & 
#ifdef NO_WADV
        &;
#else
        &      -Vwijkt(7)*divv3       &
             &      -Vwijkt(8)*divv4
#endif

     End Do
  End Do
  Return
End Subroutine advctuv
End Module advctuvmod

!======================================================================
!            End of qtcmutil.F90 file
!======================================================================
