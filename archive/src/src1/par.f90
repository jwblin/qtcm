! Generate parameters needed by QTCM model and produce vertical profiles
!c   for temperature, humidity, horizontal and vertical velocities,...
!
!  a1(p): vertical profile of perturbed temperature
!  a1p(p): a1^+
!  a1hat: a1 vertical average
!  V1(p): vertical profile of baroclinic velocity(v1); = a1^+ - <a1^+>
!          where <>denotes vertical average and ^+ is {dlnp
!  V12hat: vertical average of V1^2
!  V1s:   V1 at surface
!  Omega1(p): vertical profile of vertical velocity(omega); = {V1dp
!  b1(p): vertical humidity for preturbed humidity
!  b1hat: vertical average of b1
!  b1s:   b1 at surface, for evaporation
!  eps_1:  damping for v1
!  eps_0:  damping for v0(barotropic mode)
!  eps_01: cross damping due to v0 on v1
!  GMsr:  ref. state gross moist static energy vertically averaged;=<V1sref>
!  GMsp:  x T1 is the GMs perturbation due to T1
!  GMqr:  ref. state gross moist static energy vertically averaged;=<V1qref>
!  GMqp:  x T1 is the GMs perturbation due to q1
!
!  dp=1.e2      !data increment is 1mb=100pa
!  pt: cloud top
!  pb: cloud base
!  p0: surface pressure
!  pressure cooradinate:
!                       i=0: surface layer
!                       i=npb: cloud base
!                       i=npt: cloud top
!
Program par
  Parameter (p0=100000.,pt=15000.,pb=90000.,dp=100.,delp=p0-pt)
  Parameter (nz=p0/100.-1,npt=(p0-pt)/100.,npb=(p0-pb)/100.)
  Parameter (nvmod=2,nTmod=1)
  Real a1(0:npt),a1p(0:npt),V1(0:npt),Omega1(0:npt),b1(0:npt)
  Real V2(0:npt),Omega2(0:npt)
  Real zp(0:npt),zz(0:nz),V1z(0:nz)
  Real visc(0:npt),sref(0:npt),Tref(0:npt),qref(0:npt),qcref(0:npt)
  Real gamma(0:npt),Delta(0:npt),alphaf(0:npt)
  Real Tsat(0:npt),qsat(0:npt),ssat(0:npt)
  Real dpV1(0:npt),dpV2(0:npt)
  Real tmp(0:npt),tmp1(0:npt)
  Real V(0:nvmod-1,0:npt),Vsqhat(0:nvmod-1)
  Real Omega(0:nvmod-1,0:npt),dpV(0:nvmod-1,0:npt)
  Real Vijk(0:nvmod-1,0:nvmod-1,0:nvmod-1),tmp8(0:npt)
  Real Vwijk(0:nvmod-1,0:nvmod-1,0:nvmod-1),tmp9(0:npt)
  Real DTijk(0:nvmod-1,nTmod,nTmod)
  Real Dqijk(0:nvmod-1,nTmod,nTmod)
  Real qrefhat,bb1(0:npt),bb1hat
  Real V1zq(nz/100)
  Character(len=25) wformat ! formatting string
  !
  ! Physical contants, all in SI units 
  gravity=9.8
  Rair=287.04
  Cp=1004.
  rhoair=1.2
  Hlatent=2.43e6      !for T=30C, from Fleagle and Businger
  !
  ! model dependent parameters
  CD=1.4e-3
  Vs=10.          !m/s
  CV0=CD*Vs
  alpha=.834       !relative humidity for Betts-Miller
  !
  ! Delta: Hygrometric constant modified
  Do i=0,npt
     p=p0-i*dp
     Delta(i)=Hlatent/Cp*1.e5/p  !hsat,hsats is for p=1000mb
     sigma=(p-pt)/(p0-pt)
     alphaf(i)=0.6+0.25*sigma    !linear fit to observed RH
  !   alphaf(i)=0.95   !
  !   alphaf(i)=0.85  !
  End Do
  !
  !
  Open(12,file='Tqo3.ref',status='old')
!!$#if def QTCMPARIN
!!$  Open(21,file='qtcmpar.in',status='unknown')
!!$#endif
!!$#if def VTPARIN
!!$  Open(31,file='vtpar.in',status='unknown')
!!$#endif
#ifdef PROFOUT
  Open(23,file='prof.out',status='unknown')
#endif
#ifdef PAROUT
  Open(22,file='par.out',status='unknown')
#endif
!!$#if def V1ZOUT
!!$  Open(24,file='V1z',status='unknown')
!!$#endif
  !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  ! Tq.soundings: profiles of some warm place
  Read(12,*)Ts
  Do i=0,npt
     Read(12,*)Tref(i),qref(i)
     qref(i)=qref(i)*Hlatent/Cp    !convert humidity to Kelvin
  Enddo
  Trefhat=zbar(tref,npt)
  qrefhat=zbar(qref,npt)
  !
  !- moist adiabat, given temperature at p0
  !  Approximate pressure in Delta from current level,elsewhere middle value
  !    while temperature is from level below(small error since dp/dT small)
  Tsat(0)=Tref(0)        !Temperature at p0
  qsat(0)=Delta(0)*hsat(Tsat(0))
  qcref(0)=alphaf(0)*qsat(0)
  ! below cloud base (npb), Tsat is following dry adiabatic
  Do k=1,npb
     pp=p0-(k-.5)*dp
     Tsat(k)=Tsat(0)*(pp/p0)**(Rair/Cp)
     qsat(k)=Delta(k)*hsat(Tsat(k))
     qcref(k)=alphaf(k)*qsat(k)
  Enddo
  ! above cloud base (npb), Tsat is following moist adiabatic
  Do i=npb+1,npt
     p=p0-(i-.5)*dp
     tt=Tsat(i-1)
     Tsat(i)=tt+tt/(1.+Delta(i)*hsats(tt))                              &
          &            *(Rair/Cp+Delta(i)*hsat(tt)/tt)*(-dp)/p
     qsat(i)=Delta(i)*hsat(Tsat(i))
     qcref(i)=alphaf(i)*qsat(i)
  End Do

!  !c- adjust cloud temperature for entrainment effect -- drier upper levels
!  delT=-10.0
!  Do i=npb+1,npt
!     p=p0-(i-.5)*dp
!     Tsat(i)=Tsat(i)+delT*(pb-p)/(pb-pt)
!     qsat(i)=Delta(i)*hsat(Tsat(i))
!     qcref(i)=alphaf(i)*qsat(i)
!  End Do

  !- from the oberved reference profile (mean) in Tqo3.ref
  !       cloud top at npt from hb=ht
  Tcrefhat=zbar(Tsat,npt)
  qcrefhat=zbar(qcref,npt)
  Tcrefs=Tsat(0)
  qcrefs=qcref(0)
  !
  !- calculate a1 using Tsat
  Do i=0,npt
     p=p0-(i-1)*dp
     gamma(i)=Delta(i)*hsats(Tsat(i))
  End Do
  !
  tmp1(npb)=0.
  Do i=npb,npt-1
     p=p0-(i-.5)*dp
     tmp1(i+1)=tmp1(i)+1./(1.+.5*(gamma(i)+gamma(i+1)))*dp/p
  End Do
  Do i=npb,npt
     tmp1(i)=1./(1.+gamma(i))*Exp(-Rair/Cp*tmp1(i))
     a1(i)=tmp1(i)
  End Do
  Do k=0,npb
     p=p0-(k-.5)*dp
     a1(k)=(1./(1.+gamma(npb)))                                  &
          &         *(p/pb)**(Rair/Cp)
  Enddo
  !
  !
  !- a1hat,a1s
  a1hat=zbar(a1,npt)
  a1s=a1(0)
  !
  !- a1p(p)
  a1p(0)=0.
  Do i=0,npt-1
     pm=p0-(i+.5)*dp
     p=p0-(i+1)*dp
     a1p(i+1)=a1p(i)+.5*(a1(i)+a1(i+1))*dp/pm !medium value
!#ifdef CAP_V1
     ! mm cap the profile 
     if(p < 28000.0) then
        a1p(i+1)=a1p(i)
     endif
!#endif
  End Do
  a1phat=zbar(a1p,npt)
  !
  ! average over bndry layer for BL mean geopotential gradients:
  a1phatb=zbar(a1p,npb) 
  !
  !- V1(p),V12bar,V1s
  !        if(pt .lt. 20000.)then
  !          npt1=(p0-20000.)/100.
  !	  do i=0,npt1
  !	   V1v(i)=a1p(i)-a1phat
  !	  end do
  !          do i=npt1+1,npt
  !	   V1v(i)=V1v(npt1+npt1-i)
  !          enddo
  !        else
  !	  do i=0,npt
  !	   V1v(i)=a1p(i)-a1phat
  !	  end do
  !        endif
  Do i=0,npt
     V1(i)=a1p(i)-a1phat
     tmp(i)=V1(i)**2
  End Do
  !
  !
  V12hat=zbar(tmp,npt)
  V1s=V1(0)
  !
  !- Omega1(p)
  Omega1(0)=0.
  Do i=0,npt-1
     Omega1(i+1)=Omega1(i)+.5*(V1(i)+V1(i+1))*dp !approx. by medium value
  End Do
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc
  !- B1s=b1s
  Do i=0,npt
     bb1(i)=alphaf(i)*gamma(i)*a1(i)
  End Do
  bb1hat=zbar(bb1,npt)  ! bb1hat
  Do i=0,npt
     b1(i)=qref(i)/qref(0)
  Enddo
  b1s=b1(0)
  b1hat=zbar(b1,npt)
  !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
  !-- V2
  pvscal=.5e4           ! using the same profile for Vb as for viscosity
  pbl=8.e4-pvscal
  Do i=0,npt
     p=p0-(i-1)*dp
     tmp(i)=.6+.5*Tanh((p-pbl)/pvscal)
     tmp1(i)=tmp(i)*V1(i)
  End Do
  Vbhat=zbar(tmp,npt)
  tt1=zbar(tmp1,npt)
  Do i=0,npt
     V2(i)=tmp(i)-Vbhat-V1(i)*tt1/V12hat
  End Do
  V2s=V2(0)
  Do i=0,npt
     tmp(i)=V2(i)**2
  Enddo
  V22hat=zbar(tmp,npt)
  !
  !- Omega2(p)
  Omega2(0)=0.
  Do i=0,npt-1
     Omega2(i+1)=Omega2(i)+.5*(V2(i)+V2(i+1))*dp !approx. by medium value
  End Do
  !
  !- Interior viscosity; drops like tanh.
  ! ;assume isothermal so that rho(p)~p
  visc0=50.            !kinematic viscocity at surface; m^2/s
  Do i=0,npt
     p=p0-(i-1)*dp
     visc(i)=visc0*(.6+.5*Tanh((p-pbl)/pvscal))
  End Do
  !
  !
  !- eps_0,eps_01,eps_1,eps_10
  eps_0=gravity/delp*rhoair*CD*Vs
  eps_10=eps_0*V1s
  eps_20=eps_0*V2s
  !
  eps_01=eps_0*V1s/V12hat
  ! direct way of calculating eps_1; since dpV1=-a1/p, a simpler form is used
  Do i=0,npt-1
     p=p0-(i-1)*dp
     dpV2(i)=(V2(i)-V2(i+1))/dp
     dpV1(i)=(V1(i)-V1(i+1))/dp
     tmp1(i)=p
  End Do
  dpV1(npt)=dpV1(npt-1)       !approximate at the top
  dpV2(npt)=dpV2(npt-1)       !approximate at the top
  tmp1(npt)=p-dp
  Do i=0,npt
     tmp(i)=tmp1(i)**2*visc(i)*dpV1(i)*dpV1(i)
  Enddo
  temp=(gravity*rhoair/p0)**2*zbar(tmp,npt)
  !          write(*,*)'2nd term  new cal   =',temp/V12hat
  !
  ! ... since dpV1=-a1/p, this is used instead
  Do i=0,npt
     tmp(i)=visc(i)*a1(i)**2
  Enddo
  temp=(gravity*rhoair/p0)**2*zbar(tmp,npt)
  eps_1=eps_0*V1s**2/V12hat !  + temp/V12hat  !!changed mm
  eps_i1=temp/V12hat
  !	  write(*,*)'1st term in eps_1=',eps_0*V1s**2/V12hat
  !          write(*,*)'2nd term    =',temp/V12hat
  !
  Do i=0,npt
     tmp(i)=tmp1(i)**2*visc(i)*dpV2(i)*dpV1(i)
  Enddo
  temp=(gravity*rhoair/p0)**2*zbar(tmp,npt)
  eps_21=eps_0*V1s*V2s/V12hat+temp/V12hat
  !
  eps_02=eps_0*V2s/V22hat
  eps_12=eps_0*V1s*V2s/V22hat+temp/V22hat !temp as the same as for eps_21
  Do i=0,npt
     tmp(i)=tmp1(i)**2*visc(i)*dpV2(i)*dpV2(i)
  Enddo
  temp=(gravity*rhoair/p0)**2*zbar(tmp,npt)
  eps_2=eps_0*V2s**2/V22hat+temp/V22hat
  !
  !- ssat(p),sref(p),GMsr,GMsp,GMqr,GMqp
  ph=0.
  phref=0.
  Do i=0,npt-1
     ssat(i)=Tsat(i)+ph/Cp
     sref(i)=Tref(i)+phref/Cp
     p=p0-(i-1)*dp
     ph=ph+Rair*.5*(Tsat(i)+Tsat(i+1))*dp/p
     phref=ph+Rair*.5*(Tref(i)+Tref(i+1))*dp/p
  End Do
  ssat(npt)=Tsat(npt)+ph/Cp
  sref(npt)=Tref(npt)+phref/Cp
  GMsr=zbar2(V1,sref,npt)
  Do i=0,npt
     tmp(i)=a1(i)+a1p(i)*Rair/Cp
  Enddo
  GMsp=zbar2(V1,tmp,npt)
  GMqr=-zbar2(V1,qref,npt)
  GMqp=-zbar2(V1,b1,npt)
  !	GMqp=-zbar2(V1,b,npt)
  !- forcing GMsr and GMqr to be a fixed value because of numberical 
  !  instability
  GMsr=3.6
  GMqr=3.0
  !
  !- V0: divergence due to topographic lifting
  !       GMs0r=sref(0)+zbar(sref,npt)
  !       GMq0r=qref(0)+zbar(qref,npt)
  GMs0r=sref(npt)-sref(0)
  GMq0r=qref(npt)-qref(0)
  !
  !- PBL mode V2 (@not used)
  GMs2r=zbar2(V2,sref,npt)
  Do i=0,npt
     tmp(i)=a1(i)+a1p(i)*Rair/Cp
  Enddo
  GMs2p=zbar2(V2,tmp,npt)
  GMq2r=zbar2(V2,qref,npt)
  GMq2p=zbar2(V2,b1,npt)
  !
  Do i=0,npt
     V(0,i)=1.
  Enddo
  Vsqhat(0)=1.
  Do i=0,npt
     V(1,i)=V1(i)
  Enddo
  Vsqhat(1)=V12hat

  Do k=0,nvmod-1     !some redundance here(Vijk is somewhat symmetric)
     Do j=0,nvmod-1
        Do i=0,nvmod-1
           Do ii=0,npt
              tmp8(ii)=V(i,ii)*V(j,ii)*V(k,ii)
           Enddo
           Vijk(i,j,k)=zbar(tmp8,npt)/Vsqhat(k)
           If(Vijk(i,j,k).Lt.1.e-5) Vijk(i,j,k)=0.      !force small # zero
        End Do
     End Do
  End Do
  !
  Do i=0,npt
     p=p0-i*dp
     Omega(0,i)=p-pt
     Omega(1,i)=-Omega1(i)
     dpV(0,i)=0
     dpV(1,i)=dpV1(i)
  Enddo
  !
  Do k=0,nvmod-1     
     Do j=0,nvmod-1
        Do i=0,nvmod-1
           Do ii=0,npt
              tmp9(ii)=-V(k,ii)*Omega(i,ii)*dpV(j,ii)
           Enddo
           Vwijk(i,j,k)=zbar(tmp9,npt)/Vsqhat(k)
        End Do
     End Do
  End Do
  !
  !
  Do k=1,nTmod
     Do j=1,nTmod
        Do i=0,nvmod-1
           Do ii=0,npt
              tmp8(ii)=V(i,ii)*a1(ii)
           Enddo
           DTijk(i,j,k)=zbar(tmp8,npt)/a1hat
           Do ii=0,npt
              tmp8(ii)=V(i,ii)*b1(ii)
           Enddo
           Dqijk(i,j,k)=zbar(tmp8,npt)/b1hat
           abc=zbar(tmp8,npt)
        End Do
     End Do
  End Do
  !---
  !	Dqijk(1,1,1)=-6.5869723760904D-02
  !---
  Do kk=0,npt,50
     !           print *,kk,V(1,kk),b1(kk)
  Enddo
  !
  eps_c=1./7200.          !Betts-Miller adjustment time scale

  Cpg=Cp*delp/gravity     !conversion between column heating and flux
  !c- E, H ....
  eps_q=Cp*rhoair*CD*Vs   !evaporation        !W/m^2/K
  eps_s=Cp*rhoair*CD*Vs   !sensible heat      !W/m^2/K
  !
  !   (old)Jan(?;do annual later) Renolds SST average between 30N to 30S
  ! use a typical warm pool value
  Tsref=Ts
  !- Soil heat Capacity x Depth; assume .1m water-like soil
  slCD=4.18e3*1.e3*.1       !1cal/g/K * 1g/cm3 * .1m
  !
  qrefs=qref(0)
  Trefs=Tref(0)

  Do i=0,150,50
      p=p0-i*dp
#ifdef PROFOUT
     Write(23,*) p/100.,V1(i),a1(i),bb1(i),b1(i),-Omega1(i),     &
          &                  Tref(i),qref(i)
#endif
#ifdef VTPARIN
     Write(31,*) p,alphaf(i),Tsat(i),a1(i)
#endif
  End Do
  ! write out parameters at levels 800 - 200 mb
  Do i=200,npt,100
     p=p0-i*dp
#ifdef PROFOUT
     Write(23,*) p/100.,V1(i),a1(i),bb1(i),b1(i),-Omega1(i),     &
          &                  Tref(i),qref(i)
#endif
#ifdef VTPARIN
     Write(31,*) p,alphaf(i),Tsat(i),a1(i)
#endif
  End Do
!
#ifdef VTPARIN
  Write(31,*) a1hat
#endif
  !
  !- get V1 on z coordinate by linear interpolation
  ph=0.
  Do i=0,npt
     zp(i)=ph/gravity     !geopotential height
     p=p0-(i-1)*dp
     ph=ph+Rair*Tref(i)*dp/p
  End Do
  dz=10000/(nz+1)        !10km/nz
  Do i=0,nz
     zz(i)=(i-1)*dz
  End Do
  !
  Call lintp(npt,zp,V1,nz,zz,V1z)
!!$#if def V1ZOUT
!!$  Do i=0,nz,100
!!$     Write(24,*) V1z(i)
!!$  End Do
!!$#endif
  !ccccccccccccccccccc
  Do i=650,npt
     p=1000-i
     hhh=zp(i)*gravity+(qref(i)+Tref(i))*Cp
     hhh0=zp(50)*gravity+(Tref(50)+qref(50))*Cp
  Enddo
  !ccccccccccccccccccc
  !---------------------------
  !
  !
#ifdef PAROUT
  Do i=0,npt
     p=p0/100.-float(i)
     Write(22,'(17e13.5)') p,a1(i),V1(i),-Omega1(i)*1.e-3,b1(i)                     &
          &       ,Tsat(i),ssat(i),qsat(i),ssat(i)+qsat(i)                   &
          &       ,Tref(i),sref(i),qref(i),sref(i)+qref(i)                   &
          &       ,visc(i),visc(i)*a1(i)**2,V2(i),-Omega2(i)*1.e-3
  End Do
#endif
  !
!!$#if def QTCMPARIN
!!$  Write(21,*) a1hat,'        a1hat    dimensionless'
!!$  Write(21,*) a1s,  '        a1s'
!!$  Write(21,*) V1s,  '        V1s'
!!$  !       write(21,*) V2s,  '        V2s'
!!$  Write(21,*) b1hat,'        b1hat'
!!$  Write(21,*) b1s,  '        b1s'
!!$  Write(21,*) bb1hat,  '        bb1hat'
!!$  !       GMsr=0.9                 !using the estimation of Yu,Chou,Neelin (1997)
!!$  !       GMqr=-0.7
!!$  Write(21,*) GMsr, '        GMsr',Sqrt(GMsr*Rair/a1hat),'=c'
!!$  Write(21,*) GMsp, '        GMsp'
!!$  Write(21,*) GMqr, '        GMqr'
!!$  !       write(21,*) GMqr, '        GMqr',sqrt((GMsr+GMqr)*Rair/a1hat),'=cm'
!!$  Write(21,*) GMqp, '        GMqp'
!!$  Write(21,*) GMs0r, '        GMs0r'
!!$  Write(21,*) GMq0r, '        GMq0r'
!!$  !       write(21,*) GMs2r, '        GMs2r'
!!$  !       write(21,*) GMs2p, '        GMs2p'
!!$  !       write(21,*) GMq2r, '        GMq2r'
!!$  !       write(21,*) GMq2p, '        GMq2p'
!!$  Write(21,*) eps_c, '        eps_c   /s'
!!$  Write(21,*) CV0, '        CV0'
!!$
!!$
!!$  !
!!$  !c- allow some stress at TOA, so zonal mean of total zonal wind can be nonzero
!!$  ! HOwever, with advction this problem is gone. @Not used now
!!$  !     its contribution to V1,V2 mode is ignored for now
!!$  !ept=eps_0*.05        !arbitrarily set TOA stress 1/20 of that at surface
!!$  !eps_0=eps_0+ept
!!$  !eps_10=eps_10+ept*V1(np)
!!$  !
!!$  Write(21,*) eps_0, '        eps_0 /s ->',1./(86400*eps_0),'day'
!!$  Write(21,*) eps_10,'        eps_10 ',1./(86400*eps_10),'day'
!!$  Write(21,*) eps_01,'        eps_01',1./(86400*eps_01),'day'
!!$  Write(21,*) eps_1, '        eps_1',1./(86400*eps_1),'day'
!!$  Write(21,*) Tsref, '        Tsref     K'
!!$  Write(21,*) Cpg,   '        Cpg'
!!$  Write(21,*) Trefs, '        Trefs     K'
!!$  Write(21,*) qrefs, '        qrefs     K'
!!$  Write(21,*) Tcrefs, '       Tcrefs     K'
!!$  Write(21,*) qcrefs, '       qcrefs     K'
!!$  Write(21,*) Trefhat, '      Trefhat     K'
!!$  Write(21,*) qrefhat, '      qrefhat     K'
!!$  Write(21,*) Tcrefhat, '     Tcrefhat     K'
!!$  Write(21,*) qcrefhat, '     qcrefhat     K'
!!$  Write(21,*) Vijk,'         Vijk'
!!$  Write(21,*) Vwijk,'         Vwijk'
!!$  Write(21,*) DTijk,'        DTijk'
!!$  Write(21,*) Dqijk,'        Dqijk'
!!$  Write(21,*) (V1z(k),k=1,nz,100),'          V1(z)' 
!!$#endif

  !
  ! write fortran90 parameter module
  !
  Open(25,file='qtcmpar.f90')

  Write(25,'(a)')                                          &
       &  'Module Qtcmpar'                                 &
       & ,'  !'                                            &
       & ,'  ! QTCM mode parameter, generated by par.F90'  &
       & ,'  !'                                            &
       & ,'  Integer, Parameter  ::  &' 
  Write(25,'((a,i3,a))')                 &
       &  '       &       nvmod  = ',nvmod ,'  & !  no. v-modes' &
       & ,'       &     , nTmod  = ',nTmod ,'  & !  no. T-modes' &
       & ,'       &     , nz     = ',(nz+1)/100+4 ,'!  V1z table size'

  Write(25,'(a)')  '  REAL, Parameter     ::  &' 
  Write(25,'((a,g14.8,a))')                 &
       &  '       &       a1hat    = ',a1hat ,'  & ! hat{a_1}   [-] NZ (3.2)'  &
       & ,'       &     , a1phat   = ',a1phat,'  & ! hat{a_1^+} [-] NZ (2.16), (3.9)'  &
       & ,'       &     , a1s      = ',a1s   ,'  & ! a_{1s})    [-] NZ (5.16)' &
       & ,'       &     , V1s      = ',V1s   ,'  & ! V_{1s}     [-] NZ (4.13)' &
       & ,'       &     , V1sqhat  = ',V12hat,'  & ! hat{V_1^2} [-] NZ (4.13)' &
       & ,'       &     , b1hat    = ',b1hat ,'  & ! hat{b_1}   [-] NZ (5.17)' &
       & ,'       &     , b1s      = ',b1s   ,'  & ! b_{1s}     [-] NZ (4.32)' &
       & ,'       &     , bb1hat   = ',bb1hat,'  & ! hat{B_1}   [-] NZ (4.26)' &
       & ,'       &     , GMsr     = ',GMsr  ,'  & ! M_{sr}     [K] NZ (4.17)' &
       & ,'       &     , GMsp     = ',GMsp  ,'  & ! M_{sp}     [-] NZ (4.18)' & 
       & ,'       &     , GMqr     = ',GMqr  ,'  & ! M_{qr}     [K] NZ (5.10)' &
       & ,'       &     , GMqp     = ',GMqp  ,'  & ! M_{qp}     [-] NZ (5.10)' &
       & ,'       &     , eps_c    = ',eps_c ,'  & ! 1/tau_c [1/s] NZ (5.7)' 

  Write(25,'((a,g14.8,a,g9.2,a))')                            &  
       & '       &     , eps_i1   = ',eps_i1,'  & ! Mode 1 z-diffusion time [1/s]=>' &
       & ,1./(86400*eps_i1) ,'[day]'
  

  Write(25,'((a,g14.8,a))')                 &
       &  '       &     , Tsref    = ',Tsref   ,'  & ! ref. Ts for radiation [K]'&
       & ,'       &     , Cpg      = ',Cpg     ,'  & ! Cp*delp/gravity [J/K/m*2]'&
       & ,'       &     , Trefs    = ',Trefs   ,'  & ! surface Tref [K]' &
       & ,'       &     , qrefs    = ',qrefs   ,'  & ! surface qref [K]' &
       & ,'       &     , Tcrefs   = ',Tcrefs  ,'  & ! convect. surface Tref [K]'&
       & ,'       &     , qcrefs   = ',qcrefs  ,'  & ! convect. surface qref [K]'&
       & ,'       &     , Trefhat  = ',Trefhat ,'  & ! mean ref. temperature [K]'&
       & ,'       &     , qrefhat  = ',qrefhat ,'  & ! mean ref. humidity [K]' &
       & ,'       &     , Tcrefhat = ',Tcrefhat,'  & ! convect. mean ref. temp. [K]' &
       & ,'       &     , qcrefhat = ',qcrefhat,'    ! convect. mean ref. humidity [K]'
  Write(25,'(a)') '  !'
  Write(25,'(a)') '  ! advection coefficients'
  Write(25,'(a)') '  !'
  Write(25,'(a)') '  REAL, Parameter, Dimension(nvmod**3) ::  &' 
  Write(wformat,'(a,i1,a)') '(a,',nvmod,'(g14.8,a),i1,a,i1,a)'
  Write(wformat,'(a,i1,a)') '(a,',nvmod,'(g14.8,a),i1,a,i1,a)'
!    print*,'wformat = ',wformat
  j=0
  k=0
  Write(25,wformat) '       &    Vijkt = (/  '&
       &          ,(Vijk(i,j,k),',',i=0,nvmod-2),Vijk(nvmod-1,j,k)  &
       &                          ,'   & ! Vijk(:,',j,',',k,')'
  Do k=0,nvmod-1
     Do j=0,nvmod-1
        If ( (j==0 .And. k==0).Or.(j==nvmod-1.And.k==nvmod-1) ) Cycle ! skip first and last
        Write(25,wformat) '       &              , ',&
             &          (Vijk(i,j,k),',',i=0,nvmod-2),Vijk(nvmod-1,j,k)  &
             &                          ,'   & ! Vijk(:,',j,',',k,')'
     Enddo
  Enddo
  j=nvmod-1
  k=nvmod-1
  Write(25,wformat) '       &              , ',&
       &          (Vijk(i,j,k),',',i=0,nvmod-2),Vijk(nvmod-1,j,k)  &
       &                          ,'/) & ! Vijk(:,',j,',',k,')'

  ! Vwijk:
  j=0
  k=0
  Write(25,wformat) '       &  , Vwijkt = (/  ',&
       &         (Vwijk(i,j,k),',',i=0,nvmod-2),Vwijk(nvmod-1,j,k)  &
       &                          ,'  & ! Vwijk(:,',j,',',k,')'
  Do k=0,nvmod-1
     Do j=0,nvmod-1
        If ( (j==0 .And. k==0).Or.(j==nvmod-1.And.k==nvmod-1) ) Cycle ! skip first and last
!             &         ((Vwijk(i,j,k),','),i=0,nvmod-2),Vwijk(nvmod-1,j,k)  &
        Write(25,wformat) '       &             , ',&
             &         (Vwijk(i,j,k),',',i=0,nvmod-2),Vwijk(nvmod-1,j,k)  &
             &                          ,'    & ! Vwijk(:,',j,',',k,')'
     Enddo
  Enddo
  j=nvmod-1
  k=nvmod-1
!       &          ((Vwijk(i,j,k),','),i=0,nvmod-2),Vwijk(nvmod-1,j,k)  &
  Write(25,wformat) '       &             , ',&
       &          ((Vwijk(i,j,k)),i=0,nvmod-2),',',Vwijk(nvmod-1,j,k)  &
       &                          ,'/)    ! Vwijk(:,',j,',',k,')'

  !  write(25,'(/a)') '  REAL, parameter, dimension(0:nvmod-1,nTmod,nTmod) ::  &' 
  Write(25,'(/a)') '  REAL, parameter, dimension(nvmod*nTmod**2) ::  &' 
  ! DTijk:
  Write(25,wformat) '       &    DTijkt = (/  ',&
       &         (DTijk(i,j,k),',',i=0,nvmod-2),DTijk(nvmod-1,1,1)  &
       &                          ,'  /)  & ! DTijk(:,',1,',',1,')'
  ! Dqijk:
  Write(25,wformat) '       &  , Dqijkt = (/  ',&
       &         (Dqijk(i,j,k),',',i=0,nvmod-2),Dqijk(nvmod-1,1,1)  &
       &                          ,'  /)    ! Dqijk(:,',1,',',1,')'
  kk=0
  Do k=1,nz,100
     kk=kk+1
     V1zq(kk)=V1z(k)
  Enddo
  Write(25,'(a)') '  !'
  Write(25,'(a)') '  ! Table of V1(z) profile'
  Write(25,'(a)') '  !'
  Write(25,'(a)') '  REAL, Parameter, Dimension(2*nz) ::  V1z1d= &'
  Write(25,'(a)') '       &(/&!    z [m]     ,     V1 [-]'
  Do k=1,101,20
     Write(25,'(2(a,g14.8),a)') '       &    ',zz(k),', ',V1z(k),',   &'
  Enddo
  Do k=201,801,100
     Write(25,'(2(a,g14.8),a)') '       &    ',zz(k),', ',V1z(k),',   &'
  Enddo
  Write(25,'(4(a,g14.8),a)') '       &    ',zz(901),', ',V1z(901),'    /)'
  Write(25,'(a)') '#ifdef TOPO'
  Write(25,'(a)') '  ! Old V1(z) profile array; still used for topographic lifting'
  Write(25,'(a)') '  REAL, Parameter, Dimension(nz) ::  &'
  Write(25,'(4(a,g14.8),a)') '       & V1z = (/ ',V1z(1),(', ',V1z(k),k=21,61,20),' &'
  Write(25,'(4(a,g14.8),a)') '       & , ',V1z(81),(', ',V1z(k),k=101,301,100),'   &'
  Write(25,'(4(a,g14.8),a)') '       & , ',V1z(401),(', ',V1z(k),k=501,701,100),'   &'
  Write(25,'(2(a,g14.8),a)') '       & , ',V1z(801),', ',V1z(901),' /)'
  Write(25,'(a)') '#endif'
  Write(25,'(a)') 'End Module Qtcmpar'

  Write(25,'(/a)') 'Module T1cTableIn'
  Write(25,'(a)') '  ! Profile for nonlinear T1c lookup table'
  Write(25,'(a,i2)') '  Integer, Parameter :: np=',3+npt/100
  Write(25,'(a)') '  REAL, Parameter, Dimension(4*np) ::  table= &' 
  Write(25,'(a)') '       &(/&!  p [HPa]     ,  alpha [-]    ,  Tsat [K]    ,  a1 [-] '
  ! write out parameters at levels 1000 mb to 850 mb, at 50 mb interval
  Do i=0,150,50
     p=p0/100.-float(i)
     Write(25,'(4(a,g14.8),a)') '       &    ',p,', ',alphaf(i),',',Tsat(i),',',a1(i),',   &'
  End Do
  ! write out parameters at levels 800 - 200 mb
  Do i=200,npt-100,100
     p=p0/100.-float(i)
     Write(25,'(4(a,g14.8),a)') '       &    ',p,', ',alphaf(i),',',Tsat(i),',',a1(i),',   &'
  End Do
  p=p0/100.-float(i)
  Write(25,'(4(a,g14.8),a)') '       &    ',p,', ',alphaf(i),',',Tsat(i),',',a1(i),'    /)'
  Write(25,'(a)') 'End Module T1cTableIn'
!   Write things useful but not used in the model
     Write(25,'(a)') '                                                      '
     Write(25,'(a)') '            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     Write(25,'(a)') '            !! Things useful but not used in the model !!'
     Write(25,'(a)') '            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!'
     Write(25,'((a,g14.8,a))') &   
       & '!!  a1phatb  = ',a1phatb,' & ! ave. of {a_1^+} over subcloud layer' &
       &,'!!  GMs0r    = ',GMs0r,' & ! include divergence due to topo. lifting [K] '   &
       &,'!!  GMq0r    = ',GMq0r,' & ! include divergence due to topo. lifting [K] '   & 
       &,'!!  CV0      = ',CV0 ,'  & !  ref. drag coeff*velocity [m/s]'   
     Write(25,'((a,g14.8,a,g9.2,a))') & 
       & '!!  eps_0    = ',eps_0,' & ! [1/s]=>',1./(86400*eps_0), '[day]' &
       &,'!!  eps_10   = ',eps_10,' & ! [1/s]=>',1./(86400*eps_10), '[day]' &
       &,'!!  eps_01   = ',eps_01,' & ! [1/s]=>',1./(86400*eps_01), '[day]' &
       &,'!!  eps_1    = ',eps_1 ,' & ! [1/s]=>',1./(86400*eps_1), '[day]' 
!  Write(25,'(a)') '!! V1z For topographic lifting only' 
!  Write(25,'(a)') '!!  REAL, Parameter, Dimension(nz) ::  &'
!  Write(25,'(4(a,g14.8),a)') '!!  V1z = (/ ',V1z(1),(', ',V1z(k),k=101,301, &
!       &  100),'   &'
!  Write(25,'(4(a,g14.8),a)') '!!  & , ',V1z(401),(', ',V1z(k),k=501,701,100), ' &'
!  Write(25,'(2(a,g14.8),a)') '!!  & , ',V1z(801),', ',V1z(901), ' /)'
   Write(25,'(a)') '!!  a1p  at 1000, 950, 900, 850, 800 mb'
   Write(25,'(a)') '!!  REAL, Parameter, Dimension(5) ::  & '
   Write(25,'(3(a,g14.8),a)') '!!  & a1p = (/ ',a1p(0),(', ',a1p(k),k=49,99,50),'   &'
   Write(25,'(a,2(a,g14.8),a)') '!!  &', (', ',a1p(k),k=149,199,50),' /)'
!!   Write(25,'(5(a,g14.8),a)') '!!  & a1p = (/ ',a1p(0),(', ',a1p(k),k=49,199,50),' /)'

!
!  Write a qtcmpar.gs file
!
  open(26, file='qtcmpar.gs')
  write(26,'(a)') '* qtcmpar.gs: automatically produced by par '
  write(26,'(a)') '* '
  write(26,'(a)') '* Various QTCM parameters'
  write(26,'(a)') '* -----------------------'
  write(26,*)  
  write(26,'(a)') '* Mode 1 velocity vertical velocity structure V1(p):'
  Write(26,'(a,f10.6,a)')  '''define V1mb1000 = ',V1(0),''''
  Write(26,'(a,i3,a,f10.6,a)')  ('''define V1mb',1000-k,' = ',V1(k),'''',k=50,200,50)
  Write(26,'(a,i3,a,f10.6,a)')  ('''define V1mb',1000-k,' = ',V1(k),'''',k=300,800,100)
  write(26,'(a)') '* Surface and TOA'
  write(26,'(a)') '''define V1s=V1mb1000'''
  write(26,'(a)') '''define V1t=V1mb200'''
  write(26,*)
  write(26,'(a)') '* Vertical structure of T1, q1:'
  write(26,'(a,f10.6,a)') '''define a1hat = ',a1hat,' '''
  write(26,'(a,f10.6,a)') '''define a1s = ',a1s,' '''
  write(26,*)
  write(26,'(a)') '* Surface reference values:'
  write(26,'(a,f10.6,a)') '''define trefs = ',Trefs,' '''
  write(26,'(a,f10.6,a)') '''define qrefs = ',qrefs,' '''
  write(26,*)
  write(26,'(a)') '* Conversion factors'
  write(26,'(a)') '* '
  write(26,'(a)') '* g/kg=>K: L/C_p'
  write(26,'(a,f8.2,a)') '''define q2T = ',Hlatent/Cp*1e-3,' '''
  write(26,'(a)') '* mm/day=> W/m**2: L/(day [sec])'
  write(26,'(a,f10.6,a)') '''define rain2flux = ',Hlatent/(24.0*60.0*60.0),' '''
  write(26,'(a)') '* Column heating K/s => W/m^2: C_p*p_T/g'
  write(26,'(a,g13.7,a)') '''define Cpg = ',Cpg,' '''
  write(26,'(a)') '* ====End qtcmpar.gs==== '

End Program par

!
!
Function zbar2(a,b,n)
  !c averaging of a times b
  Real a(0:n),b(0:n)
  zbar2=0.
  Do i=0,n
     zbar2=zbar2+a(i)*b(i)
  End Do
  zbar2=zbar2/(n+1)
  !
  Return
End Function zbar2
Function zbar(a,n)
  !c averaging a(0:n)
  Real a(0:n)
  zbar=0.
  Do i=0,n
     zbar=zbar+a(i)
  End Do
  zbar=zbar/(n+1)
  !
  Return
End Function zbar
!
Real Function zbar8(a,n)
  !c averaging a(1:n)
  Real a(0:n)
  zbar8=0.d0
  Do i=0,n
     zbar8=zbar8+a(i)
  End Do
  zbar8=zbar8/(n+1)
  !
  Return
End Function zbar8
!
Function hsat(T1)
  T=T1-273.15
  !       T=T1
  esat=.6108d0*dexp(17.27d0*T/(237.3+T))  !from Shuttleworth,unit kpa (?)
  hsat=.622*esat/100.         ! unit= kg/kg, assuming Pressure=100kpa
  !
  Return
End Function hsat
!
!
Function hsats(T1)
  !-- the slope of hsat, by analytically differentiating hsat
  !       T=T1
  T=T1-273.15
  esats=2503.d0/(237.3+T)**2*dexp(17.27d0*T/(237.3+T))
  hsats=.622*esats/100.         ! unit= kg/kg, assuming Pressure=100kpa
  !
  Return
End Function hsats
!
!
Subroutine lintp(n1,x1,y1,n2,x2,y2)
  !****     linear interpolation. Interpolate y1(x1(n1)) on to y2(x2(n2)) 
  !****   if |x2|>|x1|, the result is not clear (no extrapolation at the end?)
  !   n1  # of points in x1
  !   x1  variable of input data array
  !   y1  function of input data array, x1(n1) -> y1=y1( x1(n1) )
  !   n2,x2,y2 for output
  !
  Integer n1,n2
  Real x1(0:n1),y1(0:n1),x2(0:n2),y2(0:n2)
  !
  !   write(*,*)n1,n2
  If(x2(0).Lt.x1(0).Or.x2(n2).Gt.x1(n1)) Then
     !	    write(*,*) 'terror,  |x2|>|x1|, need extropolation'
     !    stop
  End If
  !
  j=0
  Do i2=0,n2
100  j=j+1
     !   write(*,*)i2,j
     If(x2(i2).Lt.x1(j)) Then
        j=j-1
        p=( x2(i2)-x1(j) ) / ( x1(j+1)-x1(j) ) 
        y2(i2)=p*y1(j+1)+(1.-p)*y1(j)
     Else If(j.Gt.n1) Then
        y2(i2)=y1(n1)
     Else
        go to 100
     End If
  End Do

  Return
End Subroutine lintp
