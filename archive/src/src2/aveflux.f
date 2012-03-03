c Program aveflux   by H. Su  
c read QTCM output and calculate annual and monthly mean Qflux 
c assume monthly data inputs with 360 days in a year 
c aveflux.f --routine to compute monthly and annually averaged net surface flux
c             and the rate of change of surface temperature. These two terms 
c             give an approximation of ocean transport used for runs coupled 
c             with slab mixed-layer ocean model.
c 
c Note: similar to ave31.f, aveflux.f takes in monthly QTCM output and produces
c monthly and annually averaged surface net flux (fsn) and rate of change of 
c surface temperature (dts). Due to the interpolation scheme of real-time SST
c reader in QTCM non-coupled runs, monthly data for 'dts' is centered on the 
c first of each month, denoted as 0000mm01.dts. 'fsn' is still centered on the 
c 15th. 
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c Note: Please make sure the variables read by aveflux are consistent
c with the ones output from the model. 
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


c   ntmax: max # of time steps (months)
      Use Dimensions
      implicit none
      
      integer n,i,j,ntmax     
      parameter (ntmax=1200)
c baroclinic mode
      real u1(NX,NY),v1(NX,NY),T1(NX,NY),q1(NX,NY)
     2    ,advT1(NX,NY),advq1(NX,NY) 
c model physics variables associated with temperature mode 1, including:
c   convection, radiation, cloud, heat and moisture fluxes
      real    Prec(NX,NY)                    !precipitation
      real    FLWds(NX,NY),FLWus(NX,NY)      !down/up LW fluxes at surface
     2       ,OLR(NX,NY)                   !LW at top; LW column total
     3       ,S0(NX,NY)                      !incoming solar at top
     4       ,FSWds(NX,NY),FSWus(NX,NY)      !down/up SW fluxes at surface
     5       ,FSWut(NX,NY)                   !SW at top; SW absorbed by atmo.
     6       ,cl1(NX,NY)                     !cloud cover:deep cloud with anvil
     7       ,Evap(NX,NY),FTs(NX,NY)         !evap and sensible heat
     8       ,taux(NX,NY),tauy(NX,NY)        !surface stress
     9       ,GMq1(NX,NY),GMs1(NX,NY)        !gross static stability
     &       ,div0(NX,NY),div1(NX,NY)        !divergence 
     &       ,dfsT1(NX,NY),dfsq1(NX,NY)      !diffusion
c      real FLW(NX,NY),FSW(NX,NY)
c barotropic mode V0; note only vort0 and u0bar are prognostic
      real u0(NX,NY),v0(NX,NY),vort0(NX,NY),psi0(NX,NY)
c land variables
      real WD(NX,NY),Runf(NX,NY),Runs(NX,NY),Evapi(NX,NY),wet(NX,NY)
c temporary arrays
      real arr1(NX,NY),arr2(NX,NY),arr3(NX,NY),arr4(NX,NY)
      real Ts(NX,NY),STYPE(NX,NY),us(NX,NY),vs(NX,NY)
      real fsnet(NX,NY,ntmax),sst(NX,NY,ntmax),sstdt(NX,NY,ntmax)
      real Dmx,Cmx,dataint
c******************************************************************
c- Mixed-layer heat capacity
c  water heat capacity * density * Depth
      Dmx = 50.     !Depth of the mix-layer ocean;  meters
      Cmx=4.18e6*Dmx      !1cal/g/K * 1g/cm3 * Depth;  in J/K/M2
c******************************************************************
      open(10,file='qm.out',
     1      form='unformatted',status='old')
c- input data interval - 30 days
      dataint=30.*86400.
c- read model outputs
      n=0
100   read(10,end=999) u1
      read(10) v1
      read(10) T1
      read(10) q1
      read(10) u0
      read(10) v0
      read(10) vort0
      read(10) psi0
      read(10) Ts
      read(10) Prec
      read(10) cl1
      read(10) S0
      read(10) FSWds
      read(10) FSWus
      read(10) FSWut
      read(10) FLWds
      read(10) FLWus
      read(10) OLR 
      read(10) FTs
      read(10) Evap
      read(10) Runf
      read(10) WD
      read(10) STYPE
      read(10) taux
      read(10) tauy
      read(10) advT1
      read(10) advq1
      read(10) Evapi 
      read(10) GMq1 
      read(10) wet 
      read(10) Runs 
      read(10) dfsT1 
      read(10) dfsq1 
      read(10) GMs1 
      read(10) div0 
      read(10) div1 
      read(10) us 
      read(10) vs 
      n=n+1
      print *, 'read model output at month ',n
      do j=1,NY
      do i=1,NX
         fsnet(i,j,n)=FSWds(i,j)-FSWus(i,j)+FLWds(i,j)-
     &              FLWus(i,j)-Evap(i,j)-FTs(i,j)
         sst(i,j,n)=Ts(i,j)
         if(n.eq.1) then
           sstdt(i,j,1)=0.
         else
           sstdt(i,j,n)=(sst(i,j,n)-sst(i,j,n-1))/dataint
         endif
         sstdt(i,j,n)=sstdt(i,j,n)*Cmx
      enddo
      enddo
      goto 100
999   continue 
      call average(sstdt,'dts',n)
      call average(fsnet,'fsn',n)
      stop
      end

      subroutine average(f,var,ntime)
c-----calculate annual mean and monthly mean
      Use Dimensions
c      include '../src/hgrid.h'
      dimension f(nx,ny,ntime),fam(nx,ny),fmm(nx,ny,12)
      character var*3,fsmonth*20
      integer lnblnk
      open(21,file=var(1:lnblnk(var))//'am.dat',
     &     form='unformatted',status='unknown')
      open(22,file=var(1:lnblnk(var))//'mm.dat',
     &     form='unformatted',status='unknown')
      open(31,file='00000000.'//var(1:lnblnk(var)),
     &     form='formatted',status='unknown')
c initialize fmm and fam
      do j=1,NY
      do i=1,NX
         fam(i,j)=0.
         do im=1,12
            fmm(i,j,im)=0.
         enddo
      enddo
      enddo
c-----
      iyear=int(ntime/12) 
      print *, 'climatology for ', var
      print *, 'averaged over ', iyear, ' years'

      month_start=1            !assume model output starts from January
      do 200 j=1,NY
      do 200 i=1,NX
      do im=1,12
         do iy=1,iyear          
             fmm(i,j,im)=fmm(i,j,im)+f(i,j,12*(iy-1)+im)
             fam(i,j)=fam(i,j)+f(i,j,12*(iy-1)+im)
         enddo 
      enddo
200   continue

c-----output for grads display

      write(21) ((fam(i,j)/float((iyear)*12),i=1,NX),j=1,NY)
      do im=1,12
      write(22) ((fmm(i,j,im)/float(iyear),i=1,NX),j=1,NY)
      enddo

c-----output for model inputs

      do j=1,NY
      do i=1,NX
         write(31,*) fam(i,j)/float((iyear)*12)
      enddo
      enddo
      do im=1,12
         month=month_start+im-1
         if(month.gt.12) month=month-12
         if(var.eq.'dts') then
            write(fsmonth,"('0000',i2.2,'01.',a3)") month,
     &      var(1:lnblnk(var))
         else
            write(fsmonth,"('0000',i2.2,'15.',a3)") month,
     &      var(1:lnblnk(var))
         endif

         open(32,file=fsmonth,form='formatted',status='unknown')

         do j=1,NY
         do i=1,NX
            write(32,*) fmm(i,j,im)/float(iyear)
         enddo
         enddo
         close(32)
      enddo
      return
      end
