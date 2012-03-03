c This program generates various SST masks. Selected region mask=1. Otherwise,
c mask=0. 
 
      parameter(NX=64,NY=42)
      parameter(xlon_beg=0.,ylat_beg=-78.75,dlon=5.625,dlat=3.75)
      real mask_nino3(NX,NY),mask_ensopac(NX,NY)
      real mask_allobs(NX,NY),mask_allmxl(NX,NY)
      real xlon(NX),ylat(NY),stp(NX,NY)
!      open(25,file='nino3.mask',form='formatted')
      open(25,file='ensopac.mask',form='formatted')
      open(26,file='allobs.mask',form='formatted')
      open(27,file='allmxl.mask',form='formatted')
!      open(30,file='nino3mask.out',form='unformatted')
      open(30,file='ensopacmask.out',form='unformatted')
      open(41,file='../STYPE',form='formatted',status='old')

c- initializing 
c
      do j=1,NY
      do i=1,NX
         mask_nino3(i,j)=0.
         mask_ensopac(i,j)=0.
         mask_allobs(i,j)=1.
         mask_allmxl(i,j)=0.
      enddo
      enddo

c- create model grid information
c
      do i=1,NX
         xlon(i)=xlon_beg + dlon*(i-1)
      enddo 
      do j=1,NY
         ylat(j)=ylat_beg + dlat*(j-1) 
      enddo

c- set masks
c
      do j=1,NY
      do i=1,NX
         if((ylat(j).ge.(-30.)).and.(ylat(j).le.30.).and.
     &   (xlon(i).ge.180.).and.(xlon(i).le.280.)) mask_ensopac(i,j)=1.  !NINO34 SST
c         if((ylat(j).ge.(-5.)).and.(ylat(j).le.5.0).and.
c     &       (xlon(i).ge.190.).and.(xlon(i).le.240.)) mask_nino34(i,j)=1.  !NINO34 SST
c     &       (xlon(i).ge.210.).and.(xlon(i).le.270.)) mask_nino3(i,j)=1.  !NINO3 SST
      enddo
      enddo
c- apply land-sea masks
c
      do j=1,NY
      do i=1,NX
         read(41,*) stp(i,j) 
         if(stp(i,j).ne.0) then
            mask_ensopac(i,j)=0.
         endif 
        write(25,*) mask_ensopac(i,j)
        write(26,*) mask_allobs(i,j)
        write(27,*) mask_allmxl(i,j)
      enddo
      enddo

c- GrADs format
c
        write(30) mask_ensopac
      stop
      end
