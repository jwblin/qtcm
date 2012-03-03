c-- convert free-formatted data into unformatted for GrADs display
c-- dimension parameters are the same as what in /bin/README
c      parameter(NX=64,NY=32)
      parameter(NX=64,NY=42)
      real topo(NX,NY)
      open(10,file='00000115.alb',form='formatted')
      open(20,file='00000115_alb.dat',form='unformatted')
      do j=1,NY
      do i=1,NX
        read(10,*) topo(i,j)
      enddo
      enddo
     
      write(20) topo
      stop
      end 
