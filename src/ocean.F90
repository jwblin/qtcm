! File: ocean.F90
!  
! This file obtains SST by reading in observed SST, or using a mixed-layer
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
Module SSTout
  Use dimensions
  REAL Tnow(nx,ny)
End Module SSTout
!
Module Rdsst
  Use dimensions
  Use SSTout
  Real Tprior(nx,ny),Tnext(nx,ny)
  Integer date1, date2
End Module Rdsst
!
Module Mxlsst
  Use dimensions
  Use SSTout
  Real Fsn(NX,NY),Dts(NX,NY),Qfx(NX,NY)
End Module Mxlsst

#ifdef BLEND_SST
! defines SST mask for using read-in observed SST (mask=1); otherwise, using 
! mixed-layer SST (mask=0)     
Module SSTmask
    Use Dimensions
    Real, Dimension(nx,ny) :: mask
End Module SSTmask
#endif
!
!---------------------------------------------------------------------
!
      Subroutine OceanInit()
!
!  Get the initial SST data for the dates on either side of an
!  arbitrary starting date
!
  Use Input
  Use Dimensions
  Use Rdsst
  Use Calendar

  Implicit None

#ifdef BLEND_SST
  Real, Dimension(NX,NY) :: Tsobs,Tsmxl 
#endif

!
      write (*,*) "OceanInit: get the initial SST data"
!
!-  get the dates(centered on the 15th) on either side of the starting date.
! before the 15th
      if(day0.ge.1.and.day0.le.14) then
        date1 = year0*10000 + (month0-1)*100 + 15
        if(month0.eq.1)date1= (year0-1)*10000 + 12*100 + 15
        date2= year0*10000 + month0*100 + 15
        if (SSTmode.eq.'real_time')  then
      write(*,*)'Reminder: your input requires SST data before the'
      write(*,*)' model starting date! if non-existant, model will stop'
      write(*,*)' Remedy: change the starting date or get the data'
        end if
      end if
! after the 15th
      if(day0.ge.15.and.day0.le.31) then
        date1= year0*10000 + month0*100 + 15
        date2 = year0*10000 + (month0+1)*100 + 15
        if(month0.eq.12)date2= (year0+1)*10000 + 1*100 + 15
      end if
      call sstin(date1, Tprior)
      Tnow = Tprior             !Fudge the initial T ! arrays
      if (SSTmode .eq. 'perpetual') return
      call sstin(date2, Tnext)

#ifdef BLEND_SST
! for initial SST data over mixed-layer ocean, read in climatological
! SST in SSTdir   
      Tsobs=Tnow
      call bndry1('sst',SSTdir,Tsmxl)
      call blendsst(Tsobs,Tsmxl)
#endif 

      Return 
End subroutine OceanInit 
!
! ---------------------------------------------------------------------
!
      Subroutine Ocean(ndays,it)
!
  Use Dimensions
  Use SSTout 

      implicit none

      integer ndays,it

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
! ---------------------------------------------------------------------
!
      Subroutine readsst(ndays,it)
  Use Input, only : SSTmode, SSTdir
  Use Dimensions
  Use Rdsst
  Use Calendar, only : dayofmodel, dateofmodel

      implicit none

      logical TimetoReadSst
      integer julian
      integer ndays
      integer month,day1,day2,it
      REAL thistime
!
!      integer i,j
!
      if (SSTmode .eq. 'perpetual') return
!
      if (TimetoReadSst(date1, date2, dateofmodel,SSTmode)) then
         Tprior = Tnext !arrays

         date1 = date2
         date2 = date1 + 100                     !advance date2 by 1 month
         month=mod(date1-15,10000)/100           !if month=12
         if(month.eq.12)date2= date1+10000-1100   !  set to Jan of next year
         call sstin(date2, Tnext)
      endif
!
!     timeinterp requires julian dates or some other monotonic time.
      day1 = julian(date1)
      day2 = julian(date2)
      thistime = dayofmodel + ndays/2.       ! compute timemean of sst
      call TimeInterp(Tprior,Tnext,Tnow ,day1,day2,thistime,nx,ny)
!     write (*,*) 'Ocean: interpolated to t=',thistime
!    $            ,'date1,day1,date2,day2=',date1,day1,date2,day2
      Return
End Subroutine readsst 
!
!---------------------------------------------------------------------
!
Subroutine mxlocean(ndays,it)
  !  compute SST using mxlayer model
  !
  Use Surface

  Implicit None
  Integer ndays,it

  Call getQflux           !get Qflux

  Call getSfcHeat(FSWds,FSWus,FLWds,FLWus,Evap,FTs)   !get averaged sfc 
  !fluxes accumulated over previous coupled interval 
  ! temporarily stored in FSWds etc.
  !  (e.g. 1day); [cplmean] has been called in [qtcm]

  Call mxstep(ndays)      !time-stepping of the mix-layer ocean

  Return
End Subroutine mxlocean 
!
!---------------------------------------------------------------------
!
#ifdef BLEND_SST
Subroutine blendsst(Tsobs,Tsmxl)
  !  combines read-in observed SST in masked regions and 
  !  mixed-layer ocean SST in non-masked regions 
  !  
  Use Dimensions
  Use Input,   only : bnddir
  Use SSTmask
  Use SSTout

  Implicit None
  Real, Dimension(NX,NY) :: Tsobs,Tsmxl
  Integer i,j,lnblnk
  Logical, Save :: firstcall=.True.
  Character (len=305), Save :: fname

  If(firstcall) Then
    fname=bnddir(1:lnblnk(bnddir))//'/ensopac.mask'
    Open(38,file=fname,form='formatted')
    Do j=1,NY
    Do i=1,NX
       read(38,*) mask(i,j)
    End Do
    End Do
    Close(38)
    firstcall=.False. 
  Endif

! use read-in observed SST in masked (mask=1) regions and mixed-layer ocean 
! SST in non-masked regions (mask=0)

  Do j=1,NY
  Do i=1,NX  
     Tnow(i,j) = Tsobs(i,j)*mask(i,j)+Tsmxl(i,j)*(1.-mask(i,j))
  End Do
  End Do
  
  Return
End Subroutine blendsst 
#endif
!
!---------------------------------------------------------------------
!
      Subroutine getSST(T_out)
!************
!  A method for returning the Sea Surface Temperature to other models
!  In this case, it simply provides a pointer to the data
!  structure of the ocean.
!
  Use dimensions
  Use SSTout 
      implicit none

      REAL T_out(nx,ny)
!************************************************************
!
        T_out = Tnow  !array
      Return
End Subroutine getSST
!
!---------------------------------------------------------------------
!
      Logical Function TimetoReadSST(date1,date2,now,SSTmode)
      Integer date1, date2, now
      Character(Len=*) :: SSTmode
      TimetoReadSST = (now .ge. date2)
      return
      end
!
!---------------------------------------------------------------------
!
      Subroutine sstin(date, Tdummy)
!
!   Return:
!     date  -   not changed
!     Tdummy  -   the new temperatures.
!
!
!   In this version, bc file names include the date and variable name
!   in this way:
!
!      file name = yyyymmdd.var
!
!   where yyyy = year , mm = month, dd=dayofmonth, var = a variable name.
!
!   if SSTmode is seasonal, yyyy=0000
!
  use Input
  use Dimensions
  use Rdsst
  use Calendar

      Implicit None

      REAL, Dimension(nx,ny) ::  Tdummy
      Integer                :: date
      Character(Len=8)       :: SSTfile
      Character(Len=305)     :: fname
! Functions
      integer lnblnk
!
      write(SSTfile,'(i8.8)') date         !convert the date into a string
      print*,'sstin: SSTfile = ',SSTfile,', date = ',date
!
! get the climatological SST for seasonal mode.
      if (SSTmode.eq.'seasonal')  then
        SSTfile='0000'//SSTfile(5:8)   !switch the year to 0000
      end if
!
! get the perpetual SST; the sst file must be 00000000.sst
      if (SSTmode.eq.'perpetual')  then
        SSTfile='00000000'
      end if
!
! read the SST
      fname=SSTdir(1:lnblnk(SSTdir))//'/'//SSTfile//'.sst'
      open(13,file=fname,form='formatted',status='old')
      read(13,*) Tdummy                                  !in Celsius
      close(13)
!
!
      write(*,*) 'sstin: Read SST at model date', dateofmodel,          &
     &            ' from file'
      write(*,*)'        "',                                            &
     &      SSTdir(1:lnblnk(SSTdir))//'/'//SSTfile//'.sst','"'
!
      Return
End Subroutine sstin     
!
! ---------------------------------------------------------------------
!
Subroutine getQflux
  !  read in Qflux 
  !  Qflux is precalculated by aveflux.f in /src using results from
  !  a control experiment with seasonal SST 
  Use Dimensions
  Use Input,      Only : SSTmode
  Use Mxlsst

  Implicit None

  Logical firstcall
  Data firstcall/.True./
  Save firstcall

  !     
  !  get Q fluxes 
  If(SSTmode.Eq.'perpetual') Then     ! constant Qflux such as Annual Mean
     If(firstcall) Then
        Open(16,file='00000000.fsn', status='old')
        Read(16,*) qfx 
        firstcall=.False.
     Endif
  Else         !for all other SSTmodes, use seasonally varying Qflux; Fs and Ts 
     Call bndry1('fsn','./',fsn)       !  tendency are treated differently
     Call bndry2('dts','./',dts)
     qfx=fsn-dts ! arrays
  Endif

  Return
End Subroutine getQflux
!
!---------------------------------------------------------------------
!
Subroutine mxstep(intcpl)
  Use Surface
  Use Input
  Use Mxlsst
#ifdef BLEND_SST
  Use SSTmask
#endif

  Implicit None
  ! --Altered Common Variables:
  ! Mxsst 
  !
  ! --Local Variables & Functions:
  Real Rsnet,Fsnet,Dmx,Cmx
  Integer i,j,intcpl
  Data Dmx/50./     !Depth of the mix-layer ocean;  meters
  Save Dmx
  !
  !  water heat capacity * density * Depth
  Cmx=4.18e6*Dmx      !1cal/g/K * 1g/cm3 * Depth;  in J/K/M2

  Do j=1,NY
     Do i=1,NX
#ifdef BLEND_SST
       If(mask(i,j).eq.0.) Then 
#endif
        If(STYPE(i,j).Eq.0. .Or. landon.Eq.0) Then
         Rsnet= FSWds(i,j)-FSWus(i,j)+FLWds(i,j)-FLWus(i,j)  !net sfc. radiation
        ! Note: FSWds etc. temporarily stores the averages of 
        !previous coupling interval
         Fsnet= Rsnet-Evap(i,j)-FTs(i,j)   !net energy absorbed by surface
         Tnow(i,j)=Tnow(i,j)+float(intcpl)*86400.*(Fsnet-qfx(i,j))/Cmx
        End If
#ifdef BLEND_SST
       End If
#endif
     End Do
  End Do

  Return
End Subroutine mxstep
!
!---------------------------------------------------------------------
!
