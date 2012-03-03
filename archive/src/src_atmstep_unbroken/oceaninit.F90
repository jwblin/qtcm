! File: oceaninit.F90
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
!=====================================================================
