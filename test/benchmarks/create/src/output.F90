! This file contains the output subroutines for qtcm model
!
!   Add (skip) a defVar call in subroutine defOutVars to add
!   (skip) an output variable. 
!      
!   Example
!      Call defVar('u1     Mode 1 zonal velocity       [m/s]',u1)
!   For adding make sure that the variable is contained in module OutputVars
!
!   Note: The minimum value for <ntout> and <ntouti> is 1 day, so if you 
!   need to output more often than a day, you need to call the output routines 
!   from within the atmosphere-land coupling loop in [qtcm]
!
! ---------------------------------------------------------------------
!
Module OutputVars
  Use Dimensions, Only: nxy,nx,ny

  Implicit None
  ! linked list type for output variables 
  Type outvar
    Real, Pointer, Dimension(:,:) :: var          ! pointer to output variable
    Real, Dimension(nxy)          :: sum          ! storage for mean values
    Integer                       :: n            ! # entries in varsum
    Type(outvar), Pointer         :: next         ! pointer to next outvar
    Character (LEN=10)            :: name         ! name of the variable
  End Type
  Type (outvar), Pointer :: first, last
End Module OutputVars
!
! -------------------------------------------------------------------
!
Subroutine defOutVars
  ! Define Output variables using the output format dependent routine [defVar].
  ! the argument of [defVar] should be a string of the form
  ! '<variable name> <variable description> [<variable units>]'
  ! "defvar" will searches for the first blank to get the variable name
  ! for the brackets to retrieve the units.
  !
  Use Prognostic
  Use Baroclinic, Only: advT1,advq1, GMq1, div1, advwu1, advwv1, dfsT1 &
       &               , dfsq1, GMs1
  Use Barotropic, Only: vort0,psi0 , dfsu0, dfsv0, advu0,advv0  &
       & , advwu0, advwv0, div0
  Use Fluxes
  Use Surface  
  Use AuxVars
  Use Land,       Only: Evapi, WD, Runf, wet, Runs
  Use Clrad,      Only: cl1

  Implicit None
  Interface
     Subroutine defVar(longvarname,var)
       Use Dimensions, Only: nx,ny
       Real, Dimension(nx,ny), Target :: var
       Character (len=*)    :: longvarname
     End Subroutine defVar
  End Interface

  Call defVar('u1     Mode 1 zonal velocity       [m/s]', u1(1:nx,1:ny))
  Call defVar('v1     Mode 1 meridional velocity  [m/s]',v1(1:nx,1:ny))
  Call defVar('T1     Mode 1 temperature          [K]',T1(1:nx,1:ny))
  Call defVar('q1     Mode 1 humidity             [K]',q1(1:nx,1:ny))

  Call defVar('u0     Barotropic zonal velocity   [m/s]',u0(1:nx,1:ny))
  Call defVar('v0     Barotropic meridional velocity [m/s]',v0(1:nx,1:ny))

  Call defVar('vort0  Barotropic vorticity        [1/s]',vort0)
  Call defVar('psi0   Barotropic stream function  [m^2/s]',psi0(1:nx,1:ny))
  Call defVar('Ts     Surface temperature         [K]',Ts)
  Call defVar('Prec   Precipitation               [W/m^2]',Qc)

  Call defVar('cl1    Deep Cloud Amount               [-]',cl1)
  Call defVar('S0     Incoming Solar Radiation    [W/m^2]',S0)
  Call defVar('FSWds  Downward Surface SW Flux    [W/m^2]',FSWds)
  Call defVar('FSWus  Upward Surface SW Flux      [W/m^2]',FSWus)
  Call defVar('FSWut  Upward SW Flux at Top       [W/m^2]',FSWut)

  Call defVar('FLWds  Downward Surface LW Flux    [W/m^2]',FLWds)
  Call defVar('FLWus  Upward Surface LW Flux      [W/m^2]',FLWus)
  Call defVar('OLR    Upward LW Flux at Top       [W/m^2]',FLWut)
  Call defVar('FTs    Surface Sensible Heat Flux  [W/m^2]',FTs)
  Call defVar('Evap   Surface Latent Heat Flux    [W/m^2]',Evap)

  Call defVar('Runf   Run-off                     [W/m^2]',Runf)
  Call defVar('WD     Equivalent water depth      [kg/m^2]',WD)
  Call defVar('stype  Surface type                [-]',STYPE)
  Call defVar('taux   Zonal Surface Stress        [N/m^2]',taux)
  Call defVar('tauy   Meridional Surface Stress   [N/m^2]',tauy)

  Call defVar('advT1  -T1 Advection               [K/s]',advT1)
  Call defVar('advq1  -Q1 Advection               [m/s^2]',advq1)

  Call defVar('Evapi  Interception loss           [W/m^2]',Evapi)
  Call defVar('GMq1   Gross moist stratification  [K]',GMq1)
  Call defVar('wet    WD/WD0 Relative wetness     [-]',wet)
  Call defVar('Runs   Surface runoff              [W/m^2]',Runs)

  If(arr1name /= '?') Call defVar(arr1name,arr1)
  If(arr2name /= '?') Call defVar(arr2name,arr2)
  If(arr3name /= '?') Call defVar(arr3name,arr3)
  If(arr4name /= '?') Call defVar(arr4name,arr4)
  If(arr5name /= '?') Call defVar(arr5name,arr5)
  If(arr6name /= '?') Call defVar(arr6name,arr6)
  If(arr7name /= '?') Call defVar(arr7name,arr7)
  If(arr8name /= '?') Call defVar(arr8name,arr8)

  !  Call defVar('slp    Surface geopotential        [Pa]',ps)
  !  Call defVar('dpsdx  x-grad of surface geopot.   [m/s^2]',dphisdx)
  !  Call defVar('dpsdy  y-grad of surface geopot.   [m/s^2]',dphisdy(1:nx,1:ny))

  !  Call defVar('advu0   -u0 Advection              [m/s^2]',advu0)
  !  Call defVar('advv0   -v0 Advection              [m/s^2]',advv0)
  !  Call defVar('advwu0  -u0 W-Advection            [m/s^2]',advwu0)
  !  Call defVar('advwv0  -v0 W-Advection            [m/s^2]',advwv0)
  !  Call defVar('advwu1  -u1 W-Advection            [m/s^2]',advwu1)
  !  Call defVar('advwv1  -v1 W-Advection            [m/s^2]',advwv1)

     Call defVar('dfsT1   T1 Diffusion               [K/s]',dfsT1)
     Call defVar('dfsq1   q1 Diffusion               [K/s]',dfsT1)
  !  Call defVar('dfsu0   u0  Diffusion              [m/s^2]',dfsu0)
  !  Call defVar('dfsv0   v0  Diffusion              [m/s^2]',dfsv0)

     Call defVar('GMs1   Dry static stability        [K]',GMs1)
     Call defVar('div0   mode 0 divergence           [1/s]',div0)
     Call defVar('div1   mode 1 divergence           [1/s]',div1)
  !  Call defVar('top    Surface topography          [500m]',TOP)

  !  Call defVar('cdn    Surface drag coefficient    [-]',CDN)
  !  Call defVar('ub     top of mixed layer zonal velocity      [m/s]',ub)
  !  Call defVar('vb     top of mixed layer merid. velocity     [m/s]',vb)
     Call defVar('us     surface zonal velocity      [m/s]',us)
     Call defVar('vs     surface merid. velocity     [m/s]',vs)

  Return
End Subroutine DefOutVars
!
! -------------------------------------------------------------------
! -------------------------------------------------------------------
! Usually no changes should be necessary below this line
! -------------------------------------------------------------------
!
Subroutine SetOutVar(longvarname,var)
  !
  ! Adds variable var with name specification longvarname
  ! to the linked list of output variables
  !

  Use OutputVars
  Implicit None

  ! Input
  Character (len=*)    :: longvarname
  Real, Dimension(nx,ny), Target :: var

  ! Local
  Type(outvar), Pointer :: newvar
  Integer, Save :: ncall=0
  Integer :: i,i1

     
  ncall=ncall+1
  Allocate(newvar)

  ! Set the start of the linked list
  If(ncall == 1) first => newvar

  ! Set output var values
  newvar%n=0
  newvar%sum =0.0 ! array syntax
  newvar%var=>var
 

  ! use first word in longvarname as variable name
  i1=Index(longvarname,' ')-1
  If(i1 == -1) i1=Len(longvarname)
  If(i1.Ge.8) Then
     Print*,'SetOutVar: Error variable name too long!'
     Print*,'Variable: "',longvarname,'"'
     Stop
  Endif

  newvar%name=longvarname(1:i1)

  ! Link the output list
  If (Associated(last) ) last%next=>newvar
  ! Set the last entry to newvar
  last=> newvar  

  If(ncall ==1) Then
     Print*
     Print*,'Output fields:'
  Endif
  
  Print *, '"',Trim(last%name),'": ',Trim(longvarname(i1+1:Len(longvarname)))
  Return
End Subroutine SetOutVar
!
! -------------------------------------------------------------------
!
Subroutine outpAll
  !
  ! Output means and instantaneous values;
  ! called by the driver.
  !
  Use Input, Only : ntouti, ntout, ntoutr, noout, nooutr
  Use Calendar

  Implicit None

  Logical outTime

  If (outTime(ntouti,noout)) Call outpInst ! output instanteneous values 
  If (outTime(ntout ,noout)) Call outpMean ! output mean values
  If (outTime(ntoutr,nooutr)) Call out_restart  ! output restart file
!  If (dayofyear==DAYSPERYEAR) Call out_restart ! End of year restart



  Return
End Subroutine outpAll
!
! -------------------------------------------------------------------
!
Subroutine outpInst
  !
  ! Output instantaneous fields to file; 
  !
  Use OutputVars

  Implicit None
  Type(outvar), Pointer :: var

  var=>first
  Do
     Call writeI(var)
     If(Associated(var%var,last%var)) Exit
     var=> var%next
  End Do
  Return
End Subroutine outpInst
!
! -------------------------------------------------------------------
!
Subroutine outpMean
  !
  ! Output means to file 
  !
  Use OutputVars

  Implicit None
  Type(outvar), Pointer :: var

  var=>first
  Do
     Call writeM(var)
     If(Associated(var%var,last%var)) Exit
     var=> var%next
  End Do
  Return
End Subroutine outpMean
!
! -------------------------------------------------------------------
!
Subroutine varMean
  !
  ! accumulate (sum up) all variables in the linked list
  ! of output variables
  !
  Use OutputVars
  Use Input,    Only: noout, ntout, it, mt0,nastep
  Use Calendar, Only: dayofmodel

  Implicit None
  Type(outvar), Pointer :: var

  If(dayofmodel <= noout .Or. ntout==0 ) Return
  var=>first
  ! print*,'varMean: first var name = ',var%name
  Do  
  ! loop over list of output variables
     Call oacc(var)
     If(Associated(var%var,last%var)) Exit ! var was last in linked list
     var=> var%next
  End Do
  Return
End Subroutine varMean
!
! -------------------------------------------------------------------
!
Subroutine oacc(var)
  !
  ! Accumulate (sum up) for variable var 
  !
  Use OutputVars
  Implicit None
  Type(outvar) :: var

  ! count
  var%n=var%n+1
  
  Call resetVarPointer(var)
  Call sumUp(nxy,var%var(1:nx,1:ny),var%sum)

  Return
End Subroutine oacc
!
! -------------------------------------------------------------------
!
Subroutine resetVarPointer(var)
  !
  ! Reset pointer to some variable. 
  ! This is only necessary for variable whose storage location
  ! changes during model execution,
  ! i.e., the baroclinic prognostic variable which are alternating
  ! stored in two arrays in the leapfrog implementation.
  !
  Use Prognostic
  Use Barotropic, Only: psi0
  Use OutputVars
  Implicit None
  Type(outvar) :: var

  If     (var%name == 'u1') Then
    var%var=>u1
  Else If(var%name == 'v1') Then
    var%var=>v1(1:nx,1:ny)
  Else If(var%name == 'T1') Then
    var%var=>T1
  Else If(var%name == 'q1') Then
    var%var=>q1
  Else If(var%name == 'v0') Then
    var%var=>v0(1:nx,1:ny)
  Else If(var%name == 'psi0') Then
    var%var=>psi0(1:nx,1:ny)
  Endif
  Return
End Subroutine resetVarPointer
!
! -------------------------------------------------------------------
!
Subroutine sumUp(nxy,var,sum)
  !
  ! add array var to sum
  !
  Implicit None
  Integer i,nxy
  Real var(nxy), Sum(nxy)
  Do i=1,nxy  ! array syntax does not vectorize :(
     Sum(i)=Sum(i)+var(i)
  Enddo
  ! sum=sum+var ! array syntax (slower)
  Return
End Subroutine sumUp
!
! -------------------------------------------------------------------
!
Subroutine geto(varmean,var)
  !
  ! pass the average in pmean and reset psum and the counter npsum
  ! 
  Use OutputVars
  Implicit None
  Real, Dimension(nxy) :: varmean
  Type(outvar) :: var

  !- averaging; reset the counter and the array
  If(var%n /=0) Call doMean(nxy,var%n,var%sum,varmean)
  Return
End Subroutine geto
!
! -------------------------------------------------------------------
!
Subroutine doMean(nxy,nsum,sum,mean)
  !
  ! compute mean from sum of np values
  ! and reset counter and array
  !
  Implicit None
  Integer nxy,nsum,i
  Real Sum(nxy), mean(nxy)
  Real divisor
  If(nsum == 0) Then
     Stop 'doMean: nsum =0'
  Else
     divisor=1./nsum
  Endif
  Do i=1,nxy  ! array
     mean(i)=Sum(i)*divisor 
     Sum(i)=0.
  Enddo
  nsum=0
  Return
End Subroutine doMean
!
! -------------------------------------------------------------------
!
Logical Function outTime(ntout,noout)
  !
  !  Decide if the current model time is an output time
  !
  Use Calendar, Only: dayofmonth, monthofyear,monlen, dayofmodel
  Use Input,    Only: SSTmode

  Implicit None
  Integer  ntout,noout

  outTime=.False.
  If(ntout==0) Return  ! No output flag

  If(dayofmodel<=noout) Return  !no output the first noout days

  If(ntout==-30) Then     ! output at end of month after noout days
     If (SSTmode=='perpetual') then
        ntout=monlen(monthofyear)
     endif
     outTime=(dayofmonth==monlen(monthofyear))
  Else                    ! output every <ntout> days after noout days
     outTime=(Mod(dayofmodel-noout,ntout)==0)
  Endif

  Return
End Function outTime
!
! ===================================================================

#ifdef NETCDFOUT
!
! Routines for netCDF output
!

Module ncFileData

  Integer:: ncid,dims(3)
  Logical :: isinit
  
End Module ncFileData
!
! -------------------------------------------------------------------
!
Subroutine outpInit 
  !
  ! NetCDF version of outpInit
  !
  ! For netCDF output the files are open at the first output time.
  ! With isInit='.true.' All that is done by defOutVars is
  ! setting up the (linked) list of output arrays to store the
  ! accumulated values of the variables for mean value output.
  ! 
  Use ncFileData, Only: isinit
  !
  isinit=.True.
  ! defOutVars calls defVar which for isinit=.true.
  ! only sets the variable var0, which indicates the fields to average
  Call defOutVars 
  isinit=.False.
  Return
End Subroutine outpInit
!
! ---------------------------------------------------------------------
!
Subroutine ncInit(ncfnm,nlon,nlat,ismean)
  !
  ! For netCDF output
  !
  ! Called at the first output to a netCDF file.
  !
  ! Purpose:
  ! Open netCDF file and write metadata and global information 
  ! about the QTCM run to the netCDF file.
  !
  ! If  ismean =='true: the mean value netCDF file,
  ! else              : the instantaneous value netCDF file

  Use ncFileData
  Use Input
  Use Calendar
  Use Grid, Only : latt,lont 
  
  Implicit None
  Include 'netcdf.inc'
  
  Integer, Intent (IN)               :: nlat,nlon ! grid dimension
  Character (LEN=305), Intent (OUT)  :: ncfnm     ! netCDF file name
  Logical, Intent (IN)               :: ismean    ! mean file flag

  Real               :: scr(nlon)
  Integer            :: id,iret,dim_STYPE(2),dim_lat,dim_lon,dim_tim,i,lnblnk  !,ncid,dims(3)
  Character (LEN=8)  :: date
  Character (LEN=10) :: time
  Character (LEN=21) :: tstring
  Logical            :: exans
  
  iret=0
  If (runname(1:7).Eq.'runname') Then
     Call Date_and_time(date,time)
     runname=date(3:8)//'_'//time(1:4)
  Endif
  
  If(ismean) Then
     ncfnm=Trim(outdir)//'/qm_'//Trim(runname)//'.nc'
  Else
     ncfnm=Trim(outdir)//'/qi_'//Trim(runname)//'.nc'
  Endif

  Inquire (file=Trim(ncfnm),exist=exans)
  If (exans) Then 
    Print *,'Clobbering netcdf file: ',Trim(ncfnm)
    iret = nf_create(Trim(ncfnm),NF_CLOBBER,ncid)
  Else
    Print *,'NcInit: Creating netcdf file: ',Trim(ncfnm)
    iret = nf_create(Trim(ncfnm),NF_NOCLOBBER,ncid)
    If (iret/=0) Call printError('ncInit',iret)
  Endif
  !
  ! specify global attributes
  !
  i=lnblnk(title)
  iret = nf_put_att_text(ncid,NF_GLOBAL,'title',i,title)
  If (ismean) Then
     iret = nf_put_att_text(ncid,NF_GLOBAL,'Kind',11,'Mean values')
  Else
     iret = nf_put_att_text(ncid,NF_GLOBAL,'Kind',20,'Instantaneous values')
  Endif
  i=lnblnk(bnddir)
  iret = nf_put_att_text(ncid,NF_GLOBAL,'Boundary_data_directory',i,bnddir)
  i=lnblnk(outdir)
  iret = nf_put_att_text(ncid,NF_GLOBAL,'Output_directory',i,outdir)
  iret = nf_put_att_Int (ncid,NF_GLOBAL,'Landon',NF_INT,1,landon)
  i=lnblnk(SSTmode)
  iret = nf_put_att_text(ncid,NF_GLOBAL,'SSTmode',i,SSTmode)
  iret = nf_put_att_Int (ncid,NF_GLOBAL,'Mrestart',NF_INT,1,mrestart)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'dt',NF_FLOAT,1,dt)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscxu0',NF_FLOAT,1,viscxu0)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscyu0',NF_FLOAT,1,viscyu0)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscxu1',NF_FLOAT,1,viscxu1)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscyu1',NF_FLOAT,1,viscyu1)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscxT',NF_FLOAT,1,viscxT)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscyT',NF_FLOAT,1,viscyT)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscxq',NF_FLOAT,1,viscxq)
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'viscyq',NF_FLOAT,1,viscyq)
  Call Date_and_time(date)
  iret = nf_put_att_text(ncid,NF_GLOBAL,'History',30, &
       &      'Created on '//date(1:4)//'-'//date(5:6)//'-'//date(7:8)//' '//&
       &      time(1:2)//':'//time(3:4)//':'//time(5:6))
  iret = nf_put_att_text(ncid,NF_GLOBAL,'Conventions',6,'COARDS')
  iret = nf_put_att_text(ncid,NF_GLOBAL,'Source',9,'QTCM1V2.2')
  iret = nf_put_att_Real(ncid,NF_GLOBAL,'_FillValue',NF_FLOAT,1,-9999.)
  !
  ! define different dimensions
  ! 
  iret = nf_def_Dim(ncid, 'lon', nlon, dim_lon)
  iret = nf_def_Dim(ncid, 'lat', nlat, dim_lat)
  iret = nf_def_Dim(ncid, 'time', NF_UNLIMITED, dim_tim)
  !
  ! defining one-dimensional grid variables
  !
  iret = nf_def_var(ncid,'lon', NF_FLOAT, 1, dim_lon, id)
  iret = nf_put_att_text(ncid, id, 'long_name',9, 'longitude')
  iret = nf_put_att_text(ncid, id, 'units', 12, 'degrees_east')
  
  iret = nf_def_var(ncid,'lat', NF_FLOAT, 1, dim_lat, id)
  iret = nf_put_att_text(ncid, id, 'long_name', 8, 'latitude')
  iret = nf_put_att_text(ncid, id, 'units', 13, 'degrees_north')
  
  iret = nf_def_var(ncid,'time', NF_INT, 1, dim_tim, id)
  iret = nf_put_att_text(ncid,id,'long_name',10,'Model time')
  Write(tstring,'(a,i4,a,i2.2,a,i2.2)') 'days since ', yearofmodel,&
       & '-',monthofyear,'-',dayofmonth
  iret = nf_put_att_text(ncid, id, 'units', 21, tstring)
  !
  ! define 3D (lon,lat,time) fields
  !
  dims(1) = dim_lon
  dims(2) = dim_lat
  dims(3) = dim_tim
  !
  ! define output variables defOutVars is in output.F90
  ! It calls subroutine defVars below for each output variable
  !
  Call defOutVars
  !
  ! define 2D (lon,lat) fields
  !
  dim_STYPE(1) = dim_lon
  dim_STYPE(2) = dim_lat
  iret=nf_def_var(ncid,'stype',NF_INT2,2,dim_STYPE,id)
  iret=nf_put_att_text(ncid,id,'long_name',35,  &
       &' Surface Type                   [-]')
  iret=nf_put_att_text(ncid,id,'units',1,'-')
  ! define cdn
  iret=nf_def_var(ncid,'cdn',NF_FLOAT,2,dim_STYPE,id)
  iret=nf_put_att_text(ncid,id,'long_name',35,  &
       &'   Drag Coefficient (Neutral)   [-]')
  iret=nf_put_att_text(ncid,id,'units',1,'-')
  ! define top
  iret=nf_def_var(ncid,'top',NF_FLOAT,2,dim_STYPE,id)
  iret=nf_put_att_text(ncid,id,'long_name',35,  &
       &'   Topography                [500m]')
  iret=nf_put_att_text(ncid,id,'units',4,'500m')
  
  
  iret = nf_enddef(ncid)
  If (iret/=0) Call printError('ncInit',iret)
  
  iret = nf_inq_varid(ncid,'lon',id)
  If (iret/=0) Call printError('ncInit',iret)
  iret = nf_put_var_Real(ncid,id,lont)  
  If (iret/=0) Call printError('ncInit',iret)
  
  iret = nf_inq_varid(ncid,'lat',id)
  iret = nf_put_var_Real(ncid,id,latt)
  Print*,'ncInit: lat(1)=',latt(1),', latt(',nlat,') = ',latt(nlat)
  
  iret = nf_Close (ncid)
  If (iret/=0) Call printError('ncInit',iret)
  Return
End Subroutine ncInit
!
! ---------------------------------------------------------------------
!
Subroutine writeStype(ncfnm)
  !
  ! For NetCDF output
  ! 
  ! Write non-changing variables, like surface type, drag coefficient
  ! and topography only once at the first outpAll call.
  !
  Use Surface, Only: STYPE,CDN,TOP,nx,ny
  ! write stype to netcdf file as integer map
  Implicit None
  Include 'netcdf.inc'

  Integer :: start(2),Count(2),iret,ncid,id,i,j
  Integer (KIND=2) :: iSTYPE(nx,ny)
  Character (LEN=305), Intent (IN) :: ncfnm

  start     = (/1,1/)
  count     = (/nx,ny/)

  iret = nf_Open (Trim(ncfnm), NF_Write, ncid)
  If(iret>0) Print*,'writeStype: error opening file'
  iret = nf_inq_varid(ncid,'stype',id)
  !  print*,'STYPE id = ', id
  Do j=1,ny
     Do i=1,nx
        iSTYPE(i,j)=STYPE(i,j)
     Enddo
  Enddo
  ! print'(4x,64i1)',((iSTYPE(i,j),i=1,64),j=32,1,-1)
  iret = nf_put_var_int2(ncid,id,iSTYPE)
  ! also add surface drag coefficients
  iret = nf_inq_varid(ncid,'cdn',id)
  iret = nf_put_var_real(ncid,id,CDN)
  iret = nf_inq_varid(ncid,'top',id)
  iret = nf_put_var_real(ncid,id,TOP)
 
  iret = nf_Close (ncid)

  Return
End Subroutine writeStype
!
! ---------------------------------------------------------------------
!
Subroutine writeM(var)
  !
  ! NetCDF version of writeM
  ! 
  ! 
  ! Purpose: Write mean values of var to the netCDF output file
  !
  Use OutputVars
  Use Dimensions
  Use Calendar
  Use Input, Only: ntout

  Implicit None
  Include 'netcdf.inc'

 
  Type (outvar), Intent (IN)  :: var

  Character (LEN=305), Save :: ncfnm
  Integer, Save :: ncall=0, dayofmodel0
  Integer       :: iret,ncid,id,start(3),Count(3)
  Real          :: varmean(nx,ny)

  If(ncall==0) Then
     Call ncInit(ncfnm,nx,ny,.True.)       !initialize netcdf file
     Call writeStype(ncfnm) ! write surface mask
     dayofmodel0=dayofmodel
  Endif

  iret = nf_Open (Trim(ncfnm), NF_Write, ncid)
  !
  ! write new time level and increment time index
  !
  If(Associated(var%var,first%var)) Then
     Print*,'writeM: Writing mean data to  "',Trim(ncfnm),'"'
     ncall=ncall+1
     iret = nf_inq_varid(ncid,'time',id)
     iret = nf_put_var1_Int(ncid,id,ncall,dayofmodel-dayofmodel0-Abs(ntout/2))
  Endif
     !
     ! write data and zero data arrays
     !
  iret = nf_inq_varid(ncid,Trim(var%name),id)
  If (iret /= 0) Then
     !        if(ncall==1) Print*,'writeM: "',Trim(var%name)  &
     !             & , '  " not found in file "',Trim(ncfnm),'"','. Skipped'
     iret = nf_Close (ncid)
     Return
  Endif
  start     = (/1,1,ncall/)
  count     = (/nx,ny,1/)
  
  If(ncall==1) Then
     Print*,'writeM: Writing "',Trim(var%name),'".'
!     if(associated(var%next) ) then
!        print*, 'var%next is associated, var%n = ',var%n
!     else
!        print*, 'var%next is NOT associated, var%n = ',var%n 
!     endif
  Endif
  
  Call geto(varmean,var)

  
  
  iret = nf_put_vara_Real(ncid,id,start,count,varmean)
  
  If (iret /= 0) Stop 'ncwrite: error writing variable'
    iret = nf_Close (ncid)
  
  Return
End Subroutine writeM
!
! ---------------------------------------------------------------------
!
Subroutine writeI(var)
  !
  ! NetCDF version of writeI
  ! 
  ! Purpose: Write instantaneous values of var to the netCDF output file
  !
  Use OutputVars
  Use Dimensions
  Use Calendar
  
  Implicit None
  Include 'netcdf.inc'
  
  Type (outvar)   , Intent (IN)   :: var

  Character (LEN=305), Save        :: ncfni
  Integer, Save :: ncall=0, dayofmodel0, tcount=0
  Integer       :: iret,ncid,id,start(3),Count(3)

  If(ncall==0) Then
     Call ncInit(ncfni,nx,ny,.False.)  ! initialize netcdf file
     Call writeStype(ncfni)         ! write surface mask
     dayofmodel0=dayofmodel
  Endif
  iret = nf_Open (Trim(ncfni), NF_Write, ncid)
  !
  ! write new time level and increment time index
  !
  If(var%name ==first%name) Then
     Print*,'writeI: Writing inst. data to  "',Trim(ncfni),'"'
     ncall=ncall+1
     iret = nf_inq_varid(ncid,'time',id)
     If( (ncall==2) .And. (dayofmodel==dayofmodel0))  &
          & tcount=1  ! use call number as time increment
     If(tcount==1) Then
        iret = nf_put_var1_Int(ncid,id,ncall,ncall)
     Else
        iret = nf_put_var1_Int(ncid,id,ncall,dayofmodel-dayofmodel0)
     Endif
  Endif
  !
  ! write data and zero data arrays
  !
  iret = nf_inq_varid(ncid,Trim(var%name),id)
  If (iret /= 0) Then
     !        if(ncall==1) Print*,'writeI: "',Trim(var%name)     &
     !     & ,  '" not found in file "',Trim(ncfni),'"','. Skipped'
     iret = nf_Close (ncid)
     Return
  Endif
  If(ncall==1) Print*, 'writeI: Writing "',Trim(var%name),'".'
  start     = (/1,1,ncall/)
  count     = (/nx,ny,1/)
  iret = nf_put_vara_Real(ncid,id,start,count,var%var(1:nx,1:ny))  
  If (iret /= 0) Then
     Print*, 'writeI: error writing variable id = ',id,', return = ',iret
     iret = nf_Close (ncid)
     Stop
  End If
  iret = nf_Close (ncid)
  
! print*,'writeI: "',Trim(var%name),'" written.'
  Return
End Subroutine writeI
!
! ---------------------------------------------------------------------
!
Subroutine defVar(longvarname,var)
  !
  ! NetCDF version of defVar
  ! 
  ! Purpose:
  !    (i) Call by outpInit (isinit ='true')
  !        Add the variable in array "var" with description
  !        "longvarname" to the (linked) list of mean output variables
  !    (ii) At the first call of outpAll (isinit = 'false')
  !         Define the variable "var" with description "longvarname"
  !         in the netCDF file.
  ! 
  Use Dimensions, Only: nx,ny
  Use OutputVars
  Use ncFileData
  Implicit None
  Interface
     Subroutine SetOutVar(longvarname,var)
       Use Dimensions, Only: nx,ny
       Real, Dimension(nx,ny), Target :: var
       Character (len=*)    :: longvarname
     End Subroutine SetOutVar
  End Interface

  Include 'netcdf.inc'
  
  Real, Dimension(nx,ny), Target :: var
  Character (len=*)    :: longvarname
  
  Integer :: id,i1,i2,len,lnblnk,iret


  ! Don't define non-changing field again or add them to the list of
  ! mean value output variables
  If( (longvarname(1:5)=='stype') .Or. (longvarname(1:5)=='STYPE') .Or. &
   &  (longvarname(1:5)=='cdn')   .Or. (longvarname(1:5)=='CDN')   .Or. &
   &  (longvarname(1:5)=='top')   .Or. (longvarname(1:5)=='TOP') ) Return

  If(isinit) Then
     ! Add var to OutputVars list
     Call SetOutVar(longvarname,var)
     Return
  Endif
  

  ! use first word in longvarname as variable name
  i1=Index(longvarname,' ')
  If(i1.Le.1) i1=lnblnk(longvarname)
  If(i1.Ge.8) Then
     Print*,'ncInit: Error variable name too long!'
     Stop
  Endif
  
  iret=nf_def_var(ncid,longvarname(1:i1-1),NF_FLOAT,3,dims,id)
  !  print*,'defVar: name="',longvarname(1:i1-1),'", netcdf ID =',id
  If (iret/=0) Call printError('ncdefavr',iret)
  
  i2=lnblnk(longvarname)
  iret=nf_put_att_text(ncid,id,'long_name',i2-i1,longvarname(i1+1:i2))
  If (iret/=0) Call printError('ncdefavr',iret)
  
  ! get units string out of brakets in longvarname
  i1=Index(longvarname,'[')+1
  i2=Index(longvarname,']')-1
  len=i2-i1+1
  If(len.Gt.1) Then
     iret=nf_put_att_text(ncid,id,'units',len,longvarname(i1:i2))
  Else
     iret=nf_put_att_text(ncid,id,'units',7,'unknown')
  Endif
  If (iret/=0) Call printError('ncdefavr',iret)
  
  Return
End Subroutine defVar
!
! ---------------------------------------------------------------------
!
Subroutine printError(routine,iret)
  !
  ! NetCDF output error message routine
  !
  Implicit None
  Include 'netcdf.inc'
  
  Integer, Intent (IN) :: iret
  Integer              :: lnblnk
  Character (len=*)    :: routine
  Character (len=305)   :: error_string
  
  Print*,'NetCDF:  netcdf write error: no.',iret
  error_string=nf_strerror(iret)
  Print*,routine,': "',error_string(1:lnblnk(error_string)),'"'
  Stop
End Subroutine printError
!
! end of NETCDF output routines
! ---------------------------------------------------------------------
! ---------------------------------------------------------------------
#else
! .not. NETCDFOUT => Grads binary output
!
! output routines for grads binary output
!
Module DefCalls
  ! count calls to DefVar for ctl-files.
  Integer :: nDefVarCalls ! Counter of the number of output variables
  Logical :: countVars    ! Switch to turn on variable contounting with
                          ! defVar 
  Character :: defType    ! output data file ('i' or 'm')
End Module DefCalls
!
! ---------------------------------------------------------------------
!
Subroutine outpInit
  !
  ! Sequential unformatted (GrADS output) version of outpInit
  ! 
  !  Write descriptor files for instantaneous valuse and means
  !
  Use Input,    Only : ntout,ntouti

  !  write grads descriptor (ctl) files
  If( ntouti /= 0 ) Call gaCtl('i') ! instantenous 'qi...out'
  If( ntout  /= 0 ) Call gaCtl('m') ! means        'qm...out'
  Return
Contains
  !
  ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
  !
  Subroutine gaCtl(Type)
    !
    ! For sequential unformatted (GrADS output)
    !
    ! Write GrADS descriptor file 
    !     Type='m': for mean values
    !     Type='i': for instantaneous values
    !
    Use Input,    Only : noout,ntout,ntouti,lastday,outdir,runname,title
    Use Calendar, Only : monlen, daysperyear
    Use Dimensions
    Use Grid, Only  : yb
    Use DefCalls

    Implicit None
    Character (len=1), Intent(in) :: Type  ! 'i'/'m' for instantenous/mean data
    Character (len=305) :: fname
    Character(len=3), Dimension(12) ::  &
         &        moname=(/ 'JAN','FEB','MAR','APR','MAY','JUN' &
         &                 ,'JUL','AUG','SEP','OCT','NOV','DEC' /)
    Real    :: dlat
    Integer ::  syear, smonth, sday ! start date

    defType=Type
    !
    ! open descriptor file
    !
    If (runname(1:7)=='runname')  runname='.'
    fname=Trim(outdir)//'/q'//Type//'_'//Trim(runname)//'.ctl'
    Print*,'gaCtl: output file name = ',Trim(fname)
    Open(24,file=fname)
    nDefVarCalls=0
    countVars=.True.
    Call defOutVars  ! here only used for counting output vars
    countVars=.False.
    Rewind(24)
    fname='^q'//Type//'_'//Trim(runname)//'.out'

    Print*,'gaCtl: data fname = ',Trim(fname)
    Write(24,'(a,a)')  'DSET    ',Trim(fname)              &
         &           , 'TITLE   ',Trim(title)              &
         &           , 'UNDEF   ','-9.99E23'               &
#ifdef __IFC
  ! __IFC is defined for Intel's IFC compiler which only exist for little_endian
  !machines
         &           , 'OPTIONS ','sequential little_endian'
#else
         &           , 'OPTIONS ','sequential'  ! add big_endian or little_endian
#endif
    !
    ! grid information
    Write(24,'(a,i3,a,f8.3)')  &
         &      'XDEF    ',nx,' linear 0 ',(360000/nx)/1000.
    dlat=2*yb/ny
    Write(24,'(a,i3,2(a,f8.3))') 'YDEF    ',ny,' linear ',-yb+dlat/2,' ',dlat
    Write(24,'(a)') 'ZDEF      1 linear 0 1'
    !
    ! time dimension
    If(Type=='m') Then
       ! compute day and month of first output
       Call startDate(ntout,sday,smonth,syear)
       If(ntout==-30) Then
          Write(24,'(a,i5,a,a,i4.4,a)' ) 'TDEF  ' & 
               &       ,((lastday-noout)/daysperyear)*12+&
               &    Int(Mod(lastday-noout,daysperyear)/30.0+0.5),&
               &       ' linear 15',moname(smonth),syear,' 1mo'
       Else
          Write(24,'(a,i5,a,i2.2,a,i4.4,a,i3,a)' )  &
               & 'TDEF ' ,(lastday-noout)/ntout,' linear ',  &
               & sday-ntout/2,moname(smonth),syear,' ',ntout,'dy'
       Endif
    Else
       Call startDate(ntouti,sday,smonth,syear)
       If(ntouti==-30) Then
          Write(24,'(a,i5,a,i2.2,a,i4.4,a)' ) 'TDEF  ' & 
               & ,((lastday-noout)/daysperyear)*12+&
               &     Int(Mod(lastday-noout,daysperyear)/30.0+0.5),&
               &       ' linear ',monlen(smonth),moname(smonth),syear,' 1mo'
       Else
          Write(24,'(a,i5,a,i2.2,a,i4.4,a,i3,a)' )  &
               & 'TDEF ' ,(lastday-noout)/ntouti,' linear ',  &
               & sday,moname(smonth),syear,' ',ntouti,'dy'
       Endif
    Endif
    !
    ! Variables
    Write(24,'(a,i3)') 'VARS ',nDefVarCalls
    !
    ! Subroutine defOutVars in output.F90 defines the  output variables.
    ! It calls subroutine defVars below for each output variable.
    Call defOutVars
    Write(24,'(a)') 'ENDVARS'
    Close(24)
    Return
  End Subroutine gaCtl

  Subroutine startDate(dday,sday,smonth,syear)
    ! 
    ! For sequential unformatted (GrADS output)
    ! Purpose: Compute output starting date for GrADS descriptor file (.ctl)
    ! 
    Use Calendar, Only: day0,month0,year0,monlen,cummonth, daysperyear
    Use Input,    Only: noout, SSTmode
    Implicit None
    Integer, Intent(in)  :: dday ! output interval in days
    Integer, Intent(out) :: sday, smonth,syear
    Integer :: i,dnday
    !
    ! For perpetual output
    if(SSTmode == 'perpetual') then
       sday=day0
       syear=year0
       smonth=month0
       return
    endif
    ! compute day and month of first output
    !
    if(dday==-30) then ! monthly output
       dnday=monlen(month0)
    else
       dnday=Abs(dday)
    endif
    sday=cummonth(month0)+noout+dnday+day0-1
    ! year
    syear=year0
    Do i=1,100000 ! While (sday< daysperyear) ! 04/02: Alpha compiler choked.
       If (sday< daysperyear) Exit
       sday=sday-daysperyear
       syear=syear+1
    Enddo
    ! month
    smonth=12
    Do i=1,12 ! While(cummonth(smonth)<sday)
       If (cummonth(smonth)<sday) Exit
       smonth = smonth-1
    Enddo
    sday=sday-cummonth(smonth)-1
    Print*,'startDate: dd:mm:yy = ',sday,':',smonth,':',syear
    ! Print*,'cummonth(smonth) = ',cummonth(smonth)
    ! Print*,'dday = ',dday,', day0=',day0
  End Subroutine startDate
End Subroutine outpInit
!
! ---------------------------------------------------------------------
!
Subroutine defVar(longvarname,var)
  !
  ! Sequential unformatted (GrADS output) version of defVar
  !
  ! Purpose: 
  !     (i) if countVars is "true":
  !         Increment the counter of output variables.
  !     (ii) Else:
  !          Add variable in array "var" to the GrADS descriptor
  !          file and, 
  !          iff we are setting up mean value output (defType == 'm'):
  !          Add "longvarname" to the (linked) list of output variables
  !          by calling SetOutVar
  Use DefCalls, Only: countVars, nDefVarCalls,defType
  Use Dimensions
  Implicit None

  Interface
     Subroutine SetOutVar(longvarname,var)
       Use Dimensions, Only: nx,ny
       Real, Dimension(nx,ny), Target :: var
       Character (len=*)    :: longvarname
     End Subroutine SetOutVar
  End Interface

  Character (len=*), Intent(in) :: longvarname
  Real, Dimension(nx,ny), Target, Intent(in) :: var

  Integer :: i1,i2,lnblnk

  If(countVars) Then
     nDefVarCalls=nDefVarCalls+1
  Else
  ! use first word in longvarname as variable name
  i1=Index(longvarname,' ')
  If(i1.Le.1) i1=lnblnk(longvarname)
  i2=lnblnk(longvarname)
  Write(24,'(3a)') longvarname(1:i1-1),'   0 99 ',longvarname(i1+1:i2)
  If(defType=='m') Call SetOutVar(longvarname,var)
Endif

  Return
End Subroutine defVar

!
! ---------------------------------------------------------------------
!
Subroutine WriteM(var)
  !
  ! Sequential unformatted (GrADS output) version of WriteM
  ! 
  ! Purpose:  Write mean values of output variable "var"
  !
  Use Dimensions
  Use OutputVars
  Use Input, Only : outdir,runname

  Implicit None  
  Type (outvar), Intent(in)  :: var
  Real, Dimension(nx,ny)      :: varMn
  Logical, Save :: firstcall=.True.
  Character (len=305), Save :: fname

  If(firstcall) Then
     fname=Trim(outdir)//'/qm'//'_'//Trim(runname)//'.out'
     Open(22,file=fname,form='unformatted',access='sequential')
     firstcall=.False.
  Endif

  If (Trim(var%name) /= 'time') Then
     If (var%name == first%name) Then
        Print*,'WriteM: Writing mean data to '
        Print*,'       "', Trim(fname),'"'
     Endif
     Call geto(varMn,var)
     Write(22) varMn
  Endif
  Return
End Subroutine WriteM
!
! -------------------------------------------------------------------
!
Subroutine WriteI(var)
  !
  ! Sequential unformatted (GrADS output) version of WriteI
  ! 
  ! Purpose:  Write instantaneous values of output variable "var"
  !
  Use Dimensions
  Use OutputVars
  Use Input, Only : outdir,runname

  Implicit None
  Type (outvar), Intent(in)  :: var
  Logical, Save :: firstcall=.True.
  Character (len=305), Save :: fname

  If(firstcall) Then
     fname=Trim(outdir)//'/qi'//'_'//Trim(runname)//'.out'
     Open(21,file=fname,form='unformatted',access='sequential')
     firstcall=.False.
  Endif

  If (Trim(var%name) /= 'time') Then
     If (var%name == first%name) Then
        Print*,'WriteI: Writing inst. data to '
        Print*,'       "', Trim(fname),'"'
     Endif
     Write(21) var%var(1:nx,1:ny)
  Endif
  
  Return
End Subroutine WriteI

#endif
!
! ===================================================================

