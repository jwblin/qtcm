!======================================================================
! Name:
!   setbypy.F90
!
! Description:
!   Routines to set and get parameters in compiled module via Python
!
! Modification History:
! - 12 Jul 2007:  By Johnny Lin, Physics Dept., North Park University,
!   http://www.johnny-lin.com/.
!
! Notes:
! - The use of the following code:
!      Interface setitem
!        Module Procedure setitem_real, setitem_int, setitem_str
!      End Interface
!   causes a bus error on my Mac (Mac OS X 10.4, Intel).  Instead, the 
!   same functionality is included in the Python code.
!
! Copyright (c) 2007 by Johnny Lin.  For licensing, distribution 
! conditions, contact information, and additional documentation see
! the URL http://www.johnny-lin.com/py_pkgs/qtcm/doc/.
!======================================================================


Module SetbyPy
Use AuxVars
Use Barotropic
Use Calendar
Use Input
Use Prognostic
Use Qtcmpar
Use Land, only : WD, WD0
Use Surface, only : Ts, STYPE
Use Fluxes


!----------------------------------------------------------------------
! Declare module attributes that are used to pass real arrays out to
! the Python level.
!----------------------------------------------------------------------

Real, allocatable, dimension(:) :: real_rank1_array
Real, allocatable, dimension(:,:) :: real_rank2_array
Real, allocatable, dimension(:,:,:) :: real_rank3_array


!----------------------------------------------------------------------
! Declare a logical variable that is used to pass out whether or not
! a variable can be read by Python.  If is_readable is .TRUE. (the
! default set at the beginning of getitem_* and setitem_*functions), 
! Python can read the variable.  This is to prevent Python from 
! accessing pointer variables that aren't yet associated to targets.
!----------------------------------------------------------------------

Logical :: is_readable


Contains
  function getitem_real(key)
  !--------------------------------------------------------------------
  ! Get the module variable named by key, for value of real type.
  !
  ! Input:
  !   key:  String naming the variable in the module to be gotten.
  !
  ! Ouput:
  !   Returns the value of the variable of name key.  Type real.
  !   Module variable is_readable is set to .TRUE., since all
  !   module variables accessed by this function should be accessable.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
     Real :: getitem_real
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='dt') then
    getitem_real = dt
  elseif (key=='eps_c') then
    getitem_real = eps_c
  elseif (key=='u0bar') then
    getitem_real = u0bar
  elseif (key=='V1b') then
    getitem_real = V1b
  elseif (key=='viscxu0') then
    getitem_real = viscxu0
  elseif (key=='viscyu0') then
    getitem_real = viscyu0
  elseif (key=='viscxu1') then
    getitem_real = viscxu1
  elseif (key=='viscyu1') then
    getitem_real = viscyu1
  elseif (key=='visc4x') then
    getitem_real = visc4x
  elseif (key=='visc4y') then
    getitem_real = visc4y
  elseif (key=='viscxT') then
    getitem_real = viscxT
  elseif (key=='viscyT') then
    getitem_real = viscyT
  elseif (key=='viscxq') then
    getitem_real = viscxq
  elseif (key=='viscyq') then
    getitem_real = viscyq
  elseif (key=='weml') then
    getitem_real = weml
  elseif (key=='ziml') then
    getitem_real = ziml
  elseif (key=='VVsmin') then
    getitem_real = VVsmin
  else
    write(*,*) 'Error-Bad call to SetbyPy module getitem_real for ', key
    stop
  end if
  return
  end function getitem_real


  Subroutine getitem_real_array(key)
  !--------------------------------------------------------------------
  ! Get the module variable named by key, for value real array.
  !
  ! This subroutine tests module variables real_rank?_array are
  ! unallocated prior to entry into the subroutine.
  !
  ! Input:
  !   key:  String naming the variable in the module to be gotten.
  !
  ! Output:
  !   If module variable real_rank?_array can be set, sets real
  !   module variable real_rank?_array to the array of name key,
  !   and sets is_readable to .TRUE.  If output variable 
  !   real_rank?_array cannot be set (such as if the pointer variable 
  !   specified by key is not yet associated), real module variable 
  !   real_rank?_array is not set and is_readable is set to .FALSE.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
     Integer, dimension(1) :: shape_array_1d
     Integer, dimension(2) :: shape_array_2d
     Integer, dimension(3) :: shape_array_3d
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (allocated(real_rank1_array)) then
    write(*,*) 'Error-real_rank1_array should be deallocated'
    stop
  endif

  if (allocated(real_rank2_array)) then
    write(*,*) 'Error-real_rank2_array should be deallocated'
    stop
  endif

  if (allocated(real_rank3_array)) then
    write(*,*) 'Error-real_rank3_array should be deallocated'
    stop
  endif

  if (key=='Qc') then
    shape_array_2d = shape(Qc)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = Qc

  elseif (key=='FLWds') then
    shape_array_2d = shape(FLWds)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FLWds

  elseif (key=='FLWus') then
    shape_array_2d = shape(FLWus)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FLWus 

  elseif (key=='FSWds') then
    shape_array_2d = shape(FSWds)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FSWds 

  elseif (key=='FSWus') then
    shape_array_2d = shape(FSWus)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FSWus 

  elseif (key=='Evap') then
    shape_array_2d = shape(Evap)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = Evap 

  elseif (key=='FTs') then
    shape_array_2d = shape(FTs)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FTs

  elseif (key=='taux') then
    shape_array_2d = shape(taux)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = taux

  elseif (key=='tauy') then
    shape_array_2d = shape(tauy)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = tauy 

  elseif (key=='FLWut') then
    shape_array_2d = shape(FLWut)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FLWut 

  elseif (key=='FLW') then
    shape_array_2d = shape(FLW)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FLW

  elseif (key=='S0') then
    shape_array_2d = shape(S0)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = S0

  elseif (key=='FSWut') then
    shape_array_2d = shape(FSWut)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FSWut 

  elseif (key=='FSW') then
    shape_array_2d = shape(FSW)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = FSW 

  elseif (key=='u1') then
    if (associated(u1)) then
       shape_array_2d = shape(u1)
       allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
       real_rank2_array = u1
    else
       is_readable = .FALSE.
    endif

  elseif (key=='v1') then
    if (associated(v1)) then
       shape_array_2d = shape(v1)
       allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
       real_rank2_array = v1
    else
       is_readable = .FALSE.
    endif

  elseif (key=='T1') then
    if (associated(T1)) then
       shape_array_2d = shape(T1)
       allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
       real_rank2_array = T1
    else
       is_readable = .FALSE.
    endif

  elseif (key=='q1') then
    if (associated(q1)) then
       shape_array_2d = shape(q1)
       allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
       real_rank2_array = q1
    else
       is_readable = .FALSE.
    endif

  elseif (key=='u0') then
    shape_array_2d = shape(u0)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = u0

  elseif (key=='v0') then
    shape_array_2d = shape(v0)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = v0

  elseif (key=='vort0') then
    shape_array_2d = shape(vort0)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = vort0 

  elseif (key=='psi0') then
    shape_array_2d = shape(psi0)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = psi0 

  elseif (key=='rhsvort0') then
    shape_array_3d = shape(rhsvort0)
    allocate( real_rank3_array(shape_array_3d(1), shape_array_3d(2) &
       &    , shape_array_3d(3)) )
    real_rank3_array = rhsvort0 

  elseif (key=='rhsu0bar') then
    shape_array_1d = shape(rhsu0bar)
    allocate( real_rank1_array(shape_array_1d(1)) )
    real_rank1_array = rhsu0bar 

  elseif (key=='Ts') then
    shape_array_2d = shape(Ts)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = Ts 

  elseif (key=='STYPE') then
    shape_array_2d = shape(STYPE)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = STYPE

  elseif (key=='WD') then
    shape_array_2d = shape(WD)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = WD 

  elseif (key=='WD0') then
    shape_array_1d = shape(WD0)
    allocate( real_rank1_array(shape_array_1d(1)) )
    real_rank1_array = WD0

  elseif (key=='arr1') then
    shape_array_2d = shape(arr1)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr1

  elseif (key=='arr2') then
    shape_array_2d = shape(arr2)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr2 

  elseif (key=='arr3') then
    shape_array_2d = shape(arr3)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr3 

  elseif (key=='arr4') then
    shape_array_2d = shape(arr4)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr4 

  elseif (key=='arr5') then
    shape_array_2d = shape(arr5)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr5 

  elseif (key=='arr6') then
    shape_array_2d = shape(arr6)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr6 

  elseif (key=='arr7') then
    shape_array_2d = shape(arr7)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr7 

  elseif (key=='arr8') then
    shape_array_2d = shape(arr8)
    allocate( real_rank2_array(shape_array_2d(1), shape_array_2d(2)) )
    real_rank2_array = arr8 

  else
    write(*,*) 'Error-Bad call to SetbyPy module getitem_real_array for ', key
    stop
  endif
  End Subroutine getitem_real_array


  function getitem_int(key)
  !--------------------------------------------------------------------
  ! Get the module variable named key, for value of integer type.
  !
  ! Input:
  !   key:  String naming the variable in the module to be gotten.
  !
  ! Ouput:
  !   Returns the value of the variable of name key.  Type integer.
  !   Module variable is_readable is set to .TRUE., since all
  !   module variables accessed by this function should be accessable.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
     Integer :: getitem_int
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='it') then
    getitem_int = it
  elseif (key=='noout') then
    getitem_int = noout
  elseif (key=='nooutr') then
    getitem_int = nooutr
  elseif (key=='ntout') then
    getitem_int = ntout
  elseif (key=='ntouti') then
    getitem_int = ntouti
  elseif (key=='ntoutr') then
    getitem_int = ntoutr
  elseif (key=='mrestart') then
    getitem_int = mrestart
  elseif (key=='landon') then
    getitem_int = landon
  elseif (key=='nastep') then
    getitem_int = nastep
  elseif (key=='mt0') then
    getitem_int = mt0
  elseif (key=='lastday') then
    getitem_int = lastday
  elseif (key=='interval') then
    getitem_int = interval
  elseif (key=='dayofmodel') then
    getitem_int = dayofmodel
  elseif (key=='dateofmodel') then
    getitem_int = dateofmodel
  elseif (key=='yearofmodel') then
    getitem_int = yearofmodel
  elseif (key=='dayofmonth') then
    getitem_int = dayofmonth
  elseif (key=='dayofyear') then
    getitem_int = dayofyear
  elseif (key=='monthofyear') then
    getitem_int = monthofyear
  elseif (key=='year0') then
    getitem_int = year0
  elseif (key=='month0') then
    getitem_int = month0
  elseif (key=='day0') then
    getitem_int = day0
  else
      write(*,*) 'Error-Bad call to SetbyPy module getitem_int for ', key
      stop
  end if
  return
  end function getitem_int


  !@@@Subroutine getitem_str(key, value)
  function getitem_str(key)
  !--------------------------------------------------------------------
  ! Get the module variable named key, for value of character type.
  !
  ! Input:
  !   key:  String naming the variable in the module to be gotten.
  !
  ! Output:
  !   value:  The value of the variable of name key.  Type character,
  !      and of length maxitemlen=200.
  !   is_readable:  Module variable is_readable is set to .TRUE., 
  !      since all module variables accessed by this function should 
  !      be accessable.
  !
  ! Local Variables:
     Integer, Parameter :: maxitemlen=200
     Character(len=*) :: key
     Character(len=maxitemlen) :: value
     Character(len=maxitemlen) :: getitem_str !@@@
!f2py intent(in) key
!@@@!f2py intent(out) value
  !
  ! Notes:
  ! - The subroutine does a check to make sure that value has not 
  !   truncated the value of the variable named key.  The parameter 
  !   maxitemlen is a local variable.
  ! - value has to be a defined length in order for this routine
  !   to work.
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='bnddir') then
    if (len(trim(bnddir)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for bnddir'
    value = bnddir
  elseif (key=='outdir') then
    if (len(trim(outdir)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for outdir'
    value = outdir
  elseif (key=='runname') then
    if (len(trim(runname)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for runname'
    value = runname
  elseif (key=='SSTdir') then
    if (len(trim(SSTdir)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for SSTdir'
    value = SSTdir
  elseif (key=='title') then
    if (len(trim(title)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for title'
    value = title
  elseif (key=='SSTmode') then
    if (len(trim(SSTmode)) > maxitemlen) stop &
      'Error-Bad call to SetbyPy module getitem_str for SSTmode'
    value = SSTmode
  else
      write(*,*) 'Error-SetbyPy module getitem_str unknown ', key
      stop
  end if
  getitem_str = value !@@@
  !@@@End Subroutine getitem_str
  End function getitem_str


  Subroutine setitem_real(key, value)
  !--------------------------------------------------------------------
  ! Set the module variable key to value, for value of real type.
  !
  ! Input:
  !   key:    String naming the variable in the module to be changed.
  !   value:  The value to replace the module variable with.  value is
  !           scalar of type real.
  !
  ! Output:
  !   Module variable named key is set to value and set module
  !   variable is_readable to .TRUE.  If variable cannot be set (such
  !   as if the pointer variable specified by key is not yet associa-
  !   ted), is_readable is set to .FALSE. and value is not set in the
  !   model.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
     Real, intent(in) :: value
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='dt') then
    dt=value
  elseif (key=='eps_c') then
    eps_c=value
  elseif (key=='u0bar') then
    u0bar=value
  elseif (key=='V1b') then
    V1b=value
  elseif (key=='viscxu0') then
    viscxu0=value
  elseif (key=='viscyu0') then
    viscyu0=value
  elseif (key=='viscxu1') then
    viscxu1=value
  elseif (key=='viscyu1') then
    viscyu1=value
  elseif (key=='visc4x') then
    visc4x=value
  elseif (key=='visc4y') then
    visc4y=value
  elseif (key=='viscxT') then
    viscxT=value
  elseif (key=='viscyT') then
    viscyT=value
  elseif (key=='viscxq') then
    viscxq=value
  elseif (key=='viscyq') then
    viscyq=value
  elseif (key=='weml') then
    weml=value
  elseif (key=='ziml') then
    ziml=value
  elseif (key=='VVsmin') then
    VVsmin=value
  else
      write(*,*) 'Error-Bad call to SetbyPy module setitem_real'
      write(*,*) key, '  ', value
      stop
  end if
  End Subroutine setitem_real


  Subroutine setitem_real_array(key)
  !--------------------------------------------------------------------
  ! Set module variable key to value, for value a real array.
  !
  ! This subroutine will work for rank 1-3 arrays.
  !
  ! Input:
  !   key:  String naming the variable in the module to be set.
  !
  ! Output:
  !   If module variable real_rank?_array can be set, this routine
  !   assumes that another routine has already properly allocated the
  !   appropriate module variable real_rank?_array, that Python
  !   has already set that appropriate module variable to the value
  !   you want the compiled QTCM model variable to be set to, and
  !   that module variable is_readable is set to the value it should
  !   be (most of these functions can be accomplished by calling the
  !   Fortran routine getitem_real_array from Python).  Assuming all
  !   that's already done, this routine sets the compiled QTCM model
  !   module variable to the value in appropriate module variable
  !   real_rank?_array.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
  !--------------------------------------------------------------------
  if (key=='Qc') then
    Qc = real_rank2_array
  elseif (key=='FLWds') then
    FLWds = real_rank2_array
  elseif (key=='FLWus') then
    FLWus = real_rank2_array
  elseif (key=='FSWds') then
    FSWds = real_rank2_array
  elseif (key=='FSWus') then
    FSWus = real_rank2_array
  elseif (key=='Evap') then
    Evap = real_rank2_array
  elseif (key=='FTs') then
    FTs = real_rank2_array
  elseif (key=='taux') then
    taux = real_rank2_array
  elseif (key=='tauy') then
    tauy = real_rank2_array
  elseif (key=='FLWut') then
    FLWut = real_rank2_array
  elseif (key=='FLW') then
    FLW = real_rank2_array
  elseif (key=='S0') then
    S0 = real_rank2_array
  elseif (key=='FSWut') then
    FSWut = real_rank2_array
  elseif (key=='FSW') then
    FSW = real_rank2_array

  elseif (key=='u1') then
    if (associated(u1)) then
       u1 = real_rank2_array
    else
       is_readable = .FALSE.
    endif

  elseif (key=='v1') then
    if (associated(v1)) then
       v1 = real_rank2_array
    else
       is_readable = .FALSE.
    endif

  elseif (key=='T1') then
    if (associated(T1)) then
       T1 = real_rank2_array
    else
       is_readable = .FALSE.
    endif

  elseif (key=='q1') then
    if (associated(q1)) then
       q1 = real_rank2_array
    else
       is_readable = .FALSE.
    endif

  elseif (key=='u0') then
    u0 = real_rank2_array
  elseif (key=='v0') then
    v0 = real_rank2_array
  elseif (key=='vort0') then
    vort0 = real_rank2_array
  elseif (key=='psi0') then
    psi0 = real_rank2_array

  elseif (key=='rhsvort0') then
    rhsvort0 = real_rank3_array

  elseif (key=='rhsu0bar') then
    rhsu0bar = real_rank1_array

  elseif (key=='Ts') then
    Ts = real_rank2_array
  elseif (key=='STYPE') then
    STYPE = real_rank2_array
  elseif (key=='WD') then
    WD = real_rank2_array
  elseif (key=='WD0') then
    WD0 = real_rank1_array
  elseif (key=='arr1') then
    arr1 = real_rank2_array
  elseif (key=='arr2') then
    arr2 = real_rank2_array
  elseif (key=='arr3') then
    arr3 = real_rank2_array
  elseif (key=='arr4') then
    arr4 = real_rank2_array
  elseif (key=='arr5') then
    arr5 = real_rank2_array
  elseif (key=='arr6') then
    arr6 = real_rank2_array
  elseif (key=='arr7') then
    arr7 = real_rank2_array
  elseif (key=='arr8') then
    arr8 = real_rank2_array
  else
    write(*,*) 'Error-Bad call to SetbyPy setitem_real_array for ', key
    stop
  endif
  End Subroutine setitem_real_array


  Subroutine setitem_int(key, value)
  !--------------------------------------------------------------------
  ! Set the module variable key to value, for value of integer type.
  !
  ! Input:
  !   key:    String naming the variable in the module to be changed.
  !   value:  The value to replace the module variable with.  value is
  !           type integer scalar.
  !
  ! Output:
  !   Module variable named key is set to value and set module
  !   variable is_readable to .TRUE.  If variable cannot be set (such
  !   as if the pointer variable specified by key is not yet associa-
  !   ted), is_readable is set to .FALSE. and value is not set in the
  !   model.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key
     Integer, intent(in) :: value
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='it') then
    it=value
  elseif (key=='noout') then
    noout=value
  elseif (key=='nooutr') then
    nooutr=value
  elseif (key=='ntout') then
    ntout=value
  elseif (key=='ntouti') then
    ntouti=value
  elseif (key=='ntoutr') then
    ntoutr=value
  elseif (key=='mrestart') then
    mrestart=value
  elseif (key=='landon') then
    landon=value
  elseif (key=='nastep') then
    nastep=value
  elseif (key=='mt0') then
    mt0=value
  elseif (key=='lastday') then
    lastday=value
  elseif (key=='interval') then
    interval=value
  elseif (key=='dayofmodel') then
    dayofmodel=value
  elseif (key=='dateofmodel') then
    dateofmodel=value
  elseif (key=='yearofmodel') then
    yearofmodel=value
  elseif (key=='dayofmonth') then
    dayofmonth=value
  elseif (key=='dayofyear') then
    dayofyear=value
  elseif (key=='monthofyear') then
    monthofyear=value
  elseif (key=='year0') then
    year0=value
  elseif (key=='month0') then
    month0=value
  elseif (key=='day0') then
    day0=value
  else
      write(*,*) 'Error-Bad call to SetbyPy module setitem_int'
      write(*,*) key, '  ', value
      stop
  end if
  End Subroutine setitem_int


  Subroutine setitem_str(key, value)
  !--------------------------------------------------------------------
  ! Set the module variable key to value, for value of character type.
  !
  ! Input:
  !   key:    String naming the variable in the module to be changed.
  !   value:  The value to replace the module variable with.  value is
  !           type character scalar.
  !
  ! Output:
  !   Module variable named key is set to value and set module
  !   variable is_readable to .TRUE.  If variable cannot be set (such
  !   as if the pointer variable specified by key is not yet associa-
  !   ted), is_readable is set to .FALSE. and value is not set in the
  !   model.
  !
  ! Local Variables:
     Character(len=*), intent(in) :: key, value
  !--------------------------------------------------------------------
  is_readable = .TRUE.         !- Set is_readable to default

  if (key=='bnddir') then
    bnddir=value
  elseif (key=='outdir') then
    outdir=value
  elseif (key=='runname') then
    runname=value
  elseif (key=='SSTdir') then
    SSTdir=value
  elseif (key=='title') then
    title=value
  elseif (key=='SSTmode') then
    SSTmode=value
  else
      write(*,*) 'Error-Bad call to SetbyPy module setitem_str'
      write(*,*) key, '  ', value
      stop
  end if
  End Subroutine setitem_str
End Module SetbyPy
!
!======================================================================
