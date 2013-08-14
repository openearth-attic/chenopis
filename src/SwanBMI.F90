module swanbmi
  ! This module exposes a c interface and implements the bmi interface
  ! as documented at http://csdms.colorado.edu/wiki/BMI_Description

  use iso_c_binding

  implicit none

  ! Define some global constants

  ! We can't return an array of pointers with an undertimened length so
  ! we fix it to something, in this case 100.

  integer(c_int) :: MAXNAMES = 100
  integer(c_int) :: MAXDIMS = 6
  integer(c_int) :: MAXSTRINGLEN = 1024


  ! TODO get rid of this at a module level....
  CHARACTER PTYPE, PNAME *8, COMPUT *4, DTTIWR*18
  INTEGER, ALLOCATABLE :: CROSS(:)

  ! Extra local variables, these should not be here....
  integer :: inerr  ! number of the initialisation error
  integer :: ilen ! length of array, not sure which....

  real   , allocatable, target :: ac1(:,:,:), compda(:,:)
  integer, allocatable :: bgridp(:)
  real   , allocatable :: bspecs(:,:,:,:)
  real, allocatable    :: blknd(:), blkndc(:), ourqt(:)
  CHARACTER*20 NUMSTR, CHARS(1)
  INTEGER   ISTAT, IF1, IL1
  INTEGER   IOSTAT, IT0, IT, SAVITE
  INTEGER   ICOUNT
  INTEGER   IUNIT


  CHARACTER*80 MSGSTR

  logical stpnow
  external stpnow, numstr, dttiwr


contains
  subroutine windowssucks() bind(C, name="windowssucks")
    !DEC$ ATTRIBUTES DLLEXPORT :: main
    ! Somehow windows expects a main routine in the dll......
  end subroutine windowssucks

  subroutine run_model() bind(C, name="run_model")
    use swanapi, only : swan
    call swan()
  end subroutine run_model

  subroutine initialize(c_config_file) bind(C, name="initialize")
    !DEC$ ATTRIBUTES DLLEXPORT :: initialize

    ! The initialize() function accepts a string argument that
    ! gives the name (and path) of its "main input file", called
    ! a configuration file. This function should perform all tasks
    ! that are to take place before entering the model's time loop.
    use iso_c_binding, only: c_char
    use m_genarr
    use swcomm2
    use swcomm3
    use swcomm4
    use ocpcomm2
    use ocpcomm4
    use SwanGriddata
    use swanapi, only : swinit, swprep,swinco, errchk,swrbc


    ! The c name of the configuration file
    character(kind=c_char), intent(in) :: c_config_file(MAXSTRINGLEN)
    ! The fortran name of the configuration file
    character(len=strlen(c_config_file)) :: config_file



    ! Store the name
    config_file = char_array_to_string(c_config_file, strlen(c_config_file))

    ! Now we can initialize with the config_file

    ! move all of this to general init, it's now in main.... ugghh
    leverr = 0
    maxerr=1
    ITRACE=0
    INERR =0
    ISTAT =0

    call swinit(inerr)
    ! What is this.....

    IF (INERR.GT.0) RETURN
    IF (STPNOW()) RETURN

    comput = '    '

    call swread(comput)

    IF (STPNOW()) RETURN

    !       --- if last command was STOP then exit from repeat

    IF (COMPUT.EQ.'STOP') THEN
       IUNIT  = 0
       IOSTAT = 0
       FILENM = 'norm_end'
       CALL FOR (IUNIT, FILENM, 'UF', IOSTAT)
       WRITE (IUNIT, *) ' Normal end of run ', PROJNR
       return
    ENDIF

    if (numobs .gt. 0) then
       if (optg.ne.5) then
          ! structured grid
          ilen = 2*mcgrd
       else
          ! unstructured grid
          ilen = nfaces
       endif
       if (.not. allocated(cross)) allocate(cross(ilen))
    else
       if (.not.allocated(cross)) allocate(cross(0))
    endif
    if (.not.allocated(bspecs)) allocate(bspecs(mdc,msc,nbspec,2))
    if (.not.allocated(bgridp)) allocate(bgridp(6*nbgrpt))

    call swprep ( bspecs, bgridp, cross , xcgrid, ycgrid, kgrpnt, kgrbnd, spcdir, spcsig)

    ! Just terminate if we get an error....
    if (inerr > 0) then
       write(*,*) 'Initialize failed with inerr:', inerr
    end if


    IF (OPTG.EQ.5) CALL SwanPrepComp ( CROSS )
    IF (STPNOW()) RETURN
    !TIMG        CALL SWTSTO(4)

    !       --- check all possible flags and if necessary change
    !           if option is not correct

    CALL ERRCHK
    IF (STPNOW()) RETURN

    !       --- initialisation of necessary grids for depth,
    !           current, wind and friction

    IF (ALOCMP.AND.ALLOCATED(COMPDA)) DEALLOCATE(COMPDA)
    IF (.NOT.ALLOCATED(COMPDA)) THEN
       ALLOCATE(COMPDA(MCGRD,MCMVAR),STAT=ISTAT)
       ALOCMP = .FALSE.
    END IF
    IF ( ISTAT.NE.0 ) THEN
       CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
       CALL TXPBLA(CHARS(1),IF1,IL1)
       MSGSTR = &
            &         'Allocation problem: array COMPDA and return code is '//    &
            &         CHARS(1)(IF1:IL1)
       CALL MSGERR ( 4, MSGSTR )
       RETURN
    END IF

    CALL SWRBC(COMPDA)
    !       --- allocate AC1 in case of non-stationary situation or in case
    !           of using the S&L scheme

    IF ( NSTATM.EQ.1 .AND. MXITNS.GT.1 .OR. PROPSC.EQ.3 ) THEN
       IF (.NOT.ALLOCATED(AC1)) THEN
          ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
       ELSE IF (SIZE(AC1).EQ.0) THEN
          DEALLOCATE(AC1)
          ALLOCATE(AC1(MDC,MSC,MCGRD),STAT=ISTAT)
       END IF
       IF ( ISTAT.NE.0 ) THEN
          CHARS(1) = NUMSTR(ISTAT,RNAN,'(I6)')
          CALL TXPBLA(CHARS(1),IF1,IL1)
          MSGSTR =  &
               & 'Allocation problem: array AC1 and return code is '//   &
               & CHARS(1)(IF1:IL1)
          CALL MSGERR ( 4, MSGSTR )
          RETURN
       END IF
       AC1 = 0.
    ELSE
       IF(.NOT.ALLOCATED(AC1)) ALLOCATE(AC1(0,0,0))
    ENDIF

    IF (LEVERR.GT.MAXERR) THEN

       WRITE (PRINTF, 6010) LEVERR
       IF (LEVERR.LT.4) WRITE (PRINTF, 6011)
6010   FORMAT(' ** No start of computation because of error level:' &
            & ,I3)
6011   FORMAT(' ** To ignore this error, change [maxerr] with the', &
            & ' SET command')
    end IF
    !
    IF (ITEST.GE.40) THEN
       IF (NSTATC.EQ.1) THEN
          WRITE (PRINTF, '(" Type of computation: dynamic")')
       ELSE
          IF (ONED) THEN
             WRITE (PRINTF, '(" Type of computation: static 1-D")')
          ELSE
             WRITE (PRINTF, '(" Type of computation: static 2-D")')
          ENDIF
       ENDIF
    ENDIF
    !
    IF (NSTATC.EQ.1) THEN
       IT0 = 0
       IF (ICOND.EQ.1) THEN
          !
          !             --- compute default initial conditions
          !
          !TIMG              CALL SWTSTA(6)
          CALL SWINCO ( AC2   , COMPDA, XCGRID, YCGRID, &
               &                      KGRPNT, SPCDIR, SPCSIG, XYTST )
          !TIMG              CALL SWTSTO(6)
          !
          !             --- reset ICOND to prevent second computation of
          !                 initial condition
          ICOND = 0

       ENDIF
    ELSE
       IT0 = 1
    ENDIF
    IT = IT0

    !         --- synchronize nodes

    CALL SWSYNC
    IF (STPNOW()) RETURN
  end subroutine initialize

  subroutine update(delta_t) bind(C, name="update")
    !DEC$ ATTRIBUTES DLLEXPORT :: update

    use iso_c_binding, only: c_double

    use swanapi
    USE TIMECOMM
    USE OCPCOMM2
    USE OCPCOMM4
    USE SWCOMM1
    USE SWCOMM2
    USE SWCOMM3
    USE SWCOMM4
    USE OUTP_DATA
    USE M_GENARR
    USE M_PARALL
    USE SwanGriddata

    real(c_double), intent(in) :: delta_t

    ! Only use dt for non stationary.
    if (NSTATC .eq. 1) then
       DT = delta_t
    end if



    write(*,*) LEVERR, MAXERR
    IF (LEVERR.GT.MAXERR) THEN
       WRITE (PRINTF, 6030) LEVERR
       IF (LEVERR.LT.4) WRITE (PRINTF, 6011)
6030   FORMAT(' ** No continuation of computation because ', &
            &               'of error level:',I3)
6010   FORMAT(' ** No start of computation because of error level:' &
            & ,I3)
6011   FORMAT(' ** To ignore this error, change [maxerr] with the', &
            & ' SET command')

       return
    ENDIF

    !           --- synchronize nodes

    CALL SWSYNC
    IF (STPNOW()) RETURN
    !
    !           --- update boundary conditions and input fields
    !
    CALL SNEXTI ( BSPECS, BGRIDP, COMPDA, AC1   , AC2   , &
         &                    SPCSIG, SPCDIR, XCGRID, YCGRID, KGRPNT, &
         &                    XYTST , DEPTH , WLEVL , FRIC  , UXB   , &
         &                    UYB   , NPLAF , WXI   , WYI   )
    IF (STPNOW()) RETURN
    !
    !           --- synchronize nodes

    CALL SWSYNC
    IF (STPNOW()) RETURN

    IF (COMPUT.NE.'NOCO' .AND. IT.GT.0) THEN

       SAVITE = ITEST
       IF (ICOTES .GT. ITEST) ITEST = ICOTES
       !
       !             --- compute action density for current time step
       !
       !TIMG              CALL SWTSTA(8)
       IF (OPTG.NE.5) THEN
          !                structured grid
          CALL SWCOMP( AC1   , AC2   , COMPDA, SPCDIR, SPCSIG, &
               &                        XYTST , IT    , KGRPNT, XCGRID, YCGRID, &
               &                        CROSS )
       ELSE
          !                unstructured grid
          CALL SwanCompUnstruc ( AC2   , AC1   , COMPDA, &
               &                                  SPCSIG, SPCDIR, XYTST , &
               &                                  CROSS , IT    )
       ENDIF
       IF (STPNOW()) RETURN
       !
       !             --- set ICOND=4 for stationary computation, for next
       !                 (stationary) COMPUTE command
       ICOND = 4
       !
       !             --- check whether computed significant wave height at
       !                 boundary differs from prescribed value given in
       !                 boundary command values of incident Hs
       !
       IF ( BNDCHK ) THEN
          CALL HSOBND ( AC2, SPCSIG, COMPDA(1,JHSIBC), KGRPNT )
       ENDIF
       !
       ITEST = SAVITE

    ENDIF
    !
    IF ( IT.EQ.IT0 .AND. .NOT.ALLOCATED(OURQT) ) THEN
       ALLOCATE (OURQT(MAX_OUTP_REQ))
       OURQT = -9999.
    ENDIF
    !
    SAVITE = ITEST
    IF (IOUTES .GT. ITEST) ITEST = IOUTES
    !ADC
    !ADC!           --- perform output of SWAN quantities
    !ADC            CALL SwanOutput ( ITIME, IT )
    !ADC!           --- write the SWAN hot-start file, if necessary
    !ADC            IF ( WriteSwanHotStart ) THEN
    !ADC               CALL BACKUP ( AC2,SPCSIG,SPCDIR,KGRPNT,XCGRID,YCGRID )
    !ADC               WriteSwanHotStart = .FALSE.
    !ADC            ENDIF

    !           --- synchronize nodes

    CALL SWSYNC
    IF (STPNOW()) RETURN

    !           --- carry out the output requests

    !TIMG            CALL SWTSTA(9)
    CALL SWOUTP ( AC2   , SPCSIG, SPCDIR, COMPDA, XYTST , &
         & KGRPNT, XCGRID, YCGRID, OURQT )
    !TIMG            CALL SWTSTO(9)
    IF (STPNOW()) RETURN
    !
    IF (ERRPTS.GT.0) REWIND(ERRPTS)
    ITEST = SAVITE

    !           --- update time

    IF (NSTATC.EQ.1) THEN
       IF (IT.LT.MTC) THEN
          TIMCO = TIMCO + DT
          CHTIME = DTTIWR(ITMOPT, TIMCO)
          WRITE (PRINTF, 222) CHTIME, TIMCO
       ENDIF
222    FORMAT(' Time of computation ->  ',A,' in sec:', F12.0)
       IT = IT + 1
    ENDIF


  end subroutine update

  subroutine update_until(t) bind(C, name="update_until")
    use iso_c_binding, only: c_double
    real(c_double),intent(in) :: t
  end subroutine update_until

  subroutine finalize() bind(C, name="finalize")
    !DEC$ ATTRIBUTES DLLEXPORT :: finalize

  end subroutine finalize

  subroutine get_n_attributes(n) bind(C, name="get_n_attributes")
    use swcomm3
    use iso_c_binding, only: c_char, c_ptr, c_loc, c_int, c_null_char
    integer(c_int), intent(out) :: n

    n = MNUMS
  end subroutine get_n_attributes

  subroutine get_attribute_name(i, c_att_name) bind(C, name="get_attribute_name")
    use iso_c_binding, only: c_char, c_ptr, c_loc, c_int, c_null_char
    integer(c_int), intent(in) :: i
    character(kind=c_char), intent(out) :: c_att_name(MAXSTRINGLEN)
    character(len=MAXSTRINGLEN) :: name

    select case(i)
       ! From SWANCOM1
    case(1)
       name = 'DREL'
    case(2)
       name = 'DHABS'
    case(3)
       name = 'DTABS'
    case(4)
       name = 'NPNTS'
    case(5)
       name = ''
    case(6)
       name = 'CDD'
    case(7)
       name = 'CSS'
    case(8)
       name = 'NUMFRE'
    case(9)
       name = 'DIFFC'
    case(12)
       name = 'EPS2'
    case(13)
       name = 'OUTP'
    case(14)
       name = 'NITER'
    case(15)
       name = 'DHOVAL'
    case(16)
       name = 'DTOVAL'
    case(17)
       name = 'CDLIM'
    case(18)
       name = 'FROUDMAX'
    case(19)
       name = 'CFL'
    case(20)
       name = 'GRWMX'
    case(21)
       name = 'STOPC'
    case(30)
       name = 'ALFA'
    case default
       name = ''
    end select

    c_att_name = string_to_char_array(trim(name), len(trim(name)))
  end subroutine get_attribute_name

  subroutine get_attribute_type(c_att_name, c_type) bind(C, name="get_attribute_type")
    use iso_c_binding, only: c_char, c_ptr, c_null_char
    character(kind=c_char), intent(in) :: c_att_name(*)
    character(kind=c_char), intent(out) :: c_type(MAXSTRINGLEN)
    ! Use one of the following types
    ! BMI datatype        C datatype        NumPy datatype
    ! BMI_STRING          char*             S<
    ! BMI_INT             int               int16
    ! BMI_DOUBLE          double            float64

    ! The fortran name of the attribute name
    character(len=strlen(c_att_name)) :: att_name
    character(len=MAXSTRINGLEN) :: type
    ! Store the name
    att_name = char_array_to_string(c_att_name, strlen(c_att_name))

    ! Lookup the type of the variable

    type = 'BMI_INT'
    c_type = string_to_char_array(trim(type), len(trim(type)))
  end subroutine get_attribute_type

  subroutine get_double_attribute(c_att_name, value)  bind(C, name="get_double_attribute")
    use iso_c_binding, only: c_double, c_char
    character(kind=c_char), intent(in) :: c_att_name(*)
    real(c_double), intent(out)         :: value

    ! The fortran name of the attribute name
    character(len=strlen(c_att_name)) :: att_name
    ! Store the name
    att_name = char_array_to_string(c_att_name, strlen(c_att_name))

    ! Look up the value of att_name
    value = -1.0d0
  end subroutine get_double_attribute

  subroutine get_int_attribute(c_att_name, value)  bind(C, name="get_int_attribute")
    use swcomm3

    use iso_c_binding, only: c_double, c_char
    character(kind=c_char), intent(in) :: c_att_name(*)
    integer(c_int), intent(out)         :: value

    ! The fortran name of the attribute name
    character(len=strlen(c_att_name)) :: att_name
    ! Store the name
    att_name = char_array_to_string(c_att_name, strlen(c_att_name))

    ! Look up the value of att_name

    select case(att_name)
    case default
       value = -1
    case ("MDC")
       value = MDC
    case ("MSC")
       value = MSC
    case ('MCGRD')
       value = MCGRD
    end select
  end subroutine get_int_attribute

  subroutine get_string_attribute(c_att_name, c_value)  bind(C, name="get_string_attribute")
    use iso_c_binding, only: c_char
    character(kind=c_char), intent(in) :: c_att_name(*)
    character(kind=c_char), intent(out) :: c_value(MAXSTRINGLEN)

    ! The fortran name of the attribute name
    character(len=strlen(c_att_name)) :: att_name
    character(len=MAXSTRINGLEN) :: value
    ! Store the name
    att_name = char_array_to_string(c_att_name, strlen(c_att_name))

    ! Lookup the type of the variable

    value = 'some value'
    c_value = string_to_char_array(trim(value), len(trim(value)))
  end subroutine get_string_attribute


  subroutine get_n_vars(n) bind(C, name="get_n_vars")
    use iso_c_binding, only: c_char, c_ptr, c_loc, c_int, c_null_char


    integer(c_int), intent(out) :: n

    include "bmi_nvar.inc"
  end subroutine get_n_vars

  subroutine get_var_name(i, c_var_name) bind(C, name="get_var_name")
    use iso_c_binding, only: c_char, c_ptr, c_loc, c_int, c_null_char

    integer(c_int), intent(in) :: i
    character(kind=c_char), intent(out) :: c_var_name(MAXSTRINGLEN)
    character(len=MAXSTRINGLEN) :: name

    include "bmi_var_name.inc"

    c_var_name = string_to_char_array(trim(name), len(trim(name)))

  end subroutine get_var_name


  subroutine get_var_shape(c_var_name, shape) bind(C, name="get_var_shape")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_shape

    include "bmi_module.inc"
    use iso_c_binding, only: c_int, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(inout) :: shape(MAXDIMS)

    character(len=strlen(c_var_name)) :: var_name

    var_name = char_array_to_string(c_var_name, strlen(c_var_name))
    shape = (/0, 0, 0, 0, 0, 0/)

    include "bmi_shape.inc"
  end subroutine get_var_shape

  subroutine get_input_var_names(names) bind(C, name="get_input_var_names")
    use iso_c_binding, only: c_char, c_ptr
    character(kind=c_char), dimension(MAXNAMES), intent(out) :: names(*)
    !type(c_ptr), dimension(:) :: names
  end subroutine get_input_var_names


  subroutine get_var_names(names) bind(C, name="get_var_names")
    use iso_c_binding, only: c_char, c_ptr
    character(kind=c_char), dimension(MAXNAMES), intent(out) :: names(*)

    ! I can't get this to work.....

    ! http://stackoverflow.com/questions/9686532/arrays-of-strings-in-fortran-c-bridges-using-iso-c-binding
    ! The way we do it is to use a C_PTR array to point to strings. For example:

    ! CHARACTER(LEN=100), DIMENSION(numStrings), TARGET :: stringArray
    ! TYPE(C_PTR), DIMENSION(numStrings) :: stringPtrs
    ! then we set our strings in stringArray, remembering to null-terminate them such as:

    ! DO ns = 1, numStrings
    !    stringArray(ns) = "My String"//C_NULL_CHAR
    !    stringPtrs(ns) = C_LOC(stringArray(ns))
    ! END DO
    ! and pass stringPtrs to the C function.

    ! The C function has the interface:

    ! void stringFunc(int *numStrings, char **stringArray) {
    !     int i;
    !     for(i=0;i<*numStrings;++i) {
    !        printf("%s\n",stringArray[i]);
    !     }
    !  }
  end subroutine get_var_names

  subroutine get_output_var_names(names) bind(C, name="get_output_var_names")
    use iso_c_binding, only: c_char
    character(kind=c_char), dimension(MAXNAMES), intent(out) :: names(*)
  end subroutine get_output_var_names

  subroutine get_var_type(c_var_name, c_type_name)  bind(C, name="get_var_type")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_type

    use swcomm1

    character(kind=c_char), intent(in) :: c_var_name(*)
    character(kind=c_char), intent(out) :: c_type_name(MAXSTRINGLEN)
    ! Use one of the following types
    ! BMI datatype        C datatype        NumPy datatype
    ! BMI_STRING          char*             S<
    ! BMI_INT             int               int16
    ! BMI_LONG            long int          int32
    ! BMI_FLOAT           float             float32
    ! BMI_DOUBLE          double            float64


    integer(c_int) :: i
    integer :: typeid
    character(len=MAXSTRINGLEN) :: type_name
    character(len=MAXSTRINGLEN) :: var_name

    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    include 'bmi_type.inc'

    c_type_name = string_to_char_array(trim(type_name), len(trim(type_name)))

  end subroutine get_var_type

  subroutine get_var_role(long_var_name, role) bind(C, name="get_var_role")
    character(kind=c_char), intent(in) :: long_var_name(*)
    character(kind=c_char), intent(out) :: role(*)
    ! Roles:
    ! BMI_INPUT
    ! BMI_OUTPUT
    ! BMI_INPUTOUTPUT
  end subroutine get_var_role

  subroutine get_var_unit( c_var_name, c_unit )  bind(C, name="get_var_unit")
    use swcomm1
    character(kind=c_char), intent(in) :: c_var_name(*)
    character(kind=c_char), intent(out) :: c_unit(MAXSTRINGLEN)

    integer(c_int) :: i
    integer :: typeid
    character(len=MAXSTRINGLEN) :: unit

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    include "bmi_unit.inc"

    c_unit = string_to_char_array(trim(unit), len(trim(unit)))

  end subroutine get_var_unit


  subroutine get_var_rank(c_var_name, rank) bind(C, name="get_var_rank")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var_rank

    use iso_c_binding, only: c_int, c_char
    character(kind=c_char), intent(in) :: c_var_name(*)
    integer(c_int), intent(out) :: rank

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    include "bmi_rank.inc"

  end subroutine get_var_rank

  subroutine get_var(c_var_name, x) bind(C, name="get_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_var
    ! Return a pointer to the variable
    use iso_c_binding, only: c_double, c_char, c_loc


    include "bmi_module.inc"


    real(c_double), save, target :: x_0d_double_ptr
    real(c_double), allocatable, save, target :: x_1d_double_ptr(:)
    real(c_double), allocatable, save, target :: x_2d_double_ptr(:,:)
    real(c_double), allocatable, save, target :: x_3d_double_ptr(:,:,:)
    integer(c_int), save, target :: x_0d_int_ptr
    integer(c_int), allocatable, save, target :: x_1d_int_ptr(:)
    integer(c_int), allocatable, save, target :: x_2d_int_ptr(:,:)
    integer(c_int), allocatable, save, target :: x_3d_int_ptr(:,:,:)
    real(c_float), save, target  :: x_0d_float_ptr
    real(c_float), allocatable, save, target  :: x_1d_float_ptr(:)
    real(c_float), allocatable, save, target  :: x_2d_float_ptr(:,:)
    real(c_float), allocatable, save, target  :: x_3d_float_ptr(:,:,:)

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    include "bmi_get_nd.inc"

  end subroutine get_var


  subroutine get_1d_int(c_var_name, x) bind(C, name="get_1d_int")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_1d_int

    use iso_c_binding, only: c_int, c_char, c_loc
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select

  end subroutine get_1d_int

  subroutine get_1d_double(c_var_name, x) bind(C, name="get_1d_double")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_1d_double

    use iso_c_binding, only: c_double, c_char, c_loc
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select

  end subroutine get_1d_double


  subroutine get_2d_double(c_var_name, x) bind(C, name="get_2d_double")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_2d_double

    use iso_c_binding, only: c_double, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_2d_double



  subroutine get_3d_double(c_var_name, x) bind(C, name="get_3d_double")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_3d_double

    use iso_c_binding, only: c_double, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_3d_double


  subroutine get_0d_float(c_var_name, x) bind(C, name="get_0d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_0d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    include "bmi_module.inc"

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x


    real(c_float), save, target :: x0d_value

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))



    select case(var_name)
    case ("DDIR")
       x0d_value = DDIR
       x = c_loc(x0d_value)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_0d_float


  subroutine get_1d_float(c_var_name, x) bind(C, name="get_1d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_1d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_1d_float


  subroutine get_2d_float(c_var_name, x) bind(C, name="get_2d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_2d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_2d_float

  subroutine get_3d_float(c_var_name, x) bind(C, name="get_3d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_3d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call get_var(c_var_name, x)
    end select
  end subroutine get_3d_float


  subroutine set_1d_float(c_var_name, x) bind(C, name="set_1d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_1d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call set_var(c_var_name, x)
    end select
  end subroutine set_1d_float

  subroutine set_2d_float(c_var_name, x) bind(C, name="set_2d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_2d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call set_var(c_var_name, x)
    end select
  end subroutine set_2d_float

  subroutine set_3d_float(c_var_name, x) bind(C, name="set_3d_float")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_3d_float

    use iso_c_binding, only: c_float, c_char, c_loc

    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(inout) :: x

    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    select case(var_name)
    case default
       call set_var(c_var_name, x)
    end select
  end subroutine set_3d_float




  ! subroutine set_1d_double_at_index(c_var_name, index, value)  bind(C, name="set_1d_double_at_index")
  !   !DEC$ ATTRIBUTES DLLEXPORT :: set_1d_double_at_index

  !   use m_flow
  !   use network_data
  !   use m_flowparameters
  !   use iso_c_binding, only: c_double, c_char, c_int

  !   character(kind=c_char), intent(in) :: c_var_name(*)
  !   integer(c_int), intent(in)         :: index
  !   real(c_double), intent(in)         :: value

  !   character(len=strlen(c_var_name)) :: var_name

  !   var_name = char_array_to_string(c_var_name, strlen(c_var_name))

  !   select case(var_name)
  !   case("s1")
  !      s1(index + 1) = value
  !   case("ucx")
  !      ucx(index + 1) = value
  !   case("ucy")
  !      ucy(index + 1) = value
  !   case("unorm")
  !      u1(index + 1) = value
  !   case("zk")
  !      zkdropstep = value - zk(index + 1)
  !      call dropland(xz(index + 1), yz(index + 1), 1)
  !   end select

  ! end subroutine set_1d_double_at_index


  subroutine set_var(c_var_name, xptr) bind(C, name="set_var")
    !DEC$ ATTRIBUTES DLLEXPORT :: set_var
    ! Return a pointer to the variable
    use iso_c_binding, only: c_double, c_char, c_loc, c_f_pointer

    include "bmi_module.inc"
    character(kind=c_char), intent(in) :: c_var_name(*)
    type(c_ptr), intent(in) :: xptr

    real(c_double), pointer :: x_0d_double_ptr
    real(c_double), pointer :: x_1d_double_ptr(:)
    real(c_double), pointer :: x_2d_double_ptr(:,:)
    real(c_double), pointer :: x_3d_double_ptr(:,:,:)
    integer(c_int), pointer :: x_0d_int_ptr
    integer(c_int), pointer :: x_1d_int_ptr(:)
    integer(c_int), pointer :: x_2d_int_ptr(:,:)
    integer(c_int), pointer :: x_3d_int_ptr(:,:,:)
    real(c_float), pointer  :: x_0d_float_ptr
    real(c_float), pointer  :: x_1d_float_ptr(:)
    real(c_float), pointer  :: x_2d_float_ptr(:,:)
    real(c_float), pointer  :: x_3d_float_ptr(:,:,:)
    ! The fortran name of the attribute name
    character(len=strlen(c_var_name)) :: var_name
    ! Store the name
    var_name = char_array_to_string(c_var_name, strlen(c_var_name))

    include "bmi_set_nd.inc"

  end subroutine set_var

  subroutine get_long_var_name(long_var_name, name)  bind(C, name="get_long_var_name")
    character(kind=c_char), intent(in) :: long_var_name(*)
    character(kind=c_char), intent(out) :: name(*)
  end subroutine get_long_var_name

  subroutine get_time_step(step) bind(C, name="get_time_step")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_time_step

    use timecomm
    use iso_c_binding, only: c_double
    real(c_double), intent(inout) :: step
    step = DT
  end subroutine get_time_step

  subroutine get_time_units(c_unit)  bind(C, name="get_time_units")
    ! returns unit string for model time, e.g. ‘days since 1970-01-01'
    character(kind=c_char), intent(out) :: c_unit(MAXSTRINGLEN)
    character(len=MAXSTRINGLEN) :: unit

    unit = 'days since 1970-01-01 00:00'
    c_unit = string_to_char_array(trim(unit), len(trim(unit)))

  end subroutine get_time_units

  subroutine get_start_time(value) bind(C, name="get_start_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_start_time

    use timecomm
    use iso_c_binding, only: c_double
    real(c_double), intent(out) :: value
    value = 0.0d0
  end subroutine get_start_time
  subroutine get_current_time(value) bind(C, name="get_current_time")
    use timecomm
    real(c_double), intent(out) :: value


  end subroutine get_current_time

  subroutine get_end_time(value) bind(C, name="get_end_time")
    !DEC$ ATTRIBUTES DLLEXPORT :: get_end_time

    real(c_double), intent(out) :: value
  end subroutine get_end_time


  ! Utility functions, move these to interop module
  ! Make functions pure so they can be used as input arguments.
  integer(c_int) pure function strlen(char_array)
    character(c_char), intent(in) :: char_array(MAXSTRINGLEN)
    integer :: inull, i
    strlen = 0
    do i = 1, size(char_array)
       if (char_array(i) .eq. C_NULL_CHAR) then
          strlen = i-1
          exit
       end if
    end do
  end function strlen

  pure function char_array_to_string(char_array, length)
    integer(c_int), intent(in) :: length
    character(c_char),intent(in) :: char_array(length)
    character(len=length) :: char_array_to_string
    integer :: i
    do i = 1, length
       char_array_to_string(i:i) = char_array(i)
    enddo
  end function char_array_to_string

  pure function string_to_char_array(string, length)
    character(len=length), intent(in) :: string
    integer(c_int),intent(in) :: length
    character(kind=c_char,len=1) :: string_to_char_array(length+1)
    integer :: i
    do i = 1, length
       string_to_char_array(i) = string(i:i)
    enddo
    string_to_char_array(length+1) = C_NULL_CHAR
  end function string_to_char_array


!   pure function significant_wave_height()


!     EFTAIL = 1. / (PWTAIL(1) - 1.)
!     !integration over [0,inf]
!     ETOT = 0.
!     ! trapezoidal rule is applied
!     DO ID=1, MDC
!        DO IS=2,MSC
!           DS=SPCSIG(IS)-SPCSIG(IS-1)
!           EAD = 0.5*(SPCSIG(IS)*ACLOC(ID,IS)+ &
!                SPCSIG(IS-1)*ACLOC(ID,IS-1))*DS*DDIR
!           ETOT = ETOT + EAD
!        ENDDO
!        IF (MSC .GT. 3) THEN                                       10.20
!           !                contribution of tail to total energy density
!           EHFR = ACLOC(ID,MSC) * SPCSIG(MSC)                       30.72
!           ETOT = ETOT + DDIR * EHFR * SPCSIG(MSC) * EFTAIL         30.72
!        ENDIF
!     ENDDO
! !            integration over [fmin,fmax]                                 40.87
!              FMIN = PI2*OUTPAR(21)                                        40.87
!              FMAX = PI2*OUTPAR(36)                                        40.87
!              ECS  = 1.                                                    40.87
!              ETOT = SwanIntgratSpc(0. , FMIN, FMAX, SPCSIG, SPCDIR(1,1),  40.87
!      &                             WK , ECS , 0.  , 0.    , ACLOC      ,  40.87
!      &                             1  )                                   40.87
!           ENDIF                                                           40.87
!           IF (ETOT .GE. 0.) THEN                                          30.00
!             VOQ(IP,VOQR(IVTYPE)) = 4.*SQRT(ETOT)
!           ELSE
!             VOQ(IP,VOQR(IVTYPE)) = 0.                                     40.86

!   end function a

end module swanbmi
