#define ESMF_CHECK if (ESMF_LogFoundError(rcToCheck=rc, line=__LINE__, file=__FILE__)) call ESMF_Finalize(endflag=ESMF_END_ABORT)
module swan_esmf
  use esmf

  use iso_c_binding
  use swanbmi !, only : get_var_rank, get_var_shape, strlen, char_array_to_string, string_to_char_array, get_n_vars, get_var_name, initialize, update, get_var
  !use swfld_mod
  ! Import swan api
  implicit none

  integer, parameter ::  MAXSTRLEN =1024

  public swan_register

contains
  subroutine swan_register(gridcomp, rc)
    type(ESMF_GridComp) :: gridcomp
    integer, intent(out)               :: rc

    call ESMF_GridCompSetEntryPoint(gridcomp, ESMF_METHOD_INITIALIZE, userRoutine=swan_init, rc=rc)
    ESMF_CHECK
    call ESMF_GridCompSetEntryPoint(gridcomp, ESMF_METHOD_RUN, userRoutine=swan_run, rc=rc)
    ESMF_CHECK
    call ESMF_GridCompSetEntryPoint(gridcomp, ESMF_METHOD_FINALIZE, userRoutine=swan_final, rc=rc)
    ESMF_CHECK
  end subroutine swan_register


  subroutine swan_init(gridcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp)  :: gridcomp ! TODO bug in esmf? Can't sepcify intent here....
    type(ESMF_State)    :: importState
    type(ESMF_State)    :: exportState
    type(ESMF_Clock)    :: clock
    integer, intent(out)               :: rc

    type(ESMF_Grid)  :: grid
    type(ESMF_Mesh)  :: mesh
    type(ESMF_FieldBundle)  :: fieldbundle


    character(len=MAXSTRLEN) :: filename
    ! Initialize
    character(kind=c_char) :: c_filename(MAXSTRLEN)

    filename = 'INPUT'
    c_filename = string_to_char_array(trim(filename), len(trim(filename)))
    call initialize(c_filename)

    write(*,*) 'Making ESMF grid'
    ! Create the grid
    call make_swan_grid(grid=grid, rc=rc)
    !call make_swan_mesh(mesh, rc=rc)

    write(*,*) 'Making Fields'
    ! Use the same fieldbunde for import and export
    call make_swan_fieldbundle(grid=grid, fieldbundle=fieldbundle, rc=rc)
    importState = ESMF_StateCreate(name="swan import",stateintent=ESMF_STATEINTENT_IMPORT, rc=rc)
    call ESMF_StateAdd(importState, fieldbundleList=(/fieldbundle/), rc=rc)
    exportState = ESMF_StateCreate(name="swan export",stateintent=ESMF_STATEINTENT_EXPORT, rc=rc)
    call ESMF_StateAdd(exportState, fieldbundleList=(/fieldbundle/), rc=rc)
    ESMF_CHECK

    write(*,*) 'ESMF initialize complete'
  end subroutine swan_init

  subroutine swan_run(gridcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: gridcomp
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    integer, intent(out)               :: rc

    type(ESMF_TimeInterval) :: timeinterval
    real(ESMF_KIND_R8) :: esmfstep, esmft, esmfnext


    call ESMF_StateWrite(importState, 'dflowfmimport.nc', rc)


    ! Synchronize the clocks

    ! ESMF clock
    !       ESMFT -----> ESMFNEXT
    !         .             .
    !          <-ESMFSTEP-->

    ! DFLOW clock
    !     time_user --> tstop_user

    ! Get the simulation time (since tref) from esmf
    call ESMF_ClockGet(clock, currSimTime=timeinterval, rc=rc)
    call ESMF_TimeIntervalGet(timeinterval, s_r8=esmft, rc=rc)

    ! Get the timestep from esmf, the coupled timestep
    call ESMF_ClockGet(clock, timestep=timeinterval, rc=rc)
    call ESMF_TimeIntervalGet(timeinterval, s_r8=esmfstep, rc=rc)


    ! Run until we reach esmfnext
    esmfnext = (esmft+esmfstep)

    ! Assume swan runs all the way up to esmfnext....
    call update(esmfstep)

    ! State export:
    call ESMF_StateWrite(exportState, 'swanexport.nc', rc)
  end subroutine swan_run

  subroutine swan_final(gridcomp, importState, exportState, clock, rc)
    type(ESMF_GridComp) :: gridcomp
    type(ESMF_State)    :: importState, exportState
    type(ESMF_Clock)    :: clock
    integer, intent(out)               :: rc

  end subroutine swan_final


  subroutine make_swan_grid(grid, rc)
    use swcomm3
    use m_genarr
    ! use swfld_mod               ! Sea land mask

    type(ESMF_Grid), intent(inout) :: grid
    integer, intent(out)         :: rc

    real(ESMF_KIND_R8), pointer :: double2dptr(:,:)
    integer(ESMF_KIND_I4), pointer :: int2dptr(:,:)


    ! MXC    [     0]  grid points in x-direction of computational grid
    ! MYC    [     0]  grid points in y-direction of computational grid
    ! Cartesian grid
    grid = ESMF_GridCreateNoPeriDim(&
         regdecomp=(/1,1/), &
         maxIndex=(/mxc+1,myc+1/), &
         coordSys=ESMF_COORDSYS_CART, & ! Cartesian coordinates (x,y not lat,lon)
         coordDep1=(/1,2/), coordDep2=(/1,2/), & ! Use coordDep1=(/1,2/),coordDep2=(/1,2/) for curvilinear grids.
         indexflag=ESMF_INDEX_GLOBAL, & ! Needed for using XGrid
         name="swan:grid", rc=rc)
    ESMF_CHECK

    ! Add coordiantes

    call ESMF_GridAddCoord(grid, &
         staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    ESMF_CHECK


    ! Set the grid coordinates...
    ! Assuming grid is center
    call ESMF_GridGetCoord(grid, coorddim=1, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=double2dptr, rc=rc)
    ESMF_CHECK

    call ESMF_GridGetCoord(grid, coorddim=2, &
         staggerloc=ESMF_STAGGERLOC_CENTER, farrayPtr=double2dptr, rc=rc)
    ESMF_CHECK
    double2dptr(:,:) = ycgrid

    ! ! call ESMF_GridAddItem(grid, item=ESMF_GRIDITEM_MASK, &
    ! !      staggerloc=ESMF_STAGGERLOC_CENTER, rc=rc)
    ! ! ESMF_CHECK
    ! call ESMF_GridGetItem(grid,  &
    !      item=ESMF_GRIDITEM_MASK, staggerloc=ESMF_STAGGERLOC_CENTER, &
    !      fptr=int2dptr, rc=rc)
    ! ESMF_CHECK
    ! int2dptr(:,:) = mask



  end subroutine make_swan_grid


  subroutine make_swan_mesh(mesh, rc)
    use SwanGridData, only    : nverts
    use SwanGridobjects, only : cellid, &
         celltype,                      &
         cellv1,cellv2,cellv3,          &
         faceid,                        &
         facetype,                      &
         facev1,facev2,                 &
         fmarker,                       &
         gridobject,                    &
         verttype,                      &
         vertx,verty

    type(ESMF_Mesh), intent(out) :: mesh
    integer, intent(out)         :: rc



    ! Grid administration

    ! Use the same variable names as in the ESMF docs
    ! Dimensions and counting

    ! Dimension of the topology of the Mesh. (E.g. a mesh constructed
    ! of squares would have a parametric dimension of 2, whereas a
    ! Mesh constructed of cubes would have one of 3.)
    integer                      :: parametricDim = 2
    ! The number of coordinate dimensions needed to describe the
    ! locations of the nodes making up the Mesh. For a manifold, the
    ! spatial dimesion can be larger than the parametric dim (e.g. the
    ! 2D surface of a sphere in 3D space), but it can't be smaller.
    integer                      :: spatialDim = 2

    integer                      :: numNodes
    integer                      :: numQuadElems
    integer                      :: numTriElems
    integer                      :: numTotElems

    ! Variables An array containing the physical coordinates of the
    ! nodes to be created on this PET. This input consists of a 1D
    ! array the size of the number of nodes on this PET times the
    ! Mesh's spatial dimension (spatialDim). The coordinates in this
    ! array are ordered so that the coordinates for a node lie in
    ! sequence in memory. (e.g. for a Mesh with spatial dimension 2,
    ! the coordinates for node 1 are in nodeCoords(0) and
    ! nodeCoords(1), the coordinates for node 2 are in nodeCoords(2)
    ! and nodeCoords(3), etc.).
    real(ESMF_KIND_R8), allocatable :: nodeCoords(:)

    ! An array containing the global ids of the nodes to be created on
    ! this PET. This input consists of a 1D array the size of the
    ! number of nodes on this PET.
    integer, allocatable         :: nodeIds(:)
    integer, allocatable         :: nodeOwners(:)

    ! An array containing the global ids of the elements to be created
    ! on this PET. This input consists of a 1D array the size of the
    ! number of elements on this PET.
    integer, allocatable         :: elementIds(:)

    ! An array containing the types of the elements to be created on
    ! this PET. The types used must be appropriate for the parametric
    ! dimension of the Mesh. Please see Section 29.2.1 for the list of
    ! options. This input consists of a 1D array the size of the
    ! number of elements on this PET.
    integer, allocatable         :: elementTypes(:)

    ! An array containing the indexes of the sets of nodes to be
    ! connected together to form the elements to be created on this
    ! PET. The entries in this list are NOT node global ids, but
    ! rather each entry is a local index (1 based) into the list of
    ! nodes which were created on this PET by the previous
    ! ESMF_MeshAddNodes() call. In other words, an entry of 1
    ! indicates that this element contains the node described by
    ! nodeIds(1), nodeCoords(1), etc. passed into the
    ! ESMF_MeshAddNodes() call on this PET. It is also important to
    ! note that the order of the nodes in an element connectivity list
    ! matters. Please see Section 29.2.1 for diagrams illustrating the
    ! correct order of nodes in a element. This input consists of a 1D
    ! array with a total size equal to the sum of the number of nodes
    ! in each element on this PET. The number of nodes in each element
    ! is implied by its element type in elementTypes. The nodes for
    ! each element are in sequence in this array (e.g. the nodes for
    ! element 1 are elementConn(1), elementConn(2), etc.).
    integer, allocatable         :: elementConn(:) ! 4*numQuadElems+3*numTriElems

    ! For coupling only these elements are supported.
    ! Cell types
    integer :: TRI = ESMF_MESHELEMTYPE_TRI
    integer :: QUAD = ESMF_MESHELEMTYPE_QUAD

    ! iters
    integer :: i,j,k

    write(*,*) 'Step 1'
    ! Create a Mesh as a 3 step process (dims, nodes, elements)
    mesh = ESMF_MeshCreate(parametricDim, spatialDim, rc)
    ESMF_CHECK

    ! Create the nodes...

    numNodes = size(gridobject%vert_grid)
    ! Fill the indices
    allocate(nodeIds(numNodes))
    forall (i=1:numNodes:1) nodeIds(i) = i

    ! nodeCoords=(/ xk(1), yk(1), xk(2), yk(2), ....
    allocate(nodeCoords(2*numNodes))
    do i=0,(numNodes-1)
       nodeCoords((i*2)+1) = gridobject%vert_grid(i+1)%attr(VERTX)
       nodeCoords((i*2)+2) = gridobject%vert_grid(i+1)%attr(VERTY)
    end do

    ! Set all nodes owned to pet0
    allocate(nodeOwners(numNodes))
    nodeOwners=0

    write(*,*) 'Step 2'
    call ESMF_MeshAddNodes(mesh, nodeIds, nodeCoords, nodeOwners, rc=rc)
    ESMF_CHECK


    ! Let's define the elements

    ! Following example of the reference manual
    numTotElems = size(gridobject%cell_grid)
    numQuadElems = 0
    numTriElems = 0
    ! This is almost similar to the VTK data structure
    allocate(elementTypes(numTotElems))
    allocate(elementIds(numTotElems))

    do i=1,numTotElems
       select case(gridobject%cell_grid(i)%nov)
       case(3)
          elementTypes(i) = TRI
          numTriElems = numTriElems + 1
       case(4)
          elementTypes(i) = QUAD
          numQuadElems = numQuadElems + 1
       case default
          write(*,*) 'Assertion failed expecting elements of 3,4 nodes, got', gridobject%cell_grid(i)%nov, ' for ', i
       end select
    end do

    ! check...
    if (.not. (numTriElems + numQuadElems) .eq. numTotElems) rc=10
    ! A list of all nodes (without the count that vtk uses)
    allocate(elementConn(4*numQuadElems+3*numTriElems))

    j = 1
    do i=1,numTotElems
       if (.not. gridobject%cell_grid(i)%nov .eq. gridobject%cell_grid(i)%nof) then
          write(*,*) 'Assertion equal failed (nov==nof) ', gridobject%cell_grid(i)%nov , gridobject%cell_grid(i)%nof, ' for cell ', i
       end if
       select case(gridobject%cell_grid(i)%nof)
       case(3)
          do k=1,gridobject%cell_grid(i)%nof
             ! Get the first vertex of each face, that should make a triangle, assuming they're ordered.....
             elementConn(j + (k-1)) = gridobject%cell_grid(i)%face(k)%attr(FACEV1)
          end do
          j = j+3
       case(4)
          do k=1,gridobject%cell_grid(i)%nof
             ! Get the first vertex of each face, that should make a triangle, assuming they're ordered.....
             elementConn(j + (k-1)) = gridobject%cell_grid(i)%face(k)%attr(FACEV1)
          end do
          j = j+4
       case default
          write(*,*) 'Assertion failed expecting elements of 3,4 nodes, got', gridobject%cell_grid(i)%nov, ' for ', i
       end select
    end do


    ! Just the counters. (1 based)
    forall (i=1:numTotElems:1) elementIds(i) = i


    write(*,*) 'Step 3'
    ! Now we have everything, let's add them
    call ESMF_MeshAddElements(mesh, elementIds=elementIds, &
         elementTypes=elementTypes,&
         elementConn=elementConn,rc=rc)
    !call ESMF_MeshAddElements(mesh, elementIds, elementTypes, elementConn, rc)
    ESMF_CHECK

    ! Cleanup...
    deallocate(nodeIds)
    deallocate(nodeCoords)
    deallocate(nodeOwners)
    deallocate(elementIds)
    deallocate(elementTypes) !
    deallocate(elementConn)

    ! And we're done.
  end subroutine make_swan_mesh

  ! ! Also make a helper function for fields...
  subroutine make_swan_fieldbundle(mesh, grid,  fieldbundle, rc)
    ! Todo clean up
    use m_genarr
    use swcomm3
    use outp_data

    ! Then we can define the fields
    type(ESMF_Mesh), optional, intent(inout) :: mesh
    type(ESMF_Grid), optional, intent(inout) :: grid
    type(ESMF_FieldBundle), intent(out) :: fieldbundle
    integer, intent(out) :: rc

    type(ESMF_Field) :: field
    type(ESMF_TypeKind_Flag) :: typekind
    type(ESMF_MeshLoc)  :: meshloc
    type(ESMF_ArraySpec) :: arrayspec
    integer :: nelements, rank, i, j, k, idx, ix, iy, shapearr(6)
    integer :: map(8) = (/1,2,3,4,5,6,7,8/)

    integer :: nfields
    real(ESMF_KIND_R8), dimension(:), pointer :: fptr1d
    real(ESMF_KIND_R8), dimension(:,:), pointer :: fptr2d
    real(ESMF_KIND_R8), dimension(:,:,:), pointer :: fptr3d
    real(ESMF_KIND_R8), dimension(:,:,:,:), pointer :: fptr4d
    type(c_ptr) :: cptr

    character(len=MAXSTRLEN) :: name
    character(len=MAXSTRLEN) :: fieldname

    character(kind=c_char) :: c_fieldname(MAXSTRLEN)

    rc = ESMF_FAILURE

    ! Preparation
    ! TODO: check opt_zl argument....
    !call swfld_init(compda, opt_zl=0)

    ! Import states
    name = "swan"
    fieldbundle = ESMF_FieldBundleCreate(name=name, rc=rc)
    ! Lookup the number of variables using the bmi interface
    call get_n_vars(nfields)
    do i=1,nfields
       call get_var_name(i, c_fieldname)
       call get_var_rank(c_fieldname, rank)
       call get_var_shape(c_fieldname, shapearr)
       call get_var(c_fieldname, cptr)
       fieldname = char_array_to_string(c_fieldname, strlen(c_fieldname))

       ! 1 extra rank to (wetpts -> x,y)
       write(*,*) 'writing fieldname', trim(fieldname), ', shape', shapearr(1:rank), ' &', cptr
       if (rank .lt. 1) then
          write(*,*) 'incompatible rank', rank,  'skipping', trim(fieldname)
          cycle
       end if
       call ESMF_ArraySpecSet(arrayspec, typekind=ESMF_TYPEKIND_R8, rank=(rank+1), rc=rc)
       select case(rank)
       case (1)
          call c_f_pointer(cptr, fptr1d, shapearr(1:rank))
          if (.not. associated(fptr1d)) then
             write(*,*) 'No pointer to ',  trim(fieldname), 'skipping'
             cycle
          end if
          field = ESMF_FieldCreate(grid, &
               arrayspec=arrayspec, &
               gridtofieldmap=(/1,2/), &
               staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(fieldname), rc=rc)
          ESMF_CHECK
          call ESMF_FieldGet(field, farrayPtr=fptr2d, rc=rc)
          ESMF_CHECK

          ! Skip if it doesn't fit
          if (size(fptr1d,1) .ne. mcgrd) then
             write(*,*) 'skipping ', trim(fieldname), 'does not fit: ', mcgrd, '<-', size(fptr1d,1)
             cycle
          end if
          do iy = 1, myc
             do ix = 1, mxc
                ! lookup index corresponding to 2d array
                idx = kgrpnt(ix,iy)
                if (idx.gt.1) then
                   fptr2d(ix, iy) = fptr1d(idx)
                end if
             enddo
          enddo
       case (2)
          call c_f_pointer(cptr, fptr2d, shapearr(1:rank))
          field = ESMF_FieldCreate(grid, &
               arrayspec=arrayspec, &
               ungriddedLBound=(/3/), &
               ungriddedUBound=(/size(fptr2d,1)+3-1/), &
               gridtofieldmap=(/1,2/), &
               staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(fieldname), rc=rc)
          ESMF_CHECK
          call ESMF_FieldGet(field, farrayPtr=fptr3d, rc=rc)
          ESMF_CHECK
          do iy = 1, myc
             do ix = 1, mxc
                ! lookup index corresponding to 2d array
                idx = kgrpnt(ix,iy)
                if (idx.gt.1) then
                   fptr3d(ix, iy,:) = fptr2d(idx,:)
                end if
             enddo
          enddo
       case (3)
          call c_f_pointer(cptr, fptr3d, shapearr(1:rank))
          field = ESMF_FieldCreate(grid, &
               arrayspec,&
               ungriddedLBound=(/3,4/), &
               ungriddedUBound=(/size(fptr3d,1)+3-1,size(fptr3d,2)+4-1/), &
               gridtofieldmap=(/1,2/), &
               staggerloc=ESMF_STAGGERLOC_CENTER, name=trim(fieldname), rc=rc)
          ESMF_CHECK
          call ESMF_FieldGet(field, farrayPtr=fptr4d, rc=rc)
          ESMF_CHECK
          do iy = 1, myc
             do ix = 1, mxc
                ! lookup index corresponding to 2d array
                idx = kgrpnt(ix,iy)
                if (idx.gt.1) then
                   fptr4d(ix, iy,:,:) = fptr3d(:,:, idx)
                end if
             enddo
          enddo
       end select
       call ESMF_FieldBundleAdd(fieldbundle, fieldList=(/field/), rc=rc)
       ESMF_CHECK

    end do


  end subroutine make_swan_fieldbundle
end module swan_esmf
