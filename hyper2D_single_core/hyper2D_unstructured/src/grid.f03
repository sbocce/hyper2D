module grid

  implicit none

  integer :: Nele, Nnodes ! Number of elements and nodes

  real(kind=8), dimension(:,:),  allocatable :: ele_geom       ! Geometry of the cells
  real(kind=8), dimension(:,:),  allocatable :: ele_int_len    ! Interfaces length
  real(kind=8), dimension(:,:),  allocatable :: ele_int_nx     ! Interfaces normals along x
  real(kind=8), dimension(:,:),  allocatable :: ele_int_ny     ! Interfaces normals along y
  integer,      dimension(:,:),  allocatable :: ele_neigh      ! Connectivity (neighbors)
  logical,      dimension(:),    allocatable :: ele_bound_bool ! Useful for going second order
  integer,      dimension(:,:),  allocatable :: ele_nodes      ! Necessary for plotting the solution 

  real(kind=8), dimension(:,:),  allocatable :: nodes_xy ! Necessary for exporting the solution to VTK

  contains 

  ! ================================================

  subroutine load_grid_from_file

    ! This subroutine loads the grid data from a file

    implicit none
   
    character(len=512) :: dummy
    integer :: eleID, nodeID, dummyINT

    integer :: neigh_12, neigh_23, neigh_34, neigh_41
    integer :: p1, p2, p3, p4
    real(kind=8) :: A, xC, yC
    real(kind=8) :: L12, n12_x, n12_y
    real(kind=8) :: L23, n23_x, n23_y
    real(kind=8) :: L34, n34_x, n34_y
    real(kind=8) :: L41, n41_x, n41_y

    real(kind=8) :: xP, yP

    ! -------- READ NODES --------

    open(98765, file="nodes.hyp") ! Open file

    read(98765, *) Nnodes  ! Read number of elements
 
    allocate(nodes_xy(Nnodes, 2))

    do nodeID = 1, Nnodes

      read(98765, *) xP, yP

      nodes_xy(nodeID,1) = xP
      nodes_xy(nodeID,2) = yP

    end do

    close(98765)

    ! -------- READ ELEMENTS --------

    open(54321, file="mesh.hyp") ! Open file

    read(54321, *) dummy ! Read header line (useless)
    read(54321, *) Nele  ! Read number of elements

    write(*,*) "  Number of grid elements: ", Nele

    ! Allocate arrays
    allocate(ele_geom(Nele, 3))
    allocate(ele_neigh(Nele, 4))
    allocate(ele_int_len(Nele, 4))
    allocate(ele_int_nx(Nele, 4))
    allocate(ele_int_ny(Nele, 4))
    allocate(ele_nodes(Nele, 4)) 
    allocate(ele_bound_BOOL(Nele))

    ! Read line by line 
    do eleID = 1, Nele

      read(54321, *) dummyINT, A, xC, yC, L12, n12_x, n12_y, neigh_12, &
                                          L23, n23_x, n23_y, neigh_23, &
                                          L34, n34_x, n34_y, neigh_34, &
                                          L41, n41_x, n41_y, neigh_41, &
                                          p1,  p2,  p3,  p4

      ele_geom(eleID,1)  = A
      ele_geom(eleID,2)  = xC
      ele_geom(eleID,3)  = yC

      ele_int_len(eleID,1) = L12
      ele_int_len(eleID,2) = L23
      ele_int_len(eleID,3) = L34
      ele_int_len(eleID,4) = L41

      ele_int_nx(eleID,1) = n12_x ! First interface
      ele_int_nx(eleID,2) = n23_x ! second int
      ele_int_nx(eleID,3) = n34_x ! third 
      ele_int_nx(eleID,4) = n41_x ! fourth interface

      ele_int_ny(eleID,1) = n12_y ! First interface
      ele_int_ny(eleID,2) = n23_y ! second int
      ele_int_ny(eleID,3) = n34_y ! third 
      ele_int_ny(eleID,4) = n41_y ! fourth interface

      ele_neigh(eleID,1) = neigh_12
      ele_neigh(eleID,2) = neigh_23
      ele_neigh(eleID,3) = neigh_34
      ele_neigh(eleID,4) = neigh_41

      ! Set a bool to assess quickly if the element is boundarying
      ! (useful for second order)

      if ( (neigh_12 .le. 0) .or. (neigh_23 .le. 0) .or. &
           (neigh_34 .le. 0) .or. (neigh_41 .le. 0) ) then
        ele_bound_bool(eleID) = .TRUE.
      else
        ele_bound_bool(eleID) = .FALSE.
      end if 

      ! Save nodes, for plotting the solution
      ele_nodes(eleID,1) = p1
      ele_nodes(eleID,2) = p2
      ele_nodes(eleID,3) = p3
      ele_nodes(eleID,4) = p4

    end do

    close(54321)

    !! DBDBDB FOR DEBUG !! ! WRITE STUFF, FOR DEBUGGING PURPOSES 
    !! DBDBDB FOR DEBUG !! ! -------------------
    !! DBDBDB FOR DEBUG !! open(5, file="test/nodes_debug.dat", status="replace", form="formatted")
    !! DBDBDB FOR DEBUG !! do nodeID = 1, Nnodes
    !! DBDBDB FOR DEBUG !!   write(5,*)  nodes_xy(nodeID,1), nodes_xy(nodeID,2)
    !! DBDBDB FOR DEBUG !! end do
    !! DBDBDB FOR DEBUG !!
    !! DBDBDB FOR DEBUG !! close(5)
    !! DBDBDB FOR DEBUG !!
    !! DBDBDB FOR DEBUG !! ! ELEMENTS -------------------
    !! DBDBDB FOR DEBUG !! open(5, file="test/elements_debug.dat", status="replace", form="formatted")
    !! DBDBDB FOR DEBUG !!
    !! DBDBDB FOR DEBUG !! do eleID = 1, Nele
    !! DBDBDB FOR DEBUG !!   write(5,*)  ele_nodes(eleID,1), ele_nodes(eleID,2), ele_nodes(eleID,3), ele_nodes(eleID,4), &
    !! DBDBDB FOR DEBUG !!               ele_neigh(eleID,1), ele_neigh(eleID,2), ele_neigh(eleID,3), ele_neigh(eleID,4), &
    !! DBDBDB FOR DEBUG !!               ele_geom(eleID,1), ele_geom(eleID,2), ele_geom(eleID,3), &  ! A, xC, yC
    !! DBDBDB FOR DEBUG !!               ele_int_len(eleID,1),ele_int_len(eleID,2),ele_int_len(eleID,3),ele_int_len(eleID,4), &
    !! DBDBDB FOR DEBUG !!               ele_int_nx(eleID,1),ele_int_nx(eleID,2),ele_int_nx(eleID,3),ele_int_nx(eleID,4),     &
    !! DBDBDB FOR DEBUG !!               ele_int_ny(eleID,1),ele_int_ny(eleID,2),ele_int_ny(eleID,3),ele_int_ny(eleID,4)
    !! DBDBDB FOR DEBUG !! end do
    !! DBDBDB FOR DEBUG !!
    !! DBDBDB FOR DEBUG !! close(5)


  end subroutine

end module
