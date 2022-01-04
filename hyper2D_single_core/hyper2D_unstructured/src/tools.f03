module tools

  use global_module
  use grid
  use pde 

  implicit none

  contains

  ! ####################################################################
  
  subroutine export_sol_vtk(t_ID, U)

    ! ---------------------------------------------------------------------------------------
    ! This subroutine exports data in the legacy VTK format, that can be visualized using
    ! ParaView.
    ! This helps for checking that all is good.
    ! The output file is put in the "dumps" directory and is named "sol_..." where ... is the
    ! time ID "t_ID" given as an input to the subroutine.
    ! U: conserved variables
    ! 
    ! Notes: 
    !  - All cells are printed, including the two layers of ghost cells for each boundary.
    !    This helps for debugging. 
    !  - The data is exported for each cell center. In ParaView, you may visualize the data 
    !    using the option "Surface with Edges": this is misleading, and these edges do not
    !    represent the cells boundaries. Instead, each point is a cell center.
    ! ---------------------------------------------------------------------------------------
  
    integer, intent(in) :: t_ID
    real(kind=8), dimension(:,:), intent(in) :: U

    real(kind=8), dimension(:,:), allocatable :: prim ! Primitive variables
    
    integer      :: eqID, i

    character(len=512)  :: file_name  ! Name of the output file

    ! ----- Compute primitive variables on the grid ------

    allocate(prim(Neq,Nele))

    prim = 0.0 ! Init
    do i = 1, Nele
      call compute_primitive_from_conserved(U(:,i), prim(:,i))
    end do

    ! ------- Write VTK file -------

    ! Open file, and replace if existing. 11189 is just an arbitrary ID.
    write(file_name,'(A, I8.8, A)') './dumps/sol_', t_ID, '.vtk'
    open(11189, file=file_name, status="replace", form="formatted")

    ! The VTK file starts with this header
    write(11189, '(A)') "# vtk DataFile Version 2.0"
    write(11189, '(A30,I8)') "hyper2D solution at timestep ", t_ID
    write(11189, '(A)') "ASCII"

    ! Write the information about the grid: number of points etc etc
    write(11189, '(A)') "DATASET UNSTRUCTURED_GRID"
    write(11189, '(A)') " "
    write(11189, '(A,I8,A)') "POINTS ", Nnodes, " float"

    do i = 1,Nnodes
      write(11189, '(EN17.5E3,EN17.5E3,EN17.5E3)') nodes_xy(i,1), nodes_xy(i,2), 0.0 ! WRITE SOL HERE
    end do

    write(11189, '(A)') " "
    write(11189, '(A,I8,I8)') "CELLS ", Nele, Nele*5 ! For triangles, use Nele*4
    do i = 1,Nele 
      ! Note: VTK is written in C++ and uses indices starting from zero
      write(11189, '(I8, I8, I8, I8, I8)') 4, ele_nodes(i,1)-1, ele_nodes(i,2)-1, &
                                              ele_nodes(i,3)-1, ele_nodes(i,4)-1
    end do

    write(11189, '(A)') " "
    write(11189, '(A,I8)') "CELL_TYPES ", Nele
    do i = 1,Nele
      write(11189, '(I8)') 9
    end do
 
    write(11189, '(A)') " "
    write(11189, '(A,I8)') "CELL_DATA ", Nele

    do eqID = 1, Neq
      write(11189, '(A,A,A)') "SCALARS ", prim_names(eqID), "float 1"
      write(11189, '(A)')     "LOOKUP_TABLE default"

      do i = 1,Nele
        write(11189, '(EN17.5E3,A)', advance="no") prim(eqID,i), " " ! WRITE SOL HERE
      end do
 
      write(11189,*) " "

    end do

!     ! write(11189,'(A)') "METADATA"
!     ! write(11189,'(A)') "INFORMATION 0"
   
    close(11189)
  
  end subroutine

end module 
