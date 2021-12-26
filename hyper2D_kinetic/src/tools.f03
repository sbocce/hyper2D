module tools

  use global_module
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
    !  - The domain is scaled with the simulated domain size. The x and y axes will span
    !    the domain (-0.5, 0.5). This is because the physical space and velocity are usually
    !    very different, and we need scaling for the representation.
    ! ---------------------------------------------------------------------------------------
  
    integer, intent(in) :: t_ID
    real(kind=8), dimension(Neq,Nx,Ny), intent(in) :: U

    integer      :: Nx_int, Ny_int
    integer      :: i, j

    character(len=512)  :: file_name  ! Name of the output file

    real(kind=8) :: x_scale, y_scale ! for scaling the axes

    x_scale = 1.0/(x_max - x_min)
    y_scale = 1.0/(y_max - y_min)

    ! Number of internal cells
    Nx_int = Nx - 4
    Ny_int = Ny - 4

    ! ------- Write VTK file -------

    ! Open file, and replace if existing. 11189 is just an arbitrary ID.
    write(file_name,'(A, I8.8, A)') './dumps/sol_', t_ID, '.vtk'
    open(11189, file=file_name, status="replace", form="formatted")

    ! The VTK file starts with this header
    write(11189, '(A)') "# vtk DataFile Version 2.0"
    write(11189, '(A30,I8)') "hyper2D solution at timestep ", t_ID
    write(11189, '(A)') "ASCII"

    ! Write the information about the grid: number of points etc etc
    write(11189, '(A)') "DATASET STRUCTURED_POINTS"
    write(11189, '(A,I8,I8,I8)') "DIMENSIONS ", Nx_int, Ny_int, 1 ! Write only internal cells
    write(11189, '(A,F14.7,F14.7,F14.7)') "SPACING ", dx*x_scale, dy*y_scale, 0.0
    write(11189, '(A,F14.7,F14.7,F14.7)') "ORIGIN  ", (x_min+dx/2.0)*x_scale, (y_min+dy/2.0)*y_scale, 0.0

    ! Now write the fields as "point data", each point is a cell center.
    write(11189, '(A,I10)') "POINT_DATA ", Nx_int*Ny_int
  
    ! The array "prim_names" is defined in the "pde" module
    write(11189, '(A,A,A)') "SCALARS ", prim_names(1), " float 1"
    write(11189, '(A)') "LOOKUP_TABLE default"
  
    ! Print value (only internal cells, not ghost cells)
    do j = 3, Ny-2
      do i = 3, Nx-2
        write(11189, '(EN17.5E3,A)', advance="no") U(1,i,j), " " ! WRITE SOL HERE
      end do
    end do
 
    write(11189,*) " "
    write(11189,*) " "

    ! write(11189,'(A)') "METADATA"
    ! write(11189,'(A)') "INFORMATION 0"
  
    close(11189)
  
  end subroutine

end module 
