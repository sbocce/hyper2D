module tools

  use global_module
  use pde 

  implicit none

  contains

  ! ####################################################################
  
  subroutine export_sol(t_ID, t_now, U)

    ! This subroutine writes the solution on a dat file.
    ! Only internal cells (no ghost cells), are written
     
    integer,                      intent(in) :: t_ID
    real(kind=8),                 intent(in) :: t_now
    real(kind=8), dimension(:,:), intent(in) :: U

    real(kind=8), dimension(Neq,Nx) :: prim ! Primitive variables
    
    integer :: eqID, i

    character(len=512)  :: file_name  ! Name of the output file

    ! ====== Open file for writing
    write(file_name,'(A, I8.8, A)') './dumps/sol_', t_ID, '.dat'
    open(11189, file=file_name, status="replace", form="formatted")

    ! ====== Write header
    write(11189,'(A)',advance="no") '# t        x           '

    do eqID = 1, Neq
      write(11189, '(A)', advance="no") prim_names(eqID)
    end do

    write(11189, '(A)') " "

    ! ====== Write the solution 
    prim = 0.0 ! Init

    do i = 3, Nx-2 ! Only internal cells
      call compute_primitive_from_conserved(U(:,i), prim(:,i))

      write(11189, '(EN17.5E3,A,EN17.5E3,A)', advance="no") t_now, " ", x_min + dx/2.0 + dx*(i-3.0), " " ! Time and position

      write(11189, *) prim(:,i), " " ! Solution in primitive variables 

    end do

    close(11189)
  
  end subroutine

end module 
