module integration

! This module contains the routines for the time integration and the computation
! of numerical fluxes (Rusanov in this example).

  use pde
  use global_module

  implicit none
  
  real(kind=8), dimension(Neq,Nx,Ny) :: U_new ! Working variable

  contains

  ! ======================================================================== 
  
  subroutine forward_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx,Ny), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E 

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    do j = 3, Ny-2
      do i = 3, Nx-2

        call compute_fluxes_Rusanov(U, i, j, F_N, F_S, F_W, F_E)

        call introduce_flat_plate(U, i, j, F_N, F_S, F_W, F_E)

        U_new(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)
        
        ! Check that the solution did not diverge
        do eqID = 1, Neq
          if (isnan(U_new(eqID,i,j))) then 
            print*, 'Solution diverged, try with a smaller time step! Aborting.'
            stop
          end if
        end do

      end do
    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,3:Nx-2,3:Ny-2) = U_new(:,3:Nx-2,3:Ny-2)

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes_Rusanov(U, i, j, F_N, F_S, F_W, F_E)

    ! Computes Rusanov numerical fluxes for the cell (i,j)

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    real(kind=8), dimension(Neq),       intent(out) :: F_N, F_S, F_W, F_E 
    integer, intent(in) :: i, j

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R

    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_max ! wave speeds

    ! ------ North interface
    U_L = U(:,i,j)
    U_R = U(:,i,j+1)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L) 
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_N = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ South interface
    U_L = U(:,i,j-1)
    U_R = U(:,i,j)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_S = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ East interface
    U_L = U(:,i,j)
    U_R = U(:,i+1,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_E = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ West interface
    U_L = U(:,i-1,j)
    U_R = U(:,i,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_W = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

  end subroutine

  ! ==============================================

  subroutine introduce_flat_plate(U, i, j, F_N, F_S, F_W, F_E)

    ! This function introduces a flat plate, by modifying the fluxes, only for the cells 
    ! that are in contact with the plate.

    implicit none

    integer, intent(in) :: i, j
    real(kind=8), dimension(Neq), intent(inout)    :: F_N, F_S, F_W, F_E
    real(kind=8), dimension(Neq,Nx,Ny), intent(in) :: U

    real(kind=8), dimension(Neq) :: prim ! Primitive variables
    real(kind=8) :: P

    ! The indices j_plate, i_1_plate, i_2_plate identify the plate position.
    ! They are defined in global_module.f03

    if (j .EQ. j_plate) then ! Cells below the plate

      if ( (i .GE. i_1_plate) .AND. (i .LE. i_2_plate) ) then ! Cells in contact with the plate

        ! These cells are BELOW the plate. Change F_N 
        call compute_primitive_from_conserved(U(:,i,j), prim)

        P = prim(1)/M*kB*prim(4) ! P = n*kB*T = rho/M*kB*T

        ! Only pressure normal to the wall remains
        F_N(1) = 0.0
        F_N(2) = 0.0
        F_N(3) = P
        F_N(4) = 0.0

      end if

    else if (j .EQ. j_plate + 1) then

      if ( (i .GE. i_1_plate) .AND. (i .LE. i_2_plate) ) then ! Cells in contact with the plate

        ! These cells are ABOVE the plate. Change F_S
        call compute_primitive_from_conserved(U(:,i,j), prim)

        P = prim(1)/M*kB*prim(4) ! P = n*kB*T = rho/M*kB*T

        ! Only pressure normal to the wall remains
        F_S(1) = 0.0
        F_S(2) = 0.0
        F_S(3) = P
        F_S(4) = 0.0

      end if

    end if

  end subroutine

end module
