module integration

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
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E, S
    real(kind=8) :: dummy

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    do j = 3, Ny-2
      do i = 3, Nx-2

        call compute_fluxes_HLL(U, i, j, F_N, F_S, F_W, F_E)
        ! call compute_fluxes_Rusanov(U, i, j, F_N, F_S, F_W, F_E)

        call introduce_flat_plate(U, i, j, F_N, F_S, F_W, F_E)

        call compute_source_term(U(:,i,j), S, dummy, dummy, dummy)

        ! Update the solution for all terms explicitly
        U_new(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S) + dt*S
        
        ! Check that the solution did not diverge
        do eqID = 1, Neq
          if (isnan(U_new(eqID,i,j))) then 
            print*, 'Solution diverged, try with a smaller time step! Aborting.'
            print*, 'Solution that diverged: ', U_new(:,i,j)
            print*, 'in cell i,j = ', i, j
            stop
          end if
        end do

      end do
    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,3:Nx-2,3:Ny-2) = U_new(:,3:Nx-2,3:Ny-2)

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes_HLL(U, i, j, F_N, F_S, F_W, F_E)

    ! Computes Rusanov numerical fluxes for the cell (i,j)

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    real(kind=8), dimension(Neq),       intent(out) :: F_N, F_S, F_W, F_E 
    integer, intent(in) :: i, j

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R

    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    ! ------ North interface
    U_L = U(:,i,j)
    U_R = U(:,i,j+1)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)

    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_N = F_L
    else if (ws_max .lt. 0.0) then
      F_N = F_R
    else
      F_N = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

    ! Update global maximum wave speed (used for setting the time step)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ South interface
    U_L = U(:,i,j-1)
    U_R = U(:,i,j)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)

    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_S = F_L
    else if (ws_max .lt. 0.0) then
      F_S = F_R
    else
      F_S = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

    ! Update global maximum wave speed (used for setting the time step)
    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ East interface
    U_L = U(:,i,j)
    U_R = U(:,i+1,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)
    
    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_E = F_L
    else if (ws_max .lt. 0.0) then
      F_E = F_R
    else
      F_E = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

    ! Update global maximum wave speed (used for setting the time step)
    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ West interface
    U_L = U(:,i-1,j)
    U_R = U(:,i,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_min = MIN(ws_min_L, ws_min_R)
    ws_max = MAX(ws_max_L, ws_max_R)

    ! HLL fluxes
    if (ws_min .ge. 0.0) then
      F_W = F_L
    else if (ws_max .lt. 0.0) then
      F_W = F_R
    else
      F_W = (ws_min*ws_max*(U_R - U_L) + ws_max*F_L - ws_min*F_R)/(ws_max - ws_min);
    end if

    ! Update global maximum wave speed (used for setting the time step)
    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))
    ws_maxabs = MAX(ws_maxabs, ws_max)

  end subroutine

  ! ==============================================

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

    ! Update global maximum wave speed (used for setting the time step)
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ South interface
    U_L = U(:,i,j-1)
    U_R = U(:,i,j)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_S = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! Update global maximum wave speed (used for setting the time step)
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ East interface
    U_L = U(:,i,j)
    U_R = U(:,i+1,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_E = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! Update global maximum wave speed (used for setting the time step)
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ------ West interface
    U_L = U(:,i-1,j)
    U_R = U(:,i,j)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_W = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! Update global maximum wave speed (used for setting the time step)
    ws_maxabs = MAX(ws_maxabs, ws_max)

  end subroutine

  ! ==============================================


  ! ==============================================

  subroutine introduce_flat_plate(U, i, j, F_N, F_S, F_W, F_E)

    ! This function introduces a flat plate, by modifying the fluxes, only for the cells 
    ! that are in contact with the plate.

    implicit none

    integer, intent(in) :: i, j
    real(kind=8), dimension(Neq), intent(inout)    :: F_N, F_S, F_W, F_E
    real(kind=8), dimension(Neq,Nx,Ny), intent(in) :: U

    real(kind=8), dimension(Neq) :: prim_1, prim_2 ! Primitive variables
    real(kind=8) :: P_1, P_2, Pw

    ! The indices j_plate, i_1_plate, i_2_plate identify the plate position.
    ! They are defined in global_module.f03

    if (j .EQ. j_plate) then ! Cells below the plate

      if ( (i .GE. i_1_plate) .AND. (i .LE. i_2_plate) ) then ! Cells in contact with the plate

        ! -------
        ! These cells are BELOW the plate. Change F_N 
        ! -------

        ! Simple computation of the wall pressure
        ! call compute_primitive_from_conserved(U(:,i,j), prim_1)
        ! P_1 = prim_1(1)/M*kB*prim_1(4) ! P = n*kB*T = rho/M*kB*T
        ! Pw  = P_1 ! Reconstruct wall pressure

        ! Two-point extrapolation of the wall pressure (see Blazek, Computational Fluid Dynamics)
        call compute_primitive_from_conserved(U(:,i,j),   prim_1)
        call compute_primitive_from_conserved(U(:,i,j-1), prim_2)
        P_1 = prim_1(1)/M*kB*prim_1(4) ! P = n*kB*T = rho/M*kB*T
        P_2 = prim_2(1)/M*kB*prim_2(4) ! P = n*kB*T = rho/M*kB*T
        Pw  = 0.5*(3.0*P_1 - P_2) ! Reconstruct wall pressure

        ! Only pressure normal to the wall remains
        F_N(1) = 0.0
        F_N(2) = 0.0
        F_N(3) = Pw
        F_N(4) = 0.0
        F_N(5) = 0.0

      end if

    else if (j .EQ. j_plate + 1) then ! Cells ABOVE the plate

      if ( (i .GE. i_1_plate) .AND. (i .LE. i_2_plate) ) then ! Cells in contact with the plate

        ! -------
        ! These cells are ABOVE the plate. Change F_S
        ! -------

        ! Simple computation of the wall pressure
        ! call compute_primitive_from_conserved(U(:,i,j), prim_1)
        ! P_1 = prim_1(1)/M*kB*prim_1(4) ! P = n*kB*T = rho/M*kB*T
        ! Pw  = P_1 ! Reconstruct wall pressure

        ! Two-point extrapolation of the wall pressure (see Blazek, Computational Fluid Dynamics)
        call compute_primitive_from_conserved(U(:,i,j),   prim_1)
        call compute_primitive_from_conserved(U(:,i,j+1), prim_2)
        P_1 = prim_1(1)/M*kB*prim_1(4) ! P = n*kB*T = rho/M*kB*T
        P_2 = prim_2(1)/M*kB*prim_2(4) ! P = n*kB*T = rho/M*kB*T
        Pw  = 0.5*(3.0*P_1 - P_2) ! Reconstruct wall pressure

        ! Only pressure normal to the wall remains
        F_S(1) = 0.0
        F_S(2) = 0.0
        F_S(3) = Pw
        F_S(4) = 0.0
        F_S(5) = 0.0

      end if

    end if

  end subroutine

end module
