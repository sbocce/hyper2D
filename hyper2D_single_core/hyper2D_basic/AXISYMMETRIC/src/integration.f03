module integration

! This module contains the routines for the time integration and the computation
! of numerical fluxes (Rusanov or HLL in this example).

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
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E, Src
    real(kind=8) :: y_j, y_jmhalf, y_jphalf

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    do j = 3, Ny-2
      do i = 3, Nx-2

        ! call compute_fluxes_Rusanov(U, i, j, F_N, F_S, F_W, F_E)
        call compute_fluxes_HLL(U, i, j, F_N, F_S, F_W, F_E)

        ! Symmetry BC
        if (j.eq.3) then
          call compute_flux_y_SYM(U(:,i,j), F_S)
        end if

        ! Compute geometrical source term (basically, the pressure)
        call compute_source_term(U,i,j, Src)

        y_j      = y_min + dy/2.0d0 + (real(j)-3.0d0)*dy
        y_jmhalf = y_j - dy/2.0d0
        y_jphalf = y_j + dy/2.0d0

        U_new(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/(dy*y_j)*(F_N*y_jphalf - F_S*y_jmhalf) + dt*Src
        
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

  subroutine compute_fluxes_HLL(U, i, j, F_N, F_S, F_W, F_E)

    ! Computes HLL numerical fluxes for the cell (i,j)

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    real(kind=8), dimension(Neq),       intent(out) :: F_N, F_S, F_W, F_E
    integer, intent(in) :: i, j

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R

    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    ! ------ North interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i,j-1), U(:,i,j), U(:,i,j+1), U(:,i,j+2), U_L, U_R)

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

    ! ------ South interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i,j-2), U(:,i,j-1), U(:,i,j), U(:,i,j+1), U_L, U_R)

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

    ! ------ East interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i-1,j), U(:,i,j), U(:,i+1,j), U(:,i+2,j), U_L, U_R)

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

    ! ------ West interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i-2,j), U(:,i-1,j), U(:,i,j), U(:,i+1,j), U_L, U_R)

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

  end subroutine

  ! ==============================================================

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

  ! ======================================================================== 

  subroutine reconstruct_sol_interface(U_m2, U_m1, U_p1, U_p2, U_L, U_R)
  
  ! Takes the four points around an interface (minus_1, minus_2, plus_1, plus_2) and 
  ! computes the left and right states using first order or MUSCL TVD-limited reconstruction
  
    real(kind=8), dimension(Neq), intent(in)  :: U_m2, U_m1, U_p1, U_p2
    real(kind=8), dimension(Neq), intent(out) :: U_L, U_R

    ! Working variables
    integer      :: eID              ! equation ID
    real(kind=8) :: theta_L, theta_R ! ratio of consecutive gradients
    real(kind=8), dimension(Neq) :: P_m2, P_m1, P_p1, P_p2, P_L, P_R ! primitive variables 
     
    if (bool_MUSCL .EQV. .FALSE.) then
    ! ***** First order solution *****
  
      U_L = U_m1
      U_R = U_p1
  
    else
      ! ***** Second order solution *****
  
      ! Perform reconstruction in primitive variables 
      call compute_primitive_from_conserved(U_m2, P_m2)
      call compute_primitive_from_conserved(U_m1, P_m1)
      call compute_primitive_from_conserved(U_p1, P_p1)
      call compute_primitive_from_conserved(U_p2, P_p2)

      do eID = 1, Neq ! Loop on the quantities
 
        ! Compute ratio of consecutive gradients
        theta_L = (P_p1(eID) - P_m1(eID))/(P_m1(eID) - P_m2(eID) + 1.0d-15) ! Add a small tolerance
        theta_R = (P_p2(eID) - P_p1(eID))/(P_p1(eID) - P_m1(eID) + 1.0d-15) ! Add a small tolerance
      
        ! Limited linear reconstruction 
        P_L(eID) = P_m1(eID) + phi_limiter(theta_L)*(P_m1(eID) - P_m2(eID))/2.0
        P_R(eID) = P_p1(eID) - phi_limiter(1.0d0/(theta_R+1.0d-15))*(P_p2(eID) - P_p1(eID))/2.0

      end do

      ! Compute conservative variables
      call compute_conserved_from_primitive(P_L, U_L)
      call compute_conserved_from_primitive(P_R, U_R)

    end if

  end subroutine

  ! ========================================================================

  real(kind=8) function phi_limiter(theta)

    real(kind=8), intent(in) :: theta ! ratio of consecutive gradients

    phi_limiter = (theta*theta + theta)/(theta*theta + 1.0) ! van Albada symmetric limiter
    ! phi_limiter = (theta + abs(theta))/(1.0 + abs(theta)) ! van Leer limiter

  end function 



end module
