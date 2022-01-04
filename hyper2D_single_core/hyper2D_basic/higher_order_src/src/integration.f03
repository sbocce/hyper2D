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

    ! Update solution using the forward Euler integrator 
    ! INTERNAL CELLS only (4 ghost cells per side)
    do j = 5, Ny-4
      do i = 5, Nx-4

        call compute_fluxes(U, i, j, F_N, F_S, F_W, F_E)

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
    U(:,5:Nx-4,5:Ny-4) = U_new(:,5:Nx-4,5:Ny-4)

  end subroutine

  ! ======================================================================== 
  
  subroutine RK3_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx,Ny), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq, Nx, Ny) :: U1, U2
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E 

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    U1 = U ! Init
    U2 = U ! Init

    do j = 5, Ny-4 ! Internal cells only
      do i = 5, Nx-4

        call compute_fluxes(U, i, j, F_N, F_S, F_W, F_E)

        U1(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)

      end do
    end do

    call assign_BCs(U1)

    do j = 5, Ny-4 ! Internal cells only
      do i = 5, Nx-4

        call compute_fluxes(U1, i, j, F_N, F_S, F_W, F_E)

        U2(:,i,j) = 3.0d0/4.0d0*U(:,i,j) +1.0d0/4.0d0*U1(:,i,j) - &
                    1.0d0/4.0d0*dt/dx*(F_E - F_W) - 1.0d0/4.0d0*dt/dy*(F_N - F_S)

      end do
    end do

    call assign_BCs(U2)

    do j = 5, Ny-4 ! Internal cells only
      do i = 5, Nx-4

        call compute_fluxes(U2, i, j, F_N, F_S, F_W, F_E)

        U_new(:,i,j) = 1.0d0/3.0d0*U(:,i,j) + 2.0d0/3.0d0*U2(:,i,j) - & 
                       2.0d0/3.0d0*dt/dx*(F_E - F_W) - 2.0d0/3.0d0*dt/dy*(F_N - F_S)
        
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
    U(:,5:Nx-4,5:Ny-4) = U_new(:,5:Nx-4,5:Ny-4)

  end subroutine

  ! ======================================================================== 
  
  subroutine midpoint_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx,Ny), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_N, F_S, F_W, F_E 

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    U_new = U ! Init (for ghost cells)

    ! Compute first half step, and write it in U_new
    do j = 5, Ny-4
      do i = 5, Nx-4
        call compute_fluxes(U, i, j, F_N, F_S, F_W, F_E)
        U_new(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)
      end do
    end do

    call assign_BCs(U_new)

    ! Compute full step, but computing fluxes from U_new
    do j = 5, Ny-4
      do i = 5, Nx-4

        call compute_fluxes(U_new, i, j, F_N, F_S, F_W, F_E)

        U(:,i,j) = U(:,i,j) - dt/dx*(F_E - F_W) - dt/dy*(F_N - F_S)
        
        ! Check that the solution did not diverge
        do eqID = 1, Neq
          if (isnan(U(eqID,i,j))) then 
            print*, 'Solution diverged, try with a smaller time step! Aborting.'
            stop
          end if
        end do

      end do
    end do

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes(U, i, j, F_N, F_S, F_W, F_E)

    ! Select which numerical flux to use

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    real(kind=8), dimension(Neq),       intent(out) :: F_N, F_S, F_W, F_E 
    integer, intent(in) :: i, j

    if ( flux_type_Rusanov .eqv. .true. ) then 

      call compute_fluxes_Rusanov(U, i, j, F_N, F_S, F_W, F_E)

    end if 

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

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i,j-2), U(:,i,j-1), U(:,i,j), U(:,i,j+1), U(:,i,j+2), U(:,i,j+3), U_L, U_R)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L) 
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_N = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ South interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i,j-3), U(:,i,j-2), U(:,i,j-1), U(:,i,j), U(:,i,j+1), U(:,i,j+2), U_L, U_R)

    call compute_flux_ws_y(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_y(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_S = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ East interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i-2,j), U(:,i-1,j), U(:,i,j), U(:,i+1,j), U(:,i+2,j), U(:,i+3,j), U_L, U_R)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_E = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ! ------ West interface

    ! Reconstruct the solution, passing the four cells that neighbor the interface
    call reconstruct_sol_interface(U(:,i-3,j), U(:,i-2,j), U(:,i-1,j), U(:,i,j), U(:,i+1,j), U(:,i+2,j), U_L, U_R)

    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    F_W = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

  end subroutine

  ! ======================================================================== 

  subroutine reconstruct_sol_interface(U_m3, U_m2, U_m1, U_p1, U_p2, U_p3, U_L, U_R)
  
  ! Takes the four points around an interface (minus_1, minus_2, plus_1, plus_2) and 
  ! computes the left and right states using first order or MUSCL TVD-limited reconstruction
  
    real(kind=8), dimension(Neq), intent(in)  :: U_m3, U_m2, U_m1, U_p1, U_p2, U_p3
    real(kind=8), dimension(Neq), intent(out) :: U_L, U_R

    ! Working variables
    integer      :: eID              ! equation ID
    real(kind=8) :: theta_L, theta_R ! ratio of consecutive gradients
    real(kind=8), dimension(Neq) :: P_m3, P_m2, P_m1, P_p1, P_p2, P_p3, P_L, P_R ! primitive variables 
    real(kind=8) :: d0, d1, b0, b1, a0, a1, omega0, omega1, P0, P1 ! WENO 3rd order
    real(kind=8) :: b2, a2, omega2, P2                             ! additional for WENO 5th order
     
    if (space_order .EQ. 1) then  ! ***** First order solution *****
  
      U_L = U_m1
      U_R = U_p1
  
    else if (space_order .EQ. 2) then  ! ***** Second order solution - MUSCL TVD *****
  
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

    else if (space_order .eq. 3) then ! ***** Third order in space, WENO 3 *****

      ! Perform reconstruction in primitive variables 
      call compute_primitive_from_conserved(U_m2, P_m2)
      call compute_primitive_from_conserved(U_m1, P_m1)
      call compute_primitive_from_conserved(U_p1, P_p1)
      call compute_primitive_from_conserved(U_p2, P_p2)

      do eID = 1, Neq ! Loop on the quantities
 
        ! ---------- APPROXIMATE v_int^L (left state) ------------

        d0 = 2.0d0/3.0d0
        d1 = 1.0d0/3.0d0

        b0 = (P_p1(eID) - P_m1(eID))**2
        b1 = (P_m1(eID) - P_m2(eID))**2
      
        a0 = d0/(1.0d-6 + b0)**2
        a1 = d1/(1.0d-6 + b1)**2
      
        omega0 = a0/(a0+a1) ! Rescale
        omega1 = a1/(a0+a1) ! Rescale
      
        ! Approximate v for r=0 and r=1
        P0 = 0.5d0*P_m1(eID) + 0.5d0*P_p1(eID)
        P1 = -0.5d0*P_m2(eID) + 3.0d0/2.0d0*P_m1(eID)
      
        P_L(eID) = omega0*P0 + omega1*P1
      
        ! ---------- APPROXIMATE v_int^R (right state) ------------
 
        b0 = (P_p1(eID) - P_m1(eID))**2
        b1 = (P_p2(eID) - P_p1(eID))**2     

        a0 = d0/(1.0d-6 + b0)**2
        a1 = d1/(1.0d-6 + b1)**2     
 
        omega0 = a0/(a0+a1) ! Rescale
        omega1 = a1/(a0+a1) ! Rescale
      
        ! Approximate v for r=0 and r=1
        P0 = 0.5d0*P_m1(eID) + 0.5d0*P_p1(eID)
        P1 = -0.5d0*P_p2(eID) + 3.0d0/2.0d0*P_p1(eID)
      
        P_R(eID) = omega0*P0 + omega1*P1;

      end do

      ! Compute conservative variables
      call compute_conserved_from_primitive(P_L, U_L)
      call compute_conserved_from_primitive(P_R, U_R)

    else if (space_order .eq. 5) then ! ***** Fifth order in space, WENO 5 *****

      ! Perform reconstruction in primitive variables 
      call compute_primitive_from_conserved(U_m3, P_m3)
      call compute_primitive_from_conserved(U_m2, P_m2)
      call compute_primitive_from_conserved(U_m1, P_m1)
      call compute_primitive_from_conserved(U_p1, P_p1)
      call compute_primitive_from_conserved(U_p2, P_p2)
      call compute_primitive_from_conserved(U_p3, P_p3)

      do eID = 1, Neq ! Loop on the quantities

        ! ---------- APPROXIMATE v_int^L (left state) ------------

        b0 = 13.0d0/12.0d0*(P_m3(eID)-2.0*P_m2(eID) + P_m1(eID))**2 &
           + 1.0d0/4.0d0*(P_m3(eID)-4.0*P_m2(eID)+3.0*P_m1(eID))**2

        b1 = 13.0d0/12.0d0*(P_m2(eID)-2.0*P_m1(eID) + P_p1(eID))**2 &
           + 1.0d0/4.0d0*(P_m2(eID)-P_p1(eID))**2

        b2 = 13.0d0/12.0d0*(P_m1(eID)-2.0*P_p1(eID) + P_m2(eID))**2 &
           + 1.0d0/4.0d0*(3.0*P_m1(eID)-4.0*P_p1(eID)+P_p2(eID))**2

        a0 = 1.0d0/10.0d0*(1.0d0/(1.0d-6 + b0)**2)
        a1 = 6.0d0/10.0d0*(1.0d0/(1.0d-6 + b1)**2)
        a2 = 3.0d0/10.0d0*(1.0d0/(1.0d-6 + b2)**2)

        omega0 = a0/(a0+a1+a2)
        omega1 = a1/(a0+a1+a2)
        omega2 = a2/(a0+a1+a2)

        P_L(eID) = omega0*(2.0d0/6.0d0*P_m3(eID) - 7.0d0/6.0d0*P_m2(eID) + 11.0d0/6.0d0*P_m1(eID)) &
                 + omega1*(-1.0d0/6.0d0*P_m2(eID) + 5.0d0/6.0d0*P_m1(eID) + 2.0d0/6.0d0*P_p1(eID)) &
                 + omega2*(2.0d0/6.0d0*P_m1(eID) + 5.0d0/6.0d0*P_p1(eID) - 1.0d0/6.0d0*P_p2(eID))

        ! ---------- APPROXIMATE v_int^R (right state) ------------

        b0 = 13.0d0/12.0d0*(P_p1(eID)-2.0*P_p2(eID) + P_p3(eID))**2 &
           + 1.0d0/4.0d0*(3.0*P_p1(eID)-4.0*P_p2(eID)+P_p3(eID))**2

        b1 = 13.0d0/12.0d0*(P_m1(eID)-2.0*P_p1(eID) + P_p2(eID))**2 &
           + 1.0d0/4.0d0*(P_m1(eID)-P_p2(eID))**2

        b2 = 13.0d0/12.0d0*(P_m2(eID)-2.0*P_m1(eID) + P_p1(eID))**2 &
           + 1.0d0/4.0d0*(P_m2(eID)-4.0*P_m1(eID)+3.0*P_p1(eID))**2

        a0 = 1.0d0/10.0d0*(1.0d0/(1.0d-6 + b0)**2)
        a1 = 6.0d0/10.0d0*(1.0d0/(1.0d-6 + b1)**2)
        a2 = 3.0d0/10.0d0*(1.0d0/(1.0d-6 + b2)**2)

        omega0 = a0/(a0+a1+a2)
        omega1 = a1/(a0+a1+a2)
        omega2 = a2/(a0+a1+a2)

        P_R(eID) = omega2*(-1.0d0/6.0d0*P_m2(eID) + 5.0d0/6.0d0*P_m1(eID) + 2.0d0/6.0d0*P_p1(eID)) &
                 + omega1*(2.0d0/6.0d0*P_m1(eID)  + 5.0d0/6.0d0*P_p1(eID) - 1.0d0/6.0d0*P_p2(eID)) &
                 + omega0*(11.0d0/6.0d0*P_p1(eID) - 7.0d0/6.0d0*P_p2(eID) + 2.0d0/6.0d0*P_p3(eID)) 
 
      end do

      ! Compute conservative variables
      call compute_conserved_from_primitive(P_L, U_L)
      call compute_conserved_from_primitive(P_R, U_R)

    else 

      write(*,*) "Attention! Order of the space integration not implemented! Check global_module file!"
      write(*,*) "Aborting!"
      stop

    end if

  end subroutine

  ! ========================================================================

  real(kind=8) function phi_limiter(theta)

    real(kind=8), intent(in) :: theta ! ratio of consecutive gradients

    phi_limiter = (theta*theta + theta)/(theta*theta + 1.0) ! van Albada symmetric limiter
    ! phi_limiter = (theta + abs(theta))/(1.0 + abs(theta)) ! van Leer limiter

  end function 

end module
