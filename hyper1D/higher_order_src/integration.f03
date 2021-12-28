module integration

  use pde
  use global_module

  implicit none
  
  real(kind=8), dimension(Neq,Nx) :: U_new ! Working variable

  contains

  ! ======================================================================== 
  
  subroutine forward_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_W, F_E

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only
    do i = 5, Nx-4

      call compute_fluxes(U, i, F_W, F_E)

      ! Advect remaining stuff
      U_new(:,i) = U(:,i) - dt/dx*(F_E - F_W)  ! No source term

      ! Check that the solution did not diverge
      do eqID = 1, Neq
        if (isnan(U_new(eqID,i))) then 
          print*, 'Solution diverged, try with a smaller time step! Aborting.'
          print*, 'Solution that diverged: ', U_new(:,i)
          print*, 'in cell i = ', i
          stop
        end if
      end do

    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,5:Nx-4) = U_new(:,5:Nx-4)

  end subroutine

  ! ======================================================================== 
 
  subroutine midpoint_Euler_step(U, dt)

    ! This function performs one step of the Forward Euler explicit time integrator 

    implicit none
 
    real(kind=8), dimension(Neq,Nx), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq) :: F_W, F_E

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    U_new = U ! Init (for ghost cells)

    ! Compute first half step - Create U_new (it is actually "U_half")
    do i = 5, Nx-4 ! Only internal cells
      call compute_fluxes(U, i, F_W, F_E)
      U_new(:,i) = U(:,i) - dt/2.0/dx*(F_E - F_W)  ! No source term
    end do 

    call assign_BCs(U_new) ! See the pde.f03 module

    ! Now compute the second half step:
    ! use U_new to compute the fluxes (!!!), then update U.
    do i = 5, Nx-4 ! Internal cells only

      call compute_fluxes(U_new, i, F_W, F_E)

      ! Update U
      U(:,i) = U(:,i) - dt/dx*(F_E - F_W)  ! No source term

      ! Check that the solution did not diverge
      do eqID = 1, Neq
        if (isnan(U(eqID,i))) then 
          print*, 'Solution diverged, try with a smaller time step! Aborting.'
          print*, 'Solution that diverged: ', U_new(:,i)
          print*, 'in cell i = ', i
          stop
        end if
      end do

    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    ! Nope, U was already updated.

  end subroutine

  ! ======================================================================== 
 
  subroutine RK3_step(U, dt)

    ! This function performs one time integration step using the 
    ! 3rd order Runge Kutta explicit method of Gottlieb & Shu (1998)

    implicit none
 
    real(kind=8), dimension(Neq,Nx), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq,Nx) :: U1, U2
    real(kind=8), dimension(Neq)    :: F_W, F_E

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only

    U1 = U ! Init
    U2 = U ! Init

    do i = 5, Nx-4
      call compute_fluxes(U, i, F_W, F_E)
      U1(:,i) = U(:,i) - dt/dx*(F_E - F_W)  ! No source term
    end do
    call assign_BCs(U1) ! See the pde.f03 module

    do i = 5, Nx-4
      call compute_fluxes(U1, i, F_W, F_E)
      U2(:,i) = 3.0d0/4.0d0*U(:,i) + 1.0d0/4.0d0*U1(:,i) - 1.0d0/4.0d0*dt/dx*(F_E - F_W)  ! No source term
    end do
    call assign_BCs(U2) ! See the pde.f03 module
 
    do i = 5, Nx-4
      call compute_fluxes(U2, i, F_W, F_E)
      U_new(:,i) = 1.0d0/3.0d0*U(:,i) + 2.0d0/3.0d0*U2(:,i) - 2.0d0/3.0d0*dt/dx*(F_E - F_W)  ! No source term

      ! Check that the solution did not diverge
      do eqID = 1, Neq
        if (isnan(U_new(eqID,i))) then 
          print*, 'Solution diverged, try with a smaller time step! Aborting.'
          print*, 'Solution that diverged: ', U_new(:,i)
          print*, 'in cell i = ', i
          stop
        end if
      end do

    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,5:Nx-4) = U_new(:,5:Nx-4)

  end subroutine

  ! ======================================================================== 

  subroutine RK4_step(U, dt)

    ! This function performs one time integration step using the 
    ! 3rd order Runge Kutta explicit method of Gottlieb & Shu (1998)

    implicit none
 
    real(kind=8), dimension(Neq,Nx), intent(inout)  :: U

    real(kind=8), intent(in) :: dt
    integer :: i, j, eqID
    
    real(kind=8), dimension(Neq,Nx) :: U1, U2, U3
    real(kind=8), dimension(Neq)    :: F0_W, F0_E, F1_W, F1_E, F2_W, F2_E, F3_W, F3_E

    ! First, assign BCs. Wall BCs need to be updated in time, taking the value
    ! from neighboring cells.
    call assign_BCs(U) ! See the pde.f03 module

    ! Update solution using the forward Euler integrator - INTERNAL CELLS only

    U1 = U ! Init
    U2 = U ! Init
    U3 = U ! Init

print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"
print*, "WRONG! NEED TO USE TILDE OPERATOR!!!! WRONGGGGG"

    do i = 5, Nx-4
      call compute_fluxes(U, i, F0_W, F0_E)
      U1(:,i) = U(:,i) - dt/(2.0d0*dx)*(F0_E - F0_W)  ! No source term
    end do
    call assign_BCs(U1) ! See the pde.f03 module

    do i = 5, Nx-4
      call compute_fluxes(U1, i, F1_W, F1_E)
      U2(:,i) = 2.0d0/5.0d0*U(:,i) + 2.0/5.0*dt/dx*(F0_E - F0_W) &
              + 3.0d0/5.0d0*U1(:,i) - 3.0d0/5.0d0*dt/dx*(F1_E - F1_W)  ! No source term
    end do
    call assign_BCs(U2) ! See the pde.f03 module

    do i = 5, Nx-4
      call compute_fluxes(U2, i, F2_W, F2_E)
      U3(:,i) = 831.0d0/20000.0d0*U(:,i) + 1769.0d0/40000.0d0*dt/dx*(F0_E - F0_W) &
              + 4669.0d0/20000.0d0*U1(:,i) + 161.0d0/600.0d0*dt/dx*(F1_E - F1_W)  &
              + 29.0d0/40.0d0*U2(:,i) - 5.0d0/6.0d0*dt/dx*(F2_E - F2_W) ! No source term
    end do
    call assign_BCs(U3) ! See the pde.f03 module
 
    do i = 5, Nx-4
      call compute_fluxes(U3, i, F3_W, F3_E)
      U_new(:,i) = 1.0d0/3.0d0*U(:,i) + 1.0d0/3.0d0*U1(:,1) - 1.0d0/3.0d0*dt/dx*(F1_E - F1_W) &
                 + 1.0d0/3.0d0*U3(:,i) - 1.0d0/6.0d0*dt/dx*(F3_E - F3_W)  ! No source term

      ! Check that the solution did not diverge
      do eqID = 1, Neq
        if (isnan(U_new(eqID,i))) then 
          print*, 'Solution diverged, try with a smaller time step! Aborting.'
          print*, 'Solution that diverged: ', U_new(:,i)
          print*, 'in cell i = ', i
          stop
        end if
      end do

    end do

    ! Save solution (internal cells ONLY! Do not overwrite ghost cells!)
    U(:,5:Nx-4) = U_new(:,5:Nx-4)

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes(U, i, F_W, F_E)

    ! A subroutine to select which fluxes I want to use. Just for simpler
    ! implementation

    implicit none

    real(kind=8), dimension(Neq,Nx), intent(in)  :: U
    real(kind=8), dimension(Neq),    intent(out) :: F_W, F_E
    integer, intent(in) :: i

    if (flux_type_Rusanov .eqv. .true.) then

      call compute_fluxes_Rusanov(U, i, F_W, F_E)

    else if (flux_type_HLL .eqv. .true.) then

      call compute_fluxes_HLL(U, i, F_W, F_E)

    else

      write(*,*) "Attention! No flux type was selected! See global_module file!"
      write(*,*) "ABORTING!"
      stop

    end if

  end subroutine 

  ! ======================================================================== 

  subroutine compute_fluxes_HLL(U, i, F_W, F_E)

    ! Computes HLL numerical fluxes for the cell "i"
    ! 'F_W' and 'F_E' are the complete fluxes at the cell interfaces, to be then
    ! added in the time-discretized PDE.
    ! The variables 'F_L' and 'F_R' are temporary variables used for the 
    ! fluxes at the left and right of an interface.

    implicit none

    real(kind=8), dimension(Neq,Nx), intent(in)  :: U
    real(kind=8), dimension(Neq),    intent(out) :: F_W, F_E
    integer, intent(in) :: i

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_im3, U_im2, U_im1, U_i, U_ip1, U_ip2, U_ip3
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R
    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    U_im3 = U(:,i-3)
    U_im2 = U(:,i-2)
    U_im1 = U(:,i-1)
    U_i   = U(:,i)
    U_ip1 = U(:,i+1)
    U_ip2 = U(:,i+2)
    U_ip3 = U(:,i+3)

    ! ---- West interface -----

    call reconstruct_sol_interface(U_im3, U_im2, U_im1, U_i, U_ip1, U_ip1, U_L, U_R) ! UL and UR are sol across the west interface
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
    ws_max = MAX(ABS(ws_max_R), ABS(ws_min_R), ABS(ws_max_L), ABS(ws_min_L))
    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ---- East interface -----

    call reconstruct_sol_interface(U_im2, U_im1, U_i, U_ip1, U_ip2, U_ip3, U_L, U_R) ! UL and UR are sol across the east interface
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
    ws_max = MAX(ABS(ws_max_R), ABS(ws_min_R), ABS(ws_max_L), ABS(ws_min_L))
    ws_maxabs = MAX(ws_maxabs, ws_max)

  end subroutine

  ! ======================================================================== 

  subroutine compute_fluxes_Rusanov(U, i, F_W, F_E)

    ! Computes HLL numerical fluxes for the cell "i"
    ! 'F_W' and 'F_E' are the complete fluxes at the cell interfaces, to be then
    ! added in the time-discretized PDE.
    ! The variables 'F_L' and 'F_R' are temporary variables used for the 
    ! fluxes at the left and right of an interface.

    implicit none

    real(kind=8), dimension(Neq,Nx), intent(in)  :: U
    real(kind=8), dimension(Neq),    intent(out) :: F_W, F_E
    integer, intent(in) :: i

    integer :: eqID
    real(kind=8), dimension(Neq) :: U_im3, U_im2, U_im1, U_i, U_ip1, U_ip2, U_ip3
    real(kind=8), dimension(Neq) :: U_L, U_R, F_L, F_R
    real(kind=8) :: ws_min_L, ws_max_L, ws_min_R, ws_max_R, ws_min, ws_max ! wave speeds

    U_im3 = U(:,i-3)
    U_im2 = U(:,i-2)
    U_im1 = U(:,i-1)
    U_i   = U(:,i)
    U_ip1 = U(:,i+1)
    U_ip2 = U(:,i+2)
    U_ip3 = U(:,i+3)

    ! ---- West interface -----

    call reconstruct_sol_interface(U_im3, U_im2, U_im1, U_i, U_ip1, U_ip2, U_L, U_R) ! UL and UR are sol across the west interface
    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    ! West interface
    F_W    = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ws_maxabs = MAX(ws_maxabs, ws_max)

    ! ---- East interface -----

    call reconstruct_sol_interface(U_im2, U_im1, U_i, U_ip1, U_ip2, U_ip3, U_L, U_R) ! UL and UR are sol across the west interface
    call compute_flux_ws_x(U_L, F_L, ws_max_L, ws_min_L)
    call compute_flux_ws_x(U_R, F_R, ws_max_R, ws_min_R)

    ws_max = MAX(ABS(ws_max_L), ABS(ws_min_L), ABS(ws_max_R), ABS(ws_min_R))

    ! West interface
    F_E    = 0.5*(F_R + F_L) - ws_max/2.0*(U_R-U_L) ! Rusanov flux

    ws_maxabs = MAX(ws_maxabs, ws_max)

  end subroutine

  ! =========================================================================

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
     
    if (space_order .EQ. 1) then ! ============================================================================
    ! ***** First order solution *****
  
      U_L = U_m1
      U_R = U_p1
 
    else if (space_order .EQ. 2) then ! ============================================================================
      ! ***** Second order solution - MUSCL with TVD slope limiter *****
  
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
        P_R(eID) = P_p1(eID) - phi_limiter(1.0d0/(theta_R+1.0d-10))*(P_p2(eID) - P_p1(eID))/2.0

      end do

      ! Compute conservative variables
      call compute_conserved_from_primitive(P_L, U_L)
      call compute_conserved_from_primitive(P_R, U_R)

    else if (space_order .EQ. 3) then  ! ============================================================================
      ! ***** Third order solution - WENO 3 *****
  
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

    else if (space_order .EQ. 5) then  ! ============================================================================
      ! ***** Third order solution - WENO 5 *****
  
      ! Perform reconstruction in primitive variables 
      call compute_primitive_from_conserved(U_m2, P_m2)
      call compute_primitive_from_conserved(U_m1, P_m1)
      call compute_primitive_from_conserved(U_p1, P_p1)
      call compute_primitive_from_conserved(U_p2, P_p2)

      do eID = 1, Neq ! Loop on the quantities

        call compute_primitive_from_conserved(U_m3, P_m3)
        call compute_primitive_from_conserved(U_p3, P_p3)

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
