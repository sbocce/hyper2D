module pde

  use global_module
  use grid

  implicit none

  ! Adiabatic constant "gamma" (ratio of specific heats, cp/cv)
  real(kind=8), parameter :: gam = 1.4d0
  real(kind=8), parameter :: M   = 6.6337d-26 ! [kg] particles mass
 
  integer, parameter :: Neq = 4 ! Number of equations

  real(kind=8), parameter :: kB  = 1.38d-23 ! [J/K] Boltzmann constant

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'rho','UUx','UUy','TTT'/)

  ! Supported boundary conditions
  real(kind=8), dimension(Neq) :: U_inlet, U_outlet, U_wall, U_sym

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:), intent(inout) :: U

    integer      :: i
    real(kind=8) :: n0, P0

    n0   = rho0/M   ! [particles/m^3] number density
    P0   = n0*kB*T0 ! [Pa] gas pressure

    ! Initialize all cells 
    U(1,:) = rho0     ! Density
    U(2,:) = rho0*ux0 ! Momentum along x
    U(3,:) = rho0*uy0 ! Momentum along y
    U(4,:) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy

    ! Initialize all cells 
    U(1,:) = rho0     ! Density
    U(2,:) = rho0*0.0 ! Momentum along x
    U(3,:) = rho0*0.0 ! Momentum along y
    U(4,:) = rho0*(0.0**2 + 0.0**2)/2.0 + P0/(gam-1.0) ! total energy


    ! do i = 1, Nele
    !   if ( (ele_geom(i,2) .le. 0.4) ) then !  .and. (ele_geom(i,2) .le. 0.5) ) then ! xC
    !     ! if ( (ele_geom(i,3) .ge. 0.3) .and. (ele_geom(i,3) .le. 0.4) ) then ! yC

    !       U(1,i) = 5.0*rho0     ! Density
    !       U(2,i) = rho0*ux0 ! Momentum along x
    !       U(3,i) = rho0*uy0 ! Momentum along y
    !       U(4,i) = rho0*(ux0**2 + uy0**2)/2.0 + 5.0*P0/(gam-1.0) ! total energy

    !     ! end if
    !   end if
    ! end do 

    ! ----------- Also initialize BCs -----------------

    U_inlet(1) = rho0     ! Density
    U_inlet(2) = rho0*ux0 ! Momentum along x
    U_inlet(3) = rho0*uy0 ! Momentum along y
    U_inlet(4) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy

    U_outlet(1) = rho0     ! Density
    U_outlet(2) = rho0*ux0 ! Momentum along x
    U_outlet(3) = rho0*uy0 ! Momentum along y
    U_outlet(4) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy

    ! U_sym and U_wall are assigned during run time

  end subroutine 

!!!! TEST ! !!!!  ! ============================================================
!!!! TEST ! !!!! 
!!!! TEST ! !!!!  subroutine compute_wall_flux(U, nx, ny, F_dot_n, A_ele)
!!!! TEST ! !!!!
!!!! TEST ! !!!!    implicit none
!!!! TEST ! !!!!
!!!! TEST ! !!!!    real(kind=8), dimension(Neq), intent(in)  :: U
!!!! TEST ! !!!!    real(kind=8), dimension(Neq), intent(out) :: F_dot_n
!!!! TEST ! !!!!    real(kind=8), intent(in) :: nx, ny, A_ele
!!!! TEST ! !!!!
!!!! TEST ! !!!!    real(kind=8) :: rho, ux, uy, T, dummy1, dummy2
!!!! TEST ! !!!!    real(kind=8), dimension(Neq) :: U_norm, prim
!!!! TEST ! !!!!
!!!! TEST ! !!!!   ! Wave speeds
!!!! TEST ! !!!!    real(kind=8) :: ws_max, ws_min
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ! Working variables 
!!!! TEST ! !!!!    real(kind=8) :: u_abs, scal_factor
!!!! TEST ! !!!!
!!!! TEST ! !!!!    call compute_primitive_from_conserved(U, prim)
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ux = prim(2)
!!!! TEST ! !!!!    uy = prim(3)
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ! Remove the velocity components normal to the wall 
!!!! TEST ! !!!!    !!! ERROR!!! ux = ux - ux*nx
!!!! TEST ! !!!!    !!! ERROR!!! uy = uy - uy*ny
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ux = ux - (ux*nx + uy*ny)*nx
!!!! TEST ! !!!!    uy = uy - (ux*nx + uy*ny)*ny
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ! Now (ux, uy) is tangential to the wall.
!!!! TEST ! !!!!    ! Build a new state
!!!! TEST ! !!!!    prim(2) = ux
!!!! TEST ! !!!!    prim(3) = uy
!!!! TEST ! !!!!
!!!! TEST ! !!!!    call compute_conserved_from_primitive(prim, U_norm)
!!!! TEST ! !!!!
!!!! TEST ! !!!!    call compute_flux_ws(U_norm, F_dot_n, nx, ny, ws_max, ws_min)
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ! Update global maximum wave speed (used for setting the time step)
!!!! TEST ! !!!!    ws_max = MAX(abs(ws_max), abs(ws_min))
!!!! TEST ! !!!!
!!!! TEST ! !!!!    ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele))
!!!! TEST ! !!!!
!!!! TEST ! !!!!  end subroutine
!!!! TEST ! !!!!
!!!! TEST ! !!!!
!!!! TEST ! !!!!  ! ============================================================
 
  subroutine compute_wall_flux(U, nx, ny, F_dot_n, A_ele)

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: F_dot_n
    real(kind=8), intent(in) :: nx, ny, A_ele

    real(kind=8) :: rho, ux, uy, T, dummy1, dummy2
    real(kind=8), dimension(Neq) :: U_norm, prim

   ! Wave speeds
    real(kind=8) :: ws_max, ws_min

    ! Working variables 
    real(kind=8) :: u_abs, scal_factor

    call compute_primitive_from_conserved(U, prim)

    ux = prim(2)
    uy = prim(3)

    !!!! ! Remove the velocity components normal to the wall 
    !!!! ux = ux - (ux*nx + uy*ny)*nx
    !!!! uy = uy - (ux*nx + uy*ny)*ny

    ! Just say that these are zero. This will create an artificial numerical 
    ! boundary layer. But we are OK with that.
    ux = 0.0
    uy = 0.0

    ! Now (ux, uy) is tangential to the wall.
    ! Build a new state
    prim(2) = ux
    prim(3) = uy

    call compute_conserved_from_primitive(prim, U_norm)

    call compute_flux_ws(U_norm, F_dot_n, nx, ny, ws_max, ws_min)

    ! Update global maximum wave speed (used for setting the time step)
    ws_max = MAX(abs(ws_max), abs(ws_min))

    ws_over_sqrtA_maxabs = MAX(ws_over_sqrtA_maxabs, ws_max/sqrt(A_ele))

  end subroutine

  ! ============================================================

  subroutine compute_sym_state(U, nx, ny, U_sym)

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: U_sym
    real(kind=8), intent(in) :: nx, ny

    real(kind=8), dimension(Neq) :: prim

    real(kind=8) :: ux, uy, ux_norm, uy_norm

    call compute_primitive_from_conserved(U, prim)

    ux = prim(2)
    uy = prim(3)

    ! Mirror the normal velocity component (this is equivalent to removing the 
    ! normal component two times!
    !! ERROR !! ux = ux - 2.0*ux*nx
    !! ERROR !! uy = uy - 2.0*uy*ny
    ux = ux - 2.0*(ux*nx + uy*ny)*nx
    uy = uy - 2.0*(ux*nx + uy*ny)*ny

    ! Compose new state
    prim(2) = ux
    prim(3) = uy

    call compute_conserved_from_primitive(prim, U_sym)

  end subroutine 

  ! ============================================================

  subroutine compute_primitive_from_conserved(U, prim)

    ! Computes vector of primitive variables "prim" from the conserved variables "U"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: prim

    ! Working variables
    real(kind=8) :: rho, ux, uy, P

    ! Extract primitive variables
    rho = U(1)
    ux  = U(2)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
    uy  = U(3)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
    P   = (gam - 1.0)*( U(4) - rho*(ux**2.0 + uy**2.0)/2.0 )

    ! Compose array of primitive variables
    prim(1) = rho
    prim(2) = ux
    prim(3) = uy 
    prim(4) = P/(kB*rho/M) ! T = P/(n*kB) = P/(rho*Ri), with Ri=kB/M the gas constant

  end subroutine

  ! ============================================================

  subroutine compute_conserved_from_primitive(prim, U)

    ! Computes vector of conserved variables "U" from the primitive variables "prim"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: prim
    real(kind=8), dimension(Neq), intent(out) :: U

    ! Working variables
    real(kind=8) :: rho, ux, uy, T, P

    ! Extract primitive variables
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)

    P = rho/M*kB*T

    ! Compose array of conserved variables
    U(1) = rho    ! Density
    U(2) = rho*ux ! Momentum along x
    U(3) = rho*uy ! Momentum along y
    U(4) = rho*(ux**2 + uy**2)/2.0 + P/(gam-1.0) ! total energy

  end subroutine
  
  ! ============================================================

  subroutine compute_flux_ws(U, F_dot_n, nx, ny, ws_max, ws_min)

    ! Computes the convective flux along x,
    ! and also the maximum and minimum wave speeds (required by some numerical flux schemes)

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8),                 intent(in)  :: nx, ny
    real(kind=8), dimension(Neq), intent(out) :: F_dot_n
    real(kind=8),                 intent(out) :: ws_max, ws_min
  
    real(kind=8), dimension(Neq) :: prim, Fx, Fy
    real(kind=8) :: rho, ux, uy, T, P, rhoE, u_dot_n
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)

    P    = rho*kB/M*T ! Compute pressure, P = n*kB*T
    rhoE = rho*(ux*ux + uy*uy)/2.0 + P/(gam-1)

    ! Assemble fluxes Fx
    Fx(1) = rho*ux
    Fx(2) = rho*ux*ux + P
    Fx(3) = rho*ux*uy
    Fx(4) = rhoE*ux + P*ux
  
    ! Assemble fluxes Fy
    Fy(1) = rho*uy
    Fy(2) = rho*ux*uy
    Fy(3) = rho*uy*uy + P
    Fy(4) = rhoE*uy + P*uy
  
    ! Rotate in the direction of the normal
    F_dot_n = Fx*nx + Fy*ny

    ! Maximum and minimum wave speeds, normal to the interface (eigenvalues of the Euler system)
    u_dot_n  = ux*nx + uy*ny
    ws_max   = u_dot_n + sqrt(gam*P/rho)
    ws_min   = u_dot_n - sqrt(gam*P/rho)

  end subroutine

end module
