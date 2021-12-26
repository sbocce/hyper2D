module pde

  use global_module

  implicit none

  ! Adiabatic constant "gamma" (ratio of specific heats, cp/cv)
  real(kind=8), parameter :: gam = 1.66667
  real(kind=8), parameter :: M   = 6.6337d-26 ! [kg] particles mass
 
  integer, parameter :: Neq = 4 ! Number of equations

  real(kind=8), parameter :: kB  = 1.38d-23 ! [J/K] Boltzmann constant

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'rho','UUx','UUy','TTT'/)

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer      :: i, j
    real(kind=8) :: n0, P0

    n0   = rho0/M   ! [particles/m^3] number density
    P0   = n0*kB*T0 ! [Pa] gas pressure

    ! Initialize internal cells (i = "1,2" and "Nx-1, Nx" are ghost cells. Same thing for j)
    U(1, 3:Nx-2, 3:Ny-2) = rho0     ! Density
    U(2, 3:Nx-2, 3:Ny-2) = rho0*ux0 ! Momentum along x
    U(3, 3:Nx-2, 3:Ny-2) = rho0*uy0 ! Momentum along y
    U(4, 3:Nx-2, 3:Ny-2) = rho0*(ux0**2.0 + uy0**2.0)/2.0 + P0/(gam-1.0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    real(kind=8) :: n0, P0

    n0   = rho0/M   ! [particles/m^3] number density
    P0   = n0*kB*T0 ! [Pa] gas pressure

    ! ----- LEFT BC (x_min, whatever y) --------

    U(1, 1:2, :) = rho0     ! Density
    U(2, 1:2, :) = rho0*ux0 ! Momentum along x
    U(3, 1:2, :) = rho0*uy0 ! Momentum along y
    U(4, 1:2, :) = rho0*(ux0**2.0 + uy0**2.0)/2.0 + P0/(gam-1.0) ! total energy
 
    ! ----- RIGHT BC (x_max, whatever y) --------

    U(1, Nx-1:Nx, :) = rho0     ! Density
    U(2, Nx-1:Nx, :) = rho0*ux0 ! Momentum along x
    U(3, Nx-1:Nx, :) = rho0*uy0 ! Momentum along y
    U(4, Nx-1:Nx, :) = rho0*(ux0**2.0 + uy0**2.0)/2.0 + P0/(gam-1.0) ! total energy

    ! ----- BOTTOM BC (y_min, whatever x) --------

    U(1, :, 1:2) = rho0     ! Density
    U(2, :, 1:2) = rho0*ux0 ! Momentum along x
    U(3, :, 1:2) = rho0*uy0 ! Momentum along y
    U(4, :, 1:2) = rho0*(ux0**2.0 + uy0**2.0)/2.0 + P0/(gam-1.0) ! total energy
 
    ! ----- TOP BC (y_max, whatever x) --------

    U(1, :, Ny-1:Ny) = rho0     ! Density
    U(2, :, Ny-1:Ny) = rho0*ux0 ! Momentum along x
    U(3, :, Ny-1:Ny) = rho0*uy0 ! Momentum along y
    U(4, :, Ny-1:Ny) = rho0*(ux0**2.0 + uy0**2.0)/2.0 + P0/(gam-1.0) ! total energy

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

  subroutine compute_flux_ws_x(U, Fx, ws_max, ws_min)

    ! Computes the convective flux along x,
    ! and also the maximum and minimum wave speeds (required by some numerical flux schemes)

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: Fx
    real(kind=8), intent(out) :: ws_max, ws_min
  
    real(kind=8), dimension(Neq) :: prim
    real(kind=8) :: rho, ux, uy, T, P, rhoE
    
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
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = ux + sqrt(gam*P/rho)
    ws_min = ux - sqrt(gam*P/rho)

  end subroutine

  ! ============================================================

  subroutine compute_flux_ws_y(U, Fy, ws_max, ws_min)

    ! Computes the convective flux along y,
    ! and also the maximum and minimum wave speeds (required by some numerical flux schemes)

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: Fy
    real(kind=8), intent(out) :: ws_max, ws_min
  
    real(kind=8), dimension(Neq) :: prim
    real(kind=8) :: rho, ux, uy, T, P, rhoE
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)

    P    = rho*kB/M*T ! Compute pressure, P = n*kB*T
    rhoE = rho*(ux*ux + uy*uy)/2.0 + P/(gam-1)

    ! Assemble fluxes Fy
    Fy(1) = rho*uy
    Fy(2) = rho*ux*uy
    Fy(3) = rho*uy*uy + P
    Fy(4) = rhoE*uy + P*uy
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = uy + sqrt(gam*P/rho)
    ws_min = uy - sqrt(gam*P/rho)

  end subroutine


end module
