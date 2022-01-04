module pde

! This module contains all information that appears inside the PDEs.
! It specifies how many equations shall be solved, the numerical data,
! the routines to compute the fluxes from the solution, routines for
! computing primitive variables from the conserved solution U, etc.
!
! ----- EULER EQUATIONS -------

  use global_module

  implicit none

  ! Adiabatic constant "gamma" (ratio of specific heats, cp/cv)
  real(kind=8), parameter :: gam = 1.66667
  real(kind=8), parameter :: M   = 6.6337d-26 ! [kg] particles mass
 
  integer, parameter :: Neq = 4 ! Number of equations

  real(kind=8), parameter :: kB  = 1.38d-23 ! [J/K] Boltzmann constant

  ! Initial conditions in the domain
  real(kind=8), parameter :: rho0 = 1.225   ! [kg/m3]
  real(kind=8), parameter :: ux0  = 0.0     ! [m/s]
  real(kind=8), parameter :: uy0  = 0.0     ! [m/s]
  real(kind=8), parameter :: T0   = 300.0   ! [K]

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
    U(4, 3:Nx-2, 3:Ny-2) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)

  ! This subroutine assigns boundary values to the ghost cells. 
  ! In this form, we impose two jets, one at the left boundary and one at the bottom boundary.
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    real(kind=8) :: n0, P0
    
    ! Variables for imposing incoming jets as BCs
    integer      :: j_jet_start, j_jet_end, i_jet_start, i_jet_end
    real(kind=8) :: uJET, uOUTLET
    real(kind=8) :: x_jet_start, x_jet_end
    real(kind=8) :: y_jet_start, y_jet_end

    n0   = rho0/M   ! [particles/m^3] number density
    P0   = n0*kB*T0 ! [Pa] gas pressure

    ! Put jets at half of the domain +/- 1/6 of the domain
    x_jet_start = x_min + (x_max - x_min)/2.0 - (x_max - x_min)/6.0 
    x_jet_end   = x_min + (x_max - x_min)/2.0 + (x_max - x_min)/6.0

    y_jet_start = y_min + (y_max - y_min)/2.0 - (y_max - y_min)/6.0 
    y_jet_end   = y_min + (y_max - y_min)/2.0 + (y_max - y_min)/6.0

    i_jet_start = ceiling((x_jet_start - x_min)/(x_max - x_min)*Nx)
    i_jet_end   = ceiling((x_jet_end   - x_min)/(x_max - x_min)*Nx)

    j_jet_start = ceiling((y_jet_start - y_min)/(y_max - y_min)*Ny)
    j_jet_end   = ceiling((y_jet_end   - y_min)/(y_max - y_min)*Ny)

    uJET    = 600.0d0; ! [m/s] 
    uOUTLET = 1000.0d0 ! [m/s] Let stuff outflow easily (for supersonic flows)

    ! ----- LEFT BC (x_min, whatever y) --------

    U(1, 1:2, :) = rho0     ! Density
    U(2, 1:2, :) = rho0*ux0 ! Momentum along x
    U(3, 1:2, :) = rho0*uy0 ! Momentum along y
    U(4, 1:2, :) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy

    U(1, 1:2, j_jet_start:j_jet_end) = rho0     ! Density
    U(2, 1:2, j_jet_start:j_jet_end) = rho0*uJET ! Momentum along x
    U(3, 1:2, j_jet_start:j_jet_end) = rho0*0.0 ! Momentum along y
    U(4, 1:2, j_jet_start:j_jet_end) = rho0*(uJET**2 + 0.0**2)/2.0 + P0/(gam-1.0) ! total energy
 
    ! ----- RIGHT BC (x_max, whatever y) --------

    U(1, Nx-1:Nx, :) = rho0     ! Density
    U(2, Nx-1:Nx, :) = rho0*uOUTLET ! Momentum along x
    U(3, Nx-1:Nx, :) = rho0*0.0 ! Momentum along y
    U(4, Nx-1:Nx, :) = rho0*(uOUTLET**2 + 0.0**2)/2.0 + P0/(gam-1.0) ! total energy

    ! ----- BOTTOM BC (y_min, whatever x) --------

    U(1, :, 1:2) = rho0      ! Density
    U(2, :, 1:2) = rho0*ux0  ! Momentum along x
    U(3, :, 1:2) = rho0*uy0 ! Momentum along y
    U(4, :, 1:2) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy
 
    U(1, i_jet_start:i_jet_end, 1:2) = rho0      ! Density
    U(2, i_jet_start:i_jet_end, 1:2) = rho0*0.0  ! Momentum along x
    U(3, i_jet_start:i_jet_end, 1:2) = rho0*uJET ! Momentum along y
    U(4, i_jet_start:i_jet_end, 1:2) = rho0*(0.0**2 + uJET**2)/2.0 + P0/(gam-1.0) ! total energy
 
    ! ----- TOP BC (y_max, whatever x) --------

    U(1, :, Ny-1:Ny) = rho0         ! Density
    U(2, :, Ny-1:Ny) = rho0*0.0     ! Momentum along x
    U(3, :, Ny-1:Ny) = rho0*uOUTLET ! Momentum along y
    U(4, :, Ny-1:Ny) = rho0*(0.0**2 + uOUTLET**2)/2.0 + P0/(gam-1.0) ! total energy

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
