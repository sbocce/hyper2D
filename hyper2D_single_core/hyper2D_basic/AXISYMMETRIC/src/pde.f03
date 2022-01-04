module pde

! This module contains all information that appears inside the PDEs.
! It specifies how many equations shall be solved, the numerical data,
! the routines to compute the fluxes from the solution, routines for
! computing primitive variables from the conserved solution U, etc.
!
! ----- EULER EQUATIONS -------

  use global_module

  implicit none

  real(kind=8), parameter :: gam = 1.66667d0 ! Adiabatic constant (ratio of specific heats)
 
  integer, parameter :: Neq = 4 ! Number of equations

  ! Initial conditions in the domain
  real(kind=8), parameter :: rho0 = 1.0d0     ! 
  real(kind=8), parameter :: ux0  = 0.0d0     ! 
  real(kind=8), parameter :: uy0  = 0.0d0     ! 
  real(kind=8), parameter :: P0   = 1.0d0     ! 

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'rho','UUx','UUy','PPP'/)

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer      :: i, j
    
    real(kind=8) :: rhoJET, uJET, PJET, TJET

    ! Initialize all cells (i = "1,2" and "Nx-1, Nx" are ghost cells. Same thing for j)
    U(1, :, :) = rho0     ! Density
    U(2, :, :) = rho0*ux0 ! Momentum along x
    U(3, :, :) = rho0*uy0 ! Momentum along y
    U(4, :, :) = rho0*(ux0*ux0 + uy0*uy0)/2.0d0 + P0/(gam-1.0d0) ! total energy

    ! Introduce a jet at the bottom-left of the domain
    uJET   = 0.25d0
    rhoJET = 0.5d0
    PJET   = P0

    U(1, 1:floor(Nx/10.0), 1:floor(Ny/5.0)) = rhoJET      ! Density
    U(2, 1:floor(Nx/10.0), 1:floor(Ny/5.0)) = rhoJET*uJET ! Momentum along x
    U(3, 1:floor(Nx/10.0), 1:floor(Ny/5.0)) = 0.0d0       ! Momentum along y
    U(4, 1:floor(Nx/10.0), 1:floor(Ny/5.0)) = rhoJET*(uJET**2)/2.0 + PJET/(gam-1.0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)

  ! This subroutine assigns boundary values to the ghost cells. 
  ! In this form, we impose two jets, one at the left boundary and one at the bottom boundary.
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    ! Variables for imposing incoming jets as BCs
    integer      :: j_jet_start, j_jet_end, i_jet_start, i_jet_end
    real(kind=8) :: uJET, uOUTLET
    real(kind=8) :: x_jet_start, x_jet_end
    real(kind=8) :: y_jet_start, y_jet_end

!!    ! Bottom BC
!!    U(1,:,1) =  U(1,:,3)
!!    U(2,:,1) =  U(2,:,3)
!!    U(3,:,1) = -U(3,:,3) ! Opposite uy velocity
!!    U(4,:,1) =  U(4,:,3) ! Same energy
!!
!!    U(1,:,2) =  U(1,:,3)
!!    U(2,:,2) =  U(2,:,3)
!!    U(3,:,2) = -U(3,:,3) ! Opposite uy velocity
!!    U(4,:,2) =  U(4,:,3) ! Same energy

!!!     n0   = rho0/M   ! [particles/m^3] number density
!!!     P0   = n0*kB*T0 ! [Pa] gas pressure
!!! 
!!!     ! Put jets at half of the domain +/- 1/6 of the domain
!!!     x_jet_start = x_min + (x_max - x_min)/2.0 - (x_max - x_min)/6.0 
!!!     x_jet_end   = x_min + (x_max - x_min)/2.0 + (x_max - x_min)/6.0
!!! 
!!!     y_jet_start = y_min + (y_max - y_min)/2.0 - (y_max - y_min)/6.0 
!!!     y_jet_end   = y_min + (y_max - y_min)/2.0 + (y_max - y_min)/6.0
!!! 
!!!     i_jet_start = ceiling((x_jet_start - x_min)/(x_max - x_min)*Nx)
!!!     i_jet_end   = ceiling((x_jet_end   - x_min)/(x_max - x_min)*Nx)
!!! 
!!!     j_jet_start = ceiling((y_jet_start - y_min)/(y_max - y_min)*Ny)
!!!     j_jet_end   = ceiling((y_jet_end   - y_min)/(y_max - y_min)*Ny)
!!! 
!!!     uJET    = 600.0d0; ! [m/s] 
!!!     uOUTLET = 1000.0d0 ! [m/s] Let stuff outflow easily (for supersonic flows)
!!! 
!!!     ! ----- LEFT BC (x_min, whatever y) --------
!!! 
!!!     U(1, 1:2, :) = rho0     ! Density
!!!     U(2, 1:2, :) = rho0*ux0 ! Momentum along x
!!!     U(3, 1:2, :) = rho0*uy0 ! Momentum along y
!!!     U(4, 1:2, :) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy
!!! 
!!!     U(1, 1:2, j_jet_start:j_jet_end) = rho0     ! Density
!!!     U(2, 1:2, j_jet_start:j_jet_end) = rho0*uJET ! Momentum along x
!!!     U(3, 1:2, j_jet_start:j_jet_end) = rho0*0.0 ! Momentum along y
!!!     U(4, 1:2, j_jet_start:j_jet_end) = rho0*(uJET**2 + 0.0**2)/2.0 + P0/(gam-1.0) ! total energy
!!!  
!!!     ! ----- RIGHT BC (x_max, whatever y) --------
!!! 
!!!     U(1, Nx-1:Nx, :) = rho0     ! Density
!!!     U(2, Nx-1:Nx, :) = rho0*uOUTLET ! Momentum along x
!!!     U(3, Nx-1:Nx, :) = rho0*0.0 ! Momentum along y
!!!     U(4, Nx-1:Nx, :) = rho0*(uOUTLET**2 + 0.0**2)/2.0 + P0/(gam-1.0) ! total energy
!!! 
!!!     ! ----- BOTTOM BC (y_min, whatever x) --------
!!! 
!!!     U(1, :, 1:2) = rho0      ! Density
!!!     U(2, :, 1:2) = rho0*ux0  ! Momentum along x
!!!     U(3, :, 1:2) = rho0*uy0 ! Momentum along y
!!!     U(4, :, 1:2) = rho0*(ux0**2 + uy0**2)/2.0 + P0/(gam-1.0) ! total energy
!!!  
!!!     U(1, i_jet_start:i_jet_end, 1:2) = rho0      ! Density
!!!     U(2, i_jet_start:i_jet_end, 1:2) = rho0*0.0  ! Momentum along x
!!!     U(3, i_jet_start:i_jet_end, 1:2) = rho0*uJET ! Momentum along y
!!!     U(4, i_jet_start:i_jet_end, 1:2) = rho0*(0.0**2 + uJET**2)/2.0 + P0/(gam-1.0) ! total energy
!!!  
!!!     ! ----- TOP BC (y_max, whatever x) --------
!!! 
!!!     U(1, :, Ny-1:Ny) = rho0         ! Density
!!!     U(2, :, Ny-1:Ny) = rho0*0.0     ! Momentum along x
!!!     U(3, :, Ny-1:Ny) = rho0*uOUTLET ! Momentum along y
!!!     U(4, :, Ny-1:Ny) = rho0*(0.0**2 + uOUTLET**2)/2.0 + P0/(gam-1.0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine compute_conserved_from_primitive(prim, U)

    ! Computes vector of conserved variables "U" from the primitive variables "P"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: prim
    real(kind=8), dimension(Neq), intent(out) :: U

    ! Working variables
    real(kind=8) :: rho, ux, uy, P

    ! Extract primitive variables
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    P   = prim(4)

    ! Compose array of conserved variables
    U(1) = rho
    U(2) = rho*ux
    U(3) = rho*uy 
    U(4) = rho*(ux*ux + uy*uy)/2.0d0 + P/(gam-1.0d0)

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
    P   = (gam - 1.0d0)*( U(4) - rho*(ux*ux + uy*uy)/2.0d0 )

    ! Compose array of primitive variables
    prim(1) = rho
    prim(2) = ux
    prim(3) = uy 
    prim(4) = P

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
    real(kind=8) :: rho, ux, uy, P, rhoE
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    P   = prim(4)

    rhoE = rho*(ux*ux + uy*uy)/2.0d0 + P/(gam-1.0d0)

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
    real(kind=8) :: rho, ux, uy, P, rhoE
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    P   = prim(4)

    rhoE = rho*(ux*ux + uy*uy)/2.0d0 + P/(gam-1.0d0)

    ! Assemble fluxes Fy
    Fy(1) = rho*uy
    Fy(2) = rho*ux*uy
    Fy(3) = rho*uy*uy + P
    Fy(4) = rhoE*uy + P*uy
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = uy + sqrt(gam*P/rho)
    ws_min = uy - sqrt(gam*P/rho)

  end subroutine

  ! ==============================================================

  subroutine compute_flux_y_sym(U, Fy)

    ! Computes the convective flux along y at the symmetry boundary

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: Fy
  
    real(kind=8), dimension(Neq) :: prim
    real(kind=8) :: rho, ux, uy, P, rhoE
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    P   = prim(4)

    ! Impose that uy = 0.0
    uy = 0.0d0

    ! ! Compose fluxes
    ! rhoE = rho*(ux*ux + uy*uy)/2.0d0 + P/(gam-1.0d0)
    ! 
    ! ! Assemble fluxes Fy
    ! Fy(1) = rho*uy
    ! Fy(2) = rho*ux*uy
    ! Fy(3) = rho*uy*uy + P
    ! Fy(4) = rhoE*uy + P*uy
 
    Fy(1) = 0.0d0
    Fy(2) = 0.0d0
    Fy(3) = - P
    Fy(4) = 0.0d0

  end subroutine

  ! ==============================================================

  subroutine compute_source_term(U, i, j, Src)

    ! The source term for the axisymmetric equation is equal to 

    implicit none

    real(kind=8), dimension(Neq,Nx,Ny), intent(in)  :: U
    integer,                            intent(in)  :: i, j
    real(kind=8), dimension(Neq),       intent(out) :: Src
  
    real(kind=8), dimension(Neq) :: prim
    real(kind=8) :: rho, ux, uy, P, rhoE
    real(kind=8) :: y_j
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U(:,i,j), prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    P   = prim(4)

    ! Compose fluxes
    rhoE = rho*(ux*ux + uy*uy)/2.0 + P/(gam-1)

    ! Assemble source term Src = Fy/y
    ! Src(1) = rho*uy
    ! Src(2) = rho*ux*uy
    ! Src(3) = rho*uy*uy + P
    ! Src(4) = rhoE*uy + P*uy

    y_j = y_min + dy/2.0d0 + real(j-3)*dy
    ! Src = Src/y_j

    Src = 0.0d0
    Src(3) = P/y_j
  
  end subroutine


end module
