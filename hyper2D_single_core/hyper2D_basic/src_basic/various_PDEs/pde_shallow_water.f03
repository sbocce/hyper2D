module pde

! This module contains all information that appears inside the PDEs.
! It specifies how many equations shall be solved, the numerical data,
! the routines to compute the fluxes from the solution, routines for
! computing primitive variables from the conserved solution U, etc.
!
! ---------- SHALLOW WATER EQUATIONS ----------

  use global_module

  implicit none

  integer, parameter      :: Neq = 3    ! Number of equations
  real(kind=8), parameter :: g   = 9.81 ! [N/kg] Earth's gravity

  ! Initial conditions
  real(kind=8), parameter :: h0   = 1.0     ! [m]
  real(kind=8), parameter :: ux0  = 0.0     ! [m/s]
  real(kind=8), parameter :: uy0  = 0.0     ! [m/s]

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using two letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'hh','Ux','Uy'/)

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    ! Initialize internal cells (i = "1,2" and "Nx-1, Nx" are ghost cells. Same thing for j)
    U(1, 3:Nx-2, 3:Ny-2) = h0     ! Density
    U(2, 3:Nx-2, 3:Ny-2) = h0*ux0 ! Momentum along x
    U(3, 3:Nx-2, 3:Ny-2) = h0*uy0 ! Momentum along y


    U(1, floor((Nx)/3.0):floor(real(Nx)/2.0), floor((Ny)/3.0):floor(real(Ny)/2.0)) = 2.0*h0
    U(2, floor((Nx)/3.0):floor(real(Nx)/2.0), floor((Ny)/3.0):floor(real(Ny)/2.0)) = 2.0*h0*ux0
    U(3, floor((Nx)/3.0):floor(real(Nx)/2.0), floor((Ny)/3.0):floor(real(Ny)/2.0)) = 2.0*h0*uy0

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)

  ! This subroutine assigns boundary values to the ghost cells. 
  ! In this form, we impose periodic BCs, by copying into the ghost cells the values inside the 
  ! domain cells at the other side of the domain.
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    ! Periodic BCs
    U(:, 1:2, :)     = U(:, Nx-3:Nx-2, :) ! Left BC
    U(:, Nx-1:Nx, :) = U(:, 3:4, :)       ! Right BC
    U(:, :, 1:2)     = U(:, :, Ny-3:Ny-2) ! Bottom BC
    U(:, :, Ny-1:Ny) = U(:, :, 3:4)       ! Top BC

  end subroutine 

  ! ============================================================

  subroutine compute_primitive_from_conserved(U, prim)

    ! Computes vector of primitive variables "prim" from the conserved variables "U"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: prim

    ! Working variables
    real(kind=8) :: h, ux, uy

    ! Extract primitive variables
    h  = U(1)
    ux = U(2)/(h + 1.0d-25) ! Use a small tolerance, since h may be zero
    uy = U(3)/(h + 1.0d-25) ! Use a small tolerance, since h may be zero

    ! Compose array of primitive variables
    prim(1) = h
    prim(2) = ux
    prim(3) = uy 

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
    real(kind=8) :: h, ux, uy
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    h  = prim(1)
    ux = prim(2)
    uy = prim(3)

    ! Assemble fluxes Fx
    Fx(1) = h*ux
    Fx(2) = h*ux*ux + 0.5*g*h**2
    Fx(3) = h*ux*uy
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = ux + sqrt(g*h)
    ws_min = ux - sqrt(g*h)

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
    real(kind=8) :: h, ux, uy
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    h  = prim(1)
    ux = prim(2)
    uy = prim(3)

    ! Assemble fluxes Fx
    Fy(1) = h*uy
    Fy(2) = h*ux*uy
    Fy(3) = h*uy*uy + 0.5*g*h*h
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = uy + sqrt(g*h)
    ws_min = uy - sqrt(g*h)

  end subroutine

end module
