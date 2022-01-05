module pde

  use global_module

  implicit none

  integer, parameter :: Neq = 3 ! Number of equations

  real(kind=8), parameter :: gam    = 1.6667       ! adiabatic constant

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'rho','UUx','PPP'/)

  contains

  ! ============================================================

  subroutine initialize_solution(U)
 
    implicit none

    real(kind=8), dimension(:,:), intent(inout) :: U

    integer      :: i

    real(kind=8) :: rhoL, uxL, PL
    real(kind=8) :: rhoR, uxR, PR

    ! Left state
    rhoL = 2.0
    uxL  = 0.0
    PL   = 2.0

    ! Right state
    rhoR = 1.0
    uxR  = 0.0
    PR   = 1.0

    ! Initialize solution (all cells, also internal ones)
    U(1,1:floor(Nx/2.0)) = rhoL                             ! Density 
    U(2,1:floor(Nx/2.0)) = rhoL*uxL                         ! Momentum
    U(3,1:floor(Nx/2.0)) = rhoL*uxL*uxL/2.0 + PL/(gam-1.0)    ! Energy

    U(1,floor(Nx/2.0):Nx) = rhoR
    U(2,floor(Nx/2.0):Nx) = rhoR*uxR
    U(3,floor(Nx/2.0):Nx) = rhoR*uxR*uxR/2.0 + PR/(gam-1.0)

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)
 
    implicit none

    real(kind=8), dimension(:,:), intent(inout) :: U

    integer :: i

    ! ----- LEFT BC --------

    ! Do nothing!

    ! U(1, 1:2) = ...
    ! U(2, 1:2) = ...
    ! U(3, 1:2) = ...

    ! ----- RIGHT BC --------

    ! Do nothing!

    ! U(1, Nx-1:Nx) = ...
    ! U(2, Nx-1:Nx) = ...
    ! U(3, Nx-1:Nx) = ...

  end subroutine 

  ! ============================================================

  subroutine compute_conserved_from_primitive(prim, U)

    ! Computes vector of conserved variables "U" from the primitive variables "prim"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: prim
    real(kind=8), dimension(Neq), intent(out) :: U

    ! Working variables
    real(kind=8) :: rho, ux, P

    ! Extract primitive variables
    rho = prim(1)
    ux  = prim(2)
    P   = prim(3)

    ! Compose array of primitive variables
    U(1) = rho
    U(2) = rho*ux
    U(3) = rho*ux*ux/2.0d0 + P/(gam - 1.0d0)

  end subroutine
  
  ! ============================================================

  subroutine compute_primitive_from_conserved(U, prim)

    ! Computes vector of primitive variables "prim" from the conserved variables "U"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: prim

    ! Working variables
    real(kind=8) :: rho, ux, P

    ! Extract primitive variables
    rho = U(1)
    ux  = U(2)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
    P   = (gam-1.0)*(U(3) - rho*ux*ux/2.0)

    ! Compose array of primitive variables
    prim(1) = rho
    prim(2) = ux
    prim(3) = P

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
    real(kind=8) :: rho, ux, P
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    P   = prim(3)

    ! Assemble fluxes Fx
    Fx(1) = rho*ux
    Fx(2) = rho*ux*ux + P
    Fx(3) = U(3)*ux + P*ux ! note: U(3) = rho*E 
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    ws_max = ux + sqrt(gam*P/rho)
    ws_min = ux - sqrt(gam*P/rho)

  end subroutine

end module
