module pde

  use global_module

  implicit none

  real(kind=8), parameter :: M   = 5.3136d-26   ! [kg] molecular mass 
  real(kind=8), parameter :: th_vib = 2256.0d0  ! [K]  vibrational temperature
 
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

    ! Initialize internal cells (i = "1,2" and "Nx-1, Nx" are ghost cells. Same thing for j)
    U(1, 3:Nx-2, 3:Ny-2) = rho0     ! Density
    U(2, 3:Nx-2, 3:Ny-2) = rho0*ux0 ! Momentum along x
    U(3, 3:Nx-2, 3:Ny-2) = rho0*uy0 ! Momentum along y
    U(4, 3:Nx-2, 3:Ny-2) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine assign_BCs(U)
 
    implicit none

    real(kind=8), dimension(:,:,:), intent(inout) :: U

    integer :: i, j

    ! ----- LEFT BC (x_min, whatever y) --------

    U(1, 1:2, :) = rho0     ! Density
    U(2, 1:2, :) = rho0*ux0 ! Momentum along x
    U(3, 1:2, :) = rho0*uy0 ! Momentum along y
    U(4, 1:2, :) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0) ! total energy
 
    ! ----- RIGHT BC (x_max, whatever y) --------

    U(1, Nx-1:Nx, :) = rho0     ! Density
    U(2, Nx-1:Nx, :) = rho0*ux0 ! Momentum along x
    U(3, Nx-1:Nx, :) = rho0*uy0 ! Momentum along y
    U(4, Nx-1:Nx, :) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0) ! total energy

    ! ----- BOTTOM BC (y_min, whatever x) --------

    U(1, :, 1:2) = rho0     ! Density
    U(2, :, 1:2) = rho0*ux0 ! Momentum along x
    U(3, :, 1:2) = rho0*uy0 ! Momentum along y
    U(4, :, 1:2) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0) ! total energy
 
    ! ----- TOP BC (y_max, whatever x) --------

    U(1, :, Ny-1:Ny) = rho0     ! Density
    U(2, :, Ny-1:Ny) = rho0*ux0 ! Momentum along x
    U(3, :, Ny-1:Ny) = rho0*uy0 ! Momentum along y
    U(4, :, Ny-1:Ny) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0) ! total energy

  end subroutine 

  ! ============================================================

  subroutine compute_primitive_from_conserved(U, prim)

    ! Computes vector of primitive variables "prim" from the conserved variables "U"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: prim

    ! Working variables
    real(kind=8) :: rho, ux, uy, T, e_int

    ! Extract primitive variables
    rho = U(1)
    ux  = U(2)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
    uy  = U(3)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero

    e_int = U(4)/(rho + 1.0d-25) - (ux*ux + uy*uy)/2.0 ! Internal energy [J/kg]

    call compute_T_from_e_int(e_int, T)

    ! Compose array of primitive variables
    prim(1) = rho
    prim(2) = ux
    prim(3) = uy 
    prim(4) = T

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
    real(kind=8) :: rho, ux, uy, T, P, rhoE, e_int, gam
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)

    P     = rho*kB/M*T ! Compute pressure, P = n*kB*T
    e_int = compute_e_int(T) ! [J/kg]
    rhoE  = rho*(ux*ux + uy*uy)/2.0 + rho*e_int

    ! Assemble fluxes Fx
    Fx(1) = rho*ux
    Fx(2) = rho*ux*ux + P
    Fx(3) = rho*ux*uy
    Fx(4) = rhoE*ux + P*ux
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    gam = (compute_cv_int(T) + kB/M)/compute_cv_int(T)

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
    real(kind=8) :: rho, ux, uy, T, P, rhoE, e_int, gam
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)

    P     = rho*kB/M*T ! Compute pressure, P = n*kB*T
    e_int = compute_e_int(T) ! [J/kg]
    rhoE  = rho*(ux*ux + uy*uy)/2.0 + rho*e_int

    ! Assemble fluxes Fy
    Fy(1) = rho*uy
    Fy(2) = rho*ux*uy
    Fy(3) = rho*uy*uy + P
    Fy(4) = rhoE*uy + P*uy
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)
    gam = (compute_cv_int(T) + kB/M)/compute_cv_int(T)

    ws_max = uy + sqrt(gam*P/rho)
    ws_min = uy - sqrt(gam*P/rho)

  end subroutine

  ! =================================================================

  function compute_e_int(T_in)

    ! Internal energy for diatomic molecule with rotational and vibrational DOFs (RRHO model)

    implicit none

    real(kind=8), intent(in) :: T_in

    real(kind=8) :: compute_e_int
    real(kind=8) :: T

    T = abs(T_in + 0.001d0) ! [K] add a small tolerance, and make it positive (to prevent integration errors)
    compute_e_int = kB/M*(5.0/2.0*T + th_vib/(exp(th_vib/T) - 1.0)) ! [J/kg]
 
    return 

  end function

  ! =================================================================

  function compute_cv_int(T_in)

    ! Specific heat for diatomic molecule with rotational and vibrational DOFs (RRHO model)

    implicit none

    real(kind=8), intent(in) :: T_in

    real(kind=8) :: compute_cv_int
    real(kind=8) :: T
    
    T = abs(T_in + 0.001d0) ! [K] add a small tolerance, and make it positive (to prevent integration errors)
    compute_cv_int = kB/M*(5.0/2.0 + (th_vib/T)**2.0*exp(th_vib/T)/(exp(th_vib/T) - 1.0)**2.0 ) ! [J/(kg K)]
 
    return 

  end function

  ! =================================================================
 
  subroutine compute_T_from_e_int(e_int, T)

    ! This subroutine finds the temperature T that gives the internal energy e_int,
    ! using Newton iterations (cv_int is the derivative of e_int).

    implicit none

    real(kind=8), intent(in)  :: e_int
    real(kind=8), intent(out) :: T

    real(kind=8) :: error ! Error in the temperature
    real(kind=8) :: T_new

    ! ---- NEWTON METHOD ----
    ! Initialize T using the rigid rotor model without vibration

    T = 2.0/5.0*e_int*M/kB

    error = 1.0d0

    do while ( error .gt. 0.0000001d0 ) 
      T_new = T - (compute_e_int(T) - e_int)/(compute_cv_int(T) + 0.0001d0) 
      error = abs(T_new - T)
      T     = T_new
    end do

    !!! ! ---- BISECTION METHOD ----
    !!! ! This should be more stable but much heavier than Newton.
    !!! ! The problem is simple enough that we can just use Newton.
    !!! ! --------------------------
    !!! real(kind=8) :: Ta, Tb, Thalf
    !!! real(kind=8) :: e_int_a, e_int_b, e_int_half
    !!! Ta    = 100.0
    !!! Tb    = 1000000.0    
    !!! error = 100000.0
    !!! do while ( error .gt. 0.001 )
    !!!   Thalf = (Ta + Tb)/2.0
    !!!   e_int_a = compute_e_int(Ta)
    !!!   e_int_b = compute_e_int(Tb)
    !!!   e_int_half = compute_e_int(Thalf)
    !!!   error = abs(Tb - Ta)
    !!!   if (e_int_half .lt. e_int) then 
    !!!     Ta = Thalf
    !!!   else 
    !!!     Tb = Thalf
    !!!   end if
    !!! end do
    !!! T = Thalf

  end subroutine

end module
