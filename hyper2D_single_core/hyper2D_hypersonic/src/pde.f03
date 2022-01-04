module pde

  use global_module

  implicit none

  integer, parameter :: Neq = 5 ! Number of equations

  real(kind=8), parameter :: M      = 5.3136d-26   ! [kg] molecular mass 
  real(kind=8), parameter :: gam    = 1.4          ! adiabatic constant (diatomic molecule, roto-translation only)
  real(kind=8), parameter :: th_vib = 2256.0d0     ! [K]  vibrational temperature
  real(kind=8), parameter :: a_MW   = 138.0d0      ! O2, parameter for the Millikan-White model
  real(kind=8), parameter :: b_MW   = 0.03d0       ! O2, parameter for the Millikan-White model

  real(kind=8), parameter :: kB  = 1.38d-23              ! [J/K] Boltzmann constant
  real(kind=8), parameter :: PI  = 4.0d0*DATAN(1.0d0) 

  ! Name of primitive variables, used ONLY for exporting the solution to VTK file
  ! NOTE: To initialize it here, all entries need to have the same length!!!
  !       I'm using three letters for simplicity.
  character(len=20), dimension(Neq) :: prim_names = (/'rho','UUx','UUy','Ttr','Tvv'/)

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
    U(4, 3:Nx-2, 3:Ny-2) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0,T0) ! total energy
    U(5, 3:Nx-2, 3:Ny-2) = rho0*compute_e_vib(T0) ! vibrational energy

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
    U(4, 1:2, :) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0,T0) ! total energy
    U(5, 1:2, :) = rho0*compute_e_vib(T0) ! vibrational energy
 
    ! ----- RIGHT BC (x_max, whatever y) --------

    U(1, Nx-1:Nx, :) = rho0     ! Density
    U(2, Nx-1:Nx, :) = rho0*ux0 ! Momentum along x
    U(3, Nx-1:Nx, :) = rho0*uy0 ! Momentum along y
    U(4, Nx-1:Nx, :) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0,T0) ! total energy
    U(5, Nx-1:Nx, :) = rho0*compute_e_vib(T0) ! vibrational energy

    ! ----- BOTTOM BC (y_min, whatever x) --------

    U(1, :, 1:2) = rho0     ! Density
    U(2, :, 1:2) = rho0*ux0 ! Momentum along x
    U(3, :, 1:2) = rho0*uy0 ! Momentum along y
    U(4, :, 1:2) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0,T0) ! total energy
    U(5, :, 1:2) = rho0*compute_e_vib(T0) ! vibrational energy
 
    ! ----- TOP BC (y_max, whatever x) --------

    U(1, :, Ny-1:Ny) = rho0     ! Density
    U(2, :, Ny-1:Ny) = rho0*ux0 ! Momentum along x
    U(3, :, Ny-1:Ny) = rho0*uy0 ! Momentum along y
    U(4, :, Ny-1:Ny) = rho0*(ux0**2 + uy0**2)/2.0 + rho0*compute_e_int(T0,T0) ! total energy
    U(5, :, Ny-1:Ny) = rho0*compute_e_vib(T0) ! vibrational energy

  end subroutine 

  ! ============================================================

  subroutine compute_primitive_from_conserved(U, prim)

    ! Computes vector of primitive variables "prim" from the conserved variables "U"

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: prim

    ! Working variables
    real(kind=8) :: rho, ux, uy, T, Tv, e_int, e_vib, e_rt

    ! Extract primitive variables
    rho = U(1)
    ux  = U(2)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero
    uy  = U(3)/(rho + 1.0d-25) ! Use a small tolerance, since rho may be zero

    e_int = U(4)/(rho + 1.0d-25) - (ux*ux + uy*uy)/2.0 ! Internal energy [J/kg]
    e_vib = U(5)/(rho + 1.0d-25) ! Vibrational energy [J/kg]

    e_rt  = e_int - e_vib  ! Roto-translational energy [J/kg]

    T  = 2.0/5.0*e_rt*M/kB ! Roto-translational temperature
    Tv = th_vib/log(kB*th_vib/(M*e_vib)+1) ! Vibrational temperature

    ! Compose array of primitive variables
    prim(1) = rho
    prim(2) = ux
    prim(3) = uy 
    prim(4) = T
    prim(5) = Tv

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
    real(kind=8) :: rho, ux, uy, T, Tv, P, rhoE, e_rt, e_vib
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)
    Tv  = prim(5)

    P     = rho*kB/M*T ! Compute pressure, P = n*kB*T

    e_rt  = compute_e_rt(T)   ! [J/kg]
    e_vib = compute_e_vib(Tv) ! [J/kg]

    rhoE  = rho*(ux*ux + uy*uy)/2.0 + rho*(e_rt + e_vib) ! Total energy

    ! Assemble fluxes Fx
    Fx(1) = rho*ux
    Fx(2) = rho*ux*ux + P
    Fx(3) = rho*ux*uy
    Fx(4) = rhoE*ux + P*ux
    Fx(5) = rho*e_vib*ux
  
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
    real(kind=8) :: rho, ux, uy, T, Tv, P, rhoE, e_rt, e_vib
    
    ! Compute primitive variables
    call compute_primitive_from_conserved(U, prim)
    rho = prim(1)
    ux  = prim(2)
    uy  = prim(3)
    T   = prim(4)
    Tv  = prim(5)

    P     = rho*kB/M*T ! Compute pressure, P = n*kB*T

    e_rt  = compute_e_rt(T)   ! [J/kg]
    e_vib = compute_e_vib(Tv) ! [J/kg]

    rhoE  = rho*(ux*ux + uy*uy)/2.0 + rho*(e_rt + e_vib)

    ! Assemble fluxes Fy
    Fy(1) = rho*uy
    Fy(2) = rho*ux*uy
    Fy(3) = rho*uy*uy + P
    Fy(4) = rhoE*uy + P*uy
    Fy(5) = rho*e_vib*uy
  
    ! Maximum and minimum wave speeds (eigenvalues of the Euler system)

    ws_max = uy + sqrt(gam*P/rho)
    ws_min = uy - sqrt(gam*P/rho)

  end subroutine

  ! =================================================================

  function compute_e_int(T_tr, T_vib)

    ! Internal energy for diatomic molecule with rotational and vibrational DOFs (RRHO model)

    implicit none

    real(kind=8), intent(in) :: T_tr, T_vib

    real(kind=8) :: compute_e_int
 
    compute_e_int = compute_e_rt(T_tr) + compute_e_vib(T_vib) ! [J/kg]
 
    return 

  end function

  ! =================================================================

  function compute_e_rt(T_tr)

    ! Roto-vibrational energy for diatomic molecule (Rigid Rotor model)

    implicit none

    real(kind=8), intent(in) :: T_tr

    real(kind=8) :: compute_e_rt

    compute_e_rt = kB/M*5.0/2.0*T_tr ! [J/kg]
 
    return 

  end function

  ! =================================================================

  function compute_e_vib(T_vib)

    ! Vibrational energy for diatomic molecule (RRHO model)

    implicit none

    real(kind=8), intent(in) :: T_vib

    real(kind=8) :: compute_e_vib
    real(kind=8) :: T

    T = abs(T_vib + 0.00001d0) ! [K] add a small tolerance, and make it positive (to prevent integration errors)
    compute_e_vib = kB/M*th_vib/(exp(th_vib/T) - 1.0) ! [J/kg]
 
    return 

  end function

  ! =================================================================

  function compute_cv_int(T_tr, T_vib)

    ! Specific heat for diatomic molecule with rotational and vibrational DOFs (RRHO model)

    implicit none

    real(kind=8), intent(in) :: T_tr, T_vib

    real(kind=8) :: compute_cv_int

    real(kind=8) :: T_vib_tol
    
    T_vib_TOL = abs(T_vib + 0.001d0) ! [K] add a small tolerance, and make it positive (to prevent integration errors)
    compute_cv_int = kB/M*(5.0/2.0 + &
    (th_vib/T_vib_tol)**2.0*exp(th_vib/T_vib_tol)/(exp(th_vib/T_vib_tol) - 1.0)**2.0 ) ! [J/(kg K)]
 
    return 

  end function

  ! =================================================================
 
  subroutine compute_Teq_from_e_int(e_int, Teq)

    ! This subroutine finds the equilibrium temperature T that gives the internal energy e_int,
    ! using Newton iterations (cv_int is the derivative of e_int).

    implicit none

    real(kind=8), intent(in)  :: e_int
    real(kind=8), intent(out) :: Teq

    real(kind=8) :: error ! Error in the temperature
    real(kind=8) :: T_new

    ! ---- NEWTON METHOD ----
    ! Initialize Teq using the rigid rotor model without vibration

    Teq = 2.0/5.0*e_int*M/kB

    error = 1.0d0

    do while ( error .gt. 0.0000001d0 ) 
      T_new = Teq - (compute_e_int(Teq,Teq) - e_int)/(compute_cv_int(Teq,Teq) + 0.0001d0) 
      error = abs(T_new - Teq)
      Teq   = T_new
    end do

  end subroutine

  ! =======================================================

  subroutine compute_tau_vib(n, T, tau_v)

    ! Computes the vibrational relaxation time.

    implicit none

    real(kind=8), intent(in)  :: n, T  ! Input:  number density [1/m^3], and temperature [K]
    real(kind=8), intent(out) :: tau_v    ! Output: vibraitonal relaxation time [s]

    real(kind=8) :: P
    real(kind=8) :: tau_MW, tau_P, SIG

    P = n*kB*T

    ! Compute relaxation time from Millikan & White
    tau_MW = 1.0/(P*101325.0)*exp(a_MW*(T**(-1.0/3.0) - b_MW) - 18.42)

    ! Compute Park correction
    SIG    = 3.0d-21*(50000.0/T)**2 ! Estimation of the vibrational relaxation cross-section
    tau_P  = 1.0/(n*sqrt(8.0*kB*T/PI/M)*SIG) 

    ! Vibrational relaxation time
    tau_v = tau_MW + tau_P

  end subroutine

  ! =======================================================

  subroutine compute_source_term(U, S, tau_v, e_vib_eq, e_vib)

    ! Computes collisional relaxation source term

    implicit none

    real(kind=8), dimension(Neq), intent(in)  :: U
    real(kind=8), dimension(Neq), intent(out) :: S
    real(kind=8), intent(out) :: tau_v, e_vib_eq, e_vib

    real(kind=8), dimension(Neq) :: prim

    real(kind=8) :: n, T_tr, T_vib, T_eq, e_int

    ! Compute primitive variables 
    call compute_primitive_from_conserved(U, prim)

    n     = prim(1)/M ! number density n = rho/M
    T_tr  = prim(4)
    T_vib = prim(5)

    ! Vibrational and internal energy at this step
    e_vib = compute_e_vib(T_vib)
    e_int = compute_e_rt(T_tr) + e_vib

    call compute_tau_vib(n, T_tr, tau_v) ! Relaxation rate

    call compute_Teq_from_e_int(e_int, T_eq) ! Compute equilibrium temperature

    e_vib_eq = compute_e_vib(T_eq) ! Vibrational energy at equilibrium

    ! Assemble source term
    S = 0.0 ! Init

    S(5) = prim(1)*(e_vib_eq - e_vib)/tau_v

    ! print*, "tau_vib, e_vib_eq, e_vib: ", tau_v, e_vib_eq, e_vib

  end subroutine

end module
