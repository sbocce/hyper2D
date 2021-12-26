module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: t_end      = 1.0d-1 ! [s] total simulated time (from 0 to t_end)
  real(kind=8), parameter :: CFL_target = 0.25

  ! Specify domain and discretization. This version includes a flat plate in the domain,
  ! therefore the domain must contain it.

  real(kind=8), parameter :: x_min = -1.75 ! [m] needs to be NEGATIVE
  real(kind=8), parameter :: x_max = 2.75  ! [m] needs to be LARGER THAN PLATE LENGTH L_plate
  real(kind=8), parameter :: y_min = -0.7 ! [m] needs to be NEGATIVE
  real(kind=8), parameter :: y_max = 2.1  ! [m] needs to be POSITIVE

  ! Total number of cells (including 2 ghost cells for each side)
  ! (using 2 ghost cells makes it easy to reach second order accuracy)

  integer, parameter :: Nx = 100
  integer, parameter :: Ny = 100

  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-4)
  real(kind=8), parameter :: dy = (y_max-y_min)/(Ny-4)

  ! Flat plate inside the domain.
  ! The plate is placed with the leading edge at (x,y) = (0,0). 
  ! Make sure the domain is chosen correspondingly.

  real(kind=8), parameter :: L_plate = 1 ! [m] length of the flat plate

  integer, parameter :: j_plate   = floor(Ny/(y_max - y_min)*(0 - y_min)) ! find j-index of plate
  integer, parameter :: i_1_plate = floor(Nx/(x_max - x_min)*(0 - x_min)) ! starting coordinate of plate
  integer, parameter :: i_2_plate = floor(Nx/(x_max - x_min)*(L_plate - x_min)) ! starting coordinate of plate

  ! Free stream

  real(kind=8), parameter :: rho0 = 1.225d0   ! [kg/m3]
  real(kind=8), parameter :: ux0  = 0.0     ! [m/s]
  real(kind=8), parameter :: uy0  = 2867.2  ! [m/s]
  real(kind=8), parameter :: T0   = 300.0   ! [K]

  ! Utilities
  real(kind=8) :: ws_maxabs

end module
