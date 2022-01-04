module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: dt    = 5.0d-10 ! [s] time step for the simulation
  real(kind=8), parameter :: t_end = 5.0d-6  ! [s] total simulated time (from 0 to t_end)

  ! Specify domain and discretization

  real(kind=8), parameter :: x_min = -0.0 ! [m]
  real(kind=8), parameter :: x_max = 0.01  ! [m] 
  real(kind=8), parameter :: y_min = -25000.0 ! [m/s] THIS IS A VELOCITY IN THE KINETIC EQ!
  real(kind=8), parameter :: y_max = 25000.0  ! [m/s] THIS IS A VELOCITY IN THE KINETIC EQ!

  ! Total number of cells (including 2 ghost cells for each side)
  ! (using 2 ghost cells makes it easy to reach second order accuracy)

  integer, parameter :: Nx = 2*256 ! Physical space
  integer, parameter :: Ny = 2*256 ! Velocity space

  ! Periodicity and second order
  logical, parameter :: x_periodic_bool   = .True.
  logical, parameter :: second_order_bool = .True.

  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-4) ! discretization in physical space
  real(kind=8), parameter :: dy = (y_max-y_min)/(Ny-4) ! discretization in velocity space

  ! Initial 1D1V maxwellian 
  real(kind=8), parameter :: n0 = 1.0d16     ! [1/m3] number density
  real(kind=8), parameter :: u0 = 0.0d0      ! [m/s]  average velocity
  real(kind=8), parameter :: T0 = 116045.0d0  ! [K]    temperature

end module
