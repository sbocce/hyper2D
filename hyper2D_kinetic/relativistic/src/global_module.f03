module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: dt    = 1.0d-4 ! [s] time step for the simulation
  real(kind=8), parameter :: t_end = 0.5d-0  ! [s] total simulated time (from 0 to t_end)

  real(kind=8), parameter :: c = 299792458.0  ! [m/s] speed of light in vacuum

  ! Specify domain and discretization

  real(kind=8), parameter :: x_min = -c/10.0 ! [m]
  real(kind=8), parameter :: x_max = c   ! [m] 
  real(kind=8), parameter :: y_min = -0.5*9.109d-31*c     ! [kg*m/s] THIS IS A MOMENTUM
  real(kind=8), parameter :: y_max = 4.0*9.109d-31*c      ! [kg*m/s] THIS IS A MOMENTUM

  ! Total number of cells (including 2 ghost cells for each side)
  ! (using 2 ghost cells makes it easy to reach second order accuracy)

  integer, parameter :: Nx = 256   ! Physical space
  integer, parameter :: Ny = 256*4 ! Momentum space

  ! Periodicity and second order
  logical, parameter :: x_periodic_bool   = .False.
  logical, parameter :: second_order_bool = .True.

  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-4) ! discretization in physical space
  real(kind=8), parameter :: dy = (y_max-y_min)/(Ny-4) ! discretization in velocity space

!   ! Initial 1D1V maxwellian 
!   real(kind=8), parameter :: n0 = 1.0d16     ! [1/m3] number density
!   real(kind=8), parameter :: u0 = 0.0d0      ! [m/s]  average velocity
!   real(kind=8), parameter :: T0 = 116045.0d0  ! [K]    temperature

end module
