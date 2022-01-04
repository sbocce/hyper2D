module global_module

! This module contains data needed throughout the program. 
! Most data is specified as a parameter, as to make the execution faster.
! This module is loaded by most other modules.

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: dt    = 1.0d-3 ! [s] time step for the simulation
  real(kind=8), parameter :: t_end = 500.0d0 ! [s] total simulated time (from 0 to t_end)

  ! Specify domain and discretization

  real(kind=8), parameter :: x_min = 0.0d0 ! [m] 
  real(kind=8), parameter :: x_max = 3.0d0  ! [m] 
  real(kind=8), parameter :: y_min = 0.0d0  ! [m] MUST BE POSITIVE FOR AXISYMMETRIC SIMULATION
  real(kind=8), parameter :: y_max = 1.0d0  ! [m] 

  ! Total number of cells (including 2 ghost cells for each side)
  ! (using 2 ghost cells makes it easy to reach second order accuracy)

  integer, parameter :: Nx = 300
  integer, parameter :: Ny = 150

  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-4)
  real(kind=8), parameter :: dy = (y_max-y_min)/(Ny-4)

  ! Spatial accuracy
  logical, parameter :: bool_MUSCL = .TRUE.

end module
