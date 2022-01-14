module global_module

! This module contains data needed throughout the program. 
! Most data is specified as a parameter, as to make the execution faster.
! This module is loaded by most other modules.

  implicit none
  
  ! Time integration

  ! For Euler ! real(kind=8), parameter :: dt    = 1.0d-6 ! [s] time step for the simulation
  ! For Euler ! real(kind=8), parameter :: t_end = 1.0d-2 ! [s] total simulated time (from 0 to t_end)

  real(kind=8), parameter :: dt    = 1.0d-4 ! [s] time step for the simulation
  real(kind=8), parameter :: t_end = 10.0d0 ! [s] total simulated time (from 0 to t_end)

  ! Specify domain and discretization

  real(kind=8), parameter :: x_min = -1.0 ! [m] 
  real(kind=8), parameter :: x_max = 2.0  ! [m] 
  real(kind=8), parameter :: y_min = -1.0 ! [m] 
  real(kind=8), parameter :: y_max = 1.0  ! [m] 

  ! Total number of cells (including 2 ghost cells for each side)
  ! (using 2 ghost cells makes it easy to reach second order accuracy)

  integer, parameter :: Nx = 500
  integer, parameter :: Ny = 500

  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-4)
  real(kind=8), parameter :: dy = (y_max-y_min)/(Ny-4)

  ! Spacial accuracy (order)
  ! 1: uniform sol in cell, first order in space
  ! 2: MUSCL with TVD slope limiter, second order in space
  ! 3: WENO 3rd order
  ! 5: WENO 5th order
  integer, parameter :: space_order = 5

  ! Time accuracy (order)
  ! 1: forward Euler, first order in time
  ! 2: midpoint Euler, second order in time
  ! 3: explicit TVD Runge-Kutta, third order in time
  integer, parameter :: time_order = 3

  ! Space fluxes: choose ONLY one
  logical, parameter :: flux_type_Rusanov = .true.

end module
