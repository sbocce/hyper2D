module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: t_end      = 0.2d0 ! [s] total simulated time (from 0 to t_end)
  real(kind=8), parameter :: CFL_target = 0.2d0

  real(kind=8), parameter :: x_min = -0.5d0 ! [m] 
  real(kind=8), parameter :: x_max = 0.5d0  ! [m] 

  ! Total number of cells (including 4 ghost cells for each side)
  ! (using 4 ghost cells makes it easy to reach higher order accuracy)

  integer, parameter :: Nx = 200

  logical, parameter :: bool_PERIODIC = .false.

  ! Note: we use 4 ghost cells per side. This allows to reach up to 5th order
  ! (with the WENO method)
  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-8) 

  ! Order in space
  ! 1: uniform sol in cell, first order in space
  ! 2: MUSCL with TVD slope limiter, second order in space
  ! 3: WENO 3rd order
  ! 5: WENO 5th order
  integer, parameter :: space_order = 3

  ! Order in time
  ! 1: forward Euler, first order in time
  ! 2: midpoint Euler, second order in time
  ! 3: explicit TVD Runge-Kutta, third order in time
  integer, parameter :: time_order = 3

  ! Space fluxes: choose ONLY one
  logical, parameter :: flux_type_Rusanov = .false.
  logical, parameter :: flux_type_HLL     = .true.

  ! Utilities
  real(kind=8) :: ws_maxabs ! Global variable, maximum wave speed

end module
