module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: dt  = 0.0005d0 ! [s] time step
  integer,      parameter :: Nt  = 500 ! Number of time steps to be performed

  real(kind=8), parameter :: x_min = -0.5d0 ! [m] beginning of computational domain
  real(kind=8), parameter :: x_max = 0.5d0  ! [m] end of computational domain

  ! Total number of cells (including 4 ghost cells for each side)
  ! (using 4 ghost cells makes it easy to reach higher order accuracy)

  integer, parameter :: Nx = 200

  logical, parameter :: bool_PERIODIC = .false. ! Boundary conditions

  ! Note: we use here 4 ghost cells per side. This allows for a good generalization 
  ! for higher order cases.
  real(kind=8), parameter :: dx = (x_max-x_min)/(Nx-8) 

end module
