module global_module

  implicit none
  
  ! Time integration

  real(kind=8), parameter :: t_end      = 1.0d-0 ! [s] total simulated time (from 0 to t_end)
  real(kind=8), parameter :: CFL_target = 0.25

  !!!! TEST TEST TEST !!! ! Reconstruction order 
  !!!! TEST TEST TEST !!! ! integer, parameter :: reconstr_order = 0 ! 0: no reconstruction -> first order in space
  !!!! TEST TEST TEST !!! integer, parameter :: reconstr_order = 1 ! 1: linear reconstruction -> second order in space

  ! Free stream

  real(kind=8), parameter :: rho0 = 1.225   ! [kg/m3]
  real(kind=8), parameter :: ux0  = 2000.0     ! [m/s]
  real(kind=8), parameter :: uy0  = 0.0  ! [m/s]
  real(kind=8), parameter :: T0   = 300.0   ! [K]

  ! Utilities
  real(kind=8) :: ws_over_sqrtA_maxabs

end module
