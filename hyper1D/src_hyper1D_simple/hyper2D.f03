! ######################################################################
! ######################################################################
! ###########     _                           ____  ____      ##########   
! ###########    | |__  _   _ _ __   ___ _ __(___ \|  _ \     ##########
! ###########    | '_ \| | | | '_ \ / _ \ '__| __) | | | |    ##########
! ###########    | | | | |_| | |_) |  __/ |   / __/| |_| |    ##########
! ###########    |_| |_|\__, | .__/ \___|_|  /_____|____/     ##########
! ###########           |___/|_|                              ##########
! ###########                           1D version            ########## 
! ###########                                                 ########## 
! ######################################################################
! ######################################################################

program hyper2D

  use global_module   ! Simulation parameters: domain size, number of cells etc
  use pde             ! Definition of the system of equations
  use integration     ! Functions for integrating in time and numerical fluxes
  use tools           ! Output subroutines etc etc

  implicit none

  ! The solution is a 2D matrix. 
  ! The first index "eqID" represents the equation 
  ! (for Euler 1: density, 2: x-momentum, 3: total energy)
  ! the second index "i" represents the x-position.
  ! This might be counter-intuitive, but recall that Fortran represents data 
  ! in column-major order.

  real(kind=8), dimension(Neq,Nx) :: U, P

  integer      :: t_ID 
  real(kind=8) :: t_now

  write(*,*) "Initializing solution..."
  call initialize_solution(U) ! See the pde.f03 module

  write(*,*) "Assigning BCs..."
  call assign_BCs(U) ! See the pde.f03 module
  
  write(*,*) "Writing solution at time step", 0, "..."
  call export_sol(0, 0.0d0, U)

  ! $$$$$$$$$$$ Integrate in time $$$$$$$$$$$$$

  t_now = 0.0d0   ! Init
 
  do t_ID = 1, Nt

    t_now = t_now + dt

    ! ------ Integrate by dt ------
    call forward_Euler_step(U, dt)
   
    ! ----- Write solution to VTK, every ... time steps -----
    if ( mod(t_ID, 100) .EQ. 0 ) then 
      write(*,*) "Writing solution at time step", t_ID, "..."
      call export_sol(t_ID, t_now, U)
    end if

  end do

  ! Export very last timestep    
  call export_sol(t_ID, t_now, U)

end program
