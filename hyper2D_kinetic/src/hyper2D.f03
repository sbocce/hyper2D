! ######################################################################
! ######################################################################
! ###########     _                           ____  ____      ##########   
! ###########    | |__  _   _ _ __   ___ _ __(___ \|  _ \     ##########
! ###########    | '_ \| | | | '_ \ / _ \ '__| __) | | | |    ##########
! ###########    | | | | |_| | |_) |  __/ |   / __/| |_| |    ##########
! ###########    |_| |_|\__, | .__/ \___|_|  /_____|____/     ##########
! ###########           |___/|_|                              ##########
! ###########                                                 ########## 
! ######################################################################
! ######################################################################

program hyper2D

  use global_module   ! Simulation parameters: domain size, number of cells etc
  use pde             ! Definition of the system of equations
  use integration     ! Functions for integrating in time and numerical fluxes
  use tools           ! Output subroutines etc etc

  implicit none

  ! The solution is a 3D matrix. 
  ! The first index "eqID" represents the equation 
  ! (for Euler 1: density, 2: x-momentum, 3: y-momentum, 4: total energy)
  ! the second index "i" represents the x-position 
  ! the third index "j" represents the y-position.
  ! This might be counter-intuitive, but recall that Fortran represents data 
  ! in column-major order.

  real(kind=8), dimension(Neq,Nx,Ny) :: U, P

  integer      :: Nt, t_ID ! Variables for time integration
  real(kind=8) :: t_now

  real(kind=8) :: dummy1, dummy2, dummy3

  write(*,*) "Initializing solution..."
  call initialize_solution(U) ! See the pde.f03 module

  write(*,*) "Assigning BCs..."
  call assign_BCs(U) ! See the pde.f03 module
  
  write(*,*) "Writing solution at time step", 0, "..."
  call export_sol_vtk(0, U)

  ! $$$$$$$$$$$ Integrate in time $$$$$$$$$$$$$

  Nt = ceiling(t_end/dt) ! t_end and dt are defined in global_module.f90

  do t_ID = 1, Nt

    t_now = t_ID*dt
    write(*,*) 'Timestep', t_ID, 'of', Nt, ' - Current time:', t_now, '[s]'

    ! call forward_Euler_step(U, dt)
    call midpoint_Euler_step(U, dt)

    if ( mod(t_ID, 50) .EQ. 0 ) then ! Write to VTK every ... timesteps

      write(*,*) "Writing solution at time step", t_ID, "..."
      call export_sol_vtk(t_ID, U)

    end if

  end do

end program
