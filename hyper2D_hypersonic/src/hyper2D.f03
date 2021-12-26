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
  real(kind=8) :: dt, t_now, CFL_now

  write(*,*) "Initializing solution..."
  call initialize_solution(U) ! See the pde.f03 module

  write(*,*) "Assigning BCs..."
  call assign_BCs(U) ! See the pde.f03 module
  
  write(*,*) "Writing solution at time step", 0, "..."
  call export_sol_vtk(0, U)

  ! $$$$$$$$$$$ Integrate in time $$$$$$$$$$$$$

  ! Nt = ceiling(t_end/dt) ! t_end and dt are defined in global_module.f90
  ! do t_ID = 1, Nt

  t_now = 0.0d0   ! Init
  dt    = 1.0d-10 ! Very first time step
  t_ID  = 0.0     ! Init
 
  do while( t_now .le. t_end)

    ! ------ Prepare variables ------
    t_ID      = t_ID + 1
    t_now     = t_now + dt
    ws_maxabs = 0.0 ! Global variable, init to zero

    ! ------ Integrate by dt ------
    call forward_Euler_step(U, dt)

    ! ----- Write solution to VTK, every ... time steps -----
    if ( mod(t_ID, 100) .EQ. 0 ) then 
      write(*,*) "Writing solution at time step", t_ID, "..."
      call export_sol_vtk(t_ID, U)
    end if

    ! ------ Estimate current Courant number and update time step -----
    CFL_now = ws_maxabs*dt/MIN(dx,dy)
    write(*,'(A,EN15.5,A,F10.5,A,ES14.7,A)') 'Time', t_now, ' [s]. Current CFL: ', CFL_now, '. dt = ', dt, '[s]'
    dt      = dt*CFL_target/CFL_now

  end do

end program
