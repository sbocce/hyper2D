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

  integer      :: Nt, t_ID ! Variables for time integration
  real(kind=8) :: dt, t_now, CFL_now

  write(*,*) "Initializing solution..."
  call initialize_solution(U) ! See the pde.f03 module

  write(*,*) "Assigning BCs..."
  call assign_BCs(U) ! See the pde.f03 module
  
  write(*,*) "Writing solution at time step", 0, "..."
  call export_sol(0, 0.0d0, U)

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
    if (time_order .eq. 1) then
      call forward_Euler_step(U, dt)
    else if (time_order .eq. 2) then
      call midpoint_Euler_step(U, dt)
    else if (time_order .eq. 3) then
      call RK3_step(U, dt)
    else
      write(*,*) "Attention! Order of the time marching scheme unknown! Check global_module file!"
      write(*,*) "Aborting!"
      stop
    end if

    ! ----- Write solution to VTK, every ... time steps -----
    if ( mod(t_ID, 100) .EQ. 0 ) then 
      write(*,*) "Writing solution at time step", t_ID, "..."
      call export_sol(t_ID, t_now, U)
    end if

    ! ------ Estimate current Courant number and update time step -----
    CFL_now = ws_maxabs*dt/dx
    write(*,'(A,EN15.5,A,F10.5,A,ES14.7,A)') 'Time', t_now, ' [s]. Current CFL: ', CFL_now, '. dt = ', dt, '[s]'
    dt      = dt*CFL_target/CFL_now

  end do

  ! Export very last timestep    
  call export_sol(t_ID, t_now, U)

end program
